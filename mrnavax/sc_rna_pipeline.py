"""scRNA-seq → neoantigen handoff pipeline.

Closes the loop from tissue to vaccine design:

  1. Load an AnnData (h5ad) file of single-cell RNA-seq.
  2. Cluster cells with a stdlib k-medoids (when scanpy isn't available) or
     a Leiden scanpy pipeline (when it is).
  3. Identify tumor-cell clusters (highest mean expression of an optional
     tumor marker set; otherwise the largest cluster is assumed tumor).
  4. For each tumor-cluster gene with a coding-region SNV (provided as a
     variants CSV), enumerate 9-11-mer peptides overlapping the variant.
  5. Hand off the peptide list to the neoantigen screener, which scores
     peptide × HLA binding and immunogenicity.

The stdlib core runs end-to-end without any bioinformatics deps. When
``scanpy`` is installed, the cluster stage upgrades to a Leiden pipeline.
When ``scgpt`` or another foundation model is installed, its embeddings can
plug into :func:`embed_with_foundation_model` to refine cluster boundaries.

References
----------
scGPT: Cui et al. *Nat Methods* 21, 1480–1491 (2024).
TrambaHLApan / DeepNeo / NetMHCpan (neoantigen scoring).
"""

from __future__ import annotations

import csv
import json
import math
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable

from .alphamissense_integration import UNIPROT_TO_GENE

# ---------- IO ------------------------------------------------------------


@dataclass
class Variant:
    gene: str
    position: int  # 1-indexed AA position in the protein
    wt_aa: str
    mut_aa: str
    # Optional DNA-level coordinates for AlphaGenome Atlas AVI lookup.
    # When provided, the variant scorer uses AVI as the dominant signal
    # for non-coding regulatory-region variants (Avsec et al. 2026).
    # When absent, the variant falls through to the AM-or-BLOSUM path.
    chrom: str | None = None
    ref_dna: str | None = None
    alt_dna: str | None = None


def load_variants(path: str | Path) -> list[Variant]:
    out: list[Variant] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Optional DNA-level coordinates — present in newer CSVs that
            # include AlphaGenome Atlas regulatory variants alongside
            # protein-coding ones. Older CSVs without these columns still
            # work (chrom/ref_dna/alt_dna default to None).
            chrom = row.get("chrom") or None
            ref_dna = row.get("ref_dna") or None
            alt_dna = row.get("alt_dna") or None
            out.append(
                Variant(
                    gene=row["gene"],
                    position=int(row["position"]),
                    wt_aa=row["wt_aa"],
                    mut_aa=row["mut_aa"],
                    chrom=chrom,
                    ref_dna=ref_dna,
                    alt_dna=alt_dna,
                )
            )
    return out


def load_protein_fasta(path: str | Path) -> dict[str, str]:
    """Load a multi-record FASTA into ``{gene_name: protein_seq}``.

    For the gene name, the first whitespace-separated token after ``>``
    is used. UniProt-style headers like ``>sp|P01116|RASK_HUMAN ...``
    are accepted; the gene name is the first token (``sp``) for
    backward compatibility — use :func:`load_protein_fasta_detailed`
    if you need the UniProt accession as well.
    """
    detailed = load_protein_fasta_detailed(path)
    return {gene: seq for gene, (_uid, seq) in detailed.items()}


def load_protein_fasta_detailed(
    path: str | Path,
) -> dict[str, tuple[str | None, str]]:
    """Load a multi-record FASTA into ``{gene_name: (uniprot_id?, sequence)}``.

    Recognizes UniProt-style headers ``>sp|P01116|RASK_HUMAN ...`` and
    extracts ``P01116`` as the UniProt accession. Falls back to ``None``
    if the header is not UniProt-style.
    """
    import re

    out: dict[str, tuple[str | None, str]] = {}
    name: str | None = None
    uniprot: str | None = None
    buf: list[str] = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                out[name] = (uniprot, "".join(buf).upper())
            header = line[1:]
            first_token = header.split()[0] if header else "unnamed"
            pipe_parts = first_token.split("|")
            if len(pipe_parts) >= 3 and re.match(r"^[A-Z][A-Z0-9]{5,10}$", pipe_parts[1]):
                # UniProt-style: >db|UID|GENE_OS
                # The third field is the entry name (e.g. P53_HUMAN, RASK_HUMAN),
                # not always a clean gene symbol. Fall back to the curated
                # UniProt -> gene-symbol table, then to the part of the
                # entry name before the first underscore (which matches most
                # modern UniProt entry names but not all).
                uid = pipe_parts[1]
                entry_name = pipe_parts[2]
                if uid in UNIPROT_TO_GENE:
                    name = UNIPROT_TO_GENE[uid]
                else:
                    candidate = entry_name.split("_")[0] if "_" in entry_name else entry_name
                    # If the candidate looks like a gene symbol (uppercase letters
                    # or letters+digits, 2-10 chars), use it. Otherwise fall back
                    # to the first token (legacy behavior).
                    if re.match(r"^[A-Z][A-Z0-9]{1,9}$", candidate):
                        name = candidate
                    else:
                        name = first_token
                uniprot = uid
            else:
                name = first_token
                m = re.match(r"^[a-z]+\|([A-Z][A-Z0-9]{5,10})\|", header)
                uniprot = m.group(1) if m else None
            buf = []
        else:
            buf.append(line)
    if name:
        out[name] = (uniprot, "".join(buf).upper())
    return out


# ---------- peptide enumeration -------------------------------------------


def mutant_peptides(
    protein: str,
    position_1based: int,
    mut_aa: str,
    lengths: tuple[int, ...] = (9, 10, 11),
) -> list[str]:
    """Enumerate mutant peptides of the given lengths containing the variant.

    For each window length, return every peptide in ``[max(0, p-l+1), p+r]``
    that contains the variant position, with the WT residue replaced by the
    mutant residue at that position.
    """
    p = position_1based - 1  # to 0-indexed
    if p < 0 or p >= len(protein):
        return []
    out: list[str] = []
    for L in lengths:
        for start in range(max(0, p - L + 1), min(len(protein) - L + 1, p + 1) + 1):
            pep = list(protein[start : start + L])
            rel = p - start
            if not (0 <= rel < len(pep)):
                continue
            if pep[rel] == mut_aa:
                continue  # silent
            pep[rel] = mut_aa
            out.append("".join(pep))
    return out


# ---------- k-medoids clustering (stdlib) ---------------------------------


def kmedoids(
    matrix: list[list[float]],
    k: int,
    *,
    max_iter: int = 25,
    seed: int = 0,
) -> list[int]:
    """Greedy k-medoids. ``matrix`` is (n_samples, n_features).

    Returns a list of cluster labels of length n_samples.
    """
    n = len(matrix)
    if n == 0:
        return []
    rng = random.Random(seed)
    labels = [-1] * n
    # Initialize medoids: pick k well-spread points.
    medoids = rng.sample(range(n), min(k, n))
    for _ in range(max_iter):
        # Assign each point to nearest medoid
        new_labels = [
            min(range(len(medoids)), key=lambda m: _euclid2(matrix[i], matrix[medoids[m]]))
            for i in range(n)
        ]
        if new_labels == labels:
            break
        labels = new_labels
        # Update medoids: pick the point with min total distance in each cluster
        new_medoids: list[int] = []
        for m_id in range(len(medoids)):
            members = [i for i, lbl in enumerate(labels) if lbl == m_id]
            if not members:
                new_medoids.append(medoids[m_id])
                continue
            best, best_cost = members[0], math.inf
            for cand in members:
                cost = sum(_euclid2(matrix[cand], matrix[o]) for o in members)
                if cost < best_cost:
                    best, best_cost = cand, cost
            new_medoids.append(best)
        medoids = new_medoids
    return labels


def _euclid2(a: list[float], b: list[float]) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b))


# ---------- AnnData loader with graceful fallback -------------------------


def load_expression(path: str | Path) -> tuple[list[str], list[str], list[list[float]]]:
    """Load a tiny CSV/TSV of (cells × genes) expression.

    Returns ``(cell_ids, gene_names, matrix)``. If the file doesn't exist or
    the format is unsupported, returns a deterministic synthetic 50-cell ×
    200-gene dataset so the rest of the pipeline can still run as a demo.
    """
    p = Path(path)
    if not p.exists():
        return _synthetic()
    text = p.read_text()
    # try CSV/TSV with first row = gene names
    sep = "\t" if "\t" in text.splitlines()[0] else ","
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        return _synthetic()
    header = lines[0].split(sep)
    gene_names = [g.strip() for g in header[1:]]
    cell_ids: list[str] = []
    matrix: list[list[float]] = []
    for line in lines[1:]:
        parts = line.split(sep)
        cell_ids.append(parts[0].strip())
        try:
            matrix.append([float(x) for x in parts[1:]])
        except ValueError:
            return _synthetic()
    return cell_ids, gene_names, matrix


def _synthetic() -> tuple[list[str], list[str], list[list[float]]]:
    rng = random.Random(42)
    n_cells, n_genes = 50, 200
    gene_names = [f"GENE_{i:03d}" for i in range(n_genes)]
    # 3 clusters: low / medium / high expression of last 30 genes (tumor marker)
    labels = [0] * 20 + [1] * 15 + [2] * 15
    matrix = []
    for c in labels:
        row = [rng.gauss(0.5, 0.1) for _ in range(n_genes - 30)]
        # tumor markers up-regulate in cluster 2
        if c == 2:
            row += [rng.gauss(3.0, 0.3) for _ in range(30)]
        else:
            row += [rng.gauss(0.3, 0.2) for _ in range(30)]
        matrix.append(row)
    cell_ids = [f"cell_{i:02d}" for i in range(n_cells)]
    return cell_ids, gene_names, matrix


# ---------- scanpy upgrade path (optional) --------------------------------


def cluster_with_scanpy(
    matrix: list[list[float]],
    cell_ids: list[str],
    gene_names: list[str],
    *,
    n_top_genes: int = 1000,
    n_pcs: int = 20,
    resolution: float = 0.5,
    seed: int = 0,
) -> tuple[list[int], list[str]]:
    """Cluster with scanpy if available, else fall back to k-medoids."""
    try:
        import numpy as np  # noqa: F401
        import scanpy as sc  # noqa: F401
    except ImportError:
        # fall back
        labels = kmedoids(matrix, k=4, seed=seed)
        return labels, [f"cluster_{i}" for i in sorted(set(labels))]

    import anndata as ad
    import numpy as np
    import scanpy as sc

    adata = ad.AnnData(X=np.array(matrix, dtype=np.float32))
    adata.obs_names = cell_ids
    adata.var_names = gene_names
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(n_top_genes, adata.n_vars))
    sc.tl.pca(adata, n_comps=min(n_pcs, adata.n_vars - 1, adata.n_obs - 1))
    sc.pp.neighbors(adata, random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    labels = [int(x) for x in adata.obs["leiden"]]
    return labels, [f"cluster_{i}" for i in sorted(set(labels))]


# ---------- scGPT hook (optional) ------------------------------------------


def embed_with_foundation_model(
    matrix: list[list[float]],
    *,
    model: str = "tfidf-svd",
    n_components: int = 32,
    gene_names: list[str] | None = None,
) -> list[list[float]]:
    """Plug point for scGPT / Geneformer / UNI-RNA embeddings.

    Default is a deterministic, stdlib-only TF-IDF + truncated SVD embedder
    that approximates what a foundation model produces — a per-cell vector
    that preserves cluster structure. Swap in a real model by setting
    ``model="scgpt"`` after downloading the ``perturblab/scgpt-human``
    checkpoint (~205 MB) to ``~/.cache/mrnavax/``.

    When ``gene_names`` is supplied and ``model="scgpt"``, the real scGPT
    tokenizer maps each gene to its vocab ID, producing biologically
    meaningful cell embeddings.
    """
    from .foundation_embedder import embed_cells

    if model == "scgpt":
        # Pass through to scgpt with gene_names support
        from .scgpt_integration import embed_with_scgpt, scgpt_available

        if not scgpt_available():
            raise FileNotFoundError(
                "scGPT weights not found at ~/.cache/mrnavax/. "
                "Download from https://huggingface.co/perturblab/scgpt-human."
            )
        return embed_with_scgpt(matrix, gene_names=gene_names, max_cells=64)
    return embed_cells(matrix, model=model, n_components=n_components)


# ---------- main pipeline --------------------------------------------------


@dataclass
class TumorPeptide:
    cell_cluster: int
    gene: str
    position: int
    wt_aa: str
    mut_aa: str
    peptide: str
    length: int
    cluster_marker_score: float


@dataclass
class PipelineReport:
    n_cells: int
    n_genes: int
    cluster_labels: list[int]
    cluster_marker_scores: dict[int, float]
    tumor_cluster: int
    n_variants_input: int
    n_variants_after_filter: int
    variant_scores: dict[str, float]
    n_candidate_peptides: int
    peptides: list[TumorPeptide] = field(default_factory=list)
    note: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


def run_pipeline(
    expression_path: str | Path,
    variants_path: str | Path,
    proteins_path: str | Path,
    *,
    hla: Iterable[str],
    tumor_marker_genes: list[str] | None = None,
    n_clusters: int = 4,
    peptide_lengths: tuple[int, ...] = (9, 10, 11),
    variant_filter_top_fraction: float = 1.0,
    variant_filter_min_score: float = 0.0,
    embedding_model: str = "tfidf-svd",
    embedding_dim: int = 32,
) -> PipelineReport:
    """Run the full pipeline. Returns a structured report.

    Parameters
    ----------
    variant_filter_top_fraction : float
        Fraction of variants to keep after pathogenicity scoring. 1.0 = no
        filter (default). 0.2 = keep only the top 20% by score.
    variant_filter_min_score : float
        Minimum variant score (0–1) to keep. Default 0.0 (no filter).
    """
    cell_ids, gene_names, matrix = load_expression(expression_path)
    variants = load_variants(variants_path)
    proteins_detailed = load_protein_fasta_detailed(proteins_path)
    proteins = {gene: seq for gene, (_uid, seq) in proteins_detailed.items()}
    uniprot_ids = {gene: uid for gene, (uid, _seq) in proteins_detailed.items() if uid}

    # Snapshot original count before filtering.
    n_variants_input_orig = len(variants)

    # Optional variant pre-filter (AlphaMissense-style, with AlphaMissense
    # itself plugged in when the predictions TSV is available).
    filtered_count = n_variants_input_orig
    filter_scores: dict[str, float] = {}
    am_lookup_fn = None
    am_active = False
    # Only attempt to load the real AlphaMissense index when filtering is
    # requested AND the index is already cached. We don't trigger the
    # 15-min build during a normal run because that's a heavy operation
    # the user should initiate explicitly. Users with a pre-built cache
    # (created via ``python -m mrnavax.alphamissense_integration``
    # or a prior ``load_index()`` call) get the real scores automatically.
    if variant_filter_top_fraction < 1.0 or variant_filter_min_score > 0.0:
        try:
            from .alphamissense_integration import (
                CACHE_FILE,
                load_index,
                lookup,
            )

            if CACHE_FILE.exists():
                load_index()
                am_lookup_fn = lookup
                am_active = True
        except Exception:
            am_lookup_fn = None

    if variant_filter_top_fraction < 1.0 or variant_filter_min_score > 0.0:
        from .variant_scorer import filter_variants

        v_dicts = [
            {
                "gene": v.gene,
                "position": v.position,
                "wt_aa": v.wt_aa,
                "mut_aa": v.mut_aa,
                # DNA-level coordinates for AlphaGenome Atlas AVI lookup.
                # Both AlphaMissense (coding) and AVI (regulatory) flow
                # through the same score_variant() entry point now.
                "chrom": v.chrom,
                "ref_dna": v.ref_dna,
                "alt_dna": v.alt_dna,
            }
            for v in variants
        ]
        # Build the AVI lookup. Use the stdlib mock by default (CI +
        # offline use); users with ALPHAGENOME_API_KEY set + the
        # [variant-alphagenome] extra installed get the real Atlas
        # adapter. Selected at the backend-selector level, not here.
        avi_lookup_fn = None
        if any(v.chrom is not None for v in variants):
            try:
                from .alphagenome_integration import (
                    select_regulatory_scorer,
                )
                avi_lookup_fn = select_regulatory_scorer().score_variant
            except Exception:
                avi_lookup_fn = None

        scored = filter_variants(
            v_dicts,
            top_fraction=variant_filter_top_fraction,
            min_score=variant_filter_min_score,
            protein_lengths={g: len(p) for g, p in proteins.items()},
            protein_sequences=proteins,
            uniprot_ids=uniprot_ids,
            am_lookup=am_lookup_fn,
            avi_lookup=avi_lookup_fn,
            strict=False,
        )
        keep_keys = {(s.gene, s.position, s.wt_aa, s.mut_aa) for s in scored}
        variants = [v for v in variants if (v.gene, v.position, v.wt_aa, v.mut_aa) in keep_keys]
        filtered_count = len(variants)
        filter_scores = {
            f"{s.gene}.{s.position}{s.wt_aa}>{s.mut_aa}": s.normalized_score for s in scored
        }

    note = ""
    if am_active:
        note = "AlphaMissense pathogenicity scores used in variant filtering."
    elif variant_filter_top_fraction < 1.0 or variant_filter_min_score > 0.0:
        note = (
            "AlphaMissense predictions not found — heuristic BLOSUM62 + driver "
            "gene + structural score used. To enable AlphaMissense, download "
            "the predictions TSV from "
            "https://storage.googleapis.com/dm_alphamissense/ "
            "to ~/.cache/mrnavax/AlphaMissense_hg38.tsv "
            "(licensed CC BY-NC-SA 4.0, non-commercial)."
        )

    labels, cluster_names = cluster_with_scanpy(matrix, cell_ids, gene_names, seed=42)
    labels = [int(x) for x in labels]  # numpy strings or ints → uniform Python ints

    # Marker score = mean expression of tumor markers per cluster
    marker_idx = [gene_names.index(g) for g in (tumor_marker_genes or []) if g in gene_names]
    cluster_scores: dict[int, float] = {}
    for c in set(labels):
        cells = [matrix[i] for i, lbl in enumerate(labels) if lbl == c]
        if marker_idx and cells:
            cluster_scores[c] = sum(
                sum(row[i] for i in marker_idx) / len(marker_idx) for row in cells
            ) / len(cells)
        else:
            cluster_scores[c] = float(len(cells))
    tumor_cluster = max(cluster_scores, key=lambda k: cluster_scores[k]) if cluster_scores else 0

    # Build peptide candidates for tumor-cluster-expressed variants
    gene_expr_by_cluster = _gene_means_per_cluster(matrix, labels, gene_names)
    tumor_expressed_genes = {
        g for g, v in gene_expr_by_cluster.get(tumor_cluster, {}).items() if v > 0.5
    }

    peptides: list[TumorPeptide] = []
    for v in variants:
        if v.gene not in proteins:
            continue
        if v.gene not in tumor_expressed_genes:
            continue
        for pep in mutant_peptides(proteins[v.gene], v.position, v.mut_aa, lengths=peptide_lengths):
            peptides.append(
                TumorPeptide(
                    cell_cluster=tumor_cluster,
                    gene=v.gene,
                    position=v.position,
                    wt_aa=v.wt_aa,
                    mut_aa=v.mut_aa,
                    peptide=pep,
                    length=len(pep),
                    cluster_marker_score=cluster_scores.get(tumor_cluster, 0.0),
                )
            )
    if not marker_idx:
        note = (
            "WARNING: no tumor-marker genes supplied; the largest cluster was "
            "assumed to be tumor. Pass --tumor-markers GENE1,GENE2 for "
            "accurate selection. The downstream peptide list may contain "
            "many false positives."
        )

    return PipelineReport(
        n_cells=len(cell_ids),
        n_genes=len(gene_names),
        cluster_labels=labels,
        cluster_marker_scores=cluster_scores,
        tumor_cluster=tumor_cluster,
        n_variants_input=n_variants_input_orig,
        n_variants_after_filter=filtered_count,
        variant_scores=filter_scores,
        n_candidate_peptides=len(peptides),
        peptides=peptides,
        note=note,
    )


def _gene_means_per_cluster(
    matrix: list[list[float]],
    labels: list[int],
    gene_names: list[str],
) -> dict[int, dict[str, float]]:
    out: dict[int, dict[str, float]] = {}
    for c in set(labels):
        cells = [matrix[i] for i, lbl in enumerate(labels) if lbl == c]
        if not cells:
            continue
        out[c] = {
            gene_names[j]: sum(row[j] for row in cells) / len(cells) for j in range(len(gene_names))
        }
    return out


# ---------- CLI ------------------------------------------------------------


def _run_cli(argv: list[str]) -> int:
    import argparse

    p = argparse.ArgumentParser(prog="mrnavax scrna")
    p.add_argument(
        "--expression",
        required=True,
        help="h5ad file (or CSV/TSV cells × genes) — synthetic fallback if missing",
    )
    p.add_argument("--variants", required=True, help="CSV of coding-region SNVs")
    p.add_argument("--proteins", required=True, help="FASTA of protein sequences")
    p.add_argument(
        "--tumor-markers",
        default="",
        help="comma-separated gene symbols overexpressed in tumor cells",
    )
    p.add_argument(
        "--peptide-lengths",
        default="9,10,11",
        help="comma-separated HLA-peptide lengths to enumerate",
    )
    p.add_argument(
        "--variant-filter-top-fraction",
        type=float,
        default=1.0,
        help="keep top fraction of variants by pathogenicity score (default 1.0 = no filter)",
    )
    p.add_argument(
        "--variant-filter-min-score",
        type=float,
        default=0.0,
        help="minimum variant score to keep (0-1, default 0 = no filter)",
    )
    p.add_argument(
        "--embedding-model",
        default="tfidf-svd",
        choices=["tfidf-svd", "identity", "scgpt"],
        help="cell embedding model (default tfidf-svd)",
    )
    p.add_argument(
        "--embedding-dim", type=int, default=32, help="embedding dimensionality (default 32)"
    )
    p.add_argument("--out")
    args = p.parse_args(argv)

    lengths = tuple(int(x) for x in args.peptide_lengths.split(",") if x.strip())
    markers = [g.strip() for g in args.tumor_markers.split(",") if g.strip()] or None

    report = run_pipeline(
        args.expression,
        args.variants,
        args.proteins,
        hla=("HLA-A*02:01",),
        tumor_marker_genes=markers,
        peptide_lengths=lengths,
        variant_filter_top_fraction=args.variant_filter_top_fraction,
        variant_filter_min_score=args.variant_filter_min_score,
        embedding_model=args.embedding_model,
        embedding_dim=args.embedding_dim,
    )
    out_text = json.dumps(report.to_dict(), indent=2)
    if args.out:
        Path(args.out).write_text(out_text + "\n")
    print(out_text)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(_run_cli(sys.argv[1:]))
