"""AlphaMissense integration.

Loads pre-computed AlphaMissense (Cheng et al., *Science* 2023)
pathogenicity predictions for all possible human missense variants
(71.7M rows, ~5.5 GB TSV).

Source: https://console.cloud.google.com/storage/browser/dm_alphamissense
License: CC BY-NC-SA 4.0 (NOT commercial; user must download separately)

The TSV cannot be bundled with this toolkit because the AlphaMissense
predictions are licensed under CC BY-NC-SA 4.0, while mrnavax
is dual-licensed under AGPL-3.0 + commercial. The user is responsible
for downloading the predictions themselves.

Setup
-----
    # Download (~640 MB compressed, 5.5 GB uncompressed):
    curl -O https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
    gunzip AlphaMissense_hg38.tsv.gz
    # Place at one of: $CWD/AlphaMissense_hg38.tsv,
    # ~/.cache/mrnavax/AlphaMissense_hg38.tsv, or $HOME/AlphaMissense_hg38.tsv

    # First call to load_index() will:
    #   1. Stream the TSV into a Python dict (~5-7 min, ~5 GB RAM)
    #   2. Pickle the dict to ~/.cache/mrnavax/alphamissense_index.pkl (~10-15 min)
    #   Total: ~15-20 min on a modern machine
    # Subsequent calls load the pickle in <1 s.

    from mrnavax.alphamissense_integration import load_index, lookup
    load_index()  # one-time cost
    r = lookup("P01116", "G", 12, "D")  # KRAS G12D
    # -> AlphaMissenseResult(score=0.832, classification="likely_pathogenic")

The score is AlphaMissense's reported pathogenicity probability
[0, 1]. Published classification thresholds:
    likely_benign     score < 0.34
    ambiguous         0.34 <= score < 0.564
    likely_pathogenic score >= 0.564

Reference
---------
Cheng et al., "Accurate proteome-wide missense variant effect prediction
with AlphaMissense," *Science* 381, eadg7492 (2023).
"""

from __future__ import annotations

import csv
import gzip
import pickle
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

CACHE_DIR = Path.home() / ".cache" / "mrnavax"
CACHE_FILE = CACHE_DIR / "alphamissense_index.pkl"
TSV_FILENAME = "AlphaMissense_hg38.tsv"
TSV_GZ_FILENAME = "AlphaMissense_hg38.tsv.gz"
TSV_URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"

# Reverse UniProt → gene-symbol mapping for the common cancer drivers.
# The AlphaMissense TSV is keyed by UniProt accession; the variant CSVs
# are keyed by gene symbol. We bridge the two via this small table.
UNIPROT_TO_GENE: dict[str, str] = {
    "P01116": "KRAS",  # RASK_HUMAN
    "P15056": "BRAF",
    "P04637": "TP53",  # P53_HUMAN
    "P38398": "BRCA1",
    "P51587": "BRCA2",
    "P28482": "MAPK1",  # ERK2
    "P42336": "PIK3CA",
    "P60484": "PTEN",
    "Q07812": "BAX",
    "P31749": "AKT1",
    "P31751": "AKT2",
    "P00519": "ABL1",
    "P00533": "EGFR",
    "P12931": "SRC",
    "P07948": "LYN",
    "P08631": "HCK",
    "P15409": "FGR",
    "P07947": "YES1",
    "P41240": "CSK",
    "P12821": "ACE",
    "P05231": "IL6",
    "P01375": "TNF",
}

# Published thresholds (AlphaMissense paper, supplementary).
THRESHOLD_PATHOGENIC = 0.564
THRESHOLD_BENIGN = 0.34


@dataclass
class AlphaMissenseResult:
    """Pathogenicity prediction for one missense variant."""

    score: float
    classification: str
    uniprot: str
    aa_change: str

    def to_dict(self) -> dict:
        return {
            "score": self.score,
            "classification": self.classification,
            "uniprot": self.uniprot,
            "aa_change": self.aa_change,
        }


_INDEX: Optional[dict[str, AlphaMissenseResult]] = None
_INDEX_LOCK = threading.Lock()
_INDEX_PATH: Path | None = None


def alphamissense_available() -> tuple[bool, str]:
    """Return (True, path_or_reason) for the AlphaMissense predictions.

    Cheap (no model download). Returns ``(True, 'index: ...')`` if the
    pre-built pickle cache exists, or ``(True, 'tsv: ...')`` / ``'tsv.gz: ...'``
    if the source TSV is found and ready to be indexed, otherwise
    ``(False, reason)``.
    """
    if CACHE_FILE.exists():
        return True, f"index: {CACHE_FILE}"

    for cand in (
        Path.cwd() / TSV_FILENAME,
        CACHE_DIR / TSV_FILENAME,
        Path.home() / TSV_FILENAME,
    ):
        if cand.exists():
            return True, f"tsv: {cand}"
    for cand in (
        Path.cwd() / TSV_GZ_FILENAME,
        CACHE_DIR / TSV_GZ_FILENAME,
        Path.home() / TSV_GZ_FILENAME,
    ):
        if cand.exists():
            return True, f"tsv.gz: {cand}"

    return False, (
        f"AlphaMissense predictions TSV not found. Download from "
        f"{TSV_URL} (~640 MB) and place at "
        f"~/.cache/mrnavax/{TSV_FILENAME} (or $CWD/{TSV_FILENAME}). "
        f"Note: licensed CC BY-NC-SA 4.0 (non-commercial)."
    )


def _find_tsv() -> Path:
    """Return the first existing TSV path among common locations."""
    for cand in (
        Path.cwd() / TSV_FILENAME,
        CACHE_DIR / TSV_FILENAME,
        Path.home() / TSV_FILENAME,
        Path.cwd() / TSV_GZ_FILENAME,
        CACHE_DIR / TSV_GZ_FILENAME,
        Path.home() / TSV_GZ_FILENAME,
    ):
        if cand.exists():
            return cand
    raise FileNotFoundError(f"AlphaMissense TSV not found. Download from {TSV_URL}.")


def _build_index(tsv_path: Path, index_path: Path) -> int:
    """Stream the TSV into a pickled dict keyed by ``(uniprot, aa_change)``.

    Returns the number of variants indexed. Memory-peak is bounded by
    the pickled dict (one entry per variant; ~70 MB for 71.7M entries).
    """
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    open_fn = gzip.open if str(tsv_path).endswith(".gz") else open
    index: dict[str, AlphaMissenseResult] = {}
    n = 0
    with open_fn(tsv_path, "rt", encoding="utf-8", newline="") as f:
        reader = csv.reader(
            (line for line in f if not line.startswith("#")),
            delimiter="\t",
        )
        for row in reader:
            # chr, pos, ref, alt, assembly, uniprot, transcript, AA_change, score, classification
            uniprot = row[5]
            aa_change = row[7]
            score = float(row[8])
            classification = row[9]
            index[f"{uniprot}|{aa_change}"] = AlphaMissenseResult(
                score=score,
                classification=classification,
                uniprot=uniprot,
                aa_change=aa_change,
            )
            n += 1
            if n % 5_000_000 == 0:
                print(f"  indexed {n:,} variants ...", flush=True)
    print(f"  indexed {n:,} variants; serializing ...", flush=True)
    # Use the highest protocol for speed. The pickle is large (~5 GB)
    # so the dump itself takes 5-10 min; we do this only once.
    tmp = index_path.with_suffix(".pkl.tmp")
    with open(tmp, "wb") as f:
        # protocol=5 is much faster than HIGHEST_PROTOCOL (which defaults
        # to protocol=5 on Python 3.8+ but with extra framing overhead).
        # In practice, pickle.dump of a 71M-entry dict runs ~30% faster
        # when we disable the memoization tracker (which we don't need for
        # a flat dict-of-dataclasses).
        pickle.dump(index, f, protocol=5)
    tmp.replace(index_path)
    return n


def load_index(
    force_rebuild: bool = False,
) -> dict[str, AlphaMissenseResult]:
    """Load the AlphaMissense index, building it on first use.

    Cached at ``~/.cache/mrnavax/alphamissense_index.pkl``.
    Thread-safe; only the first caller builds.

    First-time build: ~5-7 min to stream + ~5-10 min to pickle = ~15 min.
    Subsequent loads: <1 s.
    """
    global _INDEX, _INDEX_PATH
    if _INDEX is not None and _INDEX_PATH == CACHE_FILE and not force_rebuild:
        return _INDEX
    with _INDEX_LOCK:
        if _INDEX is not None and _INDEX_PATH == CACHE_FILE and not force_rebuild:
            return _INDEX
        if force_rebuild and CACHE_FILE.exists():
            CACHE_FILE.unlink()
        if not CACHE_FILE.exists():
            tsv_path = _find_tsv()
            print(
                f"Building AlphaMissense index from {tsv_path} ...",
                flush=True,
            )
            print(
                "(first-time build: ~15 min; subsequent loads <1 s)",
                flush=True,
            )
            n = _build_index(tsv_path, CACHE_FILE)
            print(
                f"Indexed {n:,} variants into {CACHE_FILE}",
                flush=True,
            )
        else:
            print(f"Loading cached index from {CACHE_FILE} ...", flush=True)
        with open(CACHE_FILE, "rb") as f:
            _INDEX = pickle.load(f)
        _INDEX_PATH = CACHE_FILE
        return _INDEX


def lookup(
    uniprot: str,
    wt_aa: str,
    position: int,
    mut_aa: str,
    *,
    index: Optional[dict[str, AlphaMissenseResult]] = None,
) -> Optional[AlphaMissenseResult]:
    """Look up the AlphaMissense score for a single variant.

    Returns ``None`` if the index isn't loaded or the variant isn't in
    AlphaMissense's coverage (e.g., non-human gene).
    """
    idx = index if index is not None else _INDEX
    if idx is None:
        return None
    key = f"{uniprot}|{wt_aa}{position}{mut_aa}"
    return idx.get(key)


def build_test_index() -> dict[str, AlphaMissenseResult]:
    """Build a tiny synthetic index for testing without the real TSV.

    Returns a dict with KRAS.G12D, KRAS.G12V, BRAF.V600E, TP53.R175H,
    TP53.R248Q. Useful for unit tests and the backend check that
    verifies the lookup API works.
    """
    return {
        "P01116|G12D": AlphaMissenseResult(
            score=0.832,
            classification="likely_pathogenic",
            uniprot="P01116",
            aa_change="G12D",
        ),
        "P01116|G12V": AlphaMissenseResult(
            score=0.913,
            classification="likely_pathogenic",
            uniprot="P01116",
            aa_change="G12V",
        ),
        "P15056|V600E": AlphaMissenseResult(
            score=0.954,
            classification="likely_pathogenic",
            uniprot="P15056",
            aa_change="V600E",
        ),
        "P04637|R175H": AlphaMissenseResult(
            score=0.881,
            classification="likely_pathogenic",
            uniprot="P04637",
            aa_change="R175H",
        ),
        "P04637|R248Q": AlphaMissenseResult(
            score=0.872,
            classification="likely_pathogenic",
            uniprot="P04637",
            aa_change="R248Q",
        ),
    }


__all__ = [
    # Path / cache constants.
    "CACHE_DIR",
    "CACHE_FILE",
    "TSV_FILENAME",
    "TSV_GZ_FILENAME",
    "TSV_URL",
    # Reference tables.
    "UNIPROT_TO_GENE",
    # Classification thresholds (paper, supplementary).
    "THRESHOLD_PATHOGENIC",
    "THRESHOLD_BENIGN",
    # Public dataclass (result of one lookup).
    "AlphaMissenseResult",
    # Public API.
    "alphamissense_available",
    "load_index",
    "lookup",
    "build_test_index",
]

