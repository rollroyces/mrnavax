# `examples/` — end-to-end runnable recipes

Every example here is committed and **has been run end-to-end** against
the current release. Use these as copy-paste templates for your own
work.

## One-shot: verify your install

```bash
python examples/run_all.py
```

This single script exercises all 11 tools against the example files in
this directory and prints a one-line summary per tool. Mock mode (no
model downloads) takes ~30 seconds. With a live Atlas API key, it
also calls `predict` on BRAF V600E for a real-world smoke test.

Sample output:

```
mrnavax end-to-end demo (mock mode)

  codon              CAI=0.720, n_codons=169
  neoantigen         8 candidates, 0 binders
  trial              top-1 = NCT00000003
  lnp                target=lung, cargo=sarna
  scrna              10 cells × 9 genes, 4 clusters (tumor=3)
  manufacture        score=0.887, 7 pass / 8 total
  spatial            10 tissue modules
  construct          384 nt, CAI=0.887
  utr-design         combined=0.807
  variant-regulatory 6 variants, 6 high-impact
  predict            score=1.000, classification=high
```

## Per-tool recipes

Each tool has a copy-pasteable command. Sample outputs are committed
under `sample_outputs/`.

### 1. `codon` — Codon analysis + optimization

```bash
python -m mrnavax.cli codon --sequence examples/cas9.fasta
python -m mrnavax.cli codon --sequence examples/cas9.fasta \
    --optimize --backend lineardesign
```

Input: `cas9.fasta` (Cas9 coding sequence, 1368 nt).

### 2. `neoantigen` — Peptide × HLA binding + LM immunogenicity

```bash
python -m mrnavax.cli neoantigen \
    --variants examples/tp53_variants.csv \
    --hla HLA-A*02:01
```

Input: `tp53_variants.csv` — one **peptide** per row (the CLI expects
pre-computed MHC-binding peptides, not raw gene/variant metadata).
8 candidate peptides with their originating variant.

### 3. `trial` — Patient-to-trial matching

```bash
python -m mrnavax.cli trial \
    --patient examples/patient_summary.txt \
    --trials examples/trials.jsonl --top-k 5 \
    --matcher trialgpt-simicl
```

Inputs: free-text patient summary + JSONL of trial eligibility.

### 4. `lnp` — LNP composition advice

```bash
python -m mrnavax.cli lnp --target lung --cargo saRNA --intent "cancer vaccine"
```

No input file needed — fully synthetic.

### 5. `scrna` — scRNA-seq → tumor cluster → mutant peptides

```bash
python -m mrnavax.cli scrna \
    --expression examples/cells.csv \
    --variants examples/variants_coding.csv \
    --proteins examples/proteins.fasta \
    --tumor-markers TP53,KRAS,BRAF
```

Inputs:
- `cells.csv` — cell × gene expression (UMI)
- `variants_coding.csv` — gene/position/wt_aa/mut_aa
- `proteins.fasta` — protein reference sequences

### 6. `manufacture` — mRNA manufacturability score

```bash
python -m mrnavax.cli manufacture --cds examples/cds_gfp.json
```

Input: JSON with `amino_acid` or `cds` field.

### 7. `spatial` — Spatial transcriptomics tissue modules

```bash
python -m mrnavax.cli spatial \
    --count-file examples/st_bc2_count_matrix.tsv \
    --locations-file examples/st_bc2_locations.tsv \
    --platform ST --num-modules 10
```

Inputs: count matrix TSV + x/y locations TSV.

### 8. `construct` — Full mRNA construct assembly (5'UTR + CDS + 3'UTR + poly-A)

```bash
# Note: --sequence expects DNA, not AA. For amino-acid input use --protein
# with mrnavax.construct_designer (see docs/tools/construct.md).
python -m mrnavax.cli construct --sequence "ATGGATAAGCTACCGTATCGTATGAATAA"
```

### 9. `utr-design` — Coupled 5'UTR + CDS + 3'UTR design

```bash
python -m mrnavax.cli utr-design --protein "MVSKGEELFTGVVPILVELDGDVNGHKFS"
```

Input: amino-acid sequence (the tool reverse-translates internally).

### 10. `variant-regulatory` — AlphaGenome Atlas AVI scoring

```bash
python -m mrnavax.cli variant-regulatory --csv examples/regulatory_variants.csv
```

Input: `chrom,pos,ref,alt` CSV — `chrom` must include the `chr` prefix.

### 11. `predict` — Live AlphaGenome Atlas regulatory-variant scoring

Requires API key configured via env var, helper file, or `~/.alphagenome_key`.
See `docs/operations/atlas-activation.md`.

```bash
# Single variant
python -m mrnavax.cli predict --chrom chr7 --pos 140753336 --ref T --alt A

# Batch via CSV
python -m mrnavax.cli predict --input examples/regulatory_variants.csv
```

## Sample outputs

JSON outputs for every tool are committed under `sample_outputs/`.
Use these to compare against your own runs without burning API quota.

```
sample_outputs/
├── codon.json              (Cas9 — CAI, GC%, codon usage table)
├── neoantigen.json         (8 peptides × HLA-A*02:01)
├── trial.json              (5 ranked BRAF-V25 trial candidates)
├── lnp.json                (synthetic — no inputs)
├── scrna.json              (10 cells × 9 genes)
├── manufacture.json        (GFP CDS — 7/8 checks pass)
├── spatial.json            (10 tissue modules)
├── construct.json          (mock short CDS — 384 nt)
├── utr-design.json         (12 UTR combinations ranked)
├── variant-regulatory.json (6 cancer variants — 6/6 high)
└── (predict.json is not committed; see live Atlas)
```

## Regenerating sample outputs

```bash
for cmd in \
    "codon --sequence examples/cas9.fasta" \
    "neoantigen --variants examples/tp53_variants.csv --hla HLA-A*02:01" \
    "trial --patient examples/patient_summary.txt --trials examples/trials.jsonl --top-k 5" \
    "lnp --target lung --cargo saRNA --intent 'cancer vaccine'" \
    "scrna --expression examples/cells.csv --variants examples/variants_coding.csv --proteins examples/proteins.fasta --tumor-markers TP53,KRAS,BRAF" \
    "manufacture --cds examples/cds_gfp.json" \
    "spatial --count-file examples/st_bc2_count_matrix.tsv --locations-file examples/st_bc2_locations.tsv --platform ST --num-modules 10" \
    "construct --sequence ATGGATAAGCTACCGTATCGTATGAATAA" \
    "utr-design --protein MVSKGEELFTGVVPILVELDGDVNGHKFS" \
    "variant-regulatory --csv examples/regulatory_variants.csv"
do
    tool=$(echo "$cmd" | cut -d' ' -f1)
    MRNA_AI_FORCE_MOCK=1 python -m mrnavax.cli $cmd \
        > examples/sample_outputs/${tool}.json
done
```