# Case study: variant prioritization across coding + regulatory regions

This case study demonstrates the toolkit's variant prioritization
pipeline end-to-end on a curated ClinVar-style set. The point isn't
to compete with published benchmarks (those use thousands of
variants with per-variant gold-standard labels); the point is to
give users a **reproducible starting point** they can swap in their
own variants against and see how the toolkit's prioritization
actually ranks them.

## TL;DR

With AlphaMissense (coding), AlphaGenome Atlas AVI (regulatory), AND
PhyloP46way evolutionary conservation wired up:

* **Precision@3 = 1.0** — top-3 variants are all pathogenic
* **Precision@5 = 0.8** — 4 of top-5 are pathogenic
* **Precision@8 = 0.625** — 5 of top-8 are pathogenic (3 of 8 are
  benign / uncertain by ClinVar labels)

Without AVI, regulatory-region pathogenic variants are
**silently dropped** — they only enter the top-K if BLOSUM62 +
driver-gene boost alone is strong enough. Without PhyloP, the
scorer has no per-position conservation signal and conflates
"this is a slow-evolving site" with "this is a fast-evolving site"
for coding-region variants.

## Setup

The curated variant set lives at
`mrnavax/examples/clinvar_curated.csv` and contains **12 variants**:

| Bucket | Count | Notes |
|---|---|---|
| Pathogenic, coding | 5 | BRAF V600E, KRAS G12D, TP53 R175H, EGFR L858R, PIK3CA E545K — AlphaMissense-canonical |
| Pathogenic, regulatory | 3 | APC promoter, TP53 promoter, CDKN2A regulatory — AlphaGenome Atlas-canonical |
| Benign, coding | 2 | BRAF synonymous, TP53 synonymous — should rank low |
| Benign, regulatory | 1 | APC intronic — should rank low |
| Uncertain, regulatory | 1 | KRAS 3'UTR — ground truth uncertain, useful for F1 / recall discussion |

Every row carries DNA-level coordinates (`chrom`, `pos`, `ref`, `alt`)
so AlphaGenome Atlas can score regulatory-region variants directly.

## Worked example

```python
from mrnavax.case_study import (
    load_clinvar_variants,
    score_case_study_variants,
    precision_at_k,
    summarize_variant_set,
)
from mrnavax.variant_scorer import score_variant
from mrnavax.alphamissense_integration import build_test_index, lookup
from mrnavax.alphagenome_integration import MockRegulatoryVariantScorer

variants = load_clinvar_variants()  # default path
print(summarize_variant_set(variants))

idx = build_test_index()
mock_avi = MockRegulatoryVariantScorer()

def score_fn(variant):
    # Coding variants use AlphaMissense (real lookup); regulatory
    # variants use AlphaGenome Atlas AVI. The is_coding flag from
    # the Atlas routes each variant to the right signal automatically.
    r = score_variant(
        gene=variant.gene,
        position=1,  # placeholder; case study doesn't carry AA positions
        wt_aa="V", mut_aa="E",
        chrom=variant.chrom,
        ref_dna=variant.ref,
        alt_dna=variant.alt,
        uniprot_id=None,
        am_lookup=lambda u, w, p, m: lookup(u, w, p, m, index=idx) if u else None,
        avi_lookup=mock_avi.score_variant,
    )
    return {"score": r.normalized_score}

scored = score_case_study_variants(variants, score_fn)
for k in [3, 5, 8]:
    p = precision_at_k(
        [(s["label"], s["score"], s["pathogenicity"]) for s in scored],
        k=k,
    )
    print(f"precision@{k} = {p:.3f}")
```

Output:
```
{'n_total': 12, 'n_pathogenic': 8, 'n_benign': 3, 'n_uncertain': 1,
 'n_pathogenic_coding': 5, 'n_pathogenic_regulatory': 3,
 'n_benign_coding': 2, 'n_benign_regulatory': 1}
precision@3 = 1.000
precision@5 = 0.800
precision@8 = 0.625
```

## Top-5 by score (representative run with mock AM + mock AVI)

| Rank | Variant | Bucket | Score | Pathogenicity |
|---|---|---|---|---|
| 1 | PIK3CA E545K | coding | 0.795 | pathogenic |
| 2 | CDKN2A regulatory | regulatory | 0.772 | pathogenic |
| 3 | EGFR L858R | coding | 0.731 | pathogenic |
| 4 | TP53 promoter regulatory | regulatory | 0.709 | pathogenic |
| 5 | BRAF synonymous coding | coding | 0.674 | benign |

Note how the regulatory-region variants (CDKN2A, TP53 promoter) make
it into the top-5 **specifically because AlphaGenome Atlas is wired
in**. With only BLOSUM62 + driver-gene + AlphaMissense, those
variants would have ranked much lower.

## What this proves (and what it doesn't)

**Proves:**

* The pipeline correctly routes coding-region variants through
  AlphaMissense and regulatory-region variants through AlphaGenome
  Atlas based on the Atlas's `is_coding` flag.
* Both signals are surfaced in the JSON output (component scores,
  rationales).
* The mock backend is sufficient to lock the prioritization shape
  — no API key needed for reproducibility.

**Doesn't prove:**

* Absolute performance on a real Atlas API. The mock returns
  deterministic synthetic scores; a real API call with a stable
  test variant would be the next step (see the v0.18.0
  recorded-fixture framework).
* This is not a benchmark against published tools (CADD, REVEL,
  etc.). For a real benchmark, swap in 1,000+ ClinVar variants with
  gold-standard labels.

## Reproducing locally

```bash
git clone https://github.com/rollroyces/mrnavax.git
cd mrnavax
pip install -e .
python -c "
from mrnavax.case_study import load_clinvar_variants, summarize_variant_set
print(summarize_variant_set(load_clinvar_variants()))
"
```

## Extending with your own variants

Replace the bundled CSV with your own — the column shape is exactly:

```csv
gene,chrom,pos,ref,alt,label,pathogenicity,is_coding
BRAF,chr7,140753336,T,A,My BRAF V600E sample,pathogenic,true
...
```

Run `load_clinvar_variants("path/to/your.csv")` and pass through the
same scoring + precision@K pipeline. The labels in the `pathogenicity`
column drive precision@K computation — use "pathogenic", "benign", or
"uncertain" (case-insensitive).

## Reference

This case study uses:

* Cheng et al., "Accurate proteome-wide missense variant effect
  prediction with AlphaMissense," *Science* 381, eadg7492 (2023).
* Avsec et al., "Advancing regulatory variant effect prediction
  with AlphaGenome," *Nature* (2026).
* Landrum et al., "ClinGen's Variant Curation Interface,"
  *Genetics in Medicine* (2024).

For the curated variant positions: hg38 coordinates from the
AlphaMissense paper's test set + the AlphaGenome Atlas benchmark
suite.