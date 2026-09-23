# `construct` — full mRNA construct assembly

Compose a complete therapeutic mRNA construct (5'UTR + CDS + 3'UTR +
polyadenylation signal + poly-A tail) from a protein amino-acid
sequence. Stdlib-only — no heavy deps.

The construct designer is the natural completion of the codon
optimizer: instead of just optimizing a CDS in isolation, the
designer produces the full construct that an mRNA therapy would
actually deliver to a patient. It composes existing primitives:

* **5'UTR** — Human consensus with strong Kozak context
  (`GGGCGACGCGGTGGCGGCCGCTCATGG`).
* **CDS** — Reverse-translated from the protein AA sequence, then
  optimized via one of the [codon](../tools/codon.md) backends
  (`basic` / `ribodecode` / `lineardesign` / `ribodecode-real` /
  **`multi-objective` ← default since v0.25.0**).
* **3'UTR** — Human consensus with two ARE-stabilizing elements.
* **Polyadenylation signal** — `AAUAAA` (DNA: `AATAAA`).
* **Poly-A tail** — Default 120 nt (configurable via `--poly-a-length`).

## Usage

```bash
# Compose a full mRNA construct from a protein AA sequence
mrnavax construct --sequence "MVSKGEELFTGV"

# Use a specific CDS-optimization backend
mrnavax construct --sequence my_protein.fasta \
    --backend lineardesign

# Adjust poly-A tail length (typical: 100-150 nt for therapeutic mRNA)
mrnavax construct --sequence my_protein.fasta \
    --poly-a-length 100

# Save the full report to a JSON file
mrnavax construct --sequence my_protein.fasta \
    --out my_construct.json
```

## Output schema

```json
{
  "input_protein": "MVSKGEELFTGV",
  "construct_dna": "GGGCGACGCGGTGGCGGCCGCTCATGG...AATAAA...AAAA...",
  "utr5": "GGGCGACGCGGTGGCGGCCGCTCATGG",
  "cds": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGC...",
  "utr3": "GCCCCTGGGTTCAAGGGGCTCGAGTGAGCTGCATCTTTCTTTT...",
  "poly_a_tail": "AAAA...AAAA",
  "poly_a_signal": "AATAAA",
  "cds_length": 36,
  "construct_length": 333,
  "gc_percent": 34.23,
  "n_codons": 12,
  "cds_analysis": {
    "cai": 0.8876,
    "gc_percent": 50.0,
    "rare_codon_fraction": 0.0,
    "cpg_obs_exp": 0.92,
    ...
  },
  "optimization_notes": [
    "multi-objective score: 0.702 -> 0.805 (+0.103)"
  ]
}
```

## Why this matters

The published therapeutic-mRNA designs (Moderna mRNA-1273, BioNTech
BNT162b2) all share a common architecture: 5'UTR (Kozak context) +
CDS + 3'UTR (ARE elements) + polyadenylation signal + poly-A tail.
Before v0.26.0, mrnavax exposed each component in isolation
(`codon` for CDS, `manufacture` for QC checks). With v0.26.0, the
`construct` tool assembles the full therapeutic-mRNA architecture in
one CLI call, with the SOTA multi-objective optimizer applied to
the CDS by default.

## Reference

The 5'/3' UTR templates follow the human consensus from the published
therapeutic-mRNA designs (Moderna mRNA-1273, BioNTech BNT162b2).
The CDS optimization backends are the same five from
`mrnavax codon` (basic / ribodecode / lineardesign / ribodecode-real /
multi-objective).