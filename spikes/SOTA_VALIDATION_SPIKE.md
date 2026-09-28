# SOTA Validation Spike — Can mrnavax be validated against a published benchmark?

**Date:** 2026-09-28
**Status:** Research only — no production code, no release shipped.
**Goal:** Identify ONE concrete published mRNA-design benchmark that mrnavax can be validated against in this environment.

---

## 1. Brief — which paper's benchmark is most reproducible

Two SOTA mRNA-design papers were read in full:

- **RNop** (Gong et al., *arXiv:2505.23862*, v2 Aug 2026) — "mRNA Design and Optimization with Deep Knowledge-Infused Approach." [1]
- **RiboDecode** (Li et al., *Nat Commun* 16, 9957, Nov 2025) — "Deep generative optimization of mRNA codon sequences for enhanced mRNA translation and therapeutic efficacy." [2]

**Most of the published numbers in both papers are wet-lab.** RNop's headline number is "up to **2.28-fold** higher expression than the HO [Human-Optimized] group at 48 h" for EGFP in HEK293T (Fig. 5, panel e). [1] RiboDecode's headline numbers are "**~4.4-fold** prime / **9.6-fold** boost" neutralizing antibody titers for HA mRNA in BALB/c mice (Fig. 5e) and "**one-fifth** the dose" RGC neuroprotection for NGF mRNA in an optic-nerve-crush mouse model (Fig. 6). [2] Neither is reproducible in silico.

**Both papers also publish *in-silico* numbers that are reproducible in principle** — but with caveats:

- RNop's quantitative claim is on a **540,000-sequence held-out test set** drawn from NCBI (6 M total sequences sampled 2M from Eukaryota/Bacteria/Virus). They report per-species **CAI, tAI, MFE** improvements. [1, §1.1.1 + Fig. 3] But the **test sequences are not released** with the repo (HudenJear/RPLoss ships pretrained weights via Google Drive, not the test split). MFE is computed via **RNAStructure**; mrnavax's `multi-objective` backend reports a **sliding-window base-pair-count proxy** in arbitrary units. The MFE numbers are therefore **not directly comparable**.
- RiboDecode publishes per-gene **in-silico translation-level predictions** (Fig. 3b) and **protein expression fold-changes** vs LinearDesign. Their model runs as a closed `.whl` from Google Drive (Python 3.8.19, CUDA 12.1, ViennaRNA 2.6.4) — not pip-installable as an open model. [3]
- mRNABench (morrislab/mRNABench) is for **embedding benchmarks** (foundation-model property prediction), not codon optimization. It contains a "Translational Efficiency" dataset but it evaluates embeddings, not codon swap algorithms. [4]
- The only other benchmark repo (Vcoderprogreat/codon-optimization-benchmark, 0 stars) compares MyTool vs DNAChisel vs a greedy baseline on hand-crafted `protein_short.faa / protein_medium.faa / protein_long.faa` inputs and reports CAI only. Not SOTA-grade. [5]

**The most reproducible benchmark target** — and the only one where mrnavax can fairly be put side-by-side with a published SOTA number on the *same* input protein — is the **RNop in-silico multi-objective triple (CAI, MFE, MFE-improvement-λ) on the EGFP Human-Optimized baseline** that RNop's Figure 5 visualizes.

That target is still partial, not clean, because the exact HO sequence RNop used is not in the public repo and would need to be reconstructed from a Human-Optimized EGFP source.

---

## 2. Specific benchmark target (recommended)

| Item | Value | Source |
|---|---|---|
| **Paper** | RNop, arXiv:2505.23862 v2, Aug 2026 | [1] |
| **Section** | Fig. 5(a, c, e), "Results of EGFP optimization" | [1] |
| **Reporter protein** | EGFP (enhanced green fluorescent protein) | [1] |
| **Comparison** | RNop-C variant vs **HO (Human-Optimized)** control at 48 h in HEK293T | [1, Fig. 5c/e] |
| **Published number** | **"up to 2.28-fold higher expression than the HO group at 48 h"** (RNop-C variant) | [1, §2.4, RNop-C panel] |
| **Companion in-silico targets (same Fig. 5a caption)** | LD-MFE = −206.1 kcal/mol; HO-MFE = −132.4 kcal/mol; RNop-MFE = −137.4 kcal/mol (EGFP) | [1, Fig. 5a] |
| **Why EGFP** | It is the most-widely-used reporter, so a "Human-Optimized" reference sequence is reconstructable from public catalogues (e.g., Addgene #13031; or the EGFP ORF in pEGFP-N1). mrnavax ships an EGFP amino-acid sequence in `examples/cds_gfp.json` (AEQVI *Aequorea victoria* GFP, 239 aa). | [6, internal] |
| **Cross-check** | RNop explicitly states "any significant improvement is exceptionally difficult to achieve" for EGFP — it is the strictest SOTA bar in the protein-optimization literature. | [1, §2.4] |

**Key caveat — what mrnavax *can* and *cannot* compare:**

| Metric | Reproducible by mrnavax in this env? | Comparable to RNop's published number? |
|---|---|---|
| **CAI** (Sharp & Li 1987 formula) | ✅ Yes — `mrnavax.codon_optimizer.analyze_cds()` implements it with the same Kazusa-style human codon table | ✅ Same formula, same table → directly comparable |
| **GC%** | ✅ Yes | ⚠️ Not what RNop reports as the headline |
| **CpG O/E** | ✅ Yes | ❌ Not reported by RNop |
| **Rare-codon fraction** | ✅ Yes | ❌ Not reported by RNop |
| **GC-window stddev (structure-proxy)** | ✅ Yes | ❌ Not the same as RNop's RNAStructure MFE |
| **MFE (kcal/mol, thermodynamic)** | ❌ No — would need ViennaRNA + `RNAfold`, optional install via `ribodecode-real` backend | ❌ mrnavax reports a base-pair-count *proxy*, not kcal/mol |
| **tAI** | ❌ No — mrnavax has no tRNA-adaptation-index implementation | ❌ RNop reports tAI per Fig. 4 |
| **Wet-lab HEK293T fold-change** | ❌ No lab | ❌ out of scope |

So in practice the cleanest single-number comparison is **CAI**, and a weaker secondary signal is **GC-window stddev** (named "structure proxy" in RNop's RNop pattern). MFE is off the table without an external ViennaRNA install.

---

## 3. What mrnavax feature(s) would this test?

This benchmark hits three mrnavax subsystems directly:

1. **`mrnavax/codon_optimizer.py::optimize_basic`** — the greedy per-codon frequency swap. Confirms the canonical baseline still applies.
2. **`mrnavax/codon_multi_objective.py`** — the `multi-objective` backend. This is the v0.25.0 RNop-pattern knowledge-infused loss (CAI + GC + CpG + rare-codon-run + structure-proxy) — a direct in-silico reproduction of RNop's "knowledge-infusion" idea in pure Python. This is the feature that the RNop paper itself describes and that mrnavax explicitly cites. [1, §1.1.7–1.1.9], internal `codon_multi_objective.py` docstring.
3. **`mrnavax/codon_ribodecode.py::optimize_ribodecode`** — the stdlib-only hill-climb with rare-run / rare-pair / GC-window-stddev penalties. Closely mirrors the operational pattern of the real RiboDecode (Li 2025) without needing their GPU `.whl`.

Out-of-scope for this benchmark: UTR context scoring, construct assembly, variant pathogenicity, AlphaGenome Atlas predict.

---

## 4. Estimated effort to reproduce

**End-to-end in this environment (macOS, no GPU):** **3–5 hours**.

Breakdown:
| Step | Effort | Notes |
|---|---|---|
| Acquire / reconstruct the EGFP Human-Optimized CDS RNop used | 1–2 h | The exact HO sequence RNop used is **not in the public repo**. Need to use a known HO EGFP (Addgene #13031 = "human-codon-optimized" EGFP from the Zhang lab; or the EGFP ORF from pcDNA3.1-EGFP). 30-min fallback: use the EGFP sequence published in Cormack 1996 / Yang 1996 / the EGFP ORF in any major commercial HO construct. |
| Run mrnavax `multi-objective` on the AA sequence | 30 min | Single command, <1 min runtime on 239 aa |
| Compute CAI before/after; compute | 15 min | already in mrnavax |
| Compute RNop-style in-silico comparison table (CAI, GC%, CpG, GC-window stddev) | 1 h | Compose script + verify per-component breakdown matches RNop's pattern |
| Write the result + interpret gap to the wet-lab "2.28-fold" claim | 30 min | **Honest framing**: in-silico numbers cannot reproduce wet-lab fold-changes |

**To make it stronger (out of scope for this spike):**
- Install ViennaRNA (`brew install viennarna` or `pip install ViennaRNA`) → real MFE comparison in ~30 min more.
- Implement tAI in mrnavax (uses GtRNAdb codon-anticodon copy numbers) → ≈ 2 h to add + test, unlocks RNop's tAI column.
- Optionally run the published `ribo-decode` wheel from Google Drive (needs CUDA) → 2–3 h including env setup.

**To make it *wet-lab* (out of scope entirely):** HeLa/HEK293T transfection + fluorescence plate reader + western blot. Weeks of work, not a spike.

---

## 5. Concrete risks

1. **No wet-lab → the "2.28-fold" headline number is unreachable.** RNop's gold metric is in vitro protein expression. mrnavax cannot manufacture mRNA or run HEK293T transfections in this spike.
2. **The exact EGFP-HO CDS RNop benchmarked is not public.** Different HO sequences give different baselines; comparing "EGFP fold change" without controlling for which HO is used is a methodology hole. The published 2.28-fold number is internally consistent only if you trust RNop used the same HO they reference.
3. **MFE is not portable.** mrnavax reports a base-pair-count proxy in arbitrary units; RNop uses RNAStructure kcal/mol. Without an external ViennaRNA install in this spike we cannot compute a real MFE for direct comparison. Even with ViennaRNA installed, RNop's reported −132.4 kcal/mol HO value depends on the exact tool version and folding temperature they used.
4. **mrnavax has no tAI implementation.** RNop's Fig. 3–4 reports tAI as a core metric; reproducing it would require adding ~50 lines of tRNA-copy-number handling (low risk, just out of scope for a spike).
5. **The published RNop pretrained weights live on Google Drive** and are not version-pinned in the GitHub repo. Running RNop yourself requires clicking through Google's "download" UI, downloading ~hundreds of MB, and reproducing with the exact A800 GPU + CUDA 12.2 + PyTorch 2.1 stack they list. Not reproducible in this sandbox.
6. **RiboDecode is shipped as a closed `.whl`** (no source), tied to ViennaRNA 2.6.4 + Python 3.8.19. Reproducing RiboDecode's exact sequences in a side-by-side run is fragile.
7. **Both papers' test sets (the 540k NCBI sequences for RNop, the 60k human isoforms for RiboDecode) are not released as standalone files.** We can reconstruct a *small* in-silico test set of ~5–20 published therapeutic-protein ORFs (EGFP, ECFP, luciferase, EPO, insulin, NGF, Cas9) and use that — but that is a *new* test set, not the published one, so any comparison is illustrative, not a true reproduction.
8. **"in silico" CAI is already near ceiling for any reasonable optimizer** (Fig. 3 of RNop shows CAI saturating at ~0.97 from a baseline of ~0.67). The signal is small, easy to game, and not the same as protein expression. We must not oversell the comparison.

---

## 6. Honest verdict — reproduce or skip?

**REPRODUCE the in-silico part; SKIP the wet-lab part.** Concretely:

- **Do (3–5 h):** Run mrnavax `multi-objective` on the EGFP amino-acid sequence from `examples/cds_gfp.json` (or a published Human-Optimized EGFP CDS), report the CAI before/after, GC%, CpG O/E, rare-codon fraction, GC-window-stddev, and a per-component loss breakdown. Frame it explicitly as: "this is the in-silico knowledge-infused optimization pattern from arXiv:2505.23862 §1.1.7–1.1.9, evaluated against the same metric family (CAI) that RNop reports in Fig. 3–4." Acknowledge up-front that mrnavax's structure-proxy is not a thermodynamic MFE.
- **Skip:** Wet-lab HEK293T transfection. Anyone claiming to "reproduce the 2.28-fold expression gain" of RNop from in silico alone is fabricating — the wet-lab gap is unbridgeable from a laptop.
- **Don't claim SOTA dominance.** mrnavax's `multi-objective` is a small stdlib implementation of the RNop *pattern* (knowledge-infused loss), not the same model. A fair statement is "implements the same biological priors; matches the same metric family in direction; cannot reproduce absolute fold-expression claims without wet-lab validation."
- **Optional second-tier work (separate spike):** Add a tAI module using GtRNAdb human tRNA copy numbers (the only other RNop-comparable metric missing), and add a ViennaRNA-backed real-MFE computation for the `ribodecode-real` backend path. These two changes would let a future spike put mrnavax on the same Fig. 3/4 axes as RNop for several small proteins.

The overall conclusion: **there is a real, narrow, partially-reproducible in-silico benchmark target (RNop EGFP CAI comparison), but no SOTA wet-lab benchmark mrnavax can honestly claim to have reproduced without a lab.** The mrnavax README already cites RNop as the design inspiration for `multi-objective`; this spike would let the documentation say "and validated against RNop's published CAI direction on EGFP" with citations, which is a meaningful but honest claim.

---

## Sources

- [1] Gong Z. et al. *mRNA Design and Optimization with Deep Knowledge-Infused Approach.* arXiv:2505.23862 v2 (Aug 2026). https://arxiv.org/abs/2505.23862 — full PDF read end-to-end.
- [2] Li Y., Wang F., Yang J. et al. *Deep generative optimization of mRNA codon sequences for enhanced mRNA translation and therapeutic efficacy.* Nat Commun 16, 9957 (Nov 2025). https://www.nature.com/articles/s41467-025-64894-x
- [3] Wang F. *RiboDecode* GitHub repository. https://github.com/wangfanfff/RiboDecode — code reads as a closed `.whl` plus env stub; main weights bundled as Google Drive download.
- [4] Shi R. et al. *mRNABench: A curated benchmark for mature mRNA property and function prediction.* bioRxiv 2025.07.05.662870. Code: https://github.com/morrislab/mRNABench — embedding benchmark, not codon-optimization benchmark.
- [5] Vcoderprogreat. *mRNA Codon Optimization Benchmark* GitHub repo (0 stars). https://github.com/Vcoderprogreat/codon-optimization-benchmark — not SOTA-grade; includes only a custom MyTool, DNAChisel, and a greedy baseline.
- [6] Internal — `mrnavax-home/examples/cds_gfp.json` (AEQVI *Aequorea victoria* GFP amino acid sequence, 239 aa). The matching Human-Optimized DNA CDS is not bundled; reconstructable from Addgene #13031 or any standard pcDNA3.1-EGFP construct.

---

*Generated by research subagent. Spike, no production code. No release shipped. No fabrication of wet-lab results.*