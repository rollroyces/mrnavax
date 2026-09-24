"""mRNA manufacturability checks.

Wet-lab bridge between computational sequence design and what's
actually synthesizable at production scale. Each check returns a
:class:`CheckResult` with a ``pass_`` flag, score (0-1, higher = more
manufacturable), severity (info/warn/error), and a human-readable
explanation.

Checks implemented
------------------
1. ``poly_a_runs`` — runs of ≥5 consecutive As in the DNA template
   destabilize the plasmid during IVT. Counts occurrences and the
   longest run. Severity scales with run length.
2. ``gc_5prime_hairpin`` — GC-rich stems (≥70% GC over 30+ nt) at the
   5' UTR or CDS start block ribosome scanning. Uses a sliding
   window over the first 60 nt.
3. ``are_motif`` — AU-rich elements (ARE) in the 3' UTR trigger mRNA
   decay via TTP / HuR. Canonical nonamer ``UUAUUUAUU`` and pentamer
   ``AUUUA`` repeats. The CDS itself may contain pentamers harmlessly
   (introns mostly aren't in IVT-mRNA); we flag the nonamer explicitly
   in the 3' UTR.
4. ``kozak_strength`` — measures match to the mammalian Kozak consensus
   ``(GCC)GCC(A/G)CCATGG``. Returns a 0-1 score.
5. ``stop_context`` — termination efficiency depends on the stop codon
   identity and the +4 base. Per the literature, readthrough order is
   ``TAA < TAG < TGA``; the most efficient terminator is ``TGA-T`` or
   ``TAA-T``.
6. ``hidden_stops`` — internal in-frame stop codons (should never occur
   in a designed CDS but can sneak in via codon-table drift).
7. ``gc_window_uniformity`` — extreme local GC variation (>25% stddev
   across 30-nt windows) signals poor IVT yield and ribosomal stalling.
8. ``cpg_suppression`` — long CpG-free stretches trigger silencing in
   some contexts; extremely CpG-rich regions (>15% CpG in a window)
   cause immune activation. Returns a balanced score.

Reference
---------
- Holtkamp et al. (2006) "Modification of antigen-encoding RNA increases
  stability, translational efficacy, and T-cell stimulatory capacity of
  dendritic cells." *Blood* 108.
- Kozak (1986) "Point mutations define a sequence flanking the AUG
  initiator codon that modulates translation by eukaryotic ribosomes."
  *Cell* 44.
- Chen & Shyu (1995) "AU-rich elements: characterization and importance
  in mRNA degradation." *Trends Biochem Sci* 20.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field

# ---------------------------------------------------------------------------
# Per-check configuration
# ---------------------------------------------------------------------------

POLY_A_MIN_LEN = 5
GC_HAIRPIN_WINDOW = 30
GC_HAIRPIN_THRESHOLD = 0.70
GC_WINDOW = 30
GC_WINDOW_STDDEV_MAX = 0.15  # >15% stddev flags as warn
CPG_WINDOW = 50
CPG_LOW = 0.005  # <0.5% CpG in a window = suppressed
CPG_HIGH = 0.15  # >15% CpG in a window = too immunogenic

# Kozak strong mammalian consensus: GCCRCCATGG (R = purine)
# We score over the 9 nt preceding ATG.
KOZAK_CONSENSUS = "GCCRCCATGG"

# Stop codon readthrough efficiency (lower = stronger stop).
# Per the literature, readthrough order is TAA < TAG < TGA.
# We assign 0 (strongest) to TAA, 0.5 to TAG, 1.0 (weakest) to TGA,
# then modulate by the +4 base.
STOP_CODON_BASE_SCORE = {"TAA": 0.0, "TAG": 0.3, "TGA": 0.7}
STOP_FOLLOWER_BONUS = {"T": -0.3, "A": 0.1, "C": 0.2, "G": 0.2}

# ARE motifs (RNA alphabet). Single pentamer AUUUA is weak;
# the nonamer UUAUUUAUU is the canonical destabilizing element.
ARE_PENTAMER = "AUUUA"
ARE_NONAMER = "UUAUUUAUU"


# ---------------------------------------------------------------------------
# Check result types
# ---------------------------------------------------------------------------


@dataclass
class CheckResult:
    """Outcome of one manufacturability check."""

    name: str
    pass_: bool
    score: float  # 0..1, higher = better
    severity: str  # "info" | "warn" | "error"
    summary: str
    details: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ManufacturabilityReport:
    """Aggregated report from all manufacturability checks."""

    overall_score: float  # 0..1
    n_pass: int
    n_warn: int
    n_error: int
    checks: list[CheckResult]

    def to_dict(self) -> dict:
        return {
            "overall_score": self.overall_score,
            "n_pass": self.n_pass,
            "n_warn": self.n_warn,
            "n_error": self.n_error,
            "checks": [c.to_dict() for c in self.checks],
        }


# ---------------------------------------------------------------------------
# Individual checks
# ---------------------------------------------------------------------------


def _longest_run(s: str, char: str) -> int:
    """Length of the longest run of ``char`` in ``s``."""
    best = run = 0
    for c in s:
        if c == char:
            run += 1
            if run > best:
                best = run
        else:
            run = 0
    return best


def check_poly_a_runs(dna: str) -> CheckResult:
    """Poly-A runs ≥5 nt destabilize the DNA template.

    Returns a CheckResult that fails if any run ≥7 occurs (severe) or
    if multiple runs ≥5 occur (moderate).
    """
    dna = dna.upper().replace("U", "T")
    longest = _longest_run(dna, "A")
    # Count all runs of length >= POLY_A_MIN_LEN
    pattern = re.compile(f"A{{{POLY_A_MIN_LEN},}}")
    runs = [m.group() for m in pattern.finditer(dna)]
    n_runs = len(runs)
    if longest >= 9:
        sev = "error"
        ok = False
        score = 0.0
    elif longest >= 7:
        sev = "error"
        ok = False
        score = 0.1
    elif n_runs >= 3:
        sev = "warn"
        ok = False
        score = 0.5
    elif longest >= POLY_A_MIN_LEN:
        sev = "warn"
        ok = True  # mild
        score = 0.7
    else:
        sev = "info"
        ok = True
        score = 1.0
    return CheckResult(
        name="poly_a_runs",
        pass_=ok,
        score=score,
        severity=sev,
        summary=(f"longest poly-A run = {longest} nt; {n_runs} run(s) of length ≥{POLY_A_MIN_LEN}"),
        details={"longest_run": longest, "n_runs": n_runs, "runs": runs},
    )


def check_gc_5prime_hairpin(dna: str) -> CheckResult:
    """Sliding GC% over the first 60 nt; flag windows ≥70% GC.

    The 5' UTR's first 30-40 nt set the tone for ribosome scanning.
    Local GC-rich stems cause the 40S subunit to stall.
    """
    dna = dna.upper().replace("U", "T")
    head = dna[:60]
    if len(head) < GC_HAIRPIN_WINDOW:
        return CheckResult(
            name="gc_5prime_hairpin",
            pass_=True,
            score=1.0,
            severity="info",
            summary="sequence <60 nt; skipping 5' hairpin check",
            details={},
        )
    worst_gc = 0.0
    worst_pos = 0
    for i in range(len(head) - GC_HAIRPIN_WINDOW + 1):
        win = head[i : i + GC_HAIRPIN_WINDOW]
        gc = (win.count("G") + win.count("C")) / len(win)
        if gc > worst_gc:
            worst_gc = gc
            worst_pos = i
    if worst_gc >= 0.80:
        sev, ok, score = "error", False, 0.1
    elif worst_gc >= GC_HAIRPIN_THRESHOLD:
        sev, ok, score = "warn", False, 0.5
    else:
        sev, ok, score = "info", True, 1.0
    return CheckResult(
        name="gc_5prime_hairpin",
        pass_=ok,
        score=score,
        severity=sev,
        summary=(f"max GC in first 60 nt = {worst_gc:.0%} at position {worst_pos + 1}"),
        details={"max_gc": worst_gc, "max_position": worst_pos + 1},
    )


def check_are_motif(utr3: str) -> CheckResult:
    """AU-rich elements in the 3' UTR trigger mRNA decay.

    Counts pentamers (AUUUA) and nonamers (UUAUUUAUU) in the 3' UTR.
    Pentamers alone are weak signals; ≥2 nonamers is severe.
    """
    rna = utr3.upper().replace("T", "U")
    n_pent = len(re.findall(ARE_PENTAMER, rna))
    n_non = len(re.findall(ARE_NONAMER, rna))
    if n_non >= 2:
        sev, ok, score = "error", False, 0.1
    elif n_non == 1:
        sev, ok, score = "warn", False, 0.5
    elif n_pent >= 5:
        sev, ok, score = "warn", False, 0.6
    else:
        sev, ok, score = "info", True, 1.0
    return CheckResult(
        name="are_motif",
        pass_=ok,
        score=score,
        severity=sev,
        summary=(
            f"3' UTR: {n_pent} pentamer(s), {n_non} nonamer(s) "
            f"(pentamer = {ARE_PENTAMER}; nonamer = {ARE_NONAMER})"
        ),
        details={"n_pentamer": n_pent, "n_nonamer": n_non},
    )


def check_kozak_strength(utr5: str) -> CheckResult:
    """Score the 9 nt immediately upstream of ATG against Kozak consensus.

    Strong: ``GCCRCCATGG`` (R = purine). We score position-by-position
    with weighted matches: positions -3 (R) and +4 (G) carry more weight
    in the original Kozak paper.
    """
    rna = utr5.upper().replace("T", "U")
    # We need 9 nt upstream of ATG. Convention: ATG is in the CDS;
    # caller passes the trailing 9 nt of the 5' UTR.
    if len(rna) < 9:
        return CheckResult(
            name="kozak_strength",
            pass_=True,
            score=0.5,
            severity="info",
            summary=f"only {len(rna)} nt of 5' UTR; cannot score Kozak",
            details={},
        )
    pre = rna[-9:]
    # Position-by-position match (R = A or G)
    weights = [1, 1, 1, 1.5, 1, 1, 1.5, 1, 1]  # -3 R and +4 G heavier
    expected: list[tuple[str, ...]] = [
        ("G",),
        ("C",),
        ("C",),
        ("A", "G"),  # R
        ("C",),
        ("C",),
        ("A", "G"),  # R
        ("A",),  # A in ATG
        ("T",),
    ]
    score = 0.0
    total_w = 0.0
    for i, exp in enumerate(expected):
        if i >= len(pre):
            break
        total_w += weights[i]
        if pre[i] in exp:
            score += weights[i]
    norm = score / total_w if total_w else 0.0
    if norm >= 0.85:
        sev, ok = "info", True
    elif norm >= 0.65:
        sev, ok = "info", True
    elif norm >= 0.45:
        sev, ok = "warn", False
    else:
        sev, ok = "warn", False
    return CheckResult(
        name="kozak_strength",
        pass_=ok,
        score=norm,
        severity=sev,
        summary=(
            f"5' UTR upstream of ATG ({pre[-9:]}) matches partial Kozak consensus at {norm:.0%}"
        ),
        details={"pre_atg": pre[-9:], "score_partial": norm},
    )


def check_stop_context(dna: str) -> CheckResult:
    """Stop codon identity + +4 base. Strong = TAA or TGA, weak = TGA.

    Reads as RNA: the stop codon is the last three nt of the CDS, and
    ``+4`` is the first base of the 3' UTR.
    """
    dna = dna.upper().replace("U", "T")
    # Take the last 3 nt (CDS stop) and look at +4 if available
    if len(dna) < 3:
        return CheckResult(
            name="stop_context",
            pass_=True,
            score=1.0,
            severity="info",
            summary="sequence <3 nt; skipping stop check",
            details={},
        )
    stop = dna[-3:]
    follower = dna[3] if len(dna) > 3 else "T"  # default to T if missing
    base_score = STOP_CODON_BASE_SCORE.get(stop, 0.5)
    follower_mod = STOP_FOLLOWER_BONUS.get(follower, 0.0)
    # Combined: 0 = strong stop, 1 = leaky
    leakiness = max(0.0, min(1.0, base_score + follower_mod))
    score = 1.0 - leakiness
    if leakiness >= 0.7:
        sev, ok = "warn", False
    elif leakiness >= 0.5:
        sev, ok = "info", True
    else:
        sev, ok = "info", True
    return CheckResult(
        name="stop_context",
        pass_=ok,
        score=score,
        severity=sev,
        summary=(
            f"stop codon {stop} + follower {follower}; "
            f"leakiness={leakiness:.2f} (lower = stronger termination)"
        ),
        details={"stop": stop, "follower": follower, "leakiness": leakiness},
    )


def check_hidden_stops(cds: str) -> CheckResult:
    """Internal in-frame stop codons (should be zero)."""
    cds = cds.upper().replace("U", "T")
    codons = [cds[i : i + 3] for i in range(0, len(cds) - 2, 3)]
    stops = {"TAA", "TAG", "TGA"}
    n_internal = sum(1 for c in codons[:-1] if c in stops)
    if n_internal == 0:
        return CheckResult(
            name="hidden_stops",
            pass_=True,
            score=1.0,
            severity="info",
            summary="no internal in-frame stops",
            details={"n_internal_stops": 0},
        )
    return CheckResult(
        name="hidden_stops",
        pass_=False,
        score=0.0,
        severity="error",
        summary=f"{n_internal} internal in-frame stop codon(s) detected",
        details={"n_internal_stops": n_internal},
    )


def check_gc_window_uniformity(cds: str) -> CheckResult:
    """Local GC% variation across the CDS.

    High stddev = uneven IVT yield, more secondary structure, more
    ribosomal stalling. Threshold: stddev > 15% flags as warn.
    """
    cds = cds.upper().replace("U", "T")
    if len(cds) < GC_WINDOW:
        return CheckResult(
            name="gc_window_uniformity",
            pass_=True,
            score=1.0,
            severity="info",
            summary="sequence <30 nt; skipping GC uniformity check",
            details={},
        )
    gcs: list[float] = []
    for i in range(0, len(cds) - GC_WINDOW + 1, GC_WINDOW // 2):
        win = cds[i : i + GC_WINDOW]
        gcs.append((win.count("G") + win.count("C")) / len(win))
    n = len(gcs)
    if n < 2:
        return CheckResult(
            name="gc_window_uniformity",
            pass_=True,
            score=1.0,
            severity="info",
            summary="too few windows for stddev calculation",
            details={},
        )
    mean = sum(gcs) / n
    var = sum((x - mean) ** 2 for x in gcs) / n
    stddev = var**0.5
    if stddev > 0.25:
        sev, ok, score = "error", False, 0.2
    elif stddev > GC_WINDOW_STDDEV_MAX:
        sev, ok, score = "warn", False, 0.6
    else:
        sev, ok, score = "info", True, 1.0
    return CheckResult(
        name="gc_window_uniformity",
        pass_=ok,
        score=score,
        severity=sev,
        summary=(f"GC stddev across {n} sliding windows = {stddev:.3f} (mean = {mean:.3f})"),
        details={"stddev": stddev, "mean": mean, "n_windows": n},
    )


def check_cpg_balance(dna: str) -> CheckResult:
    """CpG density per sliding window.

    Suppressed (<0.5%) → silencing risk in some contexts.
    Excessive (>15%) → immune activation (TLR9).
    Best: a balanced profile around 1-5%.
    """
    dna = dna.upper().replace("U", "T")
    if len(dna) < CPG_WINDOW:
        return CheckResult(
            name="cpg_balance",
            pass_=True,
            score=1.0,
            severity="info",
            summary="sequence <50 nt; skipping CpG check",
            details={},
        )
    windows = []
    for i in range(0, len(dna) - CPG_WINDOW + 1, CPG_WINDOW // 2):
        win = dna[i : i + CPG_WINDOW]
        cpg = sum(1 for j in range(len(win) - 1) if win[j : j + 2] == "CG")
        windows.append(cpg / len(win))
    n = len(windows)
    if n == 0:
        return CheckResult(
            name="cpg_balance",
            pass_=True,
            score=1.0,
            severity="info",
            summary="no windows",
            details={},
        )
    mean = sum(windows) / n
    if mean < CPG_LOW:
        sev, ok, score = "warn", False, 0.6
    elif mean > CPG_HIGH:
        sev, ok, score = "warn", False, 0.5
    else:
        sev, ok, score = "info", True, 1.0
    return CheckResult(
        name="cpg_balance",
        pass_=ok,
        score=score,
        severity=sev,
        summary=f"CpG density = {mean:.3f} (target 0.5%-15%)",
        details={"mean_cpg": mean, "n_windows": n},
    )


# ---------------------------------------------------------------------------
# Aggregator
# ---------------------------------------------------------------------------


def score_manufacturability(
    cds: str,
    *,
    utr5: str = "",
    utr3: str = "",
) -> ManufacturabilityReport:
    """Run all manufacturability checks and aggregate.

    Parameters
    ----------
    cds
        Coding sequence (DNA). If it ends with a stop codon, the stop
        codon is included in the analysis.
    utr5
        5' UTR sequence (DNA). Used for Kozak scoring — only the last
        9 nt upstream of ATG are inspected.
    utr3
        3' UTR sequence (DNA). Used for AU-rich element detection.

    Returns
    -------
    ManufacturabilityReport with per-check results and an overall
    0-1 score (mean of all check scores).
    """
    cds = cds.upper().replace("U", "T")
    utr5 = (utr5 or "").upper().replace("U", "T")
    utr3 = (utr3 or "").upper().replace("U", "T")

    checks: list[CheckResult] = [
        check_poly_a_runs(cds),
        check_gc_5prime_hairpin(cds),
        check_kozak_strength(utr5) if utr5 else _skip("kozak_strength"),
        check_are_motif(utr3) if utr3 else _skip("are_motif"),
        check_stop_context(cds),
        check_hidden_stops(cds),
        check_gc_window_uniformity(cds),
        check_cpg_balance(cds),
    ]
    overall = sum(c.score for c in checks) / len(checks) if checks else 1.0
    n_pass = sum(1 for c in checks if c.pass_ and c.severity != "error")
    n_warn = sum(1 for c in checks if c.severity == "warn")
    n_error = sum(1 for c in checks if c.severity == "error")
    return ManufacturabilityReport(
        overall_score=overall,
        n_pass=n_pass,
        n_warn=n_warn,
        n_error=n_error,
        checks=checks,
    )


def _skip(name: str) -> CheckResult:
    """Build a no-op CheckResult for checks where input was not provided.

    Used by ``score_manufacturability`` to return a structured "skipped"
    result rather than skipping the check silently — keeps the report
    layout uniform across run / skip / fail.
    """
    return CheckResult(
        name=name,
        pass_=True,
        score=1.0,
        severity="info",
        summary="skipped (no input provided)",
        details={},
    )


__all__ = [
    "CheckResult",
    "ManufacturabilityReport",
    "check_poly_a_runs",
    "check_gc_5prime_hairpin",
    "check_are_motif",
    "check_kozak_strength",
    "check_stop_context",
    "check_hidden_stops",
    "check_gc_window_uniformity",
    "check_cpg_balance",
    "score_manufacturability",
]
