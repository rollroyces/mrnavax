"""Live AlphaGenome Atlas integration harness.

Roadmap item 5b: a weekly GitHub Actions cron job calls the real
Atlas API with a stable test variant, records the response, and
compares against a baseline. Catches upstream breakage in real time
(the v0.18.0 fixture catches it in 7-day lag).

Three modes:

* ``mock`` — uses the stdlib mock; always passes; no network.
  Used in unit tests + offline CI.
* ``record`` — calls the real Atlas; writes the response to a
  JSON file. Run manually when refreshing the fixture.
* ``regression`` — calls the real Atlas; compares against a
  baseline score with ±5% tolerance. Run weekly via the
  ``.github/workflows/atlas_integration.yml`` cron.

CLI::

    python -m mrnavax.live_atlas_integration --mode mock
    python -m mrnavax.live_atlas_integration --mode record \
        --output tests/fixtures/alphagenome_atlas_live.json
    python -m mrnavax.live_atlas_integration --mode regression \
        --baseline 0.72 --tolerance 0.05

The script exits 0 on pass, 1 on failure — so the workflow's
``exit 1`` semantic catches failures cleanly.

Reference
---------
Avsec et al., "Advancing regulatory variant effect prediction
with AlphaGenome," *Nature* (2026).
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Stable test variant
#
# BRAF V600E coding-region. Stable hg38 position across releases;
# cited by Avsec et al. as a benchmark variant in the Atlas paper.
# ---------------------------------------------------------------------------

DEFAULT_VARIANT: dict = {
    "chrom": "chr7",
    "pos": 140753336,
    "ref": "T",
    "alt": "A",
}


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LiveAtlasReport:
    """Result of one live Atlas call.

    Attributes
    ----------
    chrom / pos / ref / alt : str / int / str / str
        Variant coordinates queried.
    score : float
        AVI score in [0, 1] (or 0.0 if the call failed).
    classification : str
        "low" / "moderate" / "high" (or "" on failure).
    is_coding : bool
        True if Atlas reports the variant as coding.
    captured_at : str
        ISO date / datetime string for the run.
    source : str
        "mock" / "alphagenome" — which backend produced the score.
    within_tolerance : bool
        Only populated in regression mode. True if the live score
        is within ``tolerance`` of the baseline.
    delta : float
        Only populated in regression mode. ``score - baseline``.
    baseline : float | None
        The baseline score passed to ``--baseline`` (None in mock / record modes).
    tolerance : float | None
        The tolerance passed to ``--tolerance`` (None in mock / record modes).
    error : str
        Non-empty if the call failed (timeout, auth, network, etc.).
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    score: float
    classification: str
    is_coding: bool
    captured_at: str
    source: str
    within_tolerance: bool = True
    delta: float = 0.0
    baseline: float | None = None
    tolerance: float | None = None
    error: str = ""

    def to_dict(self) -> dict:
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "score": self.score,
            "classification": self.classification,
            "is_coding": self.is_coding,
            "captured_at": self.captured_at,
            "source": self.source,
            "within_tolerance": self.within_tolerance,
            "delta": self.delta,
            "baseline": self.baseline,
            "tolerance": self.tolerance,
            "error": self.error,
        }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def compare_to_baseline(
    score: float, baseline: float, tolerance: float
) -> dict:
    """Compare a live score against a baseline, returning pass/fail.

    Parameters
    ----------
    score : float
        Live AVI score.
    baseline : float
        Recorded baseline AVI score.
    tolerance : float
        Maximum absolute delta (|score - baseline|) tolerated.

    Returns
    -------
    dict with keys: within_tolerance (bool), delta (float),
    tolerance (float), baseline (float), score (float).
    """
    delta = score - baseline
    within = abs(delta) <= tolerance
    return {
        "within_tolerance": within,
        "delta": delta,
        "tolerance": tolerance,
        "baseline": baseline,
        "score": score,
    }


# ---------------------------------------------------------------------------
# Score functions per source
# ---------------------------------------------------------------------------

def _score_via_mock(chrom: str, pos: int, ref: str, alt: str, **kwargs) -> dict:
    """Score a variant via the stdlib mock. No network."""
    from .alphagenome_integration import MockRegulatoryVariantScorer

    mock = MockRegulatoryVariantScorer()
    r = mock.score_variant(chrom, pos, ref, alt)
    return {
        "score": r.score,
        "classification": r.classification,
        "is_coding": r.is_coding,
    }


def _score_via_alphagenome(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    *,
    api_key: str,
    timeout_seconds: int = 30,
) -> dict:
    """Score a variant via the AlphaGenome Atlas subprocess adapter.

    Runs the bundled ``mrnavax/_shims/alphagenome_cli.py`` shim as a
    subprocess. Stdlib only.
    """
    if not api_key:
        return {
            "score": 0.0,
            "classification": "",
            "is_coding": False,
            "error": "ALPHAGENOME_API_KEY not set",
        }
    import subprocess

    payload = json.dumps(
        {
            "api_key": api_key,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
        }
    )
    shim = Path(__file__).resolve().parent / "_shims" / "alphagenome_cli.py"
    try:
        proc = subprocess.run(
            [sys.executable, str(shim)],
            input=payload,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired:
        return {
            "score": 0.0,
            "classification": "",
            "is_coding": False,
            "error": f"timeout after {timeout_seconds}s",
        }
    except Exception as exc:  # noqa: BLE001 — adapter boundary
        return {
            "score": 0.0,
            "classification": "",
            "is_coding": False,
            "error": str(exc),
        }
    if proc.returncode != 0:
        return {
            "score": 0.0,
            "classification": "",
            "is_coding": False,
            "error": proc.stderr.strip()[:200] if proc.stderr else f"rc={proc.returncode}",
        }
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        return {
            "score": 0.0,
            "classification": "",
            "is_coding": False,
            "error": f"malformed JSON: {exc}",
        }


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_live_atlas_check(
    mode: str = "mock",
    variant: dict | None = None,
    *,
    api_key: str = "",
    baseline: float | None = None,
    tolerance: float = 0.05,
    output_path: Path | None = None,
    score_fn: Callable[..., dict] | None = None,
) -> LiveAtlasReport:
    """Run one live Atlas integration check.

    Parameters
    ----------
    mode : str
        "mock" / "record" / "regression".
    variant : dict
        {"chrom": str, "pos": int, "ref": str, "alt": str}.
        Defaults to BRAF V600E.
    api_key : str
        AlphaGenome Atlas API key. Required for ``record`` and
        ``regression`` modes; ignored in ``mock`` mode.
    baseline : float | None
        Baseline AVI score. Required for ``regression`` mode.
    tolerance : float
        Maximum absolute delta tolerated in regression mode (default 0.05).
    output_path : Path | None
        Required for ``record`` mode — where to write the captured JSON.
    score_fn : Callable | None
        Test seam: pass a stub function for unit tests. Signature
        ``() -> dict`` with keys ``score`` / ``classification`` /
        ``is_coding``. When None, the function is selected based on
        ``mode``: mock → _score_via_mock, others → _score_via_alphagenome.

    Returns
    -------
    LiveAtlasReport
    """
    from datetime import datetime, timezone

    if variant is None:
        variant = DEFAULT_VARIANT

    # Pick the score function if not provided. Use a non-shadowing name.
    if score_fn is None:
        if mode == "mock":

            def _default_score_fn() -> dict:
                return _score_via_mock(**variant)

        else:

            def _default_score_fn() -> dict:
                return _score_via_alphagenome(**variant, api_key=api_key)

        score_fn = _default_score_fn

    # Capture
    captured_at = datetime.now(timezone.utc).isoformat()
    try:
        result = score_fn()
    except Exception as exc:  # noqa: BLE001 — adapter boundary
        result = {"score": 0.0, "classification": "", "is_coding": False,
                  "error": str(exc)}

    # Build report
    report = LiveAtlasReport(
        chrom=variant["chrom"],
        pos=variant["pos"],
        ref=variant["ref"],
        alt=variant["alt"],
        score=float(result.get("score", 0.0)),
        classification=str(result.get("classification", "")),
        is_coding=bool(result.get("is_coding", False)),
        captured_at=captured_at,
        source="mock" if mode == "mock" else "alphagenome",
        baseline=baseline,
        tolerance=tolerance if mode == "regression" else None,
        error=str(result.get("error", "")),
    )

    # Regression comparison
    if mode == "regression" and baseline is not None:
        cmp = compare_to_baseline(report.score, baseline, tolerance)
        object.__setattr__(report, "within_tolerance", cmp["within_tolerance"])
        object.__setattr__(report, "delta", cmp["delta"])

    # Record mode writes the captured JSON
    if mode == "record":
        if output_path is None:
            raise ValueError("record mode requires --output path")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(
                {
                    "chrom": report.chrom,
                    "pos": report.pos,
                    "ref": report.ref,
                    "alt": report.alt,
                    "score": report.score,
                    "classification": report.classification,
                    "is_coding": report.is_coding,
                    "captured_at": report.captured_at,
                    "source": report.source,
                    "error": report.error,
                },
                f,
                indent=2,
            )
        logger.info("Recorded live Atlas response to %s", output_path)

    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="mrnavax.live_atlas_integration",
        description="Live AlphaGenome Atlas integration harness (mock / record / regression).",
    )
    parser.add_argument(
        "--mode",
        choices=["mock", "record", "regression"],
        default="mock",
        help="Execution mode (default: mock)",
    )
    parser.add_argument(
        "--variant",
        default=None,
        help="Variant as chr:pos:ref:alt (default: BRAF V600E chr7:140753336:T:A)",
    )
    parser.add_argument(
        "--baseline",
        type=float,
        default=None,
        help="Baseline AVI score (regression mode only)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.05,
        help="Maximum absolute delta tolerated in regression mode (default 0.05)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output JSON path (record mode only)",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    # Parse variant if supplied
    variant = None
    if args.variant:
        parts = args.variant.split(":")
        if len(parts) != 4:
            print(f"error: --variant must be chr:pos:ref:alt, got {args.variant!r}",
                  file=sys.stderr)
            return 2
        variant = {
            "chrom": parts[0],
            "pos": int(parts[1]),
            "ref": parts[2],
            "alt": parts[3],
        }

    # Run
    api_key = os.environ.get("ALPHAGENOME_API_KEY", "")
    if args.mode != "mock" and not api_key:
        print(
            f"error: ALPHAGENOME_API_KEY not set; required for {args.mode} mode",
            file=sys.stderr,
        )
        return 2

    report = run_live_atlas_check(
        mode=args.mode,
        variant=variant,
        api_key=api_key,
        baseline=args.baseline,
        tolerance=args.tolerance,
        output_path=args.output,
    )

    # Output
    print(json.dumps(report.to_dict(), indent=2))

    # Exit code: 0 on pass, 1 on regression failure / error
    if report.error:
        return 1
    if args.mode == "regression" and not report.within_tolerance:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
