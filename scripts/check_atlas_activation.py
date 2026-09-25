"""Verify Atlas API key activation after the user adds the GitHub secret.

Run this locally AFTER adding the ``ALPHAGENOME_API_KEY`` secret at
https://github.com/rollroyces/mrnavax/settings/secrets/actions/new to
confirm the activation worked end-to-end. The script:

  1. Calls ``_load_api_key()`` with explicit precedence (arg > env >
     helper file > ~/.alphagenome_key) and prints which path won.
  2. Imports the ``alphagenome`` package from the [variant-alphagenome]
     extra (or skips gracefully if not installed).
  3. Calls ``DnaClient.score_variant()`` against a known canonical
     variant (BRAF V600E) and reports the score + classification.
  4. Compares the result against the live fixture recorded in
     ``tests/fixtures/alphagenome_atlas_live.json`` (BRAF V600E should
     score=1.0, classification="high", is_coding=true).

If the score matches, activation is verified — the live Atlas cron
will start pulling real weekly Atlas data on its next run.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Make ``mrnavax`` importable when running this script directly from
# a fresh checkout (without `pip install -e .`).
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

REPO_ROOT = _REPO_ROOT
FIXTURE_PATH = REPO_ROOT / "tests" / "fixtures" / "alphagenome_atlas_live.json"

EXPECTED_KEY_PREFIXES = ("AIza",)  # Google AI Studio API keys start with this


def _check_helper_file() -> None:
    helper = Path("/Users/hermes/projects/alphagenome-work/.alphagenome_key")
    if helper.exists():
        n_bytes = helper.stat().st_size
        # Read the first 4 bytes without leaking the key (print length only)
        with helper.open() as f:
            first4 = f.read(4)
        print(
            f"[helper-file] OK: {helper} ({n_bytes} bytes, "
            f"prefix matches Google AI key: {first4.startswith(EXPECTED_KEY_PREFIXES)})"
        )
    else:
        print(
            f"[helper-file] NOT FOUND at {helper}; falling back to "
            f"~/.alphagenome_key or env var"
        )


def _check_alphagenome_package() -> bool:
    try:
        import alphagenome  # noqa: F401
    except ImportError:
        print(
            "[package] alphagenome not installed locally. "
            "Install with: pip install 'mrnavax[variant-alphagenome]'"
        )
        return False
    print(f"[package] alphagenome installed")
    return True


def _check_live_atlas_call(api_key: str) -> dict | None:
    """Call the alphagenome_cli shim on BRAF V600E; return result dict.

    The shim is a stdin/stdout JSON CLI script. We invoke it as a
    subprocess with the canonical BRAF V600E payload.
    """
    import subprocess

    payload = json.dumps({
        "api_key": api_key,
        "chrom": "chr7",
        "pos": 140753336,
        "ref": "T",
        "alt": "A",
    })

    # Path to the shim — import as module to honor Python path config.
    shim_path = REPO_ROOT / "mrnavax" / "_shims" / "alphagenome_cli.py"
    if not shim_path.exists():
        print(f"[live-call] shim not found at {shim_path}")
        return None

    print(f"[live-call] BRAF V600E chr7:140753336 T>A via {shim_path.name}")
    proc = subprocess.run(
        [sys.executable, str(shim_path)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if proc.returncode != 0:
        print(f"[live-call] FAILED (exit {proc.returncode}): {proc.stderr}")
        return None
    try:
        result = json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        print(f"[live-call] FAILED (bad JSON): {exc}")
        print(f"  stdout: {proc.stdout[:200]}")
        return None

    score = result.get("score", 0.0)
    classification = result.get("classification", "unknown")
    is_coding = result.get("is_coding", False)
    n_scorers = result.get("n_scorers", 0)
    print(
        f"[live-call] score={score:.3f}, classification={classification!r}, "
        f"is_coding={is_coding}, n_scorers={n_scorers}"
    )
    return {
        "score": score,
        "classification": classification,
        "is_coding": is_coding,
    }


def _compare_to_fixture(live: dict) -> bool:
    if not FIXTURE_PATH.exists():
        print(f"[fixture] no fixture at {FIXTURE_PATH} (cannot compare)")
        return False
    with FIXTURE_PATH.open() as f:
        fixture = json.load(f)
    fixture_classification = fixture.get("classification")
    fixture_score = fixture.get("score", 0.0)
    matches = (
        live["classification"] == fixture_classification
        and abs(live["score"] - fixture_score) < 0.05
    )
    print(
        f"[fixture] recorded={fixture_score:.3f}/{fixture_classification!r}, "
        f"live={live['score']:.3f}/{live['classification']!r}, "
        f"matches={matches}"
    )
    return matches


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Verify Atlas API key activation after GitHub secret setup."
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="Override the API key (skips helper-file lookup). "
        "Used for CI / testing.",
    )
    parser.add_argument(
        "--skip-live-call",
        action="store_true",
        help="Skip the actual Atlas API call (just verify key resolution).",
    )
    args = parser.parse_args(argv)

    print("=" * 60)
    print("Atlas API key activation verification")
    print("=" * 60)

    # 1. Helper file check
    _check_helper_file()

    # 2. Package check
    has_package = _check_alphagenome_package()

    # 3. Key resolution
    from mrnavax.alphagenome_integration import _load_api_key

    resolved = _load_api_key(args.api_key)
    if resolved is None:
        print("[key] FAILED to load API key from any source")
        print(
            "[key] Set it via: export ALPHAGENOME_API_KEY=... "
            "or write to ~/projects/alphagenome-work/.alphagenome_key"
        )
        return 1
    n_chars = len(resolved)
    print(f"[key] OK: resolved {n_chars}-char key (prefix matches Google AI: "
          f"{resolved.startswith(EXPECTED_KEY_PREFIXES)})")

    # 4. Live call (only if package installed AND --skip-live-call not set)
    if has_package and not args.skip_live_call:
        live = _check_live_atlas_call(resolved)
        if live is None:
            print("[result] live call failed")
            return 2
        matches = _compare_to_fixture(live)
        if matches:
            print("[result] PASS: live Atlas matches recorded fixture")
            print("[result] Activation verified — weekly Atlas cron will work")
            return 0
        print("[result] FAIL: live result does not match fixture")
        print("[result] Atlas API contract may have changed upstream")
        print("[result] Open an issue with `gh issue new` and the live output")
        return 3

    # Package not installed — just verify key resolution
    print("[result] PASS (key-resolution only; install [variant-alphagenome] for live call)")
    return 0


if __name__ == "__main__":
    sys.exit(main())