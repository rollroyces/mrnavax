"""Backend check registry — flat list of (name, callable) tuples.

Extracted from ``mrnavax/backends.py`` in v0.24.0. Family-specific
check functions live in ``mrnavax/_backends_<family>.py`` modules.
This module owns:
  * the ``CHECKS`` list (single source of truth, populated by
    ``register()`` as each family module imports + decorates)
  * the ``register()`` decorator
  * the ``_example_path()`` helper (used by codon / variant checks)
  * the ``run_all()`` CLI runner
  * the ``main()`` entry point

Stdlib-only. The check function bodies are NOT here; they live in
the family modules and import this module to register themselves.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

# Single source of truth: populated by register() decorator as each
# family module imports + decorates. Importing the family modules
# from ``mrnavax.backends`` is what populates this list.
CHECKS: list[tuple[str, Callable[[], tuple[bool, str]]]] = []


def register(name: str):
    """Decorator that appends ``(name, fn)`` to ``CHECKS``.

    Each family module uses this at module load time:

        from ._backends_registry import register
        @register("neoantigen.heuristic_A0201")
        def _check_neoantigen_heuristic() -> tuple[bool, str]:
            ...
    """

    def deco(fn: Callable[[], tuple[bool, str]]) -> Callable[[], tuple[bool, str]]:
        CHECKS.append((name, fn))
        return fn

    return deco


def _example_path(name: str) -> str:
    """Resolve an example-file path relative to the installed package.

    Works in two layouts:
    - Installed wheel: pkg_dir = site-packages/mrnavax, examples
      are bundled next to it under site-packages/mrnavax/examples.
    - Dev / repo: pkg_dir is inside the repo, examples live at the
      repo root under mrnavax/examples.

    Returns the first path that exists, falling back to the wheel
    layout (which will raise FileNotFoundError downstream if
    missing — that's the expected behavior for missing bundled
    data).
    """
    pkg_dir = Path(__file__).resolve().parent
    wheel_candidate = pkg_dir / "examples" / name
    if wheel_candidate.exists():
        return str(wheel_candidate)
    repo_candidate = pkg_dir.parent.parent / "mrnavax" / "examples" / name
    if repo_candidate.exists():
        return str(repo_candidate)
    return str(wheel_candidate)


def run_all(verbose: bool = True) -> int:
    """Run every registered check. Returns 0 if all pass, 1 otherwise.

    Prints one line per check with PASS/FAIL + the check's diagnostic
    message. ``verbose=False`` prints only the summary line.
    """
    failed: list[str] = []
    for name, fn in CHECKS:
        try:
            ok, msg = fn()
        except Exception as e:  # noqa: BLE001
            ok, msg = False, f"raised {type(e).__name__}: {e}"
        if verbose:
            status = "PASS" if ok else "FAIL"
            print(f"{status:8s} {name:50s} {msg}")
        if not ok:
            failed.append(name)
    total = len(CHECKS)
    if failed:
        print(f"\n{len(failed)}/{total} checks failed")
        for n in failed:
            print(f"  - {n}")
        return 1
    print(f"\n{total}/{total} checks passed")
    return 0


def main(argv: list[str] | None = None) -> int:
    """CLI entry point: ``python -m mrnavax.backends --check-all``."""
    argv = argv if argv is not None else sys.argv[1:]
    if "--check-all" in argv:
        return run_all(verbose=True)
    print("Usage: python -m mrnavax.backends --check-all")
    return 2


__all__ = ["CHECKS", "register", "_example_path", "run_all", "main"]
