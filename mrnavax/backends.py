"""Backend availability + integrity checks.

Used by CI to verify that the optional backend wiring is intact even
when the heavy models themselves aren't installed. Heavy model
downloads (mhcflurry, scGPT) are gated by their respective extras.

Run directly:

    python -m mrnavax.backends --check-all

Exits 0 if every check passes, 1 otherwise.

Architecture (v0.24.0)
----------------------
This module used to hold all 30 backend checks inline (1,838 lines).
The check functions are now split into family-specific modules under
``mrnavax/_backends_<family>.py`` and registered with the
``mrnavax._backends_registry`` module. ``CHECKS`` is the single
source of truth; importing this module populates the list.

Public surface preserved:

- ``CHECKS`` — flat list of (name, callable) tuples.
- ``register(name)`` — decorator to add a check.
- ``_example_path(name)`` — resolve bundled example file paths.
- ``run_all(verbose=True)`` — run every check, return 0/1.
- ``main()`` — CLI entry point (``--check-all``).
"""

from __future__ import annotations

# Importing these side-effect: registers each family's checks with
# the global CHECKS list. Import order doesn't matter for the registry
# because all checks are appended during module load.
from . import (  # noqa: F401
    _backends_codon,
    _backends_neoantigen,
    _backends_protein_lm,
    _backends_remaining,
    _backends_scrna,
    _backends_spatial,
    _backends_trial,
    _backends_variant,
)

# Re-export the public surface so existing callers
# (``from mrnavax.backends import CHECKS, run_all, register``)
# continue to work unchanged.
from ._backends_registry import (
    CHECKS,
    _example_path,
    main,
    register,
    run_all,
)

__all__ = ["CHECKS", "register", "_example_path", "run_all", "main"]


if __name__ == "__main__":
    import sys

    sys.exit(main())
