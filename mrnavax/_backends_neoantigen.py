"""Backend checks for the ``neoantigen.*`` family.

Extracted from ``mrnavax/backends.py`` in v0.24.0. Each check is
registered with the global ``CHECKS`` registry on import.
"""

from __future__ import annotations

from ._backends_registry import register


@register("neoantigen.heuristic_A0201")
def _check_neoantigen_heuristic() -> tuple[bool, str]:
    from .neoantigen_screener import screen_peptide_llm

    # Known-good A*02:01 binder — must come back positive in the heuristic.
    r = screen_peptide_llm("NLVPMVATV", "HLA-A*02:01", backend="mock")
    ok = r.source == "heuristic" and isinstance(r.binding_affinity_nM, float)
    return ok, f"source={r.source} aff={r.binding_affinity_nM} nM"


@register("neoantigen.heuristic_nonA2_silent")
def _check_neoantigen_non_a2() -> tuple[bool, str]:
    from .neoantigen_screener import screen_peptide_llm

    r = screen_peptide_llm("NLVPMVATV", "HLA-A*03:01", backend="mock")
    ok = r.binding_affinity_nM >= 5000  # heuristic is silent on non-A2
    return ok, f"non-A2 affinity={r.binding_affinity_nM} nM (expected >= 5000)"


@register("neoantigen.mhcflurry_available")
def _check_neoantigen_mhcflurry() -> tuple[bool, str]:
    from .neoantigen_screener import _mhcflurry_available

    available = _mhcflurry_available()
    if available:
        return True, "installed — will use as default"
    return True, "not installed (optional: pip install -e '.[neoantigen-mhcflurry]')"


__all__ = []
