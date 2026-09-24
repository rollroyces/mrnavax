"""Helpers for safely selecting optional real-model adapters.

The ``_safe_selector`` helper provides a uniform pattern for "try the
real adapter, fall back to a mock, swallow any ImportError or
RuntimeError". This pattern is used in three places in the codebase
(AlphaMissense, AlphaGenome Atlas, PhyloP46way) and was previously
duplicated as three near-identical functions.

New in v0.29.0 — additive refactor only. Existing callers
(``_try_am_lookup``, ``_try_avi_lookup``, ``_try_conservation_lookup``)
are preserved as thin wrappers around ``_safe_selector``; this means
no public API change.
"""

from __future__ import annotations

from typing import Any, Callable

__all__ = ["_safe_selector"]


def _safe_selector(
    import_path: str,
    factory_name: str,
    *,
    factory_args: tuple = (),
    factory_kwargs: dict | None = None,
    extra_setup: Callable[[Any], None] | None = None,
) -> tuple[Any | None, bool]:
    """Try to import ``import_path`` and call ``factory_name()`` on it.

    Returns ``(result, is_real)``:
        - On success: ``(instance_or_callable, True)``.
        - On any exception (ImportError, RuntimeError, etc.):
            ``(None, False)``.

    ``factory_args`` and ``factory_kwargs`` are passed to the factory.
    If ``extra_setup`` is supplied, it is called with the result before
    returning (e.g. for warming caches via ``load_index()``).

    The "swallow any exception" policy is intentional: the caller
    decides what to do with ``is_real == False`` (typically: fall
    back to the mock implementation). Surfacing every adapter-import
    failure would create noisy tracebacks for the common case where
    a user hasn't installed the optional extra.
    """
    if factory_kwargs is None:
        factory_kwargs = {}
    try:
        import importlib

        module = importlib.import_module(import_path, package="mrnavax")
        factory = getattr(module, factory_name)
        result = factory(*factory_args, **factory_kwargs)
        if extra_setup is not None:
            extra_setup(result)
        return result, True
    except Exception:
        return None, False
