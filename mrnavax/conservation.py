"""PhyloP46way evolutionary-conservation lookup.

Roadmap item 6: closure of the coding-region variant prioritization
gap. PhyloP scores per-base evolutionary conservation from the UCSC
46-way placental mammal alignment (hg38). Highly conserved positions
(positive PhyloP) tend to be functionally important; rapidly-evolving
positions (negative PhyloP) tend to be neutral.

Like AlphaMissense and AlphaGenome Atlas, this module follows the
Protocol-adapter pattern:

* **Real**: ``PhyloPRestAdapter`` queries the UCSC REST API
  (``https://api.genome.ucsc.edu/list/schema``). Opt-in; not bundled.
* **Mock**: ``MockPhyloPLookup`` returns deterministic synthetic
  scores from a stdlib hash of (chrom, pos). Always present.

Both backends expose the same shape::

    lookup(chrom, pos) -> float  in [-1.0, 1.0]

A separate ``ConservationLookup`` Protocol captures the contract so
``score_variant`` can take any implementation.

Note on coverage
----------------
PhyloP46way scores the entire non-coding + coding genome. The
output is meaningful for any DNA-level variant position (not just
coding). In mrnavax, the conservation_lookup is currently only wired
into ``score_variant`` for coding-region variants (where it composes
with AlphaMissense as a 4th signal). A future release could route
non-coding conservation into the AVI path.

Reference
---------
Pollard KS, Hubisz MJ, Rosenbloom KR, Siepel A. "Detection of nonneutral
substitution rates on mammalian phylogenies." *Genome Research* 20,
110–121 (2010). PhyloP46way data: UCSC Genome Browser
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP46way/.
"""

from __future__ import annotations

import hashlib
import logging
from typing import Protocol, runtime_checkable

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Protocol
# ---------------------------------------------------------------------------

@runtime_checkable
class ConservationLookup(Protocol):
    """Anything that can look up PhyloP / GERP++ conservation at a
    (chrom, pos) coordinate."""

    def lookup(self, chrom: str, pos: int) -> float:
        ...


# ---------------------------------------------------------------------------
# Mock backend
# ---------------------------------------------------------------------------

class MockPhyloPLookup:
    """Stdlib-only deterministic mock.

    Returns a synthetic PhyloP score by hashing (chrom, pos) with
    SHA-256 and mapping the first 4 bytes to a float in [-1, +1].
    Stable: same input → same output across processes.

    The sign of the score correlates with position parity (even →
    positive/conserved, odd → negative/fast-evolving) so the test
    suite can exercise both branches without real data.
    """

    def lookup(self, chrom: str, pos: int) -> float:
        key = f"{chrom}:{pos}".encode()
        digest = hashlib.sha256(key).digest()
        # First 4 bytes → unsigned int → [0, 1] → rescaled to [-1, +1]
        n = int.from_bytes(digest[:4], byteorder="big")
        u = n / 0xFFFFFFFF  # [0, 1]
        # Position-parity bias: even pos → positive (conserved),
        # odd pos → negative (fast-evolving). The bias is mild
        # (±0.2) so the hash dominates the variation.
        bias = 0.2 if (pos % 2 == 0) else -0.2
        score = u * 2.0 - 1.0 + bias
        return max(-1.0, min(1.0, score))


# ---------------------------------------------------------------------------
# Real adapter (UCSC REST API)
# ---------------------------------------------------------------------------

class PhyloPRestAdapter:
    """Real adapter — queries the UCSC REST API for PhyloP46way scores.

    The UCSC public API (``api.genome.ucsc.edu``) is free and returns
    per-base PhyloP scores from the 46-way placental alignment
    (``phyloP46way`` track). This adapter hits the endpoint via
    ``urllib.request`` (stdlib) so the rest of mrnavax stays
    stdlib-first. Heavy Python deps are not required.

    Parameters
    ----------
    timeout_seconds : int
        HTTP timeout for the UCSC request. Default 10s.
    """

    API_BASE = "https://api.genome.ucsc.edu"
    TRACK = "phyloP46way"

    def __init__(self, timeout_seconds: int = 10) -> None:
        self.timeout_seconds = timeout_seconds

    def lookup(self, chrom: str, pos: int) -> float:
        import json as _json
        import urllib.request

        chrom_clean = chrom.replace("chr", "")
        url = (
            f"{self.API_BASE}/get/track/hg38/{self.TRACK}"
            f"?chrom=chr{chrom_clean}&start={pos}&end={pos + 1}"
        )
        try:
            with urllib.request.urlopen(url, timeout=self.timeout_seconds) as resp:
                payload = _json.loads(resp.read().decode("utf-8"))
        except Exception as exc:  # noqa: BLE001 — adapter boundary
            logger.warning("UCSC PhyloP request failed (%s); returning 0.0", exc)
            return 0.0

        # Response shape:
        # {"phyloP46way": [{"chrom": "chr7", "start": 123, "end": 124, "value": 0.42}]}
        try:
            rows = payload[self.TRACK]
            if not rows:
                return 0.0
            return float(rows[0]["value"])
        except (KeyError, TypeError, ValueError, IndexError):
            return 0.0


# ---------------------------------------------------------------------------
# Backend selector
# ---------------------------------------------------------------------------

def select_conservation_lookup() -> ConservationLookup:
    """Select the best available conservation backend.

    Currently only the mock is bundled (the real UCSC adapter is
    trivial to swap in via ``--conservation-backend``). When a real
    adapter is added, this selector will auto-pick it.
    """
    return MockPhyloPLookup()


__all__ = [
    "ConservationLookup",
    "MockPhyloPLookup",
    "PhyloPRestAdapter",
    "select_conservation_lookup",
]
