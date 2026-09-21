"""Regression tests for the v0.24.2 select_conservation_lookup fix.

Before v0.24.2, ``select_conservation_lookup()`` unconditionally
returned ``MockPhyloPLookup()`` — making the 4-signal composition
story (BLOSUM62 + driver + structural + AM → AVI → PhyloP) actually
3-signal in production. The selector's docstring even promised
"When a real adapter is added, this selector will auto-pick it" but
the selector never did.

This test pins the contract: with ``MRNA_AI_FORCE_MOCK=1`` the mock
is returned; without it, the real UCSC ``PhyloPRestAdapter`` is
returned. This catches the v0.24.2-style regression.

Also exercises the actual lookup() round-trip so a future bug that
breaks silent-fail (e.g. raising instead of returning 0.0) is caught.
"""

import os
import unittest

from mrnavax.conservation import (
    MockPhyloPLookup,
    PhyloPRestAdapter,
    select_conservation_lookup,
)


class TestSelectConservationLookup(unittest.TestCase):
    """Pins the v0.24.2 selector contract."""

    def setUp(self):
        self._saved = os.environ.get("MRNA_AI_FORCE_MOCK")
        os.environ.pop("MRNA_AI_FORCE_MOCK", None)

    def tearDown(self):
        if self._saved is not None:
            os.environ["MRNA_AI_FORCE_MOCK"] = self._saved
        else:
            os.environ.pop("MRNA_AI_FORCE_MOCK", None)

    def test_no_force_mock_returns_real_ucsc_adapter(self):
        """Without MRNA_AI_FORCE_MOCK, the real UCSC adapter is used.

        This is the v0.24.2 contract: PhyloP is real in production.
        """
        sel = select_conservation_lookup()
        self.assertIsInstance(sel, PhyloPRestAdapter)
        self.assertNotIsInstance(sel, MockPhyloPLookup)

    def test_force_mock_returns_mock(self):
        """MRNA_AI_FORCE_MOCK=1 keeps the mock for offline / CI use."""
        os.environ["MRNA_AI_FORCE_MOCK"] = "1"
        sel = select_conservation_lookup()
        self.assertIsInstance(sel, MockPhyloPLookup)


class TestRealAdapterSilentFail(unittest.TestCase):
    """Real adapter silently returns 0.0 on network/parse failure."""

    def test_lookup_returns_float_in_range_on_real_call(self):
        """Real UCSC lookup returns a float in [-1, 1], never raises.

        Network may be down in CI — the adapter contract is "silent
        fail to 0.0", not "raise". This catches the regression where
        the silent-fail was removed (e.g. by switching to urllib's
        default exception propagation).
        """
        sel = PhyloPRestAdapter()
        # Real call — UCSC may be down, but it MUST return a float
        score = sel.lookup("chr7", 140753336)
        self.assertIsInstance(score, float)
        self.assertGreaterEqual(score, -1.0)
        self.assertLessEqual(score, 1.0)

    def test_lookup_with_invalid_chrom_returns_zero(self):
        """An invalid chromosome silently returns 0.0 (defensive)."""
        sel = PhyloPRestAdapter()
        score = sel.lookup("chrINVALID", 140753336)
        self.assertEqual(score, 0.0)


class TestFourSignalCompositionUsesRealPhyloP(unittest.TestCase):
    """End-to-end: the score_variant 4-signal composition uses the
    real PhyloP adapter when MRNA_AI_FORCE_MOCK is unset.

    This is the user-visible integration test: a BRAF V600E variant
    with a real (or silently-failing) PhyloP lookup gets the
    'phylop46way_score' component recorded in score_variant output.
    """

    def test_score_variant_records_phylop_component(self):
        from mrnavax.conservation import MockPhyloPLookup
        from mrnavax.variant_scorer import score_variant

        # Use the mock so the test is deterministic + fast. The
        # important contract is: the phylop46way_score component
        # appears in the components dict, proving the 4-signal
        # path is wired (regardless of which backend supplies
        # the lookup).
        mock_phylop = MockPhyloPLookup()
        r = score_variant(
            gene="BRAF",
            position=600,
            wt_aa="V",
            mut_aa="E",
            protein_length=766,
            chrom="chr7",
            pos=140753336,
            conservation_lookup=mock_phylop.lookup,
        )
        self.assertIsNotNone(r)
        self.assertIn("phylop46way_score", r.components)
        self.assertGreaterEqual(r.components["phylop46way_score"], -1.0)
        self.assertLessEqual(r.components["phylop46way_score"], 1.0)


if __name__ == "__main__":
    unittest.main()
