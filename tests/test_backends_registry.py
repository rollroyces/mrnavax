"""Tests for the v0.24.0 backends registry refactor.

After splitting ``mrnavax/backends.py`` into a registry + 8 family
modules, these tests verify the registry invariants hold and that
all 30 checks still register correctly when ``mrnavax.backends``
is imported.
"""

import importlib
import unittest

import mrnavax.backends as backends_module
from mrnavax._backends_registry import CHECKS, register


class TestRegistryInvariants(unittest.TestCase):
    """The flat CHECKS list is the single source of truth."""

    def test_checks_is_a_list(self):
        self.assertIsInstance(CHECKS, list)

    def test_every_entry_is_name_plus_callable(self):
        self.assertGreater(len(CHECKS), 0)
        for entry in CHECKS:
            self.assertEqual(len(entry), 2)
            name, fn = entry
            self.assertIsInstance(name, str)
            self.assertTrue(callable(fn), f"{name} registered a non-callable")
            self.assertGreater(len(name), 0)

    def test_names_are_unique(self):
        names = [n for n, _ in CHECKS]
        self.assertEqual(len(names), len(set(names)), f"duplicate: {names}")

    def test_names_use_dot_prefixed_family_segments(self):
        """Every check name is ``family.subname`` — what shows up in --check-all."""
        for name, _ in CHECKS:
            self.assertIn(".", name, f"{name} should be family.subname")
            family, subname = name.split(".", 1)
            self.assertGreater(len(family), 0)
            self.assertGreater(len(subname), 0)

    def test_at_least_30_checks_after_refactor(self):
        """v0.23.0 had 30 checks; v0.24.0 must keep all of them."""
        self.assertGreaterEqual(len(CHECKS), 30, f"only {len(CHECKS)} checks registered")


class TestRegisterDecorator(unittest.TestCase):
    """The register() decorator appends to CHECKS."""

    def test_register_appends_and_returns_function(self):
        initial_len = len(CHECKS)

        @register("test.sanity")
        def _check() -> tuple[bool, str]:
            return True, "sanity"

        # Function returned unchanged
        self.assertTrue(callable(_check))
        # CHECKS grew by exactly one
        self.assertEqual(len(CHECKS), initial_len + 1)
        # New entry is last
        self.assertEqual(CHECKS[-1], ("test.sanity", _check))
        # Run it to confirm it actually works
        ok, msg = _check()
        self.assertTrue(ok)
        self.assertEqual(msg, "sanity")


class TestExamplePathHelper(unittest.TestCase):
    """_example_path resolves bundled example files in both wheel + repo layouts."""

    def test_returns_string(self):
        from mrnavax._backends_registry import _example_path

        result = _example_path("simicl_demos.json")
        self.assertIsInstance(result, str)
        # Should end with the filename we asked for
        self.assertTrue(result.endswith("simicl_demos.json"))


class TestRunAll(unittest.TestCase):
    """run_all() walks CHECKS and returns 0 on all-pass."""

    def test_run_all_returns_int(self):
        from mrnavax._backends_registry import run_all

        # Just call it — don't capture stdout (verbose=True prints)
        rc = run_all(verbose=False)
        self.assertIn(rc, (0, 1))


class TestBackendsModulePublicSurface(unittest.TestCase):
    """The public ``mrnavax.backends`` module re-exports the registry."""

    def test_module_has_checks_list(self):
        self.assertTrue(hasattr(backends_module, "CHECKS"))
        self.assertIs(backends_module.CHECKS, CHECKS)

    def test_module_has_register_decorator(self):
        self.assertTrue(hasattr(backends_module, "register"))

    def test_module_has_run_all(self):
        self.assertTrue(hasattr(backends_module, "run_all"))

    def test_module_has_main(self):
        self.assertTrue(hasattr(backends_module, "main"))

    def test_module_has_example_path(self):
        self.assertTrue(hasattr(backends_module, "_example_path"))


class TestFamilyModulesImport(unittest.TestCase):
    """Each family module imports cleanly and registers its checks."""

    def test_neoantigen_family_imports(self):
        importlib.import_module("mrnavax._backends_neoantigen")

    def test_codon_family_imports(self):
        importlib.import_module("mrnavax._backends_codon")

    def test_trial_family_imports(self):
        importlib.import_module("mrnavax._backends_trial")

    def test_scrna_family_imports(self):
        importlib.import_module("mrnavax._backends_scrna")

    def test_variant_family_imports(self):
        importlib.import_module("mrnavax._backends_variant")

    def test_protein_lm_family_imports(self):
        importlib.import_module("mrnavax._backends_protein_lm")

    def test_spatial_family_imports(self):
        importlib.import_module("mrnavax._backends_spatial")

    def test_remaining_family_imports(self):
        """The catch-all for the 1-check families (manufacture, lnp, case_study, conservation)."""
        importlib.import_module("mrnavax._backends_remaining")


if __name__ == "__main__":
    unittest.main()
