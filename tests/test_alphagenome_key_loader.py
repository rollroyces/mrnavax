"""Tests for the AlphaGenome API key loader (v0.25.1).

Verifies the resolution order: explicit > env var > helper file (with
silent-skip on missing/unreadable). Locks in the v0.25.1 contract so a
future refactor that drops the helper-file path would turn these red.
"""

from __future__ import annotations

import os
import unittest
from pathlib import Path
from unittest import mock

from mrnavax.alphagenome_integration import _load_api_key


class TestLoadApiKeyPrecedence(unittest.TestCase):
    """Resolution order: explicit > env > helper files."""

    def setUp(self):
        self._saved_env = os.environ.get("ALPHAGENOME_API_KEY")
        os.environ.pop("ALPHAGENOME_API_KEY", None)

    def tearDown(self):
        if self._saved_env is not None:
            os.environ["ALPHAGENOME_API_KEY"] = self._saved_env
        else:
            os.environ.pop("ALPHAGENOME_API_KEY", None)

    def test_explicit_wins(self):
        """Caller-supplied key takes precedence over everything else."""
        os.environ["ALPHAGENOME_API_KEY"] = "env-key"
        self.assertEqual(_load_api_key("explicit-key"), "explicit-key")

    def test_env_var_when_no_explicit(self):
        os.environ["ALPHAGENOME_API_KEY"] = "env-key"
        self.assertEqual(_load_api_key(), "env-key")

    def test_helper_file_when_no_env_var(self):
        """Falls back to ~/projects/alphagenome-work/.alphagenome_key."""
        helper_path = (
            Path.home() / "projects" / "alphagenome-work" / ".alphagenome_key"
        )
        # Patch Path.home() so we don't depend on the actual file
        fake_home = Path("/tmp/_mrnavax_test_home")
        fake_helper = fake_home / "projects" / "alphagenome-work" / ".alphagenome_key"
        fake_helper.parent.mkdir(parents=True, exist_ok=True)
        fake_helper.write_text("helper-key\n")
        try:
            with mock.patch.object(Path, "home", return_value=fake_home):
                self.assertEqual(_load_api_key(), "helper-key")
        finally:
            fake_helper.unlink()
            fake_helper.parent.rmdir()
            fake_helper.parent.parent.rmdir()
            # Verify the real helper file was NOT touched
            if helper_path.exists():
                self.assertNotEqual(
                    helper_path.read_text().strip(), "helper-key",
                    "test should not have written to the real helper file",
                )

    def test_alternate_home_helper_when_first_missing(self):
        """Falls back to ~/.alphagenome_key if the projects/ one is absent."""
        fake_home = Path("/tmp/_mrnavax_test_home2")
        # First helper is absent; second one is present
        first = fake_home / "projects" / "alphagenome-work" / ".alphagenome_key"
        first.parent.mkdir(parents=True, exist_ok=True)
        if first.exists():
            first.unlink()
        first.parent.rmdir()
        first.parent.parent.rmdir()
        second = fake_home / ".alphagenome_key"
        second.write_text("home-key\n")
        try:
            with mock.patch.object(Path, "home", return_value=fake_home):
                self.assertEqual(_load_api_key(), "home-key")
        finally:
            second.unlink()

    def test_empty_string_when_nothing_found(self):
        """Returns '' when no key is available anywhere."""
        fake_home = Path("/tmp/_mrnavax_test_home3")
        fake_home.mkdir(parents=True, exist_ok=True)
        # Ensure no helper files exist in fake_home
        for p in (
            fake_home / "projects" / "alphagenome-work" / ".alphagenome_key",
            fake_home / ".alphagenome_key",
        ):
            if p.exists():
                p.unlink()
        try:
            with mock.patch.object(Path, "home", return_value=fake_home):
                self.assertEqual(_load_api_key(), "")
        finally:
            # Clean up the empty parent dirs we created
            for p in (
                fake_home / "projects" / "alphagenome-work",
                fake_home / "projects",
            ):
                if p.exists() and not any(p.iterdir()):
                    p.rmdir()

    def test_env_var_stripped(self):
        """Trailing whitespace in env var is stripped."""
        os.environ["ALPHAGENOME_API_KEY"] = "  env-key  \n"
        self.assertEqual(_load_api_key(), "env-key")

    def test_helper_file_stripped(self):
        """Trailing whitespace in helper file is stripped."""
        fake_home = Path("/tmp/_mrnavax_test_home4")
        fake_home.mkdir(parents=True, exist_ok=True)
        helper = fake_home / "projects" / "alphagenome-work" / ".alphagenome_key"
        helper.parent.mkdir(parents=True, exist_ok=True)
        helper.write_text("  helper-key  \n")
        try:
            with mock.patch.object(Path, "home", return_value=fake_home):
                self.assertEqual(_load_api_key(), "helper-key")
        finally:
            helper.unlink()
            helper.parent.rmdir()
            helper.parent.parent.rmdir()


class TestAdapterUsesLoader(unittest.TestCase):
    """AlphaGenomeCLIAdapter must accept keys from any resolution path."""

    def setUp(self):
        self._saved_env = os.environ.get("ALPHAGENOME_API_KEY")
        os.environ.pop("ALPHAGENOME_API_KEY", None)

    def tearDown(self):
        if self._saved_env is not None:
            os.environ["ALPHAGENOME_API_KEY"] = self._saved_env
        else:
            os.environ.pop("ALPHAGENOME_API_KEY", None)

    def test_adapter_loads_from_env(self):
        os.environ["ALPHAGENOME_API_KEY"] = "env-key"
        from mrnavax.alphagenome_integration import AlphaGenomeCLIAdapter
        adapter = AlphaGenomeCLIAdapter()
        self.assertEqual(adapter.api_key, "env-key")

    def test_adapter_loads_from_helper_file(self):
        fake_home = Path("/tmp/_mrnavax_test_home5")
        fake_home.mkdir(parents=True, exist_ok=True)
        helper = fake_home / "projects" / "alphagenome-work" / ".alphagenome_key"
        helper.parent.mkdir(parents=True, exist_ok=True)
        helper.write_text("helper-key")
        try:
            with mock.patch.object(Path, "home", return_value=fake_home):
                from mrnavax.alphagenome_integration import AlphaGenomeCLIAdapter
                adapter = AlphaGenomeCLIAdapter()
                self.assertEqual(adapter.api_key, "helper-key")
        finally:
            helper.unlink()
            helper.parent.rmdir()
            helper.parent.parent.rmdir()

    def test_adapter_explicit_overrides_env(self):
        os.environ["ALPHAGENOME_API_KEY"] = "env-key"
        from mrnavax.alphagenome_integration import AlphaGenomeCLIAdapter
        adapter = AlphaGenomeCLIAdapter(api_key="explicit-key")
        self.assertEqual(adapter.api_key, "explicit-key")

    def test_adapter_raises_when_no_key_anywhere(self):
        fake_home = Path("/tmp/_mrnavax_test_home6")
        fake_home.mkdir(parents=True, exist_ok=True)
        try:
            with mock.patch.object(Path, "home", return_value=fake_home):
                from mrnavax.alphagenome_integration import AlphaGenomeCLIAdapter
                with self.assertRaises(ValueError) as ctx:
                    AlphaGenomeCLIAdapter()
                # The error message must mention all three resolution paths
                msg = str(ctx.exception)
                self.assertIn("ALPHAGENOME_API_KEY", msg)
                self.assertIn(".alphagenome_key", msg)
                self.assertIn("api_key=", msg)
        finally:
            for p in (
                fake_home / "projects" / "alphagenome-work",
                fake_home / "projects",
            ):
                if p.exists() and not any(p.iterdir()):
                    p.rmdir()


class TestSelectorUsesLoader(unittest.TestCase):
    """select_regulatory_scorer returns real adapter when loader finds a key."""

    def setUp(self):
        self._saved_env = os.environ.get("ALPHAGENOME_API_KEY")
        os.environ.pop("ALPHAGENOME_API_KEY", None)

    def tearDown(self):
        if self._saved_env is not None:
            os.environ["ALPHAGENOME_API_KEY"] = self._saved_env
        else:
            os.environ.pop("ALPHAGENOME_API_KEY", None)

    def test_selector_returns_mock_when_no_key(self):
        fake_home = Path("/tmp/_mrnavax_test_home7")
        fake_home.mkdir(parents=True, exist_ok=True)
        try:
            with mock.patch.object(Path, "home", return_value=fake_home):
                from mrnavax.alphagenome_integration import (
                    MockRegulatoryVariantScorer,
                    select_regulatory_scorer,
                )
                sel = select_regulatory_scorer()
                self.assertIsInstance(sel, MockRegulatoryVariantScorer)
        finally:
            for p in (
                fake_home / "projects" / "alphagenome-work",
                fake_home / "projects",
            ):
                if p.exists() and not any(p.iterdir()):
                    p.rmdir()

    def test_selector_returns_real_when_env_set(self):
        os.environ["ALPHAGENOME_API_KEY"] = "env-key"
        from mrnavax.alphagenome_integration import (
            AlphaGenomeCLIAdapter,
            select_regulatory_scorer,
        )
        sel = select_regulatory_scorer()
        self.assertIsInstance(sel, AlphaGenomeCLIAdapter)


if __name__ == "__main__":
    unittest.main()
