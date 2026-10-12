"""Offline bounds for the clock-obstruction deposit. No network and no Lean."""
import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "checks" / "publish_clock_obstruction.py"


def load():
    spec = importlib.util.spec_from_file_location("publish_clock_obstruction", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PublishBounds(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.publisher = load()

    def test_metadata_keeps_the_already_known_boundary(self):
        metadata = self.publisher.read(self.publisher.PACKAGE / "zenodo-metadata.json")
        self.publisher.validate_metadata(metadata)

    def test_source_states_the_obstruction_and_not_a_new_shortcut(self):
        self.publisher.assert_claim_boundary(self.publisher.source_text())

    def test_only_zenodo_api_urls_are_trusted(self):
        self.assertEqual(
            self.publisher.trusted("https://zenodo.org/api/records/1"),
            "https://zenodo.org/api/records/1",
        )
        for url in (
            "http://zenodo.org/api/records/1",
            "https://zenodo.org.evil.test/api/1",
            "https://files.zenodo.org/api/1",
            "https://zenodo.org/records/1",
            "https://zenodo.org/api/records/1?token=value",
        ):
            with self.assertRaises(AssertionError):
                self.publisher.trusted(url)

    def test_metadata_rejects_a_preassigned_doi_and_a_novelty_gap(self):
        metadata = self.publisher.read(self.publisher.PACKAGE / "zenodo-metadata.json")
        broken = dict(metadata)
        broken["doi"] = "10.5281/zenodo.1"
        with self.assertRaises(AssertionError):
            self.publisher.validate_metadata(broken)
        broken = dict(metadata)
        broken["description"] = broken["description"].replace("already known", "potentially new")
        with self.assertRaises(AssertionError):
            self.publisher.validate_metadata(broken)


if __name__ == "__main__":
    unittest.main()
