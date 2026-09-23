import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class MechanismWorkflowTest(unittest.TestCase):
    def test_example_run(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temporary:
            completed = subprocess.run(
                [
                    sys.executable,
                    str(root / "workflow.py"),
                    "--input",
                    str(root / "examples" / "input.csv"),
                    "--config",
                    str(root / "config" / "workflow.json"),
                    "--output",
                    temporary,
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            output = Path(temporary)
            with (output / "results.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 6)
            self.assertGreaterEqual(float(rows[0]["composite_score"]), float(rows[-1]["composite_score"]))
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["row_count"], 6)


if __name__ == "__main__":
    unittest.main()
