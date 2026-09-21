from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from polyland.data import build_processed_tables
from polyland.table1 import table1_statistics


ROOT = Path(__file__).resolve().parents[1]


class DataLineageTests(unittest.TestCase):
    def test_processed_counts_and_corrected_ch4(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            build_processed_tables(ROOT / "data" / "raw", output)
            stats = table1_statistics(output)
            self.assertEqual(int(stats.loc["Count", ("CH4", "Ladder")]), 50)
            self.assertAlmostEqual(float(stats.loc["Mean", ("CH4", "Ladder")]), 156.1964)
            self.assertAlmostEqual(float(stats.loc["Median", ("CH4", "Ladder")]), 22.5)
            self.assertAlmostEqual(float(stats.loc["Min", ("CH4", "Ladder")]), 0.3)
            self.assertAlmostEqual(float(stats.loc["Max", ("CH4", "Ladder")]), 710.0)

    def test_manuscript_counts(self) -> None:
        stats = table1_statistics(ROOT / "data" / "processed")
        expected = {
            "O2": (489, 46),
            "N2": (486, 51),
            "H2": (311, 45),
            "CH4": (425, 50),
            "CO2": (466, 50),
        }
        for gas, (linear, ladder) in expected.items():
            self.assertEqual(int(stats.loc["Count", (gas, "Linear")]), linear)
            self.assertEqual(int(stats.loc["Count", (gas, "Ladder")]), ladder)

    def test_ladder_n2_and_ch4_are_not_duplicates(self) -> None:
        n2 = pd.read_csv(ROOT / "data" / "processed" / "final_ladder_data_N2.csv")
        ch4 = pd.read_csv(ROOT / "data" / "processed" / "final_ladder_data_CH4.csv")
        paired = n2[["PID", "SMILES", "N2"]].merge(
            ch4[["PID", "SMILES", "CH4"]], on=["PID", "SMILES"], how="inner"
        )
        self.assertLess((paired["N2"] == paired["CH4"]).mean(), 0.25)


if __name__ == "__main__":
    unittest.main()

