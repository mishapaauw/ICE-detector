"""End-to-end test: full pipeline on data/test must reproduce the expected candidates.

Slow (PPanGGOLiN + MacSyFinder + nhmmscan). Run inside the ICE-detector env:
    pytest -m slow
Needs the reference data in db/ (see db/README.md).
"""
import subprocess
import sys

import pandas as pd
import pytest

from conftest import ROOT

EXPECTED = ROOT / "tests" / "expected" / "test_candidates.tsv"
COLS = ["region_id", "category", "mobility_genes", "evidence_count"]


@pytest.mark.slow
def test_pipeline_reproduces_expected_candidates(tmp_path):
    if not (ROOT / "db" / "macsy_models" / "CONJScan").exists() or \
            not (ROOT / "db" / "orit" / "oriT.hmm").exists():
        pytest.skip("reference data missing in db/ (see db/README.md)")

    subprocess.run([sys.executable, ROOT / "scripts" / "run_all.py", ROOT / "data" / "test",
                    "-o", tmp_path, "--cpu", "4"], check=True)

    got = pd.read_csv(tmp_path / "final" / "ice_candidates.tsv", sep="\t", usecols=COLS)
    exp = pd.read_csv(EXPECTED, sep="\t", comment="#")
    norm = lambda d: d[COLS].fillna("").astype(str).sort_values("region_id").reset_index(drop=True)
    pd.testing.assert_frame_equal(norm(got), norm(exp))
