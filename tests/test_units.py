"""Fast unit tests for the pipeline scripts (no external tools needed).

Run inside the ICE-detector env:  pytest
"""
import numpy as np
import pandas as pd
import pytest

import build_table
import extract_regions
import gc_content
import make_fasta_list
import summarize_orit


# --- make_fasta_list ---------------------------------------------------------

def write_fasta(path, records):
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))


def test_fasta_list_marks_contigs_circular(tmp_path):
    write_fasta(tmp_path / "g1.fasta", [("c1", "ACGT"), ("c2", "ACGT")])
    write_fasta(tmp_path / "g2.fasta", [("c3", "ACGT")])
    rows = make_fasta_list.build_list(tmp_path)
    assert [r[0] for r in rows] == ["g1", "g2"]
    assert rows[0][2:] == ["c1", "c2"]
    assert make_fasta_list.build_list(tmp_path, circular=False)[0][2:] == []


def test_fasta_list_rejects_duplicate_contig_names(tmp_path):
    write_fasta(tmp_path / "g1.fasta", [("contig_1", "ACGT")])
    write_fasta(tmp_path / "g2.fasta", [("contig_1", "ACGT")])
    with pytest.raises(SystemExit, match="unique"):
        make_fasta_list.build_list(tmp_path)


# --- extract_regions -----------------------------------------------------------

def test_parse_and_cut_coordinates():
    seq = "ABCDEFGHIJ"
    assert extract_regions.cut_sequence(seq, extract_regions.parse_coordinates("2..4")) == "BCD"
    # region wrapping the origin of a circular contig
    segs = extract_regions.parse_coordinates("9..10,1..2")
    assert extract_regions.cut_sequence(seq, segs) == "IJAB"


# --- summarize_orit ------------------------------------------------------------

def test_orit_region_mapping():
    regions = [{"region_id": "R1", "contig": "c1", "coordinates": "100..200"}]
    segs = summarize_orit.region_segments(regions)
    assert summarize_orit.find_region(segs, "c1", 190, 250) == "R1"   # overlaps
    assert summarize_orit.find_region(segs, "c1", 201, 250) == ""     # just outside
    assert summarize_orit.find_region(segs, "c2", 150, 160) == ""     # other contig


def test_orit_coverage_filter(tmp_path):
    """A strong but short hit must fail the profile-coverage filter."""
    hmm = tmp_path / "oriT.hmm"
    hmm.write_text("HMMER3/f\nNAME  P1\nLENG  100\n//\n")
    tbl_dir = tmp_path / "nhmmscan"
    tbl_dir.mkdir()
    # profile, acc, contig, acc, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to,
    # modlen, strand, evalue, score, bias, description
    (tbl_dir / "g1.tbl").write_text(
        "# comment\n"
        "P1 - c1 - 1 80 150 229 150 229 100 + 1e-10 50 0 -\n"   # 80% coverage: pass
        "P1 - c1 - 1 20 500 519 500 519 100 + 1e-10 50 0 -\n"   # 20% coverage: fail
        "P1 - c1 - 1 90 900 989 900 989 100 + 0.5 10 0 -\n")    # bad E-value: fail
    regions = tmp_path / "regions.tsv"
    pd.DataFrame([{"region_id": "R1", "contig": "c1", "coordinates": "100..1000"}]).to_csv(
        regions, sep="\t", index=False)
    import sys
    argv = sys.argv
    sys.argv = ["summarize_orit.py", "--regions", str(regions), "--nhmmscan", str(tbl_dir),
                "--hmm", str(hmm), "-o", str(tmp_path / "out")]
    try:
        summarize_orit.main()
    finally:
        sys.argv = argv
    hits = pd.read_csv(tmp_path / "out" / "orit_hits.tsv", sep="\t")
    assert hits["pass"].tolist() == [True, False, False]
    orit = pd.read_csv(tmp_path / "out" / "orit.tsv", sep="\t")
    assert orit.loc[0, "n_orit_hits"] == 1


# --- gc_content -------------------------------------------------------------------

def test_flank_gc_linear_and_circular():
    seq = "GGGGG" + "AAAAA" + "CCCCC"     # region = positions 6..10
    gc, acgt = gc_content.counts(seq)
    up, down, both = gc_content.flank_gc(gc, acgt, 6, 10, 5, circular=False)
    assert (up, down, both) == (1.0, 1.0, 1.0)
    # flank larger than the contig: linear truncates, circular wraps
    up, _, _ = gc_content.flank_gc(gc, acgt, 1, 5, 3, circular=False)
    assert up is None
    up, _, _ = gc_content.flank_gc(gc, acgt, 1, 5, 3, circular=True)
    assert up == 1.0                      # wraps to the last 3 bases (CCC)


def test_window_gc():
    gc, acgt = gc_content.counts("GGAA" * 10)
    win = gc_content.window_gc(gc, acgt, 4, 4, circular=False)
    assert np.allclose(win, 0.5)


# --- build_table --------------------------------------------------------------------

@pytest.mark.parametrize("t4ss,relaxase,orit,dconj,expected", [
    (True, True, True, False, ("ICE", "high")),
    (True, True, False, False, ("ICE", "medium")),
    (True, False, False, False, ("ICE", "low")),
    (False, False, False, True, ("decayed_ICE", "")),
    (False, True, True, False, ("mobilizable", "")),
    (False, False, True, False, ("none", "")),
])
def test_classify(t4ss, relaxase, orit, dconj, expected):
    row = {"t4ss": t4ss, "relaxase": relaxase, "orit": orit, "dconj_system": dconj}
    assert build_table.classify(row) == expected


# --- direct_repeats -------------------------------------------------------------------

import argparse
import random

import direct_repeats


def random_dna(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_maximal_matches_extends_to_full_length():
    a = "TTTT" + "ACGTGCATGCAAGT" + "CCCC"
    b = "GGGG" + "ACGTGCATGCAAGT" + "AAAA"
    assert direct_repeats.maximal_matches(a, b, 10) == [(14, 4, 4)]


def test_low_complexity():
    assert direct_repeats.low_complexity("AT" * 10)
    assert not direct_repeats.low_complexity("ACGTGCATGCAAGTCCGA")


def test_find_dr_in_flanks_of_inserted_element():
    """Simulated integration: att copy before and after a 5 kb element."""
    att = "GCTCAGTTGGTAGAGCGCAACC"                       # 22 bp
    left, element, right = random_dna(10000, 1), random_dna(5000, 2), random_dna(10000, 3)
    genome = left + att + element + att + right
    start = len(left) + len(att) + 1                     # element (region) start, 1-based
    end = start + len(element) - 1
    args = argparse.Namespace(outside=3000, inside=3000, min_len=18, max_copies=2)
    res = direct_repeats.find_dr(genome, [genome], start, end, args, circular=False)
    assert res["dr"] and res["dr_length"] >= len(att)
    assert res["dr_left_offset"] < 0 and res["dr_right_offset"] < 0   # both copies outside
    assert res["element_start"] <= start - len(att) + 1


def test_find_dr_rejects_genome_wide_repeat():
    """A repeat with a third copy elsewhere (e.g. REP element) is not a DR."""
    att = "GCTCAGTTGGTAGAGCGCAACC"
    left, element, right = random_dna(10000, 1), random_dna(5000, 2), random_dna(10000, 3)
    genome = left + att + element + att + right + att
    start = len(left) + len(att) + 1
    end = start + len(element) - 1
    args = argparse.Namespace(outside=3000, inside=3000, min_len=18, max_copies=2)
    assert not direct_repeats.find_dr(genome, [genome], start, end, args, circular=False)["dr"]


def test_find_dr_wraps_circular_contig():
    """Right att copy split across the origin: found only if the contig is circular."""
    att = "GCTCAGTTGGTAGAGCGCAACC"
    left, element, right = random_dna(10000, 1), random_dna(5000, 2), random_dna(1000, 3)
    linear = left + att + element + att + right
    cut = len(left) + len(att) + len(element) + 10      # inside the second att copy
    genome = linear[cut:] + linear[:cut]                 # rotate the circular contig
    start = genome.find(element) + 1
    end = start + len(element) - 1
    args = argparse.Namespace(outside=3000, inside=3000, min_len=18, max_copies=2)
    assert direct_repeats.find_dr(genome, [genome], start, end, args, circular=True)["dr"]
    assert not direct_repeats.find_dr(genome, [genome], start, end, args, circular=False)["dr"]
