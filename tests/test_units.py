"""Fast unit tests for the pipeline scripts (no external tools needed).

Run inside the ICE-detector env:  pytest
"""
import numpy as np
import pandas as pd
import pytest

import build_table
import common
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


def test_extract_chromosomes_keeps_largest_contig_only(tmp_path):
    inp, out = tmp_path / "in", tmp_path / "chrom"
    inp.mkdir()
    write_fasta(inp / "g1.fna", [("p1", "ACGT"), ("chr1", "ACGTACGTAC\nGTAC"), ("p2", "AC")])
    write_fasta(inp / "g2.fasta", [("chr2", "AAAA"), ("p3", "CCCC")])  # tie: first wins
    before = (inp / "g1.fna").read_text()
    report = make_fasta_list.extract_chromosomes(inp, out)
    assert (inp / "g1.fna").read_text() == before                 # input untouched
    assert (out / "g1.fasta").read_text() == ">chr1\nACGTACGTAC\nGTAC\n"
    assert [r["kept_contig"] for r in report] == ["chr1", "chr2"]
    assert report[0]["kept_length"] == 14 and report[0]["assembly_length"] == 20
    assert report[0]["dropped_contigs"] == "p1:4,p2:2"
    assert [r[2:] for r in make_fasta_list.build_list(out)] == [["chr1"], ["chr2"]]


def test_extract_chromosomes_allows_same_plasmid_names_across_genomes(tmp_path):
    inp = tmp_path / "in"
    inp.mkdir()
    write_fasta(inp / "g1.fasta", [("chr1", "ACGTACGT"), ("plasmid1", "ACG")])
    write_fasta(inp / "g2.fasta", [("chr2", "ACGTACGT"), ("plasmid1", "ACG")])
    make_fasta_list.extract_chromosomes(inp, tmp_path / "chrom")
    assert len(make_fasta_list.build_list(tmp_path / "chrom")) == 2


# --- extract_regions -----------------------------------------------------------

def test_parse_and_cut_coordinates():
    seq = "ABCDEFGHIJ"
    assert extract_regions.cut_sequence(seq, common.parse_coordinates("2..4")) == "BCD"
    # region wrapping the origin of a circular contig
    segs = common.parse_coordinates("9..10,1..2")
    assert segs == [(9, 10), (1, 2)]
    assert extract_regions.cut_sequence(seq, segs) == "IJAB"


# --- common ---------------------------------------------------------------------

def test_tsv_roundtrip_and_comments(tmp_path):
    path = tmp_path / "x.tsv"
    common.write_tsv(path, [{"a": 1, "b": "x"}, {"a": 2}])        # missing b -> ""
    assert path.read_text() == "a\tb\n1\tx\n2\t\n"
    path.write_text("# note\na\tb\n1\t2\n")
    assert common.read_tsv(path, comment="#") == [{"a": "1", "b": "2"}]


def test_read_fasta_list(tmp_path):
    path = tmp_path / "list.tsv"
    path.write_text("g1\t/x/g1.fasta\tc1\tc2\ng2\t/x/g2.fasta\n\n")
    genomes, circular = common.read_fasta_list(path)
    assert genomes == {"g1": "/x/g1.fasta", "g2": "/x/g2.fasta"}
    assert circular == {"c1", "c2"}


def test_region_segments_wrapping_region():
    segs = common.region_segments([{"region_id": "R", "contig": "c", "coordinates": "90..100,1..10"}])
    assert common.find_region(segs, "c", 5, 6) == "R"
    assert common.find_region(segs, "c", 50, 60) == ""


# --- summarize_orit ------------------------------------------------------------

def test_orit_region_mapping():
    regions = [{"region_id": "R1", "contig": "c1", "coordinates": "100..200"}]
    segs = common.region_segments(regions)
    assert common.find_region(segs, "c1", 190, 250) == "R1"   # overlaps
    assert common.find_region(segs, "c1", 201, 250) == ""     # just outside
    assert common.find_region(segs, "c2", 150, 160) == ""     # other contig


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

@pytest.mark.parametrize("t4ss,relaxase,orit,dconj,expected,genes", [
    (True, True, True, False, "ICE", "T4SS+MOB+oriT"),
    (True, True, True, True, "ICE", "T4SS+MOB+oriT"),       # complete system wins over dCONJ
    (False, True, True, True, "dICE", "dCONJ+MOB+oriT"),
    (False, True, True, False, "IME", "MOB+oriT"),
    (False, False, True, False, "iOriT", "oriT"),
    (True, True, False, False, "ICE_noOriT", "T4SS+MOB"),
    (False, True, False, True, "dICE_noOriT", "dCONJ+MOB"),
    (False, True, False, False, "IME_noOriT", "MOB"),
    (True, False, True, False, "other", "T4SS+oriT"),
    (False, False, False, True, "other", "dCONJ"),
    (False, False, False, False, "none", "-"),
])
def test_classify(t4ss, relaxase, orit, dconj, expected, genes):
    row = {"t4ss": t4ss, "relaxase": relaxase, "orit": orit, "dconj_system": dconj}
    assert build_table.classify(row) == expected
    assert build_table.mobility_genes(row) == genes


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
    """A high-copy repeat (REP-like, 6 copies) at both ends is not a DR, even if
    chance bases next to the two flanking copies make the match longer."""
    rep = "GCTCAGTTGGTAGAGCGCAACC"
    left, element, right = random_dna(10000, 1), random_dna(5000, 2), random_dna(10000, 3)
    # force a 2 bp chance extension: same 2 bases before both flanking copies
    left, element = left[:-2] + "GA", element[:-2] + "GA"
    elsewhere = "".join(random_dna(500, 10 + i) + rep for i in range(4))
    genome = left + rep + element + rep + right + elsewhere
    start = len(left) + len(rep) + 1
    end = start + len(element) - 1
    args = argparse.Namespace(outside=3000, inside=3000, min_len=18, max_copies=2)
    assert not direct_repeats.find_dr(genome, [genome], start, end, args, circular=False)["dr"]


def test_find_dr_allows_one_near_copy_elsewhere():
    """A real att whose 2bp-trimmed core occurs once more elsewhere is still a DR."""
    att = "GCTCAGTTGGTAGAGCGCAACC"
    left, element, right = random_dna(10000, 1), random_dna(5000, 2), random_dna(10000, 3)
    genome = left + att + element + att + right + att[2:]   # partial third copy
    start = len(left) + len(att) + 1
    end = start + len(element) - 1
    args = argparse.Namespace(outside=3000, inside=3000, min_len=18, max_copies=2)
    assert direct_repeats.find_dr(genome, [genome], start, end, args, circular=False)["dr"]


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


def brute_copies(seq, genome):
    """Slow reference: count every trimmed variant on both strands position by position."""
    from Bio.Seq import reverse_complement as rc

    def count(v):
        strands = {v, rc(v)}
        return sum(genome[p:p + len(v)] in strands for p in range(len(genome) - len(v) + 1))

    n = len(seq)
    trimmed = [seq[i:n - j] for i in range(3) for j in range(3) if (i, j) != (0, 0)]
    return count(seq), max(count(v) for v in trimmed)


@pytest.mark.parametrize("seq", [
    "GCTCAGTTGGTAGAGCGCAACC",     # ordinary DR
    "GGGCCCGGGCCCGGGCCC",         # reverse-palindromic, trimmed versions are not
    "AAGCTTGCATGCAAGCTT",
])
def test_dr_copies_matches_brute_force(seq):
    from Bio.Seq import reverse_complement as rc
    genome = (random_dna(3000, 4) + seq + random_dna(2000, 5) + rc(seq) + random_dna(2000, 6)
              + seq[2:] + random_dna(1000, 7) + seq[:-1] + random_dna(500, 8))
    assert direct_repeats.dr_copies(seq, [genome]) == brute_copies(seq, genome)


# --- summarize_integrase ----------------------------------------------------------------

import summarize_integrase
import trna_proximity


def domtbl_line(gene, profile, acc, evalue, ievalue):
    # target acc tlen query acc qlen E score bias # of c-E i-E dscore dbias
    # hmm_from hmm_to ali_from ali_to env_from env_to acc description
    return (f"{gene} - 300 {profile} {acc} 170 {evalue} 50.0 0.1 1 1 {ievalue} {ievalue} "
            f"50.0 0.1 1 160 10 170 5 175 0.95 -\n")


def test_integrase_roles_regions_and_end_distance(tmp_path):
    tbl = tmp_path / "g1.domtbl"
    tbl.write_text("# header\n"
                   + domtbl_line("gA", "Phage_integrase", "PF00589.28", 1e-30, 1e-30)
                   + domtbl_line("gA", "Phage_integrase", "PF00589.28", 1e-30, 1e-5)  # 2nd domain
                   + domtbl_line("gB", "rve", "PF00665.34", 1e-20, 1e-20)
                   + domtbl_line("gC", "Resolvase", "PF00239.28", 1e-20, 1e-20)
                   + domtbl_line("gX", "Phage_integrase", "PF00589.28", 1e-30, 1e-30))  # no region
    hits = summarize_integrase.read_domtbl(tbl)
    assert len(hits) == 4                                   # best domain per gene x profile
    regions = [{"region_id": "R1", "coordinates": "1000..1999", "length": "1000"},
               {"region_id": "R2", "coordinates": "5000..5999", "length": "1000"}]
    genes = [{"region_id": "R1", "gene": "gA", "start": "1500", "stop": "1800"},  # 199 bp to end
             {"region_id": "R2", "gene": "gB", "start": "5000", "stop": "5300"},
             {"region_id": "R2", "gene": "gC", "start": "5010", "stop": "5300"}]  # 10 bp to start
    rows = {r["region_id"]: r for r in summarize_integrase.summarize(regions, genes, hits)}
    assert rows["R1"]["integrase"] and rows["R1"]["integrase_types"] == "Tyr"
    assert rows["R1"]["integrase_end_dist"] == 199 and not rows["R1"]["dde_transposase"]
    assert rows["R2"]["integrase_types"] == "Ser" and rows["R2"]["dde_transposase"]
    assert rows["R2"]["integrase_end_dist"] == 10
    assert {h["gene"]: h["region_id"] for h in hits}["gX"] == ""


def test_integrase_dde_alone_is_not_integrase():
    hits = [{"gene": "g", "accession": "PF00665"}]
    rows = summarize_integrase.summarize(
        [{"region_id": "R", "coordinates": "1..100", "length": "100"}],
        [{"region_id": "R", "gene": "g", "start": "1", "stop": "90"}], hits)
    assert not rows[0]["integrase"] and rows[0]["dde_transposase"]


def test_end_distance_in_region_wrapping_origin():
    segs = [(9001, 10000), (1, 500)]                       # 1500 bp, wraps the origin
    assert summarize_integrase.end_distance({"start": "400", "stop": "450"}, segs, 1500) == 50
    assert summarize_integrase.end_distance({"start": "9100", "stop": "9200"}, segs, 1500) == 99


# --- trna_proximity ---------------------------------------------------------------------

TRNA = {"type": "tRNA", "name": "tRNA-Gly(gcc)", "start": 1000, "end": 1075, "strand": "+"}


def test_trna_near_region_end():
    res = trna_proximity.nearest_trna([TRNA], 1300, 50000, 100000, False, 1000)
    assert res["trna_adjacent"] and res["trna_end"] == "left" and res["trna_dist"] == 225
    far = trna_proximity.nearest_trna([TRNA], 3000, 50000, 100000, False, 1000)
    assert not far["trna_adjacent"]


def test_trna_overlapping_dr_counts_even_if_far():
    dr = {"dr": "True", "dr_left_start": "1060", "dr_left_end": "1074",
          "dr_right_start": "60000", "dr_right_end": "60014"}
    res = trna_proximity.nearest_trna([TRNA], 3000, 50000, 100000, False, 1000, dr)
    assert res["trna_adjacent"] and res["trna_dr_overlap"]


def test_trna_distance_wraps_circular_contig():
    t = {**TRNA, "start": 99900, "end": 99975}
    assert trna_proximity.nearest_trna([t], 100, 5000, 100000, True, 1000)["trna_adjacent"]
    assert not trna_proximity.nearest_trna([t], 100, 5000, 100000, False, 1000)["trna_adjacent"]


# --- build_table end to end on tiny tables ---------------------------------------------

def test_build_table_counts_seven_evidence_items(tmp_path):
    res = tmp_path / "run"
    def tsv(rel, rows):
        path = res / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    rid = ["R1", "R2"]
    tsv("regions/regions.tsv", [{"region_id": r, "spot_id": "spot_1", "genome": "g", "contig": "c",
                                 "start": 1, "end": 100, "length": 100, "n_orfs": 10,
                                 "coordinates": "1..100", "contig_border": False,
                                 "whole_contig": False} for r in rid])
    tsv("conjugation/conjugation.tsv", [
        {"region_id": "R1", "conj_system": True, "conj_systems": "x", "dconj_system": False,
         "dconj_systems": "", "relaxase": True, "relaxase_types": "MOBH", "t4cp": True,
         "virb4": True, "n_conj_hits": 20},
        {"region_id": "R2", "conj_system": False, "conj_systems": "", "dconj_system": False,
         "dconj_systems": "", "relaxase": False, "relaxase_types": "", "t4cp": False,
         "virb4": False, "n_conj_hits": 0}])
    tsv("orit/orit.tsv", [{"region_id": r, "orit": r == "R1", "n_orit_hits": 1,
                           "orit_profiles": "", "orit_best_evalue": ""} for r in rid])
    tsv("gc/gc.tsv", [{"region_id": r, "region_gc": .5, "genome_gc": .6, "gc_diff": -.1,
                       "gc_zscore": -3 if r == "R1" else 0, "flank_up_gc": .6,
                       "flank_down_gc": .6, "flank_gc": .6, "gc_diff_flank": -.1} for r in rid])
    tsv("dr/dr.tsv", [{**dict.fromkeys(direct_repeats.FIELDS), "region_id": r, "dr": r == "R1"}
                      for r in rid])
    tsv("integrase/integrase.tsv", [{"region_id": r, "integrase": True, "integrase_types": "Tyr",
                                     "n_integrases": 1, "integrase_end_dist": 0,
                                     "dde_transposase": False} for r in rid])
    tsv("trna/trna.tsv", [{"region_id": r, "trna_adjacent": r == "R1", "trna_name": "",
                           "trna_end": "", "trna_dist": None, "trna_dr_overlap": False}
                          for r in rid])
    tsv("ppanggolin/spots.tsv", [{"spot_id": "spot_1", "rgp_id": r} for r in rid])
    tsv("ppanggolin/regions_of_genomic_plasticity.tsv", [{"region": r, "genome": "g"} for r in rid])
    (res / "regions/fasta").mkdir()
    for r in rid:
        (res / "regions/fasta" / f"{r}.fasta").write_text(f">{r}\nACGT\n")

    import sys
    argv = sys.argv
    sys.argv = ["build_table.py", "--results", str(res)]
    try:
        build_table.main()
    finally:
        sys.argv = argv
    out = pd.read_csv(res / "final/all_regions.tsv", sep="\t").set_index("region_id")
    assert out.loc["R1", "evidence_count"] == 7 and out.loc["R1", "category"] == "ICE"
    assert out.loc["R2", "evidence_count"] == 1 and out.loc["R2", "category"] == "none"
    cand = pd.read_csv(res / "final/ice_candidates.tsv", sep="\t")
    assert cand["region_id"].tolist() == ["R1"]   # R2: integrase only, no mobility marker
