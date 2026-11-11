#!/usr/bin/env python
"""Merge all evidence into the final ICE table (milestone 5).

Per region (RGP with >= --min-orfs genes), evidence is:
  t4ss       conjugation system (CONJScan T4SS model, >= 50% of its genes in region)
  relaxase   MOB relaxase gene in the region
  orit       oriT hit passing the E-value and profile-coverage filters
  gc         atypical GC: |gc_zscore| >= --gc-z (default 2)
  dr         direct repeats (putative attL/attR) found at both region ends
evidence_count = number of these five that are True.

category (from oriT, relaxase and conjugation system):
  ICE    oriT + relaxase + complete conjugation system (T4SS model)
  dICE   oriT + relaxase + decayed conjugation system (dCONJ model), no complete one
  IME    oriT + relaxase, no conjugation system
  iOriT  oriT only: no relaxase, no conjugation system
  ICE_noOriT, dICE_noOriT, IME_noOriT
         as ICE / dICE / IME, but no oriT found (the oriT profiles cover only
         known oriT families, so many real elements will land here)
  other  any other combination of mobility genes (see mobility_genes, e.g.
         "T4SS" = conjugation system without relaxase)
  none   no oriT, relaxase or conjugation system
mobility_genes: which of T4SS, dCONJ, MOB (relaxase), oriT were found, joined by "+".

Outputs (in --outdir):
  all_regions.tsv     every region with all evidence columns
  ice_candidates.tsv  regions with evidence_count >= --min-evidence (default 2)
                      AND at least one mobility marker (t4ss, relaxase or orit),
                      so ICEs and IME candidates, but not e.g. prophages that
                      only have atypical GC + DRs. Sorted by category,
                      and evidence
  ice_fasta/          FASTA of each candidate region (copied from regions/fasta)
element_start/element_end: outer edges of the DR pair, when found; a refined
estimate of the element boundaries (region start/end come from PPanGGOLiN).
  versions.txt        settings + versions.txt of every upstream step
"""
import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

from common import write_versions

CATEGORY_ORDER = {"ICE": 0, "dICE": 1, "IME": 2, "iOriT": 3, "ICE_noOriT": 4,
                  "dICE_noOriT": 5, "IME_noOriT": 6, "other": 7, "none": 8}


def mobility_genes(row):
    dconj = row["dconj_system"] and not row["t4ss"]
    found = [("T4SS", row["t4ss"]), ("dCONJ", dconj),
             ("MOB", row["relaxase"]), ("oriT", row["orit"])]
    return "+".join(name for name, present in found if present) or "-"


def classify(row):
    t4ss, relaxase, orit = row["t4ss"], row["relaxase"], row["orit"]
    dconj = row["dconj_system"] and not t4ss
    if relaxase:
        base = "ICE" if t4ss else "dICE" if dconj else "IME"
        return base if orit else f"{base}_noOriT"
    if orit and not (t4ss or dconj):
        return "iOriT"
    if t4ss or dconj or relaxase:
        return "other"
    return "none"


def spot_genome_counts(spots_tsv, rgp_tsv):
    spots = pd.read_csv(spots_tsv, sep="\t")
    rgps = pd.read_csv(rgp_tsv, sep="\t", usecols=["region", "genome"])
    m = spots.merge(rgps, left_on="rgp_id", right_on="region")
    return m.groupby("spot_id")["genome"].nunique().rename("spot_n_genomes")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True,
                    help="run folder containing ppanggolin/, regions/, conjugation/, orit/, gc/")
    ap.add_argument("-o", "--outdir", help="default: <results>/final")
    ap.add_argument("--min-evidence", type=int, default=2,
                    help="min evidence_count for ice_candidates.tsv (default: 2)")
    ap.add_argument("--gc-z", type=float, default=2.0,
                    help="|GC z-score| at or above which GC counts as atypical (default: 2)")
    args = ap.parse_args()

    res = Path(args.results)
    out = Path(args.outdir) if args.outdir else res / "final"
    out.mkdir(parents=True, exist_ok=True)

    regions = pd.read_csv(res / "regions" / "regions.tsv", sep="\t")
    conj = pd.read_csv(res / "conjugation" / "conjugation.tsv", sep="\t")
    orit = pd.read_csv(res / "orit" / "orit.tsv", sep="\t")
    gc = pd.read_csv(res / "gc" / "gc.tsv", sep="\t")
    dr = pd.read_csv(res / "dr" / "dr.tsv", sep="\t")

    df = (regions.merge(conj, on="region_id", how="left", validate="1:1")
                 .merge(orit, on="region_id", how="left", validate="1:1")
                 .merge(gc, on="region_id", how="left", validate="1:1")
                 .merge(dr, on="region_id", how="left", validate="1:1"))
    if len(df) != len(regions):
        sys.exit("Region tables do not match; re-run upstream steps.")

    df = df.rename(columns={"conj_system": "t4ss"})
    for col in ("t4ss", "dconj_system", "relaxase", "orit", "t4cp", "virb4", "dr"):
        df[col] = df[col].fillna(False).astype(bool)
    df["gc_atypical"] = df["gc_zscore"].abs() >= args.gc_z
    df["evidence_count"] = df[["t4ss", "relaxase", "orit", "gc_atypical", "dr"]].sum(axis=1)
    for col in ("dr_length", "dr_genome_copies", "dr_left_offset", "dr_right_offset",
                "dr_left_start", "dr_left_end", "dr_right_start", "dr_right_end",
                "element_start", "element_end"):
        df[col] = df[col].astype("Int64")
    df["category"] = df.apply(classify, axis=1)
    df["mobility_genes"] = df.apply(mobility_genes, axis=1)
    df = df.join(spot_genome_counts(res / "ppanggolin" / "spots.tsv",
                                    res / "ppanggolin" / "regions_of_genomic_plasticity.tsv"),
                 on="spot_id")
    df["spot_n_genomes"] = df["spot_n_genomes"].fillna(1).astype(int)

    cols = ["region_id", "genome", "contig", "start", "end", "length", "n_orfs",
            "spot_id", "spot_n_genomes", "category", "mobility_genes", "evidence_count",
            "t4ss", "relaxase", "orit", "gc_atypical", "dr",
            "element_start", "element_end",
            "conj_systems", "dconj_systems", "relaxase_types", "t4cp", "virb4",
            "orit_profiles", "orit_best_evalue", "region_gc", "genome_gc", "gc_diff",
            "gc_zscore", "flank_up_gc", "flank_down_gc", "flank_gc", "gc_diff_flank",
            "dr_length", "dr_seq", "dr_genome_copies", "dr_left_start", "dr_left_end",
            "dr_right_start", "dr_right_end", "dr_left_offset", "dr_right_offset",
            "coordinates", "contig_border"]
    df = df[cols].sort_values(
        ["category", "evidence_count", "genome", "start"],
        key=lambda s: s.map(CATEGORY_ORDER) if s.name == "category" else s,
        ascending=[True, False, True, True])
    df.to_csv(out / "all_regions.tsv", sep="\t", index=False)

    mobility = df["t4ss"] | df["relaxase"] | df["orit"]
    ice = df[(df["evidence_count"] >= args.min_evidence) & mobility]
    ice.to_csv(out / "ice_candidates.tsv", sep="\t", index=False)
    fasta_out = out / "ice_fasta"
    if fasta_out.exists():
        shutil.rmtree(fasta_out)
    fasta_out.mkdir()
    for rid in ice["region_id"]:
        shutil.copy(res / "regions" / "fasta" / f"{rid}.fasta", fasta_out)

    write_versions(out / "versions.txt", pandas=pd.__version__, gc_z=args.gc_z,
                   min_evidence=args.min_evidence)
    with open(out / "versions.txt", "a") as fh:
        for step in ("ppanggolin", "regions", "conjugation", "orit", "gc", "dr"):
            v = res / step / "versions.txt"
            if v.exists():
                fh.write(f"\n## {step}\n{v.read_text()}")
        f = res / "orit" / "orit_filters.txt"
        if f.exists():
            fh.write(f"\n## orit filters\n{f.read_text()}")

    print(f"{len(df)} regions: " + ", ".join(
        f"{k}={v}" for k, v in df["category"].value_counts().items()))
    other = df[df["category"] == "other"]
    print("other: " + (", ".join(
        f"{k}={v}" for k, v in other["mobility_genes"].value_counts().items()) or "none"))
    print(f"candidates (evidence >= {args.min_evidence} + mobility marker): {len(ice)} (" + ", ".join(
        f"{k}={v}" for k, v in ice["category"].value_counts().items()) + ")")
    print(f"-> {out / 'ice_candidates.tsv'}, {out / 'all_regions.tsv'}, {fasta_out}/")


if __name__ == "__main__":
    main()
