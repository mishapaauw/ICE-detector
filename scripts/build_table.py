#!/usr/bin/env python
"""Merge all evidence into the final ICE table (milestone 5).

Per region (RGP with >= --min-orfs genes), evidence is:
  t4ss       conjugation system (CONJScan T4SS model, >= 50% of its genes in region)
  relaxase   MOB relaxase gene in the region
  orit       oriT hit passing the E-value and profile-coverage filters
  gc         atypical GC: |gc_zscore| >= --gc-z (default 2)
evidence_count = number of these four that are True.

category:
  ICE           t4ss (Cury et al. 2017: spots encoding a conjugative system)
  decayed_ICE   no t4ss, but a decayed conjugation system (dCONJ model)
  mobilizable   no conjugation system, but a relaxase (possible IME)
  none          none of the above
confidence (ICE only):
  high    t4ss + relaxase + orit
  medium  t4ss + relaxase or orit
  low     t4ss only

Outputs (in --outdir):
  all_regions.tsv     every region with all evidence columns
  ice_candidates.tsv  regions with evidence_count >= --min-evidence (default 2),
                      any category (ICEs and e.g. IMEs), sorted by category,
                      confidence and evidence
  ice_fasta/          FASTA of each candidate region (copied from regions/fasta)
  versions.txt        settings + versions.txt of every upstream step
"""
import argparse
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

CONF_ORDER = {"high": 0, "medium": 1, "low": 2, "": 3}


def classify(row):
    if row["t4ss"]:
        n = row["relaxase"] + row["orit"]
        return "ICE", {2: "high", 1: "medium", 0: "low"}[n]
    if row["dconj_system"]:
        return "decayed_ICE", ""
    if row["relaxase"]:
        return "mobilizable", ""
    return "none", ""


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

    df = (regions.merge(conj, on="region_id", how="left", validate="1:1")
                 .merge(orit, on="region_id", how="left", validate="1:1")
                 .merge(gc, on="region_id", how="left", validate="1:1"))
    if len(df) != len(regions):
        sys.exit("Region tables do not match; re-run upstream steps.")

    df = df.rename(columns={"conj_system": "t4ss"})
    for col in ("t4ss", "dconj_system", "relaxase", "orit", "t4cp", "virb4"):
        df[col] = df[col].fillna(False).astype(bool)
    df["gc_atypical"] = df["gc_zscore"].abs() >= args.gc_z
    df["evidence_count"] = df[["t4ss", "relaxase", "orit", "gc_atypical"]].sum(axis=1)
    df[["category", "confidence"]] = df.apply(lambda r: pd.Series(classify(r)), axis=1)
    df = df.join(spot_genome_counts(res / "ppanggolin" / "spots.tsv",
                                    res / "ppanggolin" / "regions_of_genomic_plasticity.tsv"),
                 on="spot_id")
    df["spot_n_genomes"] = df["spot_n_genomes"].fillna(1).astype(int)

    cols = ["region_id", "genome", "contig", "start", "end", "length", "n_orfs",
            "spot_id", "spot_n_genomes", "category", "confidence", "evidence_count",
            "t4ss", "relaxase", "orit", "gc_atypical",
            "conj_systems", "dconj_systems", "relaxase_types", "t4cp", "virb4",
            "orit_profiles", "orit_best_evalue", "region_gc", "genome_gc", "gc_diff",
            "gc_zscore", "flank_up_gc", "flank_down_gc", "flank_gc", "gc_diff_flank",
            "coordinates", "contig_border"]
    df = df[cols].sort_values(
        ["category", "confidence", "evidence_count", "genome", "start"],
        key=lambda s: s.map(CONF_ORDER) if s.name == "confidence"
        else (s.map({"ICE": 0, "decayed_ICE": 1, "mobilizable": 2, "none": 3})
              if s.name == "category" else s),
        ascending=[True, True, False, True, True])
    df.to_csv(out / "all_regions.tsv", sep="\t", index=False)

    ice = df[df["evidence_count"] >= args.min_evidence]
    ice.to_csv(out / "ice_candidates.tsv", sep="\t", index=False)
    fasta_out = out / "ice_fasta"
    if fasta_out.exists():
        shutil.rmtree(fasta_out)
    fasta_out.mkdir()
    for rid in ice["region_id"]:
        shutil.copy(res / "regions" / "fasta" / f"{rid}.fasta", fasta_out)

    with open(out / "versions.txt", "w") as fh:
        fh.write(f"date: {datetime.now(timezone.utc):%Y-%m-%dT%H:%M:%SZ}\n"
                 f"python: {sys.version.split()[0]}\npandas: {pd.__version__}\n"
                 f"gc_z: {args.gc_z}\nmin_evidence: {args.min_evidence}\n")
        for step in ("ppanggolin", "regions", "conjugation", "orit", "gc"):
            v = res / step / "versions.txt"
            if v.exists():
                fh.write(f"\n## {step}\n{v.read_text()}")
        f = res / "orit" / "orit_filters.txt"
        if f.exists():
            fh.write(f"\n## orit filters\n{f.read_text()}")

    print(f"{len(df)} regions: " + ", ".join(
        f"{k}={v}" for k, v in df["category"].value_counts().items()))
    ices = df[df["category"] == "ICE"]
    print("ICE confidence: " + (", ".join(
        f"{k}={v}" for k, v in ices["confidence"].value_counts().items()) or "none"))
    print(f"candidates (evidence >= {args.min_evidence}): {len(ice)} (" + ", ".join(
        f"{k}={v}" for k, v in ice["category"].value_counts().items()) + ")")
    print(f"-> {out / 'ice_candidates.tsv'}, {out / 'all_regions.tsv'}, {fasta_out}/")


if __name__ == "__main__":
    main()
