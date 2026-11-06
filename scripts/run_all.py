#!/usr/bin/env python
"""Run the ICE-detector pipeline on a folder of genome assemblies.

Run inside the ICE-detector conda env:
    python scripts/run_all.py data/test -o results/test --cpu 8 --min-orfs 10

Steps (each writes its own subfolder with versions.txt):
    1. genomes.fasta.list  (make_fasta_list.py)
    2. ppanggolin/          (run_ppanggolin.sh)
    3. regions/             (extract_regions.py)
    4. conjugation/         (detect_conjugation.sh + summarize_conjugation.py)
Later milestones add oriT, GC and the final table.
"""
import argparse
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent


def run(cmd, log):
    print(f"[run_all] {' '.join(map(str, cmd))}  (log: {log})", flush=True)
    with open(log, "w") as fh:
        result = subprocess.run([str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        sys.exit(f"[run_all] failed (exit {result.returncode}); see {log}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_dir", help="folder with one FASTA per genome")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--min-orfs", type=int, default=10)
    ap.add_argument("--linear", action="store_true", help="treat contigs as linear")
    ap.add_argument("--skip-ppanggolin", action="store_true",
                    help="reuse an existing <outdir>/ppanggolin")
    ap.add_argument("--skip-macsyfinder", action="store_true",
                    help="reuse existing <outdir>/conjugation/macsyfinder results")
    ap.add_argument("--models-dir", default=str(SCRIPTS.parent / "db" / "macsy_models"),
                    help="folder containing the CONJScan model package")
    args = ap.parse_args()

    out = Path(args.outdir)
    logs = out / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    fasta_list = out / "genomes.fasta.list"
    py = sys.executable

    cmd = [py, SCRIPTS / "make_fasta_list.py", args.input_dir, "-o", fasta_list]
    if args.linear:
        cmd.append("--linear")
    run(cmd, logs / "make_fasta_list.log")

    if not args.skip_ppanggolin:
        run(["bash", SCRIPTS / "run_ppanggolin.sh", fasta_list, out / "ppanggolin", args.cpu],
            logs / "ppanggolin.log")

    run([py, SCRIPTS / "extract_regions.py", "--ppanggolin", out / "ppanggolin",
         "--fasta-list", fasta_list, "-o", out / "regions", "--min-orfs", args.min_orfs],
        logs / "extract_regions.log")

    conj = out / "conjugation"
    if not args.skip_macsyfinder:
        run(["bash", SCRIPTS / "detect_conjugation.sh", out / "ppanggolin", conj, args.cpu,
             args.models_dir, "linear" if args.linear else "circular"],
            logs / "detect_conjugation.log")
    run([py, SCRIPTS / "summarize_conjugation.py", "--regions", out / "regions" / "regions.tsv",
         "--tables", out / "ppanggolin" / "table", "--macsyfinder", conj / "macsyfinder",
         "-o", conj], logs / "summarize_conjugation.log")

    print(f"[run_all] done: {out}/regions/regions.tsv, {conj}/conjugation.tsv")


if __name__ == "__main__":
    main()
