#!/usr/bin/env python
"""Run the ICE-detector pipeline on a folder of genome assemblies.

Run inside the ICE-detector conda env:
    python scripts/run_all.py data/test -o results/test --cpu 8 --min-orfs 10

Only the largest contig of each genome (the chromosome) is analysed; plasmids
and other contigs are dropped. The input FASTA files are only read: temporary
copies with just the largest contig are written to <outdir>/tmp_chromosomes/
and deleted when the run ends (also when it fails). input_contigs.tsv records
what was kept and dropped.

Steps (each writes its own subfolder with versions.txt):
    1. genomes.fasta.list  (make_fasta_list.py) + input_contigs.tsv
    2. ppanggolin/          (run_ppanggolin.sh)
    3. regions/             (extract_regions.py)
    4. conjugation/         (detect_conjugation.sh + summarize_conjugation.py)
    4b. integrase/          (detect_integrase.sh + summarize_integrase.py)
    5. orit/                (detect_orit.sh + summarize_orit.py)
    6. gc/                  (gc_content.py)
    6b. dr/                 (direct_repeats.py)
    6c. trna/               (predict_trnas.py before step 3; trna_proximity.py)
    7. final/               (build_table.py): ice_candidates.tsv, all_regions.tsv, ice_fasta/
Reference data must be in db/ first (see db/README.md).
"""
import argparse
import shutil
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
    ap.add_argument("--min-orfs", type=int, default=10,
                    help="keep regions with at least this many genes (default: 10)")
    ap.add_argument("--merge-gap", type=int, default=3,
                    help="merge neighbouring RGPs at most this many genes apart "
                         "(default: 3; 0 = no merging)")
    ap.add_argument("--linear", action="store_true", help="treat contigs as linear")
    ap.add_argument("--skip-ppanggolin", action="store_true",
                    help="reuse an existing <outdir>/ppanggolin")
    ap.add_argument("--skip-macsyfinder", action="store_true",
                    help="reuse existing <outdir>/conjugation/macsyfinder results")
    ap.add_argument("--models-dir", default=str(SCRIPTS.parent / "db" / "macsy_models"),
                    help="folder containing the CONJScan model package")
    ap.add_argument("--skip-nhmmscan", action="store_true",
                    help="reuse existing <outdir>/orit/nhmmscan results")
    ap.add_argument("--skip-hmmsearch", action="store_true",
                    help="reuse existing <outdir>/integrase/hmmsearch results")
    ap.add_argument("--integrase-hmm",
                    default=str(SCRIPTS.parent / "db" / "pfam" / "integrase.hmm"))
    ap.add_argument("--trna-dist", type=int, default=1000,
                    help="max bp between a tRNA/tmRNA and a region end (default: 1000)")
    ap.add_argument("--orit-hmm", default=str(SCRIPTS.parent / "db" / "orit" / "oriT.hmm"))
    ap.add_argument("--orit-min-cov", type=float, default=0.7,
                    help="min fraction of the oriT profile a hit must cover (default: 0.7)")
    ap.add_argument("--orit-max-evalue", type=float, default=0.01,
                    help="max E-value for oriT hits (default: 0.01)")
    ap.add_argument("--min-evidence", type=int, default=2,
                    help="min evidence count (of 7) for ice_candidates.tsv (default: 2)")
    ap.add_argument("--dr-outside", type=int, default=5000,
                    help="bp searched for direct repeats outside each region end (default: 5000)")
    ap.add_argument("--dr-inside", type=int, default=5000,
                    help="bp searched for direct repeats inside each region end (default: 5000)")
    ap.add_argument("--dr-min-len", type=int, default=15,
                    help="minimum direct repeat length (default: 15)")
    ap.add_argument("--dr-max-copies", type=int, default=2,
                    help="max genome-wide copies of a direct repeat (default: 2)")
    ap.add_argument("--gc-flank", type=int, default=35000,
                    help="bp up- and downstream used for flank GC (default: 35000)")
    ap.add_argument("--gc-z", type=float, default=2.0,
                    help="|GC z-score| counted as atypical GC (default: 2)")
    args = ap.parse_args()
    preflight(args)

    out = Path(args.outdir)
    logs = out / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    fasta_list = out / "genomes.fasta.list"
    chrom_dir = out / "tmp_chromosomes"
    try:
        run_steps(args, out, logs, fasta_list, chrom_dir)
    finally:
        shutil.rmtree(chrom_dir, ignore_errors=True)  # temporary chromosome copies

    run([sys.executable, SCRIPTS / "build_table.py", "--results", out, "--gc-z", args.gc_z,
         "--min-evidence", args.min_evidence],
        logs / "build_table.log")

    print((logs / "build_table.log").read_text().strip())
    print(f"[run_all] done: {out}/final/ice_candidates.tsv")


def preflight(args):
    """Fail early (before PPanGGOLiN) when inputs or reference data are missing."""
    problems = []
    if not Path(args.input_dir).is_dir():
        problems.append(f"input folder not found: {args.input_dir}")
    if not args.skip_macsyfinder and not (Path(args.models_dir) / "CONJScan").is_dir():
        problems.append(f"CONJScan models not found in {args.models_dir} "
                        "(msf_data install CONJScan==2.1.0 -t db/macsy_models)")
    if not args.skip_nhmmscan and not Path(args.orit_hmm).is_file():
        problems.append(f"oriT HMMs not found: {args.orit_hmm} "
                        "(python scripts/download_orit_hmms.py -o db/orit)")
    if not args.skip_hmmsearch and not Path(args.integrase_hmm).is_file():
        problems.append(f"integrase HMMs not found: {args.integrase_hmm} "
                        "(python scripts/download_pfam_hmms.py -o db/pfam)")
    if problems:
        sys.exit("[run_all] cannot start:\n  " + "\n  ".join(problems))


def run_steps(args, out, logs, fasta_list, chrom_dir):
    """Steps 1-6c; they read the genomes through fasta_list (copies in chrom_dir)."""
    py = sys.executable
    cmd = [py, SCRIPTS / "make_fasta_list.py", args.input_dir, "-o", fasta_list,
           "--chromosome-dir", chrom_dir, "--report", out / "input_contigs.tsv"]
    if args.linear:
        cmd.append("--linear")
    run(cmd, logs / "make_fasta_list.log")

    if not args.skip_ppanggolin:
        run(["bash", SCRIPTS / "run_ppanggolin.sh", fasta_list, out / "ppanggolin", args.cpu],
            logs / "ppanggolin.log")

    run([py, SCRIPTS / "predict_trnas.py", "--fasta-list", fasta_list, "-o", out / "trna"],
        logs / "predict_trnas.log")

    run([py, SCRIPTS / "extract_regions.py", "--ppanggolin", out / "ppanggolin",
         "--fasta-list", fasta_list, "-o", out / "regions", "--min-orfs", args.min_orfs,
         "--merge-gap", args.merge_gap, "--trnas", out / "trna" / "trnas.tsv"],
        logs / "extract_regions.log")

    conj = out / "conjugation"
    if not args.skip_macsyfinder:
        run(["bash", SCRIPTS / "detect_conjugation.sh", out / "ppanggolin", conj, args.cpu,
             args.models_dir, "linear" if args.linear else "circular"],
            logs / "detect_conjugation.log")
    run([py, SCRIPTS / "summarize_conjugation.py", "--regions", out / "regions" / "regions.tsv",
         "--region-genes", out / "regions" / "region_genes.tsv",
         "--tables", out / "ppanggolin" / "table", "--macsyfinder", conj / "macsyfinder",
         "-o", conj], logs / "summarize_conjugation.log")

    integ = out / "integrase"
    if not args.skip_hmmsearch:
        run(["bash", SCRIPTS / "detect_integrase.sh", conj / "proteins", integ, args.cpu,
             args.integrase_hmm], logs / "detect_integrase.log")
    run([py, SCRIPTS / "summarize_integrase.py", "--regions", out / "regions" / "regions.tsv",
         "--region-genes", out / "regions" / "region_genes.tsv",
         "--hmmsearch", integ / "hmmsearch", "-o", integ], logs / "summarize_integrase.log")

    orit = out / "orit"
    if not args.skip_nhmmscan:
        run(["bash", SCRIPTS / "detect_orit.sh", fasta_list, orit, args.cpu, args.orit_hmm],
            logs / "detect_orit.log")
    run([py, SCRIPTS / "summarize_orit.py", "--regions", out / "regions" / "regions.tsv",
         "--nhmmscan", orit / "nhmmscan", "--hmm", args.orit_hmm, "-o", orit,
         "--min-cov", args.orit_min_cov, "--max-evalue", args.orit_max_evalue],
        logs / "summarize_orit.log")

    run([py, SCRIPTS / "gc_content.py", "--regions", out / "regions" / "regions.tsv",
         "--fasta-list", fasta_list, "-o", out / "gc", "--flank", args.gc_flank],
        logs / "gc_content.log")

    run([py, SCRIPTS / "direct_repeats.py", "--regions", out / "regions" / "regions.tsv",
         "--fasta-list", fasta_list, "-o", out / "dr", "--outside", args.dr_outside,
         "--inside", args.dr_inside, "--min-len", args.dr_min_len,
         "--max-copies", args.dr_max_copies], logs / "direct_repeats.log")

    run([py, SCRIPTS / "trna_proximity.py", "--regions", out / "regions" / "regions.tsv",
         "--trnas", out / "trna" / "trnas.tsv", "--dr", out / "dr" / "dr.tsv", "-o", out / "trna",
         "--dist", args.trna_dist], logs / "trna_proximity.log")


if __name__ == "__main__":
    main()
