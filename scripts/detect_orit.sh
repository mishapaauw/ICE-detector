#!/usr/bin/env bash
# Search oriT HMM profiles (DNA) against whole genomes with nhmmscan.
# Hits are reported loosely (E <= REPORT_E); summarize_orit.py applies the
# real E-value and profile-coverage cutoffs, so those can change without rerun.
# Usage: detect_orit.sh <genomes.fasta.list> <outdir> [cpu] [oriT.hmm]
set -euo pipefail

FASTA_LIST=$1
OUTDIR=$2
CPU=${3:-4}
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
HMM=${4:-$ROOT/db/orit/oriT.hmm}
REPORT_E=1

if [ ! -s "$HMM" ]; then
  echo "oriT HMMs not found: $HMM. Run: python scripts/download_orit_hmms.py -o db/orit" >&2
  exit 1
fi
[ -s "$HMM.h3m" ] || hmmpress -f "$HMM" > /dev/null

mkdir -p "$OUTDIR/nhmmscan"
{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "hmmer: $(nhmmscan -h | awk '/^# HMMER/{print $3; exit}')"
  echo "orit_hmm: $HMM ($(grep -c '^HMMER3' "$HMM") profiles)"
  echo "orit_source: $(head -1 "$(dirname "$HMM")/manifest.tsv" 2>/dev/null | sed 's/^# //')"
  echo "nhmmscan report E-value: $REPORT_E (filtering happens in summarize_orit.py)"
  echo "filters (summarize_orit.py): see orit_filters.txt"
  echo "cpu: $CPU"
  echo "seed: nhmmscan default (--seed 42)"
} > "$OUTDIR/versions.txt"

while IFS=$'\t' read -r genome path _; do
  [ -z "$genome" ] && continue
  echo "nhmmscan: $genome"
  nhmmscan --cpu "$CPU" -E "$REPORT_E" --seed 42 --noali \
    --tblout "$OUTDIR/nhmmscan/$genome.tbl" -o /dev/null "$HMM" "$path"
done < "$FASTA_LIST"
echo "oriT search finished: $OUTDIR"
