#!/usr/bin/env bash
# Search integrase / recombinase Pfam profiles in all proteins of each genome.
# Uses the per-genome protein FASTA made by detect_conjugation.sh (PPanGGOLiN gene IDs).
# Cutoff: Pfam gathering thresholds (--cut_ga), as curated per family by Pfam.
# Usage: detect_integrase.sh <proteins_dir> <outdir> [cpu] [integrase.hmm]
set -euo pipefail

PROT=$1
OUTDIR=$2
CPU=${3:-4}
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
HMM=${4:-$ROOT/db/pfam/integrase.hmm}

if [ ! -s "$HMM" ]; then
  echo "Pfam integrase HMMs not found: $HMM. Run: python scripts/download_pfam_hmms.py -o db/pfam" >&2
  exit 1
fi
ls "$PROT"/*.faa > /dev/null 2>&1 || { echo "No protein FASTA (*.faa) in $PROT" >&2; exit 1; }

mkdir -p "$OUTDIR/hmmsearch"
{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "hmmer: $(hmmsearch -h | awk '/^# HMMER/{print $3; exit}')"
  echo "integrase_hmm: $HMM ($(grep -c '^HMMER3' "$HMM") profiles)"
  echo "pfam: $(head -1 "$(dirname "$HMM")/manifest.tsv" 2>/dev/null | sed 's/^# //')"
  echo "cutoff: --cut_ga (Pfam gathering thresholds)"
  echo "cpu: $CPU"
  echo "seed: hmmsearch default (--seed 42)"
} > "$OUTDIR/versions.txt"

for faa in "$PROT"/*.faa; do
  genome=$(basename "$faa" .faa)
  echo "hmmsearch: $genome"
  hmmsearch --cut_ga --cpu "$CPU" --seed 42 --noali \
    --domtblout "$OUTDIR/hmmsearch/$genome.domtbl" -o /dev/null "$HMM" "$faa"
done
echo "Integrase search finished: $OUTDIR"
