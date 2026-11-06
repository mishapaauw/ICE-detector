#!/usr/bin/env bash
# Run PPanGGOLiN on a fasta list and write RGP/spot tables.
# Usage: run_ppanggolin.sh <genomes.fasta.list> <outdir> [cpu]
set -euo pipefail

FASTA_LIST=$1
OUTDIR=$2
CPU=${3:-4}
CONFIG="$(cd "$(dirname "$0")/.." && pwd)/config/ppanggolin.yaml"

mkdir -p "$OUTDIR"
cp "$CONFIG" "$OUTDIR/ppanggolin_config.yaml"
{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "ppanggolin: $(ppanggolin --version 2>&1 | awk '{print $2}')"
  echo "mmseqs2: $(mmseqs version 2>&1)"
  echo "seed (partition): $(awk '/^partition:/{f=1} f&&/seed:/{print $2; exit}' "$CONFIG")"
  echo "cpu: $CPU"
} > "$OUTDIR/versions.txt"

ppanggolin all --fasta "$FASTA_LIST" -o "$OUTDIR" --config "$CONFIG" \
  -c "$CPU" -f --log "$OUTDIR/ppanggolin_all.log"

# Explicitly (re)write the tables we depend on, so later steps don't rely on
# what `all` happens to write by default.
ppanggolin write_pangenome -p "$OUTDIR/pangenome.h5" -o "$OUTDIR" -f \
  --regions --spots --borders --stats --log "$OUTDIR/write_pangenome.log"
ppanggolin write_genomes -p "$OUTDIR/pangenome.h5" -o "$OUTDIR" -f \
  --table --log "$OUTDIR/write_genomes.log"

echo "PPanGGOLiN finished: $OUTDIR"
