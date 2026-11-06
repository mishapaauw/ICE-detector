#!/usr/bin/env bash
# Detect conjugation systems and relaxases with MacSyFinder + CONJScan.
# Runs on whole genomes (ordered replicons) so systems are not cut at region
# borders; summarize_conjugation.py maps hits onto regions afterwards.
# Usage: detect_conjugation.sh <ppanggolin_dir> <outdir> [cpu] [models_dir] [topology]
set -euo pipefail

PPG=$1
OUTDIR=$2
CPU=${3:-4}
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODELS_DIR=${4:-$ROOT/db/macsy_models}
TOPOLOGY=${5:-circular}

if [ ! -d "$MODELS_DIR/CONJScan" ]; then
  echo "CONJScan models not found in $MODELS_DIR." >&2
  echo "Install with: msf_data install CONJScan==2.1.0 -t $MODELS_DIR" >&2
  exit 1
fi

mkdir -p "$OUTDIR/proteins" "$OUTDIR/macsyfinder"
{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "macsyfinder: $(macsyfinder --version 2>&1 | awk '/MacSyFinder/{print $2; exit}')"
  echo "hmmer: $(hmmsearch -h | awk '/^# HMMER/{print $3; exit}')"
  echo "CONJScan: $(awk '/^vers:/{print $2}' "$MODELS_DIR/CONJScan/metadata.yml")"
  echo "models: CONJScan/Chromosome all"
  echo "replicon_topology: $TOPOLOGY"
  echo "cpu: $CPU"
  echo "seed: not applicable (MacSyFinder is deterministic)"
} > "$OUTDIR/versions.txt"

ppanggolin fasta -p "$PPG/pangenome.h5" -o "$OUTDIR/proteins" --genes all -f \
  --log "$OUTDIR/ppanggolin_fasta.log"
python "$ROOT/scripts/make_proteins.py" --genes "$OUTDIR/proteins/all_genes.fna" \
  --tables "$PPG/table" -o "$OUTDIR/proteins"
rm "$OUTDIR/proteins/all_genes.fna"

for faa in "$OUTDIR"/proteins/*.faa; do
  genome=$(basename "$faa" .faa)
  echo "MacSyFinder: $genome"
  macsyfinder --db-type ordered_replicon --replicon-topology "$TOPOLOGY" \
    --sequence-db "$faa" --models-dir "$MODELS_DIR" \
    --models CONJScan/Chromosome all \
    -o "$OUTDIR/macsyfinder/$genome" -w "$CPU" --mute --force
done
echo "Conjugation detection finished: $OUTDIR"
