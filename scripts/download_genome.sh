#!/usr/bin/env bash
set -euo pipefail

ACCESSION="$1"
OUTDIR="data/raw"

mkdir -p "$OUTDIR"

echo "Pobieranie $ACCESSION z NCBI..."
datasets download genome accession "$ACCESSION" \
  --include genome \
  --filename "${OUTDIR}/ncbi_dataset.zip"

echo "Rozpakowywanie..."
unzip -o "${OUTDIR}/ncbi_dataset.zip" -d "${OUTDIR}/tmp"

find "${OUTDIR}/tmp" -name "*.fna" -exec mv {} "${OUTDIR}/genome.fasta" \;

rm -rf "${OUTDIR}/tmp" "${OUTDIR}/ncbi_dataset.zip"

echo "Gotowe. Liczba sekwencji w assembly: $(grep -c '^>' ${OUTDIR}/genome.fasta)"