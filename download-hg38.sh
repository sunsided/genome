#!/usr/bin/env bash

set -euo pipefail

BASE="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips"
OUT="hg38"

mkdir -p "$OUT"
cd "$OUT"

# Big file: segmented + resumable.
aria2c -x 8 -s 8 -c \
  "$BASE/hg38.fa.gz"

aria2c -x 1 -s 1 -c --allow-overwrite=true \
  "$BASE/hg38.chrom.sizes"

aria2c -x 1 -s 1 -c --allow-overwrite=true \
  "$BASE/hg38.chromAlias.txt"

aria2c -x 1 -s 1 -c --allow-overwrite=true \
  "$BASE/md5sum.txt"

md5sum -c md5sum.txt --ignore-missing

echo "Unzipping reference genome"
gunzip -k hg38.fa.gz

