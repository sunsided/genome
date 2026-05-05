#!/usr/bin/env bash

set -euo pipefail

BASE="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips"
OUT="hg38"

mkdir -p "$OUT"
cd "$OUT"

# Big file: segmented + resumable.
aria2c -x 8 -s 8 -c \
  "$BASE/latest/hg38.fa.gz"

# Small files: no segmentation, no resume weirdness.
aria2c -x 1 -s 1 --allow-overwrite=true \
  "$BASE/hg38.chrom.sizes" \
  "$BASE/hg38.chromAlias.txt" \
  "$BASE/md5sum.txt"

# Verify downloaded files that are present in md5sum.txt.
md5sum -c md5sum.txt --ignore-missing

