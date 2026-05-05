#!/usr/bin/env bash

set -euo pipefail

mkdir -p hg38
cd hg38

aria2c -x 8 -s 8 -c \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.txt \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt

