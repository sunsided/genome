# 🧬 Human Genome Playground

<div align="center">
  <img src="https://raw.githubusercontent.com/sunsided/genome/refs/heads/main/.readme/viz.png" alt="Cypher crate hero picture" />
</div>

> Output of the `genome-viz` TUI, run with `task viz`. Mitochondrial chromosome in the hg38/GRCh38 FASTA file.


## `genome-inspect` output

```text
file: hg38/hg38.fa.gz

summary
-------
contigs:              455
total bases:          3209286105
canonical chromosomes: 25

base composition
----------------
A:     898285419 (27.99%)
C:     623727342 (19.44%)
G:     626335137 (19.52%)
T:     900967885 (28.07%)
N:     159970322 (4.98%)
other: 0 (0.000000%)

inferred metadata
-----------------
assembly:             unknown, likely hg38/GRCh38 if downloaded from UCSC hg38
provider:             unknown from FASTA alone, likely UCSC if file is hg38.fa.gz
contig style:         UCSC-style chr-prefixed names
mitochondrial contig: yes
contains alt/patch-like contigs: yes
contains random/unplaced contigs: yes
soft masked:          no / not detected
lowercase bases:      0

note
----
A FASTA file alone does not reliably prove the assembly/provider.
Keep the download URL, md5sum.txt, and chromAlias.txt next to the FASTA.

chromAlias
----------
file: hg38/hg38.chromAlias.txt
alias rows: 455
```

## Reference genome data

This project may use the UCSC hg38 / GRCh38 human reference genome downloaded
from the UCSC Genome Browser download server:

https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/

The UCSC Genome Browser group states that no license is required for the raw
data files and database tables used by the Genome Browser, and that these data
are freely available for public and commercial use, except for specific files
such as liftOver chain files.

Please cite the UCSC Genome Browser when using these data in published work.
