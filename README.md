JunctionJuror
=============
Given a junction.bed file and .gff genome file, this script will identify the
amount of alternative splicing arising from differential splice-site usage (i.e.
it ignores intron retention).

Arguments
---------
**Help** (`-h`, `--help`). Display help.

**Junction.bed list** (`-j`, `--junction JUNCTION.BED.LIST`). The input junction.bed
files should be referenced in a space-delimited list, with each line comprising
a path to a junction.bed file, followed by a condition. The current version of
JunctionJuror analyses single conditions, and hence this file should only
contain reference to one condition.

**Reference gene models** (`-g`, `--ref_gff GFF_FILE`). The gff parser is designed
for use with gff files from [EuPathDB](http://eupathdb.org). It only parses
features marked as `CDS` (and `tRNA` and `rRNA`), since UTR information is not
available for all genes.

**Output file** (`-o`, `--output OUTPUT_PATH`). JunctionJuror outputs the list of
genes that are alternatively spliced. If the `-e` flag has been used, it will
output the list of multi-exon genes instead. 

**Threshold for confirmation** (`-t`, `--threshold NUMBER`). Junctions must be
present in NUMBER replicates to be accepted as real. This defaults to two.

**Only list multi-exon genes** (`-e`, `--multi-exons`). Instead of reporting
alternatively spliced genes, JunctionJuror can also list multi-exon genes, i.e.
those containing at least one confirmed junction.
