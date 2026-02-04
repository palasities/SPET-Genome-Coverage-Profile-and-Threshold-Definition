# SPET Genome Coverage Profile and Threshold Definition

This section describes the workflow used to define the genome coverage profile and calculate the breadth of coverage at different depth thresholds for SPET (Single Primer Enrichment Technology) data.

1. Quality Filtering for "Clean" Read. We filter the BAM files to keep only primary, mapped reads with a mapping quality (MAPQ) >= 20.

Standard Filtering

```bash
samtools view -b -F 0x904 -q 20 118-If-Qi123-P2QiIF-F03.q10.bam > 118-If-Qi123-P2QiIF-F03.q20_NO_unmapped_sec_suppl.bam
# 0x904 = unmapped(0x4) + secondary(0x100) + supplementary(0x800)
```
Strict Filtering (Properly Paired)

```bash
samtools view -b -F 0x904 -q 20 -f 0x2 118-If-Qi123-P2QiIF-F03.q10.bam > 118-If-Qi123-P2QiIF-F03.q20_properlyPaired.bam
# 0x2 = properly paired
```

Read Count Statistics

```bash
samtools view -c
```
- 118-If-Qi123-P2QiIF-F03.q10.bam   					        419879
- 118-If-Qi123-P2QiIF-F03.q20_NO_unmapped_sec_suppl.bam		    417611
- 118-If-Qi123-P2QiIF-F03.q20_properlyPaired.bam			    415601

2. Genome Coverage Calculation

Using bedtools genomecov to generate a BedGraph file.

```bash
# Basic command
bedtools genomecov -ibam 118-If-Qi123-P2QiIF-F03.q20_properlyPaired.bam -bga > 118-If-Qi123-P2QiIF-F03.q20_properlyPaired_cov.bedgraph

# More robust approach using a genome size file
cut -f1,2 ref.fasta.fai > genome.txt

bedtools genomecov -ibam ../118-If-Qi123-P2QiIF-F03.q20_properlyPaired.bam -bga -g genome.txt > 118-If-Qi123-P2QiIF-F03.q20_properlyPaired_cov.bedgraph

```

Note: When providing a BAM file as input, bedtools extracts chromosome sizes directly from the BAM header. Using the -g flag alongside a BAM might trigger a warning as the tool prioritizes the header info.

3. Breadth of Coverage at Different Thresholds

We calculated the total number of base pairs (bp) covered at >=1X, >=5X, and >=10X depth. A 10X threshold is generally considered adequate for reliable variant calling.

```bash
# Coverage >= 10X
awk '$4>=10{cov+=($3-$2)} END{print cov}' 118-If-Qi123-P2QiIF-F03.q20_properlyPaired_cov.bedgraph > bp_10X_SPET.txt

# Coverage >= 5X
awk '$4>=5{cov+=($3-$2)} END{print cov}' 118-If-Qi123-P2QiIF-F03.q20_properlyPaired_cov.bedgraph > bp_5X_SPET.txt

# Coverage >= 1X
awk '$4>=1{cov+=($3-$2)} END{print cov}' 118-If-Qi123-P2QiIF-F03.q20_properlyPaired_cov.bedgraph > bp_1X_SPET.txt
```

Based on the total genome size, the percentages of captured regions are:

- 10X Depth:  0.16%
- 5X Depth:   0.35%
- 1X Depth:   1.09%