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

3. Breadth of Coverage at Different Thresholds and Number of Regions

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

-  1X Depth:   1.09%
-  5X Depth:   0.35%
- 10X Depth:   0.16%

Number of genomic regions captured at different coverage thresholds

```bash
awk '$4>10 {print $1,$2,$3}' OFS="\t" 118-If-Qi123-P2QiIF-F03.q20_properlyPaired_cov.bedgraph \
| bedtools merge \
| wc -l > n.regions_10X_bedtools.txt  ##for 10X
```
-  >1X Depth: 33868
-  >5X Depth: 21602
- >10X Depth: 13136

For each depth threshold (>1×, >5×, >10×), genomic regions were defined as continuous intervals where coverage met or exceeded the threshold. Adjacent bases meeting the criterion were merged into single regions (See the -d parameter for merging nearby continuous regions)

```bash
for T in 1 5 10; do   awk -v T=$T '$4>T {print $1,$2,$3}' OFS="\t" ../118-If-Qi123-P2QiIF-F03.q20_properlyPaired_cov.bedgraph   | bedtools merge -i -   > regions_ge${T}x.bed;    n=$(wc -l < regions_ge${T}x.bed);   bp=$(awk '{sum+=($3-$2)} END{print sum}' regions_ge${T}x.bed);    mean=$(awk '{sum+=($3-$2); n++} END{print sum/n}' regions_ge${T}x.bed)
printf ">%dx\tregions: %d\tbases: %d\tmean_len: %.1f\n" "$T" "$n" "$bp" "$mean"; done
```
- >1x     regions: 33868  bases: 6107706  mean_len: 180.3
- >5x     regions: 21602  bases: 2357203  mean_len: 109.1
- >10x    regions: 13136  bases: 1202682  mean_len: 91.6


Last command, define T (coverage threshold) and D (-d parameter maximum gap between regions):

```
IN="118-If-Qi123-P2QiIF-F03.q20_properlyPaired_cov.bedgraph"
OUTDIR="bed_files/with_gaps_fusion/beds_cov"
mkdir -p "$OUTDIR"

for T in 1 5 10 15 20; do         ## umbrales de cobertura
  for D in 10 20 50; do    ## gaps máximos para fusionar

    out="$OUTDIR/$(basename "${IN%.bedgraph}").gt${T}.merge_d${D}.bed"

    awk -v T="$T" '$4 > T {print $1, $2, $3}' OFS="\t" "$IN" \
      | LC_ALL=C sort -k1,1 -k2,2n \
      | bedtools merge -i - -d "$D" \
      > "$out"

    n=$(wc -l < "$out")

    # bases totales
    bp=$(awk '{sum+=($3-$2)} END{print sum+0}' "$out")

    if [ "$n" -gt 0 ]; then
      read mean min max < <(
        awk '
          {len=$3-$2; sum+=len}
          NR==1 {min=len; max=len}
          {if (len<min) min=len; if (len>max) max=len}
          END {print sum/NR, min, max}
        ' "$out"
      )
    else
      mean=0; min=0; max=0
    fi

    printf ">%dx (merge_d=%d)\tregions: %d\tbases: %d\tmean_len: %.1f\tmin_len: %d\tmax_len: %d\tbed: %s\n" \
      "$T" "$D" "$n" "$bp" "$mean" "$min" "$max" "$out"

  done
done
```

Something like:

- >1x (merge_d=10)        regions: 29678  bases: 6122727  mean_len: 206.3 min_len: 1      max_len: 1377   bed: bed_files/with_gaps_fusion/beds_cov/118-If-Qi123-P2QiIF-F03.q20_properlyPaired_cov.gt1.merge_d10.bed
- >1x (merge_d=20)        regions: 28215  bases: 6144966  mean_len: 217.8 min_len: 1      max_len: 1377   bed: bed_files/with_gaps_fusion/beds_cov/118-If-Qi123-P2QiIF-F03.q20_properlyPaired_cov.gt1.merge_d20.bed

......
......
......

- >20x (merge_d=50)       regions: 5250   bases: 508886   mean_len: 96.9  min_len: 1      max_len: 1081   bed: bed_files/with_gaps_fusion/beds_cov/118-If-Qi123-P2QiIF-F03.q20_properlyPaired_cov.gt20.merge_d50.bed



