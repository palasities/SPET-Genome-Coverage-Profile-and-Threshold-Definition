#!/usr/bin/env bash

IN="/mnt2/fscratch/users/forescent_001_upm/jpallares/ICIFOR/Sana_RadSeq/bam_files/subset_n15_analisis_SPET/intersection_minCov_gt0.bedgraph"
OUTDIR="/mnt2/fscratch/users/forescent_001_upm/jpallares/ICIFOR/Sana_RadSeq/bam_files/subset_n15_analisis_SPET/3_results_coverage_and_fusionGaps"
mkdir -p "$OUTDIR"

for T in 5 10 15 20; do         ## umbrales de cobertura
  for D in 10 20 30 40 50; do    ## gaps máximos para fusionar

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
