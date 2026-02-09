
#########1#############

#!/usr/bin/env bash

bam_dir="/mnt2/fscratch/users/forescent_001_upm/jpallares/ICIFOR/Sana_RadSeq/bam_files/original"
Properly_paired="/mnt2/fscratch/users/forescent_001_upm/jpallares/ICIFOR/Sana_RadSeq/bam_files/original/properly_paired"

for bam in "$bam_dir"/*.bam; do
    base=$(basename "$bam" .bam)

    samtools view -b -F 0x904 -q 20 -f 0x2 "$bam" \
      > "${base}.q20_properlyPaired.bam"
done

#########2##########

module load bedtools


#!/usr/bin/env bash

PP_bam_dir="/mnt2/fscratch/users/forescent_001_upm/jpallares/ICIFOR/Sana_RadSeq/bam_files/original/properly_paired"
OUT='/mnt2/fscratch/users/forescent_001_upm/jpallares/ICIFOR/Sana_RadSeq/bam_files/subset_n15_analisis_SPET/1_bedGraphs'

for bam in "$PP_bam_dir"/*bam; do
    base=$(basename "$bam" .bam)
    
bedtools genomecov -ibam "$bam" -bga > "$OUT/${base}_cov.bedgraph"

done


########3##########

intersect directo entre bedGraphs no te alinea bien los cortes (cada archivo tiene su propia segmentación).
	  Find overlapping intervals in various ways.
unionbedg sí te parte el genoma en intervalos compatibles entre todos y te permite aplicar un AND exacto.
          Combines coverage intervals from multiple BEDGRAPH files.


#!/usr/bin/env bash
set -euo pipefail

##intersection_minCov_gt0.bedgraph contiene:
##solo regiones presentes en todas las muestras
##4ª columna = mínima cobertura entre todas


BDG_DIR="/mnt2/fscratch/users/forescent_001_upm/jpallares/ICIFOR/Sana_RadSeq/bam_files/subset_n15_analisis_SPET/1_bedGraphs"
OUTDIR="/mnt2/fscratch/users/forescent_001_upm/jpallares/ICIFOR/Sana_RadSeq/bam_files/subset_n15_analisis_SPET/2_intersection"

mkdir -p "$OUTDIR"

# 1) lista de bedGraphs
ls "$BDG_DIR"/*_cov.bedgraph > "$OUTDIR/bedgraphs.list"

# 2) ordena in-place (imprescindible)
while read -r f; do
  LC_ALL=C sort -k1,1 -k2,2n "$f" -o "$f"
done < "$OUTDIR/bedgraphs.list"

# 3) intersección: unionbedg + min across samples
bedtools unionbedg -i $(cat "$OUTDIR/bedgraphs.list") \
| awk 'BEGIN{OFS="\t"}
  {
    # columnas: chr start end cov1 cov2 ... covN
    min=$4; ok=1
    for(i=4;i<=NF;i++){
      if($i<=0){ok=0; break}  ##Si alguna muestra tiene 0, ese intervalo se descarta
      if($i<min) min=$i
    }
    if(ok) print $1,$2,$3,min
  }' \
| awk 'BEGIN{OFS="\t"}
  NR==1{c=$1;s=$2;e=$3;v=$4;next}
  ($1==c && $2==e && $4==v){e=$3;next}
  {print c,s,e,v; c=$1;s=$2;e=$3;v=$4}
  END{print c,s,e,v}' \                ##Regiones donde TODAS las muestras tienen cobertura > 0
> "$OUTDIR/intersection_minCov_gt0.bedgraph"  


########4##########

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
