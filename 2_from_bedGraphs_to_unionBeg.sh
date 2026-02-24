#!/usr/bin/env bash
# Leave only one comment symbol on selected options
# Those with two commets will be ignored:
# The name to show in queue lists for this job:
#SBATCH -J bedgraphs_to_unionbedg_GCoverage.sh

# Number of desired cpus (can be in any node):
#SBATCH --ntasks=1

# Number of desired cpus (all in same node):
#SBATCH --cpus-per-task=10

# Amount of RAM needed for this job:
#SBATCH --mem=100gb

# The time the job will be running:
#SBATCH --time=168:00:00

# To use GPUs you have to request them:
##SBATCH --gres=gpu:1

# If you need nodes with special features leave only one # in the desired SBATCH constraint line. cal is selected by default:
# * to request any machine without GPU - DEFAULT
#SBATCH --constraint=cal
# * to request only the machines with 128 cores and 1800GB of usable RAM
##SBATCH --constraint=bigmem
# * to request only the machines with 128 cores and 450GB of usable RAM (
##SBATCH --constraint=sr
# * to request only the machines with 52 cores and 187GB of usable RAM (
##SBATCH --constraint=sd

# Set output and error files
#SBATCH --error=job.bedgraphs_to_unionbedgGCoverage.err
#SBATCH --output=job.bedgraphs_to_unionbedgGCoverage.out

# Leave one comment in following line to make an array job. Then N jobs will be launched. In each one SLURM_ARRAY_TASK_ID will take one value from 1 to 100
##SBATCH --array=1-100

module load bedtools

##############
##intersection_minCov_gt0.bedgraph contiene:
##solo regiones presentes en todas las muestras
##4ª columna = mínima cobertura entre todas

BED='/mnt/home/users/forescent_001_upm/jpallares/fscratch/ICIFOR/Sana_RadSeq/bam_files/subset_n15_analisis_SPET/2A_3A_without_lowcoverage_samples'
OUTDIR="/mnt/home/users/forescent_001_upm/jpallares/fscratch/ICIFOR/Sana_RadSeq/bam_files/subset_n15_analisis_SPET/2A_3A_without_lowcoverage_samples/2A_intersection"


mkdir -p "$OUTDIR"

# 2) ordena in-place (imprescindible)
while read -r f; do
  LC_ALL=C sort -k1,1 -k2,2n "$f" -o "$f"
done < "$BED/bed_list_goodCoverage.list"

# 3) intersección: unionbedg + min across simples
 ##guarda el archivo y sigue
bedtools unionbedg -i $(cat "$BED/bed_list_goodCoverage.list") \
| tee "$OUTDIR/unionbedg_raw_GCoverage.tsv" \
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
  END{print c,s,e,v}' \
> "$OUTDIR/intersection_minCov_gt0_GCoverage.bedgraph"
