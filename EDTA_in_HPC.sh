#!/usr/bin/env bash
# Leave only one comment symbol on selected options
# Those with two commets will be ignored:
# The name to show in queue lists for this job:
#SBATCH -J EDTA_Qilex_localscratch.sh

# Number of desired cpus (can be in any node):
#SBATCH --ntasks=1

# Number of desired cpus (all in same node):
#SBATCH --cpus-per-task=50

# Amount of RAM needed for this job:
#SBATCH --mem=400gb

# The time the job will be running:
#SBATCH --time=150:00:00

# To use GPUs you have to request them:
##SBATCH --gres=gpu:1

# If you need nodes with special features leave only one # in the desired SBATCH constraint line. cal is selected by default:
# * to request any machine without GPU - DEFAULT
#SBATCH --constraint=cal
# * to request only the machines with 128 cores and 1800GB of usable RAM
##SBATCH --constraint=bigmem
# * to request only the machines with 128 cores and.4_new50GB of usable RAM (
##SBATCH --constraint=sr
# * to request only the machines with 52 cores and 187GB of usable RAM (
##SBATCH --constraint=sd

# Set output and error files
#SBATCH --error=job.EDTA_localscratch.err
#SBATCH --output=job.EDTA_localscratch.out

# Leave one comment in following line to make an array job. Then N jobs will be launched. In each one SLURM_ARRAY_TASK_ID will take one value from 1 to 100
##SBATCH --array=1-100

##modules

module load singularity/3.7.2

# ==============================
# RUTAS
# ==============================

OUTDIR="ICIFOR/AdaptiveSampling/1_EDTA_reps"
path_EDTA="HOME/modules/EDTA/EDTA"
fastita="subset_Quercus/0_parseo_NCBI_SANA/Qilex_chr_only_UPPER.fa"

# ==============================
# WORKDIR EN SCRATCH LOCAL
# ==============================

WORKDIR="${TMPDIR%/}/EDTA_${SLURM_JOB_ID:-manual}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "======================================"
echo "Working directory: $WORKDIR"
echo "======================================"

# ==============================
# COPIAR GENOMA A SCRATCH LOCAL
# ==============================

cp -v "$fastita" .

# ==============================
# EJECUTAR EDTA EN LOCAL
# ==============================

echo "Starting EDTA..."
singularity exec \
  --bind /mnt/home/users/forescent_001_upm \
  --bind "${TMPDIR%/}" \
  "$path_EDTA/EDTA.sif" \
  EDTA.pl \
    --genome "$(basename "$fastita")" \
    --step all \
    --species others \
    --sensitive 1 \
    --threads 50

echo "EDTA finished."

# ==============================
# INSPECCIÓN ANTES DE COPIAR
# ==============================

echo "======================================"
echo "Listing outputs before copy:"
echo "======================================"

ls -lh
echo
echo "Directory sizes:"
du -sh *
echo

echo "Key EDTA result files:"
ls -lh *.TElib.fa *.gff3 *.sum *.out 2>/dev/null || true
echo "======================================"

# ==============================
# COPIAR RESULTADOS A FSCRATCH
# ==============================

FINALDIR="$OUTDIR/EDTA_results_${SLURM_JOB_ID:-manual}"
mkdir -p "$FINALDIR"

echo "Copying results to $FINALDIR"

rsync -av \
  --exclude 'round-*' \
  --exclude '*.stk' \
  --exclude '*.nhr' --exclude '*.nin' --exclude '*.nsq' \
  . \
  "$FINALDIR/"

echo "======================================"
echo "Results copied successfully."
echo "Scratch directory preserved at:"
echo "$WORKDIR"
echo "======================================"
