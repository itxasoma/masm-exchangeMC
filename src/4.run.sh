#!/bin/bash
#SBATCH -J exchMC_part4
#SBATCH -o exchMC_part4_%A_%a.out
#SBATCH -e exchMC_part4_%A_%a.err
#SBATCH --mail-user=imunozal74@alumnes.ub.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --array=1-5

set -euo pipefail

BASE_SEED=200000
SEED=$((BASE_SEED + SLURM_ARRAY_TASK_ID))
TAG=$(printf "part4_s%02d" "${SLURM_ARRAY_TASK_ID}")

mkdir -p ../inputs/part4_samples
mkdir -p ../results/part4

INFILE="../inputs/part4_samples/${TAG}.in"

awk -v newseed="${SEED}" '
  NR==9 { print newseed; next }
  { print }
' ../inputs/part4.in > "${INFILE}"

echo "Running sample ${SLURM_ARRAY_TASK_ID} with seed ${SEED}"
./exchange_mc "${INFILE}"

mv "timeseries_${TAG}.dat" "../results/part4/"
mv "swap_stats_${TAG}.dat" "../results/part4/"