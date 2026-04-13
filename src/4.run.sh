#!/bin/bash
#SBATCH -J exchMC_part4
#SBATCH -o exchMC_part4_%a.out
#SBATCH -e exchMC_part4_%a.err
#SBATCH --mail-user=imunozal74@alumnes.ub.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04-00:00:00
#SBATCH --array=1-5

set -e

BASE_SEED=200000
SEED=$((BASE_SEED + SLURM_ARRAY_TASK_ID))
TAG=$(printf "part4_s%04d" "${SLURM_ARRAY_TASK_ID}")
INFILE="../inputs/part4_samples/${TAG}.in"

mkdir -p ../inputs/part4_samples ../results/part4 ../logs/part4

awk -v newseed="${SEED}" 'NR==9{$0=newseed} {print}' \
    ../inputs/part4.in > "${INFILE}"

echo "Sample ${SLURM_ARRAY_TASK_ID}/5  seed=${SEED}"
./exchange_mc "${INFILE}" > "../logs/part4/${TAG}.log" 2>&1

mv "timeseries_${TAG}.dat" "../results/part4/"
mv "swap_stats_${TAG}.dat" "../results/part4/"