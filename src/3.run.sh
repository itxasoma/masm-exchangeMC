#!/bin/bash
#SBATCH -J exchMC_part3
#SBATCH -o exchMC_part3.out
#SBATCH -e exchMC_part3.err
#SBATCH --mail-user=imunozal74@alumnes.ub.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04-00:00:00

set -e

NSAMPLES=1000
BASE_SEED=123456

TEMPLATE="../inputs/part3.in"

make clean && make

mkdir -p ../inputs/part3_samples
mkdir -p ../results/part3
mkdir -p ../logs/part3

for s in $(seq 1 ${NSAMPLES}); do
  SEED=$((BASE_SEED + s))
  TAG=$(printf "part3_s%04d" "${s}")
  INFILE="../inputs/part3_samples/${TAG}.in"
  LOGFILE="../logs/part3/${TAG}.log"

  awk -v newseed="${SEED}" 'NR==9{$0=newseed} {print}' "${TEMPLATE}" > "${INFILE}"

  echo "Launching sample ${s}/${NSAMPLES}  seed=${SEED}  ->  ${INFILE}"
  ./exchange_mc "${INFILE}" > "${LOGFILE}" 2>&1

  mv "histogram_${TAG}.dat" "../results/part3/"
  mv "swap_stats_${TAG}.dat" "../results/part3/"
done

echo "All runs finished."
