#!/bin/bash
#SBATCH -J exchMC_part6
#SBATCH -o exchMC_part6_%A_%a.out
#SBATCH -e exchMC_part6_%A_%a.err
#SBATCH --mail-user=imunozal74@alumnes.ub.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --array=1-5

set -euo pipefail

TAG4=$(printf "part4_s%02d" "${SLURM_ARRAY_TASK_ID}")
TAG6=$(printf "part6_s%02d" "${SLURM_ARRAY_TASK_ID}")

SRC_IN="../inputs/part4_samples/${TAG4}.in"
if [ ! -f "${SRC_IN}" ]; then
  SRC_IN="../inputs/part4.in"
fi

EXE=./exchange_mc_part6

mkdir -p ../inputs/part6_samples
mkdir -p ../results/part6

INFILE="../inputs/part6_samples/${TAG6}.in"

# Claude helper to create the input file for part 6 from the sample input files of part 4
awk '
  NR==1  { d=$0 }
  NR==2  { L=$0 }
  NR==8  { nMCS=$0 }
  NR==9  { seed=$0 }
  NR==10 { bond=$0 }
  END {
    print d
    print L
    print "0.2"
    print "0.2"
    print "1"
    print "1"
    print "100"
    print nMCS
    print seed
    print bond
    print "linear"
  }
' "${SRC_IN}" > "${INFILE}"

echo "Running Metropolis sample ${SLURM_ARRAY_TASK_ID}"
echo "Input file: ${INFILE}"

"${EXE}" "${INFILE}"

mv "timeseries_${TAG6}.dat" "../results/part6/"
