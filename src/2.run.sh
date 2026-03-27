#!/bin/bash
#SBATCH --job-name=part2_pt
#SBATCH --output=part2_%j.out
#SBATCH --error=part2_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -e

echo "Job started on $(date)"

module purge
module load gcc

mkdir -p ../results ../figures

make clean
make

make run2
make binning

echo "Job finished on $(date)"
