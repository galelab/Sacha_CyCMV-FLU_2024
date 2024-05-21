#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 5800-12
#SBATCH --partition HoldingPen
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
PATH_FASTQ="./trimgalore_results/"
samples=("$PATH_FASTQ"*.fq.gz)

mkdir trimgalore_fastqc_results
for sample in  ${samples[*]}
do
	srun -c 8 fastqc "$sample" --noextract -t 8 -o ./trimgalore_fastqc_results/
	wait
done
