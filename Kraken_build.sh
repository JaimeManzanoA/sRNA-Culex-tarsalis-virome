#!/bin/bash
#SBATCH -J Kraken2_build
#SBATCH -o Kraken2_build.output
#SBATCH -e Kraken2_build.err
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jmm9413@psu.edu
#SBATCH --mail-type=ALL
#SBATCH -t 48:00:00
#SBATCH --partition=open 


module purge
module load anaconda
conda activate virome

cd ~/scratch/Virome2

make Kraken2_build