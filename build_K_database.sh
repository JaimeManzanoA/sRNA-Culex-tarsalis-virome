#!/bin/bash
#SBATCH -J K_database
#SBATCH -o K_database.output
#SBATCH -e K_database.err
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jmm9413@psu.edu
#SBATCH --mail-type=ALL
#SBATCH -t 48:00:00
#SBATCH --partition=open 


module purge
module load anaconda
conda activate virome

cd ~/scratch/Virome2

kraken2-build --build --threads 32 --db refs/bact_database