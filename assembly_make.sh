#!/bin/bash
#SBATCH -J assembly
#SBATCH -o assembly-%j.output
#SBATCH -e assembly-%j.err
#SBATCH -N 1
#SBATCH -n 45
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jmm9413@psu.edu
#SBATCH --mail-type=ALL
#SBATCH -t 48:00:00
#SBATCH --partition=open 


module purge
module load anaconda
conda activate Trinity

cd ~/storage/work/jmm9413/Virome/trim2

make assembly