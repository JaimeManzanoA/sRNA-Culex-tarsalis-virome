#!/bin/bash
#SBATCH -J bacteria_download
#SBATCH -o bacteria_download.output
#SBATCH -e bacteria_download.err
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jmm9413@psu.edu
#SBATCH --mail-type=ALL
#SBATCH -t 48:00:00
#SBATCH --partition=open 


module purge
module load anaconda
conda activate virome

cd ~/scratch/Virome2

make bacteria_download