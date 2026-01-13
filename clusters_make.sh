#!/bin/bash
#SBATCH -J cluster
#SBATCH -o cluster-%j.output
#SBATCH -e cluster-%j.err
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jmm9413@psu.edu
#SBATCH --mail-type=ALL
#SBATCH  -t 48:00:00


module purge
module load anaconda
conda activate virome

cd /storage/group/jlr54/default/Jaime_Manzano/Virome2/NA_ba

make cluster_all