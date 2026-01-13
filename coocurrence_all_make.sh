#!/bin/bash
#SBATCH -J coocurrence_all
#SBATCH -o coocurrence_all-%j.output
#SBATCH -e coocurrence_all-%j.err
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jmm9413@psu.edu
#SBATCH --mail-type=ALL
#SBATCH  -t 24:00:00


module purge
module load anaconda
conda activate bioinfo

cd /storage/group/jlr54/default/Jaime_Manzano/Virome2


make coocurrence_all ACC=reference_viruses
