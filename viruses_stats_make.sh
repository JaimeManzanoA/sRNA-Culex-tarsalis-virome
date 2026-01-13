#!/bin/bash
#SBATCH -J viruses_stats
#SBATCH -o viruses_stats-%j.output
#SBATCH -e viruses_stats-%j.err
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jmm9413@psu.edu
#SBATCH --mail-type=ALL
#SBATCH  -t 48:00:00


module purge
module load anaconda
conda activate /storage/home/jmm9413/micromamba/envs/stats

cd /storage/group/jlr54/default/Jaime_Manzano/Virome2/stats

make -f ../Makefile stats_virus SRR=Culex_Iflavi-like_virus_1
make -f ../Makefile stats_virus SRR=Culex_Iflavi-like_virus_4
make -f ../Makefile stats_virus SRR=Hubei_mosquito_virus_4
make -f ../Makefile stats_virus SRR=Partitivirus_like_virus_4
make -f ../Makefile stats_virus SRR=Wuhan_Mosquito_Virus_6

