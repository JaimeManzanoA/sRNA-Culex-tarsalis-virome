#!/bin/bash
#SBATCH -J viruses
#SBATCH -o viruses-%j.output
#SBATCH -e viruses-%j.err
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jmm9413@psu.edu
#SBATCH --mail-type=ALL
#SBATCH  -t 48:00:00


module purge
module load anaconda
conda activate bioinfo

cd /storage/group/jlr54/default/Jaime_Manzano/Virome2

make align_Bunya_like ACC=Culex_Bunya_like_virus
make align_Bunyavirus_2 ACC=Culex_Bunyavirus_2
make align_Iflavi_like_1 ACC=Culex_Iflavi-like_virus_1
make align_Iflavi_like_4 ACC=Culex_Iflavi-like_virus_4
make align_Hubei ACC=Hubei_mosquito_virus_4
make align_Partitivirus ACC=Partitivirus_like_virus_4
make align_Wuhan ACC=Wuhan_Mosquito_Virus_6