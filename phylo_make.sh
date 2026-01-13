#!/bin/bash
#SBATCH -J phylo
#SBATCH -o phylo-%j.output
#SBATCH -e phylo-%j.err
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jmm9413@psu.edu
#SBATCH --mail-type=ALL
#SBATCH  -t 48:00:00


module purge
module load anaconda
conda activate virome

cd /storage/group/jlr54/default/Jaime_Manzano/Virome2

make tree VIR=Bunya_like_MA MOD=K80+I
make tree VIR=Bunyavirus_2_MA MOD=K80+I
make tree VIR=Hubei_MA MOD=TIM2ef+G4
make tree VIR=Iflavi_4_MA MOD=TIM2+I
make tree VIR=Marma_MA MOD=TIM2ef+G4
make tree VIR=Narnaviridae_MA MOD=TPM3uf+G4
make tree VIR=Partitivirus_MA MOD=TIM2+I
make tree VIR=Tombusviridae_MA MOD=K80+G4
make tree VIR=Wuhan_MA MOD=TPM1+G4
