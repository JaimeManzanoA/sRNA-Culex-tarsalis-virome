#!/bin/bash
#SBATCH -J align
#SBATCH -o align-%j.output
#SBATCH -e align-%j.err
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jmm9413@psu.edu
#SBATCH --mail-type=ALL
#SBATCH  -t 48:00:00


module purge
module load anaconda
conda activate bioinfo

cd /storage/work/jmm9413/Virome

make align SRR=KWNR_S148_L002_R1_001
make align SRR=P10_S8_L002_R1_001
make align SRR=P11_S13_L002_R1_001
make align SRR=P12_S9_L002_R1_001
make align SRR=P13_S14_L002_R1_001
make align SRR=P14_S10_L002_R1_001
make align SRR=P15_S15_L002_R1_001
make align SRR=P16_S16_L002_R1_001
make align SRR=P17_S17_L002_R1_001
make align SRR=P1_S45_L002_R1_001
make align SRR=P2_S46_L002_R1_001
make align SRR=P3_S47_L002_R1_001
make align SRR=P4_S4_L002_R1_001
make align SRR=P5_S11_L002_R1_001
make align SRR=P6_S5_L002_R1_001
make align SRR=P7_S6_L002_R1_001
make align SRR=P8_S7_L002_R1_001
make align SRR=P9_S12_L002_R1_001