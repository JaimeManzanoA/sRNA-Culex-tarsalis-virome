###############Run this in cluster

salloc -N 1 -n 2 --mem-per-cpu=5gb -t 4:00:00


#for a bit more intense work
salloc -N 1 -n 4 --mem-per-cpu=8gb -t 6:00:00

ssh -Y jmm9413@submit.hpc.psu.edu

#####################################


Uploaded the raw files to the default folder in jlr54 group roar collab allocation. All the data is inside Sultan's folder.

Copied some files from his folder to my scratch folder

cp -r /storage/group/jlr54/default/Sultan_Asad_2022/Wild_tarsalis_sRNA_sequencing_data/Novagene/1735-1771 /scratch/jmm9413/Virome


mv -r /scratch/jmm9413/reads /scratch/jmm9413/Virome

Removed P1, 2, and 3 because only had one read each.

cat > samples.txt << EOF
KWNR_S0_L001
P10_S8_L002
P11_S13_L002
P12_S9_L002
P13_S14_L002
P14_S10_L002
P15_S15_L002
P16_S16_L002
P17_S17_L002
P4_S4_L002
P5_S11_L002
P6_S5_L002
P7_S6_L002
P8_S7_L002
P9_S12_L002
EOF

######will run make parallel
#then created report with fastqc and multiqc
cat samples.txt | parallel --lb make fastqc SRR={}

make report

#then trimmed the reads for quality reasons (standard parameters) and selected 21-30 nt reads.
cat samples.txt | parallel --lb make trim SRR={}


make trim SRR=KWNR_S0_L001
make trim SRR=P10_S8_L002
make trim SRR=P11_S13_L002
make trim SRR=P12_S9_L002
make trim SRR=P13_S14_L002
make trim SRR=P14_S10_L002
make trim SRR=P15_S15_L002
make trim SRR=P16_S16_L002
make trim SRR=P17_S17_L002
make trim SRR=P4_S4_L002
make trim SRR=P5_S11_L002
make trim SRR=P6_S5_L002
make trim SRR=P7_S6_L002
make trim SRR=P8_S7_L002
make trim SRR=P9_S12_L002

#command used for readin the files which are compressed
zcat KWNR_S0_L001_R1_001_val_1.fq.gz | head -10


######Once I get the trimmed reads, I will run the fastqc and multiqc again to check the quality of the reads after trimming.

cat samples | parallel --lb make qc_trim SRR={}


#############################For long reads only.



### Moving the reference genome 

cp -r /storage/home/jmm9413/work/Blasting/refs/refs_culex_tarsalis /scratch/jmm9413/Virome

### Aligning the genome
sbatch align1_make



	Trinity --seqType fq \
            --left KWNR_S0_L001_R1_001_val_1.fq.gz \
            --right KWNR_S0_L001_R2_001_val_2.fq.gz \
            --output /scratch/jmm9413/Virome/contigs/combined_trinity_out






####################################
##########################################
####################################
Here for the small reads SE50

cp -r /storage/group/jlr54/default/Sultan_Asad_2022/Wild_tarsalis_sRNA_sequencing_data/Novagene/1735-1771-trimmed /storage/home/jmm9413/scratch/Virome2

or

cp -r /storage/group/jlr54/default/Sultan_Asad_2022/Wild_tarsalis_sRNA_sequencing_data/Novagene/1735-1771-trimmed /storage/work/jmm9413/Virome



samples:

KWNR_S148
P10_S8
P11_S13
P12_S9
P13_S14
P14_S10
P15_S15
P16_S16
P17_S17
P1_S45
P2_S46
P3_S47
P4_S4
P5_S11
P6_S5
P7_S6
P8_S7
P9_S12




#######Qc process of the raw reads (trimmed by novogene)
cat samples.txt | parallel --lb make fastqc SRR={}

#multiqc
make report

########initial trimming
cat samples.txt | parallel --lb make trim SRR={}

make trim SRR=KWNR_S148
make trim SRR=P10_S8        
make trim SRR=P11_S13
make trim SRR=P12_S9
make trim SRR=P13_S14        
make trim SRR=P14_S10        
make trim SRR=P15_S15
make trim SRR=P16_S16
make trim SRR=P17_S17
make trim SRR=P1_S45
make trim SRR=P2_S46
make trim SRR=P3_S47
make trim SRR=P4_S4        
make trim SRR=P5_S11
make trim SRR=P6_S5
make trim SRR=P7_S6
make trim SRR=P8_S7
make trim SRR=P9_S12
ls


make trim2 SRR=KWNR_S148
make trim2 SRR=P10_S8
make trim2 SRR=P11_S13
make trim2 SRR=P12_S9
make trim2 SRR=P13_S14
make trim2 SRR=P14_S10
make trim2 SRR=P15_S15
make trim2 SRR=P16_S16
make trim2 SRR=P17_S17
make trim2 SRR=P1_S45
make trim2 SRR=P2_S46
make trim2 SRR=P3_S47
make trim2 SRR=P4_S4
make trim2 SRR=P5_S11
make trim2 SRR=P6_S5
make trim2 SRR=P7_S6
make trim2 SRR=P8_S7
make trim2 SRR=P9_S12
ls

######qc of trimmed reads

cat samples.txt | parallel --lb make qc_trim SRR={}


######move references .fa .gff and .gtf from Virome to Virome2 (alternativelly can take them from the rasgonlab folder)

cp /scratch/jmm9413/Virome/refs/Culex_tarsalis.fa /scratch/jmm9413/Virome2/refs
cp /scratch/jmm9413/Virome/refs/Culex_tarsalis.gff /scratch/jmm9413/Virome2/refs
cp /scratch/jmm9413/Virome/refs/Culex_tarsalis.gff3 /scratch/jmm9413/Virome2/refs
cp /scratch/jmm9413/Virome/refs/Culex_tarsalis.gtf /scratch/jmm9413/Virome2/refs

cp -r /storage/home/jmm9413/work/Blasting/refs/refs_culex_tarsalis /scratch/jmm9413/Virome2

#############Alignment to Culex tarsalis genome (to extract ummaped reads)

cat samples.txt | parallel --lb make align SRR={}

#####################Now trim the ummaped reads to select the 20-30 bp reads

cat samples.txt | parallel --lb make trim2 SRR={}

#####If I have any problem with sbatch scripts, run the dos2unix command
Example:
dos2unix assembly_make.sh

####Decided to use Kraken2 for mapping the reads against bacteria genome.

##############################################################Do not forget to create the kraken directory
####Need to follow this youtube video https://www.youtube.com/watch?v=rlirpGw8O20

mkdir KRAKEN2_DIR
./install_kraken2.sh KRAKEN2_DIR

####Used the script for download_bac in the makefile


####Then ran thew kraken2 command from makefile

cat samples.txt | parallel --lb make Kraken2 SRR={}



conda remove --name virome kraken2



./install_kraken2.sh /storage/home/jmm9413/scratch/kraken2



kraken2-build --standard --db bact_database

############Used the script 
sbatch build_K_database.sh

kraken2-build --build --db BNAME


kraken2-build --build --threads 24 --db bact_database
########I ran the command without the threads option.... Hope for the best .__. 



##################New attempt. Running a script with this information. 

kraken2-build --download-taxonomy --db refs/bact_database
kraken2-build --download-library bacteria --db refs/bact_database
kraken2-build --build --db refs/bact_database


###########It worked, but it was out of memory, so it gave an error at the end. 

Updated the number of threads used (cores).

kraken2-build --build --threads 62 --db refs/bact_database



##############################finally database assembled!!!

Now runned the Kraken2 script as sbatch





####################Trying to run velvet

make velvet SRR=KWNR_S148_L002_R1_001          ####It worked. 

######Now let's run it for all the samples 


cat samples.txt | parallel --lb make velvet SRR={}



##################################################### Now lets blast. First set up blast reference database for viruses (contigs >50bp but <200)

###########create database
make a blastdb folder and run the following command

mkdir blastdb
update_blastdb.pl ref_viruses_rep_genomes --decompress

########Then collect the contigs from each folder that contain less than 200 bp and more than 50 bp, 
and run blast in all of them

awk '/^>/ {if (length(seq) >= 50 && length(seq) < 200) {print header "\n" seq} header = $0; seq = ""} !/^>/ {seq = seq $0} END {if (length(seq) >= 50 && length(seq) < 200) {print header "\n" seq}}' contigs.fa > contigs50.fa

awk '/^>/ {if (length(seq) >= 200) {print header "\n" seq} header = $0; seq = ""} !/^>/ {seq = seq $0} END {if (length(seq) >= 200) {print header "\n" seq}}' contigs.fa > contigs200.fa

################ I included the previous scripts into a makefile section called split

####################Run split for all contigs

cat samples.txt | parallel --lb make split SRR={}


#####################################################################Blasting sequences that have less than 200bp   (ENVIRONMENT )
######################blast them looking to get species name and gene name

export BLASTDB=$BLASTDB:/scratch/jmm9413/Virome2/blastdb
blastn -db ref_viruses_rep_genomes \
 -query contigs50.fa \
 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames sblastnames" \
 -out contigs50_blastn.out

#######Made a section in the makefile called blast50, so I can run it in all the samples at the same time. 

#############MUST run both lines of code at the same time!!!!  because I was having problems with adding the first line to the makefile
export BLASTDB=$$BLASTDB:/scratch/jmm9413/Virome2/blastdb
cat samples.txt | parallel --lb make blast50 SRR={}

#####Check the variaty of viruses in the output files
cat velvet_output/*/contigs50_blastn.out | cut -13 | sort | head -10

###############Now check which contigs200.fa have more than 0kb

The ones with contigs, I used them for blast in the NCBI website

I will use them for manual blast against everything and use the contigs that match viral genomes for further steps.

########################################Blast results

velvet_output/KWNR_S148_L002_R1_001/contigs200.fa  ###No virus
velvet_output/P14_S10_L002_R1_001/contigs200.fa   ####No virus
velvet_output/P15_S15_L002_R1_001/contigs200.fa    ###### NODE_28_length_279_cov_35.992832  ###Culex Iflavi-like virus 4
velvet_output/P17_S17_L002_R1_001/contigs200.fa    ###### NODE_19_length_328_cov_110.512192 ###Culex Iflavi-like virus 4   ###All contigs of sample
                                                   ###### NODE_36_length_242_cov_54.735538  ###Hubei mosquito virus 4
velvet_output/P1_S45_L002_R1_001/contigs200.fa  #########   
#No virus:    NODE_85_length_271_cov_33.025829    ###NODE_326_length_214_cov_44.869160

velvet_output/P2_S46_L002_R1_001/contigs200.fa  ### Virus
velvet_output/P4_S4_L002_R1_001/contigs200.fa   ### Virus
velvet_output/P6_S5_L002_R1_001/contigs200.fa   ### NO virus ###
velvet_output/P7_S6_L002_R1_001/contigs200.fa   ### NO virus

######################Continue to next step

velvet_output/P15_S15_L002_R1_001/contigs200.fa
velvet_output/P17_S17_L002_R1_001/contigs200.fa
velvet_output/P1_S45_L002_R1_001/contigs200.fa
velvet_output/P2_S46_L002_R1_001/contigs200.fa
velvet_output/P4_S4_L002_R1_001/contigs200.fa

######## I now use an R script to get a combined hit table with the first result, and only keeping viruses.

update_blast_results_with_taxonomy("P15_S15_L002_R1_001.csv", "P15_S15_L002_R1_001.csv_VBlast.csv")
update_blast_results_with_taxonomy("P17_S17_L002_R1_001.csv", "P17_S17_L002_R1_001.csv_VBlast.csv")
update_blast_results_with_taxonomy("P1_S45_L002_R1_001.csv", "P1_S45_L002_R1_001.csv_VBlast.csv")
update_blast_results_with_taxonomy("P2_S46_L002_R1_001.csv", "P2_S46_L002_R1_001.csv_VBlast.csv")
update_blast_results_with_taxonomy("P4_S4_L002_R1_001.csv", "P4_S4_L002_R1_001.csv_VBlast.csv")

I will save those files in a folder called blast200_results. Then I will extract the selected contigs using
the following command in the Virome2 directory.

cut -d ',' -f 1 blast200_results/P15_S15_L002_R1_001_VBlast.csv | tail -n +2 | sed 's/"//g' | \
grep -w -A1 -f - velvet_output/P15_S15_L002_R1_001/contigs200.fa | grep -v '^--$' > velvet_output/P15_S15_L002_R1_001/selected_contigs200.fa

cut -d ',' -f 1 blast200_results/P17_S17_L002_R1_001_VBlast.csv | tail -n +2 | sed 's/"//g' | \
grep -w -A1 -f - velvet_output/P17_S17_L002_R1_001/contigs200.fa | grep -v '^--$' > velvet_output/P17_S17_L002_R1_001/selected_contigs200.fa

cut -d ',' -f 1 blast200_results/P1_S45_L002_R1_001_VBlast.csv | tail -n +2 | sed 's/"//g' | \
grep -w -A1 -f - velvet_output/P1_S45_L002_R1_001/contigs200.fa | grep -v '^--$' > velvet_output/P1_S45_L002_R1_001/selected_contigs200.fa

cut -d ',' -f 1 blast200_results/P2_S46_L002_R1_001_VBlast.csv | tail -n +2 | sed 's/"//g' | \
grep -w -A1 -f - velvet_output/P2_S46_L002_R1_001/contigs200.fa | grep -v '^--$' > velvet_output/P2_S46_L002_R1_001/selected_contigs200.fa

cut -d ',' -f 1 blast200_results/P4_S4_L002_R1_001_VBlast.csv | tail -n +2 | sed 's/"//g' | \
grep -w -A1 -f - velvet_output/P4_S4_L002_R1_001/contigs200.fa | grep -v '^--$' > velvet_output/P4_S4_L002_R1_001/selected_contigs200.fa

velvet_output/P4_S4_L002_R1_001/contigs50_blastn.out

##########Now I will make the samples200.txt file to run the next steps

cat > samples200.txt << EOF
P15_S15_L002_R1_001
P17_S17_L002_R1_001
P1_S45_L002_R1_001
P2_S46_L002_R1_001
P4_S4_L002_R1_001
EOF

######Firs try using one sample
make profile200 SRR=P15_S15_L002_R1_001

###########It worked!!!
#########Now run the makefile for all the samples
cat samples200.txt | parallel --lb make profile200 SRR={}

#################For the following steps I installed the stats environment from biostars handbook. As well as the reshape package of R.


#############got the the plotGeralDistributionPerBaseByReads.pl
##########from here https://github.com/ericgdp/sRNA-virome/blob/main/plotGeralDistributionPerBaseByReads.pl


perl plotGeralDistributionPerBaseByReads.pl -sam ../sam200/P15_S15_L002_R1_001.sam -s 15 -e 30 -p P15 [--plot]
perl plotGeralDistributionPerBaseByReads.pl -sam ../sam200/P17_S17_L002_R1_001.sam -s 15 -e 30 -p P17 [--plot]
perl plotGeralDistributionPerBaseByReads.pl -sam ../sam200/P1_S45_L002_R1_001.sam -s 15 -e 30 -p P1 [--plot]
perl plotGeralDistributionPerBaseByReads.pl -sam ../sam200/P2_S46_L002_R1_001.sam -s 15 -e 30 -p P2 [--plot]
perl plotGeralDistributionPerBaseByReads.pl -sam ../sam200/P4_S4_L002_R1_001.sam -s 15 -e 30 -p P4  [--plot]



######Need to install this first :)
conda install bioconda::perl-bioperl

###################To activate stats environment I should use

conda activate /storage/home/jmm9413/micromamba/envs/stats


#######################IMPORTANT

Have to add this to all perl scripts just before use " Bio::SeqIO "

# Add this line to include the path where BioPerl modules are located
use lib '/storage/home/jmm9413/micromamba/envs/stats/lib/perl5/site_perl/5.22.0';


samtools view -H input.sam > header.sam  # Extract header to a separate file
samtools view -F 4 input.sam | awk -v OFS="\t" '{print $0 > $3".sam"}'  # Split SAM by reference sequence

#Modified looks like this in the makefile!!!
    mkdir -p sam50/$(SRR)
	samtools view -H sam50/$(SRR).sam header.sam  # Extract header to a separate file
	samtools view -F 4 input.sam | awk -v OFS="\t" '{print $$0 > sam50/$(SRR)/$$3".sam"}'  # Split SAM by reference sequence


#####################################################
#######################################################  IMPORTANTTTT

I need to move everything to my advisors folder.

cp -r /storage/home/jmm9413/scratch/Virome2 /storage/group/jlr54/default/Jaime_Manzano

cp -r /storage/home/jmm9413/scratch/CT_Gerard /storage/group/jlr54/default/Jaime_Manzano


I can use this command to go to the folder

/storage/group/jlr54/default/Jaime_Manzano

##################testing script
make split_sam200 SRR=P15_S15_L002_R1_001
samtools view -F 4 sam200/P15_S15_L002_R1_001.sam | awk -v OFS="\t" '{print $0 > "sam200/P15_S15_L002_R1_001/" $3 ".sam"}'


#######Now running the split in all the samples

cat samples200.txt | parallel --lb make split_sam200 SRR={}


#####Toca arreglar el script. Tambien hacer un P15_200sam.txt para hacer un parallel

perl plotGeralDistributionPerBaseByReads.pl -sam ../sam200/P15_S15_L002_R1_001.sam -s 15 -e 30 -p P15 [--plot]

make -f ../makefile

ls ../sam200/P15_S15_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P15_S15_L002_R1_001 CONT200={}
###Now run it for all the samples 
ls ../sam200/P17_S17_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P17_S17_L002_R1_001 CONT200={}
ls ../sam200/P1_S45_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P1_S45_L002_R1_001 CONT200={}

##########Funciono!!!!!!!

ls ../sam200/P15_S15_L002_R1_001/






########################################3 Redoing the contig assembly using velvet 15k and metaspades 15,17k (kmer) 

Run the following scripts for making the assembly
(virome env)
make velvet SRR=KWNR_S148_L002_R1_001
make spades SRR=KWNR_S148_L002_R1_001
(cdhit env)
I created a new environment for cd-hit (cdhit)
make combine_contigs SRR=KWNR_S148_L002_R1_001

It worked!!!

##########################New contig assembly workflow using the makefile.

######Use virome env!!!!!
cat samples.txt | parallel --lb make velvet SRR={}
cat samples.txt | parallel --lb make spades SRR={}


######Use cdhit env
cat samples.txt | parallel --lb make combine_contigs SRR={}

###Then keep the contigs of 200 or more using bioinfo env

cat samples.txt | parallel --lb make split SRR={}



########Now make a list of samples which have more than 0 contigs of 200bp or more

find . -maxdepth 1 -type f -name "*.fasta" -size +0c | sed -e 's|^\./||' -e 's|_200\.fasta$||' > samples200.txt

$ cat samples200.txt 
KWNR_S148_L002_R1_001
P14_S10_L002_R1_001
P15_S15_L002_R1_001
P17_S17_L002_R1_001
P1_S45_L002_R1_001
P2_S46_L002_R1_001
P3_S47_L002_R1_001
P4_S4_L002_R1_001
P5_S11_L002_R1_001
P6_S5_L002_R1_001
P7_S6_L002_R1_001
P8_S7_L002_R1_001

####Saved some stats of the contigs in contig_stats.csv
seqkit stats *.fasta > contig_stats.csv


#############Now I have to blast them and continue the other steps.

In P3, P6 all of them are unknown. So, there is no hit table.

#continue the code with all the other samples from samples200.txt 

update_blast_results_with_taxonomy("KWNR.csv", "KWNR_S148_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P14.csv", "P14_S10_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P15.csv", "P15_S15_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P17.csv", "P17_S17_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P1.csv", "P1_S45_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P2.csv", "P2_S46_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P4.csv", "P4_S4_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P5.csv", "P5_S11_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P7.csv", "P7_S6_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P8.csv", "P8_S7_L002_R1_001_VBlast.csv")

P5 and P7 have only one value in the blast file, so it becomes challenging for the code to process. So, I modified the code!! works now. The problem was that I was using header=TRUE, when it was really false.

#################################################
#################################################I added new lines to the R code, I was able to retrieve a summary table with the name of known
viruses.

################### I moved the content of the folder to a new one callled renamed_contigs

Now align the reads to the contigs

cat samples200.txt | parallel --lb make profile200 SRR={}

####And split the contigs!
cat samples200.txt | parallel --lb make split_sam200 SRR={}



make split_bam21 SRR=KWNR_S148_L002_R1_001
cat samples200.txt | parallel --lb make split_bam21 SRR={}

##############Lets just try to use GVI for visualization




###########For now lets get stats 

ls ../sam200/KWNR_S148_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=KWNR_S148_L002_R1_001 CONT200={}
## 

KWNR_S148_L002_R1_001_renamed.fasta  P1_S45_L002_R1_001_renamed.fasta  P5_S11_L002_R1_001_renamed.fasta
P14_S10_L002_R1_001_renamed.fasta    P2_S46_L002_R1_001_renamed.fasta  P6_S5_L002_R1_001_renamed.fasta
P15_S15_L002_R1_001_renamed.fasta    P3_S47_L002_R1_001_renamed.fasta  P7_S6_L002_R1_001_renamed.fasta
P17_S17_L002_R1_001_renamed.fasta    P4_S4_L002_R1_001_renamed.fasta   P8_S7_L002_R1_001_renamed.fasta

##Fill the code with the rest of the samples 

#############Activate stats env
conda activate /storage/home/jmm9413/micromamba/envs/stats

ls ../sam200/KWNR_S148_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=KWNR_S148_L002_R1_001 CONT200={}
ls ../sam200/P1_S45_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P1_S45_L002_R1_001 CONT200={}
ls ../sam200/P2_S46_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P2_S46_L002_R1_001 CONT200={}
ls ../sam200/P3_S47_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P3_S47_L002_R1_001 CONT200={}
ls ../sam200/P4_S4_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P4_S4_L002_R1_001 CONT200={}
ls ../sam200/P5_S11_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P5_S11_L002_R1_001 CONT200={}
ls ../sam200/P6_S5_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P6_S5_L002_R1_001 CONT200={}
ls ../sam200/P7_S6_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P7_S6_L002_R1_001 CONT200={}
ls ../sam200/P8_S7_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P8_S7_L002_R1_001 CONT200={}
ls ../sam200/P14_S10_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P14_S10_L002_R1_001 CONT200={}
ls ../sam200/P15_S15_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P15_S15_L002_R1_001 CONT200={}
ls ../sam200/P17_S17_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P17_S17_L002_R1_001 CONT200={}



##############Just found out about a mistake in the contigs naming. I need to re run everything since the sam elaboration

copy profiles in pdf into a profiles folder

mkdir -p profiles
cp /storage/group/jlr54/default/Jaime_Manzano/Virome2/stats/sam200/*/*distribution_by_reads.pdf /storage/group/jlr54/default/Jaime_Manzano/Virome2/profiles/



#############################I should consider all reads from the bam file, not only 21nt. So I am running again the bam conversion

cat samples200.txt | parallel --lb make split_bam SRR={}

#####################################Progress!!!!!!!######################

Finished the manual curation, Now I am using cdhit (remember to use the environment dedicated for it) to remove redundancy using this code in the virome directory

cd-hit -i curated_contigs.fasta -o non_redundant_contigs.fasta -c 0.90 -aS 0.90

#########################################

I got a curated list of contigs. Now I need to map all my reads to them 

cat samples.txt | parallel --lb make coocurrence SRR={}



#######################Now I need to count the mapped reads. But I need to create an annotation file from my contigs for that

seqkit fx2tab -l -n -i non_redundant_contigs.fasta > contig_lengths.txt

awk 'BEGIN {OFS="\t"} {print $1, ".", "contig", 1, $2, ".", "+", ".", "gene_id \"" $1 "\";"}' contig_lengths.txt > non_redundant_contigs.gtf

######This one use it with stats environment
featureCounts -F GTF -t contig -a non_redundant_contigs.gtf -o counts.txt co_sam/co_bam/KWNR_S148_L002_R1_001.bam co_sam/co_bam/P14_S10_L002_R1_001.bam co_sam/co_bam/P15_S15_L002_R1_001.bam co_sam/co_bam/P17_S17_L002_R1_001.bam co_sam/co_bam/P1_S45_L002_R1_001.bam co_sam/co_bam/P2_S46_L002_R1_001.bam co_sam/co_bam/P3_S47_L002_R1_001.bam co_sam/co_bam/P4_S4_L002_R1_001.bam co_sam/co_bam/P5_S11_L002_R1_001.bam co_sam/co_bam/P6_S5_L002_R1_001.bam co_sam/co_bam/P7_S6_L002_R1_001.bam co_sam/co_bam/P8_S7_L002_R1_001.bam

####################Updated version
featureCounts -F GTF -t contig -a non_redundant_contigs.gtf -o counts.txt co_sam/co_bam/KWNR_S148_L002_R1_001.bam co_sam/co_bam/P14_S10_L002_R1_001.bam co_sam/co_bam/P15_S15_L002_R1_001.bam co_sam/co_bam/P17_S17_L002_R1_001.bam co_sam/co_bam/P1_S45_L002_R1_001.bam co_sam/co_bam/P2_S46_L002_R1_001.bam co_sam/co_bam/P3_S47_L002_R1_001.bam co_sam/co_bam/P4_S4_L002_R1_001.bam co_sam/co_bam/P5_S11_L002_R1_001.bam co_sam/co_bam/P6_S5_L002_R1_001.bam co_sam/co_bam/P7_S6_L002_R1_001.bam co_sam/co_bam/P8_S7_L002_R1_001.bam co_sam/co_bam/P10_S8_L002_R1_001.bam co_sam/co_bam/P11_S13_L002_R1_001.bam co_sam/co_bam/P12_S9_L002_R1_001.bam co_sam/co_bam/P13_S14_L002_R1_001.bam co_sam/co_bam/P16_S16_L002_R1_001.bam co_sam/co_bam/P9_S12_L002_R1_001.bam


####check chatgpt for RPKM code. Now in a R script in the RPKM folder. 

#####################Now for contig elongation. These are the samples per contig cluster.

Cluster 1: P15, P16, P17, P2, P3, P14, P17.
Cluster 2: P1, P4, P15, P17, P2, P13, P3, P14, P7.
Cluster 3: P1, P4, P13.                       #######Very high hopes on this one!!
Cluster 4: all, but KWNR and P10.
Cluster 5: P1, P4, P16, P15, P17, P2, P6, P13, P8, P9, P3, P14, P7.
Cluster 6: P1, P4, P17, P2, P6, P11, P12, P3, P14, P7.
Cluster 7: P1, P4, P2, P6, P11, P13, P8, P9, P10**, P3
Cluster 8: All.
Cluster 9: P1, P4, P17, P2, P6, P11, P13, P8, P9, P10**, P12, P3, P14, P7.
Cluster 10: P1, P2, P3. 
Cluster 11: All but KWNR. 



Now I need to combine the output files from clusters and use cd-hit on them.

Then remove contigs shorter than 300.

Use this in virome2 folder
cat extended_contigs/*/contigs.fasta > final_contigs/contigs.fasta

then use cd-hit

cd-hit-est -i final_contigs/contigs.fasta -o final_contigs/contigs2.fasta -c 0.95 -n 10

###Finally, reduce the number of the contigs to only the ones with >200 bp
####Use bioinfo environment 

seqtk seq -L 300 final_contigs/contigs2.fasta > final_contigs/contigs_final.fasta
or
seqtk seq -L 200 final_contigs/contigs2.fasta > final_contigs/contigs_final.fasta

seqtk seq -L 200 final_contigs/contigs2.fasta | seqtk seq -U 300 > final_contigs/contigs_final200.fasta

seqtk seq -L 200 contigs2.fasta | seqtk seq -U 300 > contigs_final200.fasta

Then sort them by size 

seqkit sort -l -r final_contigs/contigs_final.fasta > final_contigs/contigs_final_sorted.fasta

Then blast using the final set of contigs.




###########################Getting the reference genome for each of the found viruses
IMPORTANT
Use it in the virome2 folder

Wuhan Mosquito Virus 6
Partitivirus-like Culex mosquito virus
Culex Iflavi-like virus 4
NA
Culex pipiens pallens
Hubei mosquito virus 4
Marma virus
Culex Bunyavirus 2
Culex quinquefasciatus
Culex Iflavi-like virus 1
Culex Bunya-like virus
Tombusviridae sp.


# Download and save the file with wget, redirecting the output to a new file
###Use this in the virome2 folder. They will automatically go to refs 

####No reference for Marma virus   :(

wget -O refs/Wuhan_Mosquito_Virus_6.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/031/525/425/GCA_031525425.1_ASM3152542v1/GCA_031525425.1_ASM3152542v1_genomic.fna.gz
wget -O refs/Partitivirus_like_virus_4.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/031/175/225/GCA_031175225.1_ASM3117522v1/GCA_031175225.1_ASM3117522v1_genomic.fna.gz
wget -O refs/Culex_Iflavi-like_virus_4.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/133/705/GCF_004133705.1_ASM413370v1/GCF_004133705.1_ASM413370v1_genomic.fna.gz
wget -O refs/Hubei_mosquito_virus_4.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/923/215/GCF_001923215.1_ViralProj358929/GCF_001923215.1_ViralProj358929_genomic.fna.gz
wget -O refs/Culex_Bunyavirus_2.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/031/175/155/GCA_031175155.1_ASM3117515v1/GCA_031175155.1_ASM3117515v1_genomic.fna.gz
wget -O refs/Culex_Iflavi-like_virus_1.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/133/665/GCF_004133665.1_ASM413366v1/GCF_004133665.1_ASM413366v1_genomic.fna.gz
wget -O refs/Culex_Bunya_like_virus.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/031/175/145/GCA_031175145.1_ASM3117514v1/GCA_031175145.1_ASM3117514v1_genomic.fna.gz


#################Now I have to map all the reads to each of them.

##first I create the fasta file with all the viruses, and delete the names of the segments of each virus (if any, manually; did for wuhan)

cat Culex_Bunya_like_virus.fa Culex_Bunyavirus_2.fa Culex_Iflavi-like_virus_1.fa Culex_Iflavi-like_virus_4.fa Hubei_mosquito_virus_4.fa Partitivirus_like_virus_4.fa Wuhan_Mosquito_Virus_6.fa > reference_viruses.fa

Now added a section to the makefile called "coocurrence_final". And made a script for that called coocurrence_final_make

bowtie2-build refs/reference_viruses.fa refs/reference_viruses.fa
cat samples.txt | parallel --lb make coocurrence_final SRR={}


#########Then count the reads

seqkit fx2tab -l -n -i refs/reference_viruses.fa > refs/refs_lengths.txt

awk 'BEGIN {OFS="\t"} {print $1, ".", "contig", 1, $2, ".", "+", ".", "gene_id \"" $1 "\";"}' refs/refs_lengths.txt > refs/reference_viruses.gtf

featureCounts -F GTF -t contig -a refs/reference_viruses.gtf -o refs/counts.txt co_sam2/co_bam/KWNR_S148_L002_R1_001.bam co_sam2/co_bam/P14_S10_L002_R1_001.bam co_sam2/co_bam/P15_S15_L002_R1_001.bam co_sam2/co_bam/P17_S17_L002_R1_001.bam co_sam2/co_bam/P1_S45_L002_R1_001.bam co_sam2/co_bam/P2_S46_L002_R1_001.bam co_sam2/co_bam/P3_S47_L002_R1_001.bam co_sam2/co_bam/P4_S4_L002_R1_001.bam co_sam2/co_bam/P5_S11_L002_R1_001.bam co_sam2/co_bam/P6_S5_L002_R1_001.bam co_sam2/co_bam/P7_S6_L002_R1_001.bam co_sam2/co_bam/P8_S7_L002_R1_001.bam co_sam2/co_bam/P10_S8_L002_R1_001.bam co_sam2/co_bam/P11_S13_L002_R1_001.bam co_sam2/co_bam/P12_S9_L002_R1_001.bam co_sam2/co_bam/P13_S14_L002_R1_001.bam co_sam2/co_bam/P16_S16_L002_R1_001.bam co_sam2/co_bam/P9_S12_L002_R1_001.bam


I already got the RPKM2 of them.

Now I need to align only the libraries where the virus was present to the viral genomes to do the profiling of the viruses.


>MH188002.1 Culex Bunya-like virus strain CBunVL/Fresno, complete genome
All but KWNR

>MH188052.1 Culex Bunyavirus 2 strain CBunV2/Fresno, complete genome
All but KWNR

>NC_040646.1 Culex Iflavi-like virus 1 strain CIVL1/Sutter, complete genome
P1, P4

>NC_040832.1 Culex Iflavi-like virus 4 strain CIVL4-Sonoma, complete genome
P1, P4, P15, P17, P2 

>NC_032231.1 Hubei mosquito virus 4 strain 3
P1, P4, P15, P17, P2, P3, P7 

>MH188050.1 Partitivirus-like Culex mosquito virus strain CPLV/Lake, complete genome
P1, P4, P11, P2, P3, P6, P8

>MF176248.1 Wuhan Mosquito Virus 6 strain mos172gb65361 segment 1 to 6, complete sequence
P1, P4, P17, P2, P3, P6, P8, P7, P16, P14, P5 


###########I first need to prepare an index for each virus genome fasta file 

Culex_Bunya_like_virus.fa     Culex_Iflavi-like_virus_4.fa  Wuhan_Mosquito_Virus_6.fa  
Culex_Bunyavirus_2.fa         Hubei_mosquito_virus_4.fa     
Culex_Iflavi-like_virus_1.fa  Partitivirus_like_virus_4.fa  


cat viruses.txt | parallel --lb make index ACC={}



make index ACC=Culex_Bunya_like_virus
make index ACC=Culex_Bunyavirus_2
make index ACC=Culex_Iflavi-like_virus_1
make index ACC=Culex_Iflavi-like_virus_4
make index ACC=Hubei_mosquito_virus_4
make index ACC=Partitivirus_like_virus_4
make index ACC=Wuhan_Mosquito_Virus_6


Now lets prepare the code for the alignment.


make align_Bunya_like ACC=Culex_Bunya_like_virus
make align_Bunyavirus_2 ACC=Culex_Bunyavirus_2
make align_Iflavi_like_1 ACC=Culex_Iflavi-like_virus_1
make align_Iflavi_like_4 ACC=Culex_Iflavi-like_virus_4
make align_Hubei ACC=Hubei_mosquito_virus_4
make align_Partitivirus ACC=Partitivirus_like_virus_4
make align_Wuhan ACC=Wuhan_Mosquito_Virus_6

##################Activate stats env
conda activate /storage/home/jmm9413/micromamba/envs/stats

make -f ../Makefile stats_virus SRR=Culex_Bunya_like_virus
make -f ../Makefile stats_virus SRR=Culex_Bunyavirus_2

######## It takes too much time. I am running the leftover ones in a script bash 

make -f ../Makefile stats_virus SRR=Culex_Iflavi-like_virus_1
make -f ../Makefile stats_virus SRR=Culex_Iflavi-like_virus_4
make -f ../Makefile stats_virus SRR=Hubei_mosquito_virus_4
make -f ../Makefile stats_virus SRR=Partitivirus_like_virus_4
make -f ../Makefile stats_virus SRR=Wuhan_Mosquito_Virus_6


Now convert the sam files to bam files (sorted and indexed)

cat viruses.txt | parallel --lb make viruses_bam SRR={}


##################################Contigs that are NA, but could contain virus RDRP

>final_contig012_387_NA_NA
ttattgaggaagtgcacctccggttctcaaagaaaattgaaggctctaccatctcttttg
gaggagaagagaagaaaaacaagatcagcgaaactgtggagtctctgttctgccatggaa
agcacacgctccggatgcaggggacagaagacgcgaccaaatggaatgaaacgctctcag
cggctttgtttggaatggtgcacaagacgatgctggacgaccccactcgagtcaagtttg
ggctcccgaagatgacagagcaagaaaggatctacctccgcctgtgcatggcttcccatt
tcatcctggcgataaagcg



wuhan mosquito virus = Orthomyxoviridae (segmented -RNA)
Hubei mosquito virus 4 isolate = unclasified riboviria
(iflaviviridae) = Picornavirales (singlke +RNA)
Narnavirus = Narnaviridae (single +RNA)
Partitivirus-like Culex mosquito virus = Partitiviridae(dsRNA)



##########################################Phylogeny
Longest contigs in the RdRP file: 


>final_contig008_590_Partitivirus-like Culex mosquito virus
>final_contig003_532_Wuhan Mosquito Virus 6 isolate
>final_contig002_387_MAG: Hubei mosquito virus 4
>final_contig005_328_MAG: Marma virus isolate
>final_contig007_596_MAG: Culex Iflavi-like virus 4

from the file with 200-300bps contig
>final200_contig065_239_Narnaviridae environmental sample clone sfra.cpip_contig1516 RNA-dependent RNA polymerase gene, complete cds_virus
>final200_contig018_232_MAG: Tombusviridae sp. isolate MEMO007_6B72p genomic sequence_virus

For bunyavirus 2 and culex bunya-like virus we used the longest contig of nucleoprotein

>final_contig005_462_MAG: Culex Bunyavirus 2 isolate CMS002_040a_COAV nucleoprotein gene, complete cds_virus
>final_contig015_594_Culex Bunya-like virus strain CBunVL/Fresno, complete genome_virus



For Narnavirus (200-300bp contig), we did find RdRP genes. The longest will be integrated to the RdRP fasta file. Althought none of 
them gave match to a conserved domain in NCBI conserved domain thing. (could it be that it is outdated??, or it is not a conserved area in the genes???)

For bunyavirus2 and bunya-like-virus the analyses will proceed using the longest contig of the nucleoprotein segment because we did not find
RdRP in the contigs. 

The phylogenetic analyses proceed as explained in the paper. I got the blast hits for each selected contig of each virus, then used 
multiple alignment website (in paper), the did the analyses in the terminal with the tool xxxx.

To run modeltest, I need to change the file type of the multiple alignment ".aln-fasta" to ".fasta" 

The code was placed in the makefile under model_testing function. VIR= will be use for calling the files name.

####################################Model testing!!!

make model_fit VIR=Bunya_like_MA                                   =   K80+I
make model_fit VIR=Bunyavirus_2_MA                                 =   K80+I
make model_fit VIR=Hubei_MA                                        =   TIM2ef+G4
make model_fit VIR=Iflavi_4_MA                                     =   TIM2+I
make model_fit VIR=Marma_MA                                        =   TIM2ef+G4
make model_fit VIR=Narnaviridae_MA                                 =   TPM3uf+G4
make model_fit VIR=Partitivirus_MA                                 =   TIM2+I
make model_fit VIR=Tombusviridae_MA                                =   K80+G4
make model_fit VIR=Wuhan_MA                                        =   TPM1+G4

#######################################Tree building!!!

make tree VIR=Bunya_like_MA MOD=K80+I
make tree VIR=Bunyavirus_2_MA MOD=K80+I
make tree VIR=Hubei_MA MOD=TIM2ef+G4
make tree VIR=Iflavi_4_MA MOD=TIM2+I
make tree VIR=Marma_MA MOD=TIM2ef+G4
make tree VIR=Narnaviridae_MA MOD=TPM3uf+G4
make tree VIR=Partitivirus_MA MOD=TIM2+I
make tree VIR=Tombusviridae_MA MOD=K80+G4
make tree VIR=Wuhan_MA MOD=TPM1+G4






contig for narnavirus phylogeny


##################################
########################################################################
######################################################## New virome alignment

I added new genomes (Tombusviridae, Marma virus) to the refs/reference_viruses.fa file, so I need to build the bowtie index again 

bowtie2-build refs/reference_viruses.fa refs/reference_viruses.fa

awk 'BEGIN {OFS="\t"} !/^@/ {if ($1 in seen) {if ($5 > seen[$1]) {seen[$1] = $5; line[$1] = $0}} else {seen[$1] = $5; line[$1] = $0}} END {for (r in line) print line[r]}' viruses/all.sorted.sam > viruses/best_alignments.sam

or

samtools view -b -F 256 input.bam > primary_alignments.bam

####I prefer the second option better!!!

######This one works for SAM files. It removes secondary, and chimeric alignments 

samtools view -h -F 256 -F 2048 viruses/all.sam > viruses/all_filtered.sam


$ samtools view -c -f 256 all.sam
138054

$ samtools view -c -f 256 all_filtered.sam
0

#########Now I need to process the sam into bam file 

samtools view -bS viruses/all_filtered.sam > viruses/all_filtered.bam
samtools sort viruses/all_filtered.bam -o viruses/all_filtered.bam
samtools index viruses/all_filtered.bam

mkdir -p split
awk '{if($1 ~ /^@/){header=header $0 "\n"} else {print header $0 > "split/"$3".sam"}}' all_filtered.sam

##################Use this version only. Fixed!!

awk '{if(NF && $1 ~ /^@/){header=header $0 "\n"} else if(NF) {file="split/"$3".sam"; if(!(file in files)){print header > file; files[file]=1} print $0 >> file}}' all_filtered.sam

##############It only keeps mapped reads. Deleted the file called "*.sam" 

############################ I also need to prepare a gtf file for the reference_viruses, and ref_lenghts.




seqkit fx2tab -l -n -i refs/reference_viruses.fa > refs/refs_lengths.txt

awk 'BEGIN {OFS="\t"} {print $1, ".", "contig", 1, $2, ".", "+", ".", "gene_id \"" $1 "\";"}' refs/refs_lengths.txt > refs/reference_viruses.gtf


#Change this
featureCounts -F GTF -t contig -a refs/reference_viruses.gtf -o refs/counts.txt co_sam2/co_bam/KWNR_S148_L002_R1_001.bam co_sam2/co_bam/P14_S10_L002_R1_001.bam co_sam2/co_bam/P15_S15_L002_R1_001.bam co_sam2/co_bam/P17_S17_L002_R1_001.bam co_sam2/co_bam/P1_S45_L002_R1_001.bam co_sam2/co_bam/P2_S46_L002_R1_001.bam co_sam2/co_bam/P3_S47_L002_R1_001.bam co_sam2/co_bam/P4_S4_L002_R1_001.bam co_sam2/co_bam/P5_S11_L002_R1_001.bam co_sam2/co_bam/P6_S5_L002_R1_001.bam co_sam2/co_bam/P7_S6_L002_R1_001.bam co_sam2/co_bam/P8_S7_L002_R1_001.bam co_sam2/co_bam/P10_S8_L002_R1_001.bam co_sam2/co_bam/P11_S13_L002_R1_001.bam co_sam2/co_bam/P12_S9_L002_R1_001.bam co_sam2/co_bam/P13_S14_L002_R1_001.bam co_sam2/co_bam/P16_S16_L002_R1_001.bam co_sam2/co_bam/P9_S12_L002_R1_001.bam


################################# I will keep this on pending for now (Not sure why i did that :,v  ). I reran the coocurrence2_make.sh script and now I have to wait 

Then prepared the gtf file as previously explained and ran 
featureCounts -F GTF -t contig -a refs/reference_viruses.gtf -o refs/counts.txt co_sam2/co_bam/KWNR_S148_L002_R1_001.bam co_sam2/co_bam/P14_S10_L002_R1_001.bam co_sam2/co_bam/P15_S15_L002_R1_001.bam co_sam2/co_bam/P17_S17_L002_R1_001.bam co_sam2/co_bam/P1_S45_L002_R1_001.bam co_sam2/co_bam/P2_S46_L002_R1_001.bam co_sam2/co_bam/P3_S47_L002_R1_001.bam co_sam2/co_bam/P4_S4_L002_R1_001.bam co_sam2/co_bam/P5_S11_L002_R1_001.bam co_sam2/co_bam/P6_S5_L002_R1_001.bam co_sam2/co_bam/P7_S6_L002_R1_001.bam co_sam2/co_bam/P8_S7_L002_R1_001.bam co_sam2/co_bam/P10_S8_L002_R1_001.bam co_sam2/co_bam/P11_S13_L002_R1_001.bam co_sam2/co_bam/P12_S9_L002_R1_001.bam co_sam2/co_bam/P13_S14_L002_R1_001.bam co_sam2/co_bam/P16_S16_L002_R1_001.bam co_sam2/co_bam/P9_S12_L002_R1_001.bam


################List of viral genomes, I should change the name of the sam files accordingly.

MH188002.1,Culex_Bunya_like_virus
MH188052.1,Culex_Bunyavirus_2
NC_040646.1,Culex_Iflavi-like_virus_1
NC_040832.1,Culex_Iflavi-like_virus_4     (Sonoma)
NC_032231.1,Hubei_mosquito_virus_4
MH188050.1,Partitivirus-like_Culex_mosquito_virus    ####Change name in poster and graphs!!!
MF176248.1,Wuhan_Mosquito_Virus_6
PP076491.1,Tombusviridae                          (Consider assuming that this is just HMV4)
MW434901.1,Marma_virus
KP642120.1,Narnaviridae_RdRP

I used this for changing names in the RPKM graph, and the sub-sam files from the general alignment
Check line 434 for splitting sam files :)


##########################################consider using this if the other code fails for the general alignment 

USE THIS ONE ONLY FOR SPLITING SAM FILESSSSSSSSS
it requires creating a bam file and then indexing it first.
#################
##############

for contig in $(samtools view -H all_filtered.bam | grep '^@SQ' | awk '{print $2}' | cut -d ':' -f 2); do samtools view -h all_filtered.bam $contig > split/${contig}.sam; done



#############################################################

I ran the viruses_stats2_make.sh script and got the profiles from the viruses in general (at some point I want to contrast them)

Now I will attempt to run the R script in the bam files. So first I need to create the bam files and index them.

###I updated the viruses.txt file to add all the viruses. So Now I can do

cat viruses.txt | parallel --lb make viruses_bam2 SRR={}

Now I have the bam files and their index. I can proceed to make the coverage files. 


#####################################
##########################################################

Exploring contigs 100-200 bps 


in Virome2/final_contigs folder

seqtk seq -L 100 contigs2.fasta | awk '/^>/ {header=$0; next} {if(length($0) <= 199) print header "\n" $0}' > contigs_final100.fasta
################## Now sort them
seqkit sort -l -r contigs_final100.fasta > contigs_final100_sorted.fasta

################Blasted the contigs_final100_sorted.fasta and saved the csv file as contigs100.csv in the contigs folder at one drive


Moving the original reads to my reads folder again 

cp -r /storage/group/jlr54/default/Sultan_Asad_2022/Wild_tarsalis_sRNA_sequencing_data/Novagene/1735-1771-trimmed /storage/group/jlr54/default/Jaime_Manzano/Virome2

cp -r /storage/home/jmm9413/work/Virome/trim /storage/group/jlr54/default/Jaime_Manzano/Virome2

#####Checking maped reads
samtools flagstat P4_S4_L002_R1_001.bam > alignment_stats_P4.txt              Total  64319507     mapped 20631486
samtools flagstat P4_S4_L002_R1_001_unmapped.bam > alignment_stats_P4U.txt          Total 43688021 (remaining)

seqkit stats P4_S4_L002_R1_001_unmapped.fq.gz


samtools fastq P4_S4_L002_R1_001_unmapped.bam -o ../NA_Cu2/P4_S4_L002_R1_001_unmapped_2.fq

##########################
################UPDATE January.

Re aligned the reads to CT genome using the code "create" since the previous code was duplicating the reads in the new unmapped files (Unmapped reads were duplicated)
Which may not affect contig assembly, so contig assembly will not be repeated. However, the rest of the steps should be done again.
Also changed the name of directory NA_Cu to NA_Cu_OLD, and NA_Cu2 to NA_Cu. In order to continue the processing. 

Ran Kaken2.sh 

Now I need to map all my reads to contigs

cat samples.txt | parallel --lb make coocurrence SRR={}

After some consideration, decided to run the stats virus after removing Iflavi 1, bunya-like and tombusviridae. I need to compare profiles before and after. I removed them because the first 5 blast hits were not consistent
and actually some of them fell into existent viruses in the list. 

####I removed the viruses, updated the reference viruses file, indexed them and re ran the analysis using "coocurrence_all"

First need to make the fasta file has the same line length 

seqtk seq -l 60 refs/reference_viruses.fa > refs/reference_viruses_fixed.fa

follow from line 922 once I do the alignments with coocurrence_all

###########Basically this: 

samtools view -h -F 256 -F 2048 viruses/all.sam > viruses/all_filtered.sam

samtools view -bS viruses/all_filtered.sam > viruses/all_filtered.bam
samtools sort viruses/all_filtered.bam -o viruses/all_filtered.bam
samtools index viruses/all_filtered.bam


#Then split the bam files using this

for contig in $(samtools view -H all_filtered.bam | grep '^@SQ' | awk '{print $2}' | cut -d ':' -f 2); do samtools view -h all_filtered.bam $contig > split/${contig}.sam; done

I have to change the name of the viruses in the files for getting the graphs 

I ran the viruses_stats2_make.sh script and got the profiles from the viruses in general (at some point I want to contrast them)

Convert sam to bam for coverage plots
cat viruses.txt | parallel --lb make viruses_bam2 SRR={}


This to get stats of the fq.gz files or fastq.gz
seqkit stats *.fq.gz > stats_summary.txt
seqkit stats *.fastq.gz > stats_summary.txt

Used this to know the amount of reads in each bam files in the "bam" folder (To know the number of reads mapped to Culex)
for file in *1.bam; do samtools flagstat "$file" > "alignment_stats_$(basename "$file" .bam).txt"; done

To know how many reads were between 20 and 30bps 
seqkit seq -m 20 -M 30 *.fq.gz | seqkit stats > filtered_stats.txt

######################################Redid everything again due to the lack of filtering of reads withouth adapters.
Need to run coocurrence again. 


samtools view -h -F 256 -F 2048 viruses/all.sam > viruses/all_filtered.sam

#########Now I need to process the sam into bam file 

samtools view -bS viruses/all_filtered.sam > viruses/all_filtered.bam
samtools sort viruses/all_filtered.bam -o viruses/all_filtered.bam
samtools index viruses/all_filtered.bam
cd viruses

mkdir -p split
for contig in $(samtools view -H all_filtered.bam | grep '^@SQ' | awk '{print $2}' | cut -d ':' -f 2); do samtools view -h all_filtered.bam $contig > split/${contig}.sam; done

Converted sam into file in order to get the coverage files 

cat viruses.txt | parallel --lb make viruses_bam2 SRR={}

featureCounts -F GTF -t contig -a refs/reference_viruses.gtf -o co_sam2/counts.txt co_sam2/co_bam/KWNR_S148_L002_R1_001.bam co_sam2/co_bam/P14_S10_L002_R1_001.bam co_sam2/co_bam/P15_S15_L002_R1_001.bam co_sam2/co_bam/P17_S17_L002_R1_001.bam co_sam2/co_bam/P1_S45_L002_R1_001.bam co_sam2/co_bam/P2_S46_L002_R1_001.bam co_sam2/co_bam/P3_S47_L002_R1_001.bam co_sam2/co_bam/P4_S4_L002_R1_001.bam co_sam2/co_bam/P5_S11_L002_R1_001.bam co_sam2/co_bam/P6_S5_L002_R1_001.bam co_sam2/co_bam/P7_S6_L002_R1_001.bam co_sam2/co_bam/P8_S7_L002_R1_001.bam co_sam2/co_bam/P10_S8_L002_R1_001.bam co_sam2/co_bam/P11_S13_L002_R1_001.bam co_sam2/co_bam/P12_S9_L002_R1_001.bam co_sam2/co_bam/P13_S14_L002_R1_001.bam co_sam2/co_bam/P16_S16_L002_R1_001.bam co_sam2/co_bam/P9_S12_L002_R1_001.bam



#####Re did contig assembly 

cat samples.txt | parallel --lb make velvet SRR={}
cat samples.txt | parallel --lb make spades SRR={}

#Those gave arror, so I am repeating them. They only ran with -k 15
make spades SRR=P9_S12_L002_R1_001
make spades SRR=P13_S14_L002_R1_001
make spades SRR=P11_S13_L002_R1_001

#Now run cdhit and split (To get the 200bp contigs in the contigs folder)
cat samples.txt | parallel --lb make combine_contigs SRR={}
cat samples.txt | parallel --lb make split SRR={}

Now get the stats of the contigs in each of the folders (Combined_output for total contigs, contigs folder for contigs longer than 200)

##In "contigs" folder
seqkit stats *.fasta > contigs_summary.txt

###In combined output
seqkit stats *1_contigs.fasta > contigs_summary.txt


############################################SIDE NOTE#########################
Called the IRS service (8008293903, I was directed to this line, called around 1:10pm in Feb 4th 2025)

1005136765 ID transfered me to an account manager 1004985163.
############################################END OF SIDE NOTE#########################


In P6, the 8 contigs resulted in unknown. So I couldnt download the csv file, But I need to highlight that we got 8 unknown contigs on it.

###############NEED to do the QC of the reads from the trim folder, to double check!!!!!
cat samples.txt | parallel --lb make qc_trim SRR={}

###And this one just once for the multiqc 
make report_trim



####################Blasted the contigs, now put the files into renamed_contigs

cat samples200.txt | parallel --lb make profile200 SRR={}

cat samples200.txt | parallel --lb make split_sam200 SRR={}

#############Activate stats env
conda activate /storage/home/jmm9413/micromamba/envs/stats

ls ../sam200/KWNR_S148_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=KWNR_S148_L002_R1_001 CONT200={}
ls ../sam200/P1_S45_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P1_S45_L002_R1_001 CONT200={}
ls ../sam200/P2_S46_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P2_S46_L002_R1_001 CONT200={}
ls ../sam200/P3_S47_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P3_S47_L002_R1_001 CONT200={}
ls ../sam200/P4_S4_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P4_S4_L002_R1_001 CONT200={}
ls ../sam200/P5_S11_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P5_S11_L002_R1_001 CONT200={}
ls ../sam200/P6_S5_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P6_S5_L002_R1_001 CONT200={}
ls ../sam200/P7_S6_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P7_S6_L002_R1_001 CONT200={}
ls ../sam200/P14_S10_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P14_S10_L002_R1_001 CONT200={}
ls ../sam200/P15_S15_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P15_S15_L002_R1_001 CONT200={}
ls ../sam200/P17_S17_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P17_S17_L002_R1_001 CONT200={}


####Move the contigs into a single folder to facilitate curation. 

Trusted contigs 

##KWNR =
None 
##P1 =
P1_contig006_978_MAG
P1_contig012_422_MAG
P1_contig013_306_MAG
P1_contig014_283_Culex
P1_contig016_276_MAG
P1_contig017_251_MAG
P1_contig019_244_MAG
P1_contig020_242_Culex
##P2 =
P2_contig001_434_Partitivirus-like
P2_contig002_477_NA_NA
P2_contig003_328_MAG
P2_contig004_205_MAG
P2_contig005_239_Culex
P2_contig006_1013_NA_NA
P2_contig007_537_NA_NA
P2_contig010_416_Partitivirus-like
P2_contig011_387_NA_NA
P2_contig012_379_MAG
P2_contig014_354_Partitivirus-like
P2_contig015_335_Partitivirus-like
P2_contig017_284_NA_NA
P2_contig019_250_NA_NA
P2_contig021_239_NA_NA
##P3 = 
P3_contig008_268_NA_NA
##P4 =
P4_contig003_268_NA_NA
##P6 = 
P6_contig003_750_NA_NA
##P7 = 
P7_contig001_419_NA_NA
##P17 = 
P17_contig001_449_Culex

###################Created the curated_contigs.fasta with the selected contigs. Then I have to filter for repetitive ones. Use cdhit environment in the Virome2 folder
cd-hit -i curated_contigs.fasta -o non_redundant_contigs.fasta -c 0.90 -aS 0.90

cat samples.txt | parallel --lb make coocurrence SRR={}

featureCounts -F GTF -t contig -a non_redundant_contigs.gtf -o co_sam/co_bam/counts.txt co_sam/co_bam/KWNR_S148_L002_R1_001.bam co_sam/co_bam/P14_S10_L002_R1_001.bam co_sam/co_bam/P15_S15_L002_R1_001.bam co_sam/co_bam/P17_S17_L002_R1_001.bam co_sam/co_bam/P1_S45_L002_R1_001.bam co_sam/co_bam/P2_S46_L002_R1_001.bam co_sam/co_bam/P3_S47_L002_R1_001.bam co_sam/co_bam/P4_S4_L002_R1_001.bam co_sam/co_bam/P5_S11_L002_R1_001.bam co_sam/co_bam/P6_S5_L002_R1_001.bam co_sam/co_bam/P7_S6_L002_R1_001.bam co_sam/co_bam/P8_S7_L002_R1_001.bam co_sam/co_bam/P10_S8_L002_R1_001.bam co_sam/co_bam/P11_S13_L002_R1_001.bam co_sam/co_bam/P12_S9_L002_R1_001.bam co_sam/co_bam/P13_S14_L002_R1_001.bam co_sam/co_bam/P16_S16_L002_R1_001.bam co_sam/co_bam/P9_S12_L002_R1_001.bam


For the new contig assembly (extension)

Cluster1: P1, P4, P3, P2, P6, P13, P15, P17, P7.
Cluster2: P1, P4, P13
Cluster3: P1, P4, P3, P2, P6, P8, P5, P16, P7, P14
Cluster4: ALL
Cluster5: P1, P4, P11, P3, P2, P6, P13, P9, P8
Cluster6: P1, P4, P11, P3, P2, P6, P17, P7, P14, P12
Cluster7: P1, P4, P11, P3, P2, P6, P7, P14
Cluster8: P1, P4, P11, P3, P2, P6, P13, P9, P5, P7, P14, P12

#run the sbatch script 
clusters_make.sh

Use this in virome2 folder
cat extended_contigs/*/contigs.fasta > final_contigs/contigs.fasta

then use cd-hit

cd-hit-est -i final_contigs/contigs.fasta -o final_contigs/contigs2.fasta -c 0.95 -n 10

###Finally, reduce the number of the contigs to only the ones with >200 bp

####Use bioinfo environment 
seqtk seq -L 200 final_contigs/contigs2.fasta > final_contigs/contigs_final.fasta

Then sort them by size 

seqkit sort -l -r final_contigs/contigs_final.fasta > final_contigs/contigs_final_sorted.fasta

Then blast using the final set of contigs.

I got 145 contigs and got most of the viruses back. (all of them)



MH188050

MH188002.1,Culex_Bunya_like_virus            ######Not considered for now
MH188052.1,Culex_Bunyavirus_2             ###Same
NC_040646.1,Culex_Iflavi-like_virus_1            #############Not considered
NC_040716.1,Culex_Iflavi-like_virus_4         ####Updated to Fresno 
NC_032231.1,Hubei_mosquito_virus_4                    ###########Same
MH188050.1,Partitivirus-like_Culex_mosquito_virus    ####Change name in poster and graphs!!! ####Updated
MF176248.1,Wuhan_Mosquito_Virus_6                      ###########Same
PP076491.1,Tombusviridae                          (Consider assuming that this is just HMV4), ###########it is HMV4
MW434901.1,Marma_virus                               ##############Same
MK628543.1,Culex_narnavirus_1                             #########Updated Culex narnavirus 1 isolate CxNV1



MH188052.1,Culex_Bunyavirus_2
NC_040716.1,Culex_Iflavi-like_virus_4                       same as MH188011.1 (Genbank)
NC_032231.1,Hubei_mosquito_virus_4                          same as KX883008.1 (Genbank)
MH188050.1,Partitivirus-like_Culex_mosquito_virus
MF176248.1,Wuhan_Mosquito_Virus_6
MW434901.1,Marma_virus
MK628543.1,Culex_narnavirus_1


##################Use this to make a summary of contigs 

seqkit stats *_renamed.fasta > fasta_summary.txt



for file in *_renamed.fasta; do
    virus_count=$(seqkit seq "$file" | awk -F'_' '{if ($NF == "virus") print $0}' | wc -l)
    no_virus_count=$(seqkit seq "$file" | awk -F'_' '{if ($NF == "no virus") print $0}' | wc -l)
    unknown_count=$(seqkit seq "$file" | awk '{if ($1 ~ /_NA_NA$/) print $0}' | wc -l)
    echo -e "$file\tVirus: $virus_count\tNo_Virus: $no_virus_count\tUnknown: $unknown_count"
done > fasta_classification_summary.txt



seqkit seq -m 20 -M 30 *.fastq.gz | seqkit stats > filtered_stats.txt

seqkit stats *.fastq.gz > stats_summary.txt