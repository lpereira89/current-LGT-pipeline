#!/bin/bash
#$ -j y
#$ -l mem=16G
#$ -l rmem=16G
#$ -l h_rt=96:00:00

#####################################################################################
#	Script Name: 	Initial_blastn_filter.sh
#	Description: 	Generate database and run individual scripts for each genome to blast CDS
#	Author:		LTDunning
#	Last updated:	15/02/2022, LPG
#####################################################################################

source /usr/local/extras/Genomics/.bashrc

#### Directories and input files

wd=/mnt/fastdata/bo1lpg/pangenome-pipeline/
# input_files_tab is a list of databases to be used in the LGT identification
# at this step, only genomes! Not transcriptomic data nor short-read sequencing
# the format should be DB_NAME<\t>IDENTIFIER<\t>PATH_TO_FILE
input_files_tab=${wd}/ctrl_files/genomes-n67-screen.txt
DB_Directory=${wd}/BlastDB-combined
genomes=${wd}/ctrl_files/list_genomes.txt
results=${wd}/results_01_initial_blastn_filter

#### Scripts
loop=${wd}/nested_scripts/PG_01b_blastn_loop.sh

#### Step 1: generate a file containing CDS from all genomes in input_files_tab
# add identifier to fasta sequences and copy to a new file in blast directory
mkdir results_01_initial_blastn_filter
cat ${input_files_tab} | while read line ; do fasta_DB=$(echo "$line" | cut -f 3); identifier=$(echo "$line" | cut -f 2) ; DB_name=$(echo "$line" | cut -f 1) ; \
	cat ${fasta_DB} | cut -f 1 -d ' ' | sed 's/>/>'${identifier}'/g' >> ${DB_Directory}/Combined.fa ; done

#### Step 2: make a blastDB - all genomes in one unique database
makeblastdb -in ${DB_Directory}/Combined.fa -dbtype nucl

#### Step 3: run blastn script for each genome
cat ${genomes} | cut -f 1 | while read line ; do cat ${loop} | sed 's/XXXX/'$line'/g' > ${results}/$line.sh ; done
cd  ${results}
source /usr/local/extras/Genomics/.bashrc
ls *sh | while read line ; do qsub "$line" ; done
