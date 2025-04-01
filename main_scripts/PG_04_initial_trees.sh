#!/bin/bash
#$ -l h_rt=96:00:00
#$ -t 1-21680
#$ -l mem=8G
#$ -l rmem=8G
#$ -j y

#####################################################################################
#	Script Name: 	initial_trees.sh
#	Description: 	use SMS to generate ML trees with best model and 100 bootstrap reps
#	Author:		LTDunning
#	Last updated:	04/05/2022, LPG
#####################################################################################

source /usr/local/extras/Genomics/.bashrc

i=$(expr $SGE_TASK_ID)

#### Directories and input files
wd=/mnt/fastdata/bo1lpg/pangenome-pipeline
fasta_list=${wd}/results_03_alignments_filtered/10_species_or_more/less_than_200seqs.txt
fasta_dir=${wd}/results_03_alignments_filtered/10_species_or_more/less_than_200seqs
outdir=${wd}/results_04_initial_trees/less_than_200seqs
scriptname=PG_04_initial_trees.sh

#### Scripts
# 'make' from this file /shared/dunning_lab/Shared/scripts/programs/sms-1.8.1.zip
# !!! sms fails to construct the tree if seq names are too long - first step to clean names
sms=/shared/dunning_lab/Shared/programs/sms-1.8.1/./sms.sh
fasta_to_phylip=${wd}/nested_scripts/Fasta2Phylip.pl

#### Step 2: create directories and convert fasta to phylip
head -$i ${fasta_list} | tail -1 | while read line ; do mkdir -p ${outdir}/individual/"$line" ; mkdir -p ${outdir}/combined ; mkdir -p ${outdir}/logs ;  cd ${outdir}/individual/"$line" ; perl ${fasta_to_phylip} ${fasta_dir}/"$line" "$line" ; done

#### Step 3: run sms and make a copy of the tree
head -$i ${fasta_list} | tail -1 | while read line ; do cd ${outdir}/individual/"$line" ; ${sms} -i $line -d nt -t -b 100; cp *phyml_tree.txt ${outdir}/combined ; done

#### Step 4: clean up by moving log files
mv ${scriptname}.o*.$i ${outdir}/logs
mv ${scriptname}.e*.$i ${outdir}/logs

#### Step 5: check that none of the tree files are empty
cd ${outdir}/combined
ls | while read line; do echo ${line} >> ../check_empty.txt; \
  [ -s ${line} ] && echo "full" >> ../check_empty.txt;
  [ -s ${line} ] || echo "empty" >> ../check_empty.txt; done
