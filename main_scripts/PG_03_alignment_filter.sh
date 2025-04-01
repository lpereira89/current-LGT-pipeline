#!/bin/bash
#$ -j y
#$ -l mem=8G
#$ -l rmem=8G
#$ -l h_rt=96:00:00


#####################################################################################
#	Script Name: 	alignment_filter.sh
#	Description: 	classify alignments depending on # of taxa and seq
#	Author:		LTDunning
#	Last updated:	22/02/2022, LPG
#####################################################################################

#### Directories and input files
wd=/mnt/parscratch/users/bo1lpg/NEOF-project/phylobased-LGT
input=${wd}/results_02_blastn_to_aln

#### Parameters
rmname='evm.model'
# Generic genome ID we attached to each CDS ID in script PG_00 - Zea_mays in maize, for example
ID="Zuloagaea_bulbosa"

#### step 1: create directories
cd ${wd}
mkdir -p results_03_alignments_filtered
cd results_03_alignments_filtered

mkdir 10_species_or_more
mkdir less_than_10
mkdir -p 10_species_or_more/less_than_200seqs
mkdir -p 10_species_or_more/more_than_200seqs

#### step 2: copy and clean alignments. Remove irrelevant information from seq ID and conflicting characters
# !!! sms fails to construct the tree if seq names are too long - first step to clean names
# the portion to be removed has to be adjusted - in grasses, current genome database, this parameter works
cp ${input}/*/Fasta_mafft_alignments/* .

ls | grep ${ID} | while read file; do sed -i "s/${rmname}//g" ${file}; done
ls | grep ${ID} | while read file; do sed -i "s/:/_/g" ${file}; done
ls | grep ${ID} | while read file; do sed -i "s/\//_/g" ${file}; done

#### step 3: check number of species and move alignments with <10 taxa
ls * | while read line ; do grep ">" "$line" | cut -f 1,2 -d '_' | sort | uniq | wc -l | sed 's/^/'$line'\t/g' >> number_sp.txt ; done
cat number_sp.txt | awk '$2 < 10' | cut -f 1 | while read line ; do mv "$line" less_than_10 ; done
rm number_sp.txt

#### step 4: move alignments with >200 seq
grep -c ">" * | sed 's/:/\t/g' | awk '$2 >=200' | cut -f 1 | while read line ; do mv "$line" 10_species_or_more/more_than_200seqs/ ; done

#### step 5: move rest of the alignments
mv ${ID}* 10_species_or_more/less_than_200seqs/

#### step 6: create a list of files to make trees for
cd 10_species_or_more
ls less_than_200seqs > less_than_200seqs.txt
