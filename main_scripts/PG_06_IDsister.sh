#!/bin/bash
#$ -j y
#$ -l h_rt=24:00:00

#####################################################################################
#	Script Name: 	IDsister.sh
#	Description: 	Edit names in the tree file and determine the sister relationships
#	Author:		LTDunning
#	Last updated:	22/02/2022, LPG
#####################################################################################

source /usr/local/extras/Genomics/.bashrc

#### Directories and input files
wd=/mnt/fastdata/bo1lpg/pangenome-pipeline
mid_trees=${wd}/results_05_rooted_trees
species_list=${wd}/ctrl_files/species_list

#### Scripts
idsis=${wd}/nested_scripts/id_sister.pl

#### Parameters
min_bootstrap=50
target_group="andropogoneae"

#### step 1: run id_sister.pl perl script to look at sister relationships of target sequences
# perl id_sister.pl folder_name species_list_name
# the name of the file and the name of the target sequence have to be identical
# the name has to be 'Genus_species_geneID'
mkdir ${wd}/results_06_IDsister
cd ${wd}/results_06_IDsister
perl ${idsis} ${mid_trees} ${species_list}

#### step 2: filter results based on sister1 and sister2 relationships (i.e. they need to be the same)
cat results_sister | awk '$5 == $8' | grep -v "mixed" | grep -v ${target_group} | cut -f 1 > trees_samesister.txt
cat trees_samesister.txt | while read line; do grep "$line" results_sister >> results_sister_samesister ; done
mkdir midpoint_root_samesister
cat trees_samesister.txt | while read line ; do cp ${mid_trees}/"$line"* midpoint_root_samesister ; done

#### step 3: filter results based on sister1 and sister2 relationships (i.e. they need to be the same)
cat results_sister_samesister | awk -v min="$min_bootstrap" '$2 >=min' | cut -f 1,2 | grep -v "NA" | cut -f 1 > intermediate_1
cat results_sister_samesister | awk -v min="$min_bootstrap" '$3 >=min' | cut -f 1,3 | grep -v "NA" | cut -f 1 > intermediate_2
cat intermediate_1 intermediate_2 | sort | uniq > trees_samesister_min${min_bootstrap}bootstrap.txt
rm intermediate_1 intermediate_2
mkdir trees_samesister_min${min_bootstrap}bootstrap
cat trees_samesister_min${min_bootstrap}bootstrap.txt | while read line ; do cp ${mid_trees}/"$line"* trees_samesister_min${min_bootstrap}bootstrap ; done

~
