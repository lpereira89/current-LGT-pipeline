#!/bin/bash
#$ -j y
#$ -l mem=8G
#$ -l rmem=8G
#$ -l h_rt=96:00:00

#####################################################################################
#	Script Name: 	Blastn_to_align_step2.sh
#	Description: 	blast query sequences against multiple databases and generate mafft alingements of the blast matches for each blastDB
#	Author:		LTDunning
#	Last updated:	04/05/2022, LPG
#####################################################################################

source /usr/local/extras/Genomics/.bashrc

#### Directories and input files
wd=/mnt/fastdata/bo1lpg/pangenome-pipeline
BlastDB=${wd}/ctrl_files/BlastDBs_n67.txt
BlastDB_location=${wd}/BlastDB-separate
query_dir=${wd}/results_01_initial_blastn_filter
CDS=${wd}/CDS
good_aln=${wd}/results_04_initial_trees/less_than_200seqs/trees_samesister_min50bootstrap.txt

#### Scripts
consensus=${wd}/nested_scripts/consensus.pl

#### Parameters
min_blast_aln_length=300

#### Step 1: add pan-genomes to blastDB
# !!!! until this step, BlastDB-separate contains one database for each of the genomes from BlastDBs
# Here we add all the accessions from the pangenome being analized
# For example, the 27 maize genomes will be each an independent blast database

ls ${CDS} > ${wd}/ctrl_files/BlastDBs_PG.txt
BlastDB_PG=${wd}/ctrl_files/BlastDBs_PG.txt
cat ${ctrl_files}/list_genomes.txt | while read line ; do cp ${CDS}/${line}/${line}_final.cds.fa ${BlastDB_location} ; done
cd ${BlastDB_location}
# Add DB_ at the beggining of all CDS to be able to differentiate from the original query
cat ${ctrl_files}/list_genomes.txt | while read line ; do sed -i "s/>/>DB_/g" ${line}_final.cds.fa ; cd ../ ; done
cat ${BlastDB_PG} | while read line ; do makeblastdb -in "$line" -dbtype nucl ; done

#### Step 2: add transcriptomes to blastDB
# Currently, copying databases already in a shared directory.
# Adjust if more databases need to be included
blastdb_trans=${wd}/ctrl_files/blastDBS_extra_transcriptomes.txt
blastdb_trans_dir=/shared/christin_lab1/shared/Luke/BLASTDBs/BLAST
cat ${blastdb_trans} | while read line ; do cp ${blastdb_trans_dir}/"$line"* . ; done

#### Step 3: list of all databases
cat ${BlastDB_genome} ${BlastDB_PG} ${blastdb_trans} > ${wd}/ctrl_files/BlastDBs_ALL.txt
BlastDB=${wd}/ctrl_files/BlastDBs_ALL.txt

#### Step 4: generate query sequence file and make one file for each query sequence
cd ${wd}
mkdir results_07_blastn_to_aln
cd results_07_blastn_to_aln
mkdir query
cd query
cat ${CDS}/*/*_final.cds.fa > to_fish
cat ${good_aln} | while read line ; do grep "$line" -A 1 to_fish > "$line" ; done
rm to_fish
ls > ../query.txt

#### Step 5: blastn each query sequence against all the blastDBs outputting sequence id, the length of alignment, and aligned part of the subject sequence.
mkdir ../Blastn_results
cat ${BlastDB} | while read line ; do mkdir ../Blastn_results/"$line" ; done
cat ../query.txt | while read line ; do cat ${BlastDB} | while read line2 ; do blastn -query "$line" -db ${BlastDB_location}/"$line2" -outfmt '6 sseqid length sseq' \
    >> ../Blastn_results/"$line2"/"$line"_"$line2" ; done ; done

#### Step 6: process blastn results to make a fasta file of all matchs > ${min_blast_aln_length}
mkdir ../Blastn_over_${min_blast_aln_length}
cat ../query.txt | while read line ; do cat ${BlastDB} | while read line2 ; do cat ../Blastn_results/"$line2"/"$line"_"$line2" |  \
    awk -F "\t" 'NF {a[$1]+=$2}  END {for(i in a)print i, a[i]}' | awk -v min="$min_blast_aln_length" '$2 >=min'  |    cut -f 1 -d ' ' |  \
    while read line3 ; do grep "^$line3\s"  ../Blastn_results/"$line2"/"$line"_"$line2" | cut -f 1,3 | sed 's/^/>/g' \
    | sed 's/\s/\n/g' | sed '/>/!s/-//g' >> ../Blastn_over_${min_blast_aln_length}/"$line";  done ; done ; done


#### Step 7: report number of blastn matches above threshold length to compare with final results later
ls ../Blastn_over_${min_blast_aln_length} | while read line ; do grep ">" ../Blastn_over_${min_blast_aln_length}/"$line" | sort | uniq | wc -l | sed 's/^/'$line'\t/g' >> ../blastn_matches_per_gene.txt ; done

#### Step 8: align the blastn fragments to the original gene sequence using MAFFT & unwrap alignment
source /usr/local/extras/Genomics/apps/anaconda_python/etc/profile.d/conda.sh
conda activate mafft
mkdir ../ALN1-mafft
unset MAFFT_BINARIES
cat ../query.txt  | while read line ; do mafft --addfragments ../Blastn_over_${min_blast_aln_length}/"$line" "$line" > ../ALN1-mafft/"$line" ; done
cat ../query.txt | while read line ; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ../ALN1-mafft/"$line" > ../ALN1-mafft/"$line"_unwrap ; done
cat ../query.txt | while read line ; do tail -n +2 ../ALN1-mafft/"$line"_unwrap > ../ALN1-mafft/"$line" ; rm ../ALN1-mafft/"$line"_unwrap ; done
source /usr/local/extras/Genomics/.bashrc

#### Step 9: identify blastn matches represented by only one sequence fragment and move these to their own folder
mkdir ../ALN2-unique
cat ../query.txt | while read line ; do grep ">" ../ALN1-mafft/"$line" | cut -f 2 -d ">" | sort | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | \
    grep "\s1$" | cut -f 1 -d ' '  | while read line2 ; do grep "$line2$" -A 1 ../ALN1-mafft/"$line" >> ../ALN2-unique/"$line" ; done; done

#### Step 10: identify blastn matches represented by more than one sequence fragment and move these to their own folder
mkdir ../ALN3-duplicated
cat ../query.txt | while read line ; do grep ">" ../ALN1-mafft/"$line" | cut -f 2 -d ">" | sort | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | \
grep -v "\s1$" | cut -f 1 -d ' '  | while read line2 ; do grep "$line2$" -A 1 ../ALN1-mafft/"$line" >> ../ALN3-duplicated/"$line"_"$line2" ; done; done

#### Step 11: generate a single consensus sequence for blastn matches represented by more than one sequence
mkdir ../ALN4-consensus
ls ../ALN3-duplicated | while read line ; do perl ${consensus} -in ../ALN3-duplicated/"$line"  -out ../ALN4-consensus/"$line" -iupac ; done
ls ../ALN4-consensus  | while read line ; do sed -i '/>/c\>'$line'' ../ALN4-consensus/"$line" ; done
cat ../query.txt | while read line ; do sed -i 's/>'$line'_/>/g' ../ALN4-consensus/"$line"_*  ; done

#### Step 12: merge the consensus and unique sequences into a single alignment
mkdir ../Fasta_mafft_alignments
cat ../query.txt | while read line ; do cat ../ALN2-unique/"$line" ../ALN4-consensus/"$line"_* > ../Fasta_mafft_alignments/"$line" ; done
cat ../query.txt | while read line ; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ../Fasta_mafft_alignments/"$line" > ../Fasta_mafft_alignments/"$line"_unwrap ; done
cat ../query.txt | while read line ; do tail -n +2 ../Fasta_mafft_alignments/"$line"_unwrap > ../Fasta_mafft_alignments/"$line" ; rm ../Fasta_mafft_alignments/"$line"_unwrap ; done

#### Step 13: check correct number of sequences are in the alignment
ls find -name '../Fasta_mafft_alignments/*' -size 0 -delete
ls ../Fasta_mafft_alignments | while read line ; do grep ">" ../Fasta_mafft_alignments/"$line" | sort | uniq | wc -l | sed 's/^/'$line'\t/g' >> ../sequences_in_alignment.txt ; done
paste ../sequences_in_alignment.txt ../blastn_matches_per_gene.txt | awk 'BEGIN { OFS = "\t" } NR == 0 { $5 = "diff." } NR >= 0 { $5 = $2 - ($4+1) } 1' | cut -f 1,5 | grep -v "\s0$" > ../Check_for_errors.txt
cat ../sequences_in_alignment.txt | grep -v "\s0" > ../sequences_in_alignment_no0.txt
cat ../sequences_in_alignment_no0.txt ../blastn_matches_per_gene.txt | cut -f 1 | sort | uniq -u > ../missing_alignments.txt
cp ../blastn_matches_per_gene.txt ../blastn_matches_per_gene_no-missing-aln.txt
cat ../missing_alignments.txt | while read line ; do sed -i '/'$line'/d' ../blastn_matches_per_gene_no-missing-aln.txt ; done
paste ../sequences_in_alignment_no0.txt ../blastn_matches_per_gene_no-missing-aln.txt | awk 'BEGIN { OFS = "\t" } NR == 0 { $5 = "diff." } NR >= 0 { $5 = $2 - ($4+1) } 1' | cut -f 1,5 | grep -v "\s0$" > ../Check_for_errors2.txt

#### Step 14: remove intermediate files
cd ../
rm -r ALN* Blastn_results
