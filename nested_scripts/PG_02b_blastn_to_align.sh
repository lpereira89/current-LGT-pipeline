#!/bin/bash
#$ -j y
#$ -l mem=16G
#$ -l rmem=16G
#$ -l h_rt=96:00:00

#####################################################################################
#	Script Name: 	blastn_to_align.sh
#	Description: 	blast query sequences against multiple databases and generate mafft alignments of the blast matches for each blastDB
#	Author:		LTDunning
#	Last updated:	15/02/2021, LPG
#####################################################################################

module load Anaconda3/2022.05

#### Directories and input files
# Genome is defined by script PG_02a_call_blast_to_align.sh - sed function substitutes XXXX by each genome ID
Genome=XXXX
wd=/mnt/fastdata/bo1lpg/pangenome-pipeline/
BlastDB=${wd}/ctrl_files/BlastDBs_n67.txt
# Generate one database per genome (CDS, not genomic DNA) in the folder BlastDB_location - DBName as in BlastDB file
BlastDB_location=${wd}/BlastDB-separate
query_dir=${wd}/results_01_initial_blastn_filter
CDS=${wd}/CDS

#### Script
consensus=${wd}/nested_scripts/consensus.pl

#### Parameters
min_blast_aln_length=300

#### Step 1: create directories
cd ${wd}
mkdir -p results_02_blastn_to_aln/${Genome}
cd results_02_blastn_to_aln/${Genome}

#### Step 2: generate query sequence file and make 1 file for each query sequence
# filter only hits > than min_blast_aln_length
mkdir query
cd query
cat ${query_dir}/${Genome}_4_Blastn_results_tophit_nonNative.txt | awk -v min="$min_blast_aln_length" '$4 >=min' | cut -f 1 \
   | while read line ; do grep "$line" -A 1 ${CDS}/${Genome}/${Genome}_final.cds.fa | cut -f 1 -d ' ' >> intermediate_A ; done
grep ">" intermediate_A | sed 's/>//g' | while read line ; do grep "$line" -A 1 intermediate_A > "$line" ; done
rm intermediate_A
ls > ../query.txt

#### Step 3: blastn each query sequence against all the blastDBs
# outputting sequence id, the length of alignment, and aligned part of the subject sequence.
mkdir ../Blastn_results
source activate blast
cat ${BlastDB} | while read line ; do mkdir ../Blastn_results/"$line" ; done
cat ../query.txt | while read line ; do cat ${BlastDB} | while read line2 ; \
    do blastn -query "$line" -db ${BlastDB_location}/"$line2" -outfmt '6 sseqid length sseq' \
    > ../TEMP ; cat ../TEMP | cut -f 1 | head -n 1 | while read line3 ; do grep "$line3" ../TEMP \
    >> ../Blastn_results/"$line2"/"$line"_"$line2" ; done ; done ; rm ../TEMP ; done

#### Step 4: process blastn results to make a fasta file of all matchs > ${min_blast_aln_length}
mkdir ../Blastn_over_${min_blast_aln_length}
cat ../query.txt | while read line ; do cat ${BlastDB} | while read line2 ; do cat ../Blastn_results/"$line2"/"$line"_"$line2" |  \
    awk -F "\t" 'NF {a[$1]+=$2}  END {for(i in a)print i, a[i]}' | awk -v min="$min_blast_aln_length" '$2 >=min'  |    cut -f 1 -d ' ' |  \
    while read line3 ; do grep "^$line3\s"  ../Blastn_results/"$line2"/"$line"_"$line2" | cut -f 1,3 | sed 's/^/>/g' \
    | sed 's/\s/\n/g' | sed '/>/!s/-//g' >> ../Blastn_over_${min_blast_aln_length}/"$line";  done ; done ; done

#### Step 5: report number of blastn matches above threshold length to compare with final results later
ls ../Blastn_over_${min_blast_aln_length} | while read line ; do grep ">" ../Blastn_over_${min_blast_aln_length}/"$line" \
    | sort | uniq | wc -l | sed 's/^/'$line'\t/g' >> ../blastn_matches_per_gene.txt ; done

#### Step 6: align the blastn fragments to the original gene sequence using MAFFT & unwrap alingment
conda activate mafft 
# --addfragments argument doesn't work in all mafft versions, newest gave me problems (currently using version 7.453)
mkdir ../ALN1-mafft
unset MAFFT_BINARIES
cat ../query.txt  | while read line ; do mafft --addfragments ../Blastn_over_${min_blast_aln_length}/"$line" "$line" > ../ALN1-mafft/"$line" ; done
cat ../query.txt | while read line ; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' \
    < ../ALN1-mafft/"$line" > ../ALN1-mafft/"$line"_unwrap ; done
cat ../query.txt | while read line ; do tail -n +2 ../ALN1-mafft/"$line"_unwrap > ../ALN1-mafft/"$line" ; rm ../ALN1-mafft/"$line"_unwrap ; done
source /usr/local/extras/Genomics/.bashrc

#### Step 7: identify blastn matches represented by only one sequence fragment and move these to their own folder
mkdir ../ALN2-unique
cat ../query.txt | while read line ; do grep ">" ../ALN1-mafft/"$line" | cut -f 2 -d ">" | sort \
    | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | grep "\s1$" | cut -f 1 -d ' ' \
    | while read line2 ; do grep "$line2$" -A 1 ../ALN1-mafft/"$line" >> ../ALN2-unique/"$line" ; done; done

#### Step 8: identify blastn matches represented by more than one sequence fragment and move these to their own folder
mkdir ../ALN3-duplicated
cat ../query.txt | while read line ; do grep ">" ../ALN1-mafft/"$line" | cut -f 2 -d ">" | sort \
    | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | grep -v "\s1$" | cut -f 1 -d ' '  \
    | while read line2 ; do grep "$line2$" -A 1 ../ALN1-mafft/"$line" >> ../ALN3-duplicated/"$line"_"$line2" ; done; done

#### Step 9: generate a single consensus sequence for blastn matches represented by more than one sequence
mkdir ../ALN4-consensus
module load BioPerl/1.7.8-GCCcore-12.2.0
ls ../ALN3-duplicated | while read line ; do perl ${consensus} -in ../ALN3-duplicated/"$line"  -out ../ALN4-consensus/"$line" -iupac ; done
ls ../ALN4-consensus  | while read line ; do sed -i '/>/c\>'$line'' ../ALN4-consensus/"$line" ; done
cat ../query.txt | while read line ; do sed -i 's/>'$line'_/>/g' ../ALN4-consensus/"$line"_*  ; done

#### Step 10: merge the consensus and unique sequences into a single alignment
mkdir ../Fasta_mafft_alignments
cat ../query.txt | while read line ; do cat ../ALN2-unique/"$line" ../ALN4-consensus/"$line"_* > ../Fasta_mafft_alignments/"$line" ; done
cat ../query.txt | while read line ; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < \
    ../Fasta_mafft_alignments/"$line" > ../Fasta_mafft_alignments/"$line"_unwrap ; done
cat ../query.txt | while read line ; do tail -n +2 ../Fasta_mafft_alignments/"$line"_unwrap > ../Fasta_mafft_alignments/"$line" ; \
    rm ../Fasta_mafft_alignments/"$line"_unwrap ; done

#### Step 11: check whether the number of sequences in the alignment is correct - the 'check' files must be manually inspected
ls find -name '../Fasta_mafft_alignments/*' -size 0 -delete
ls ../Fasta_mafft_alignments | while read line ; do grep ">" ../Fasta_mafft_alignments/"$line" | sort | uniq | wc -l | sed 's/^/'$line'\t/g' >> ../sequences_in_alignment.txt ; done
paste ../sequences_in_alignment.txt ../blastn_matches_per_gene.txt | awk 'BEGIN { OFS = "\t" } NR == 0 { $5 = "diff." } NR >= 0 { $5 = $2 - ($4+1) } 1' \
    | cut -f 1,5 | grep -v "\s0$" > ../Check_for_errors.txt
cat ../sequences_in_alignment.txt | grep -v "\s0" > ../sequences_in_alignment_no0.txt
cat ../sequences_in_alignment_no0.txt ../blastn_matches_per_gene.txt | cut -f 1 | sort | uniq -u > ../missing_alignments.txt
cp ../blastn_matches_per_gene.txt ../blastn_matches_per_gene_no-missing-aln.txt
cat ../missing_alignments.txt | while read line ; do sed -i '/'$line'/d' ../blastn_matches_per_gene_no-missing-aln.txt ; done
paste ../sequences_in_alignment_no0.txt ../blastn_matches_per_gene_no-missing-aln.txt \
    | awk 'BEGIN { OFS = "\t" } NR == 0 { $5 = "diff." } NR >= 0 { $5 = $2 - ($4+1) } 1' | cut -f 1,5 | grep -v "\s0$" > ../Check_for_errors2.txt

#### Step 12: remove intermediate files
cd ../
rm -r ALN* Blastn_results
