#!/bin/bash

wdir=$(pwd);
if [[ ! -d "$wdir" ]]; then $wdir="$PWD"; fi
. "$wdir/setup.sh"

echo " #### GenomeFASTA : $GenomeFASTA"
echo " #### GenomeGTF : $GenomeGTF"
echo " #### Cutadapt : $Cutadapt"
echo " #### SeqPrep : $SeqPrep"
echo " #### BWA : $BWA"
echo " #### Bedtools : $BedTools"
echo " #### SAMTools : $SAMTools"

cd $wdir
python simulator.py

fastq1="/Users/hqyone/PycharmProjects/eccDNA_Explorer/test_data/test_1.fastq"
fastq2="/Users/hqyone/PycharmProjects/eccDNA_Explorer/test_data/test_2.fastq"

cut_fastq1="/Users/hqyone/PycharmProjects/eccDNA_Explorer/test_data/cut_test_1.fastq"
cut_fastq2="/Users/hqyone/PycharmProjects/eccDNA_Explorer/test_data/cut_test_2.fastq"

s1_adapter_5=TTTTTTT
s1_adapter_3=GGGGGGG
s2_adapter_5=`echo $s1_adapter_3 | tr ACGTacgt TGCAtgca |rev`
s2_adapter_3=`echo $s1_adapter_5 | tr ACGTacgt TGCAtgca |rev`

#$Cutadapt -g $adapter_5 $fastq1 > $cut_fastq1 2> report1.txt
$Cutadapt --pair-adapters \
-g $s1_adapter_5 \
-a $s1_adapter_3 \
-G $s2_adapter_5 \
-A $s2_adapter_3 \
-o $cut_fastq1 \
-p $cut_fastq2 \
$fastq1 $fastq2

$BWA mem -B 8 -a $GenomeFASTA $cut_fastq1 $cut_fastq2 > aln-pe.sam
#awk '/^*@SQ/ || ($7=="=" && ($9>100 && $9<100000)){print $0;}' aln-pe.sam>temp.sam
#awk '/^*@SQ/ || ($3=="chr1" && $4<169301 && $4>167301){print $0;}' aln-pe.sam>temp.sam
# awk '/^*@SQ/ || ($1~/seq_9\|/ && $3=="chr1" && $4<170301 && $4>166301){print $0;}' aln-pe.sam>temp.sam
awk '($6~/[0-9]+[HS][0-9]+M/ || $6~/[0-9]+M[0-9]+[HS]/) && $4~/^16/ {print $0}' aln-pe.sam> temp.sam
$SAMTools view -bS aln-pe.sam > aln-pe.bam
$SAMTools sort -o aln-pe_sort.bam aln-pe.bam
rm aln-pe.bam
$SAMTools view -u -f 1 -F 12 aln-pe_sort.bam >  aln-pe_sort_map_map.bam
