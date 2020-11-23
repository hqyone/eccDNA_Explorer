#!/bin/bash
###
# @Author: Quanyuan(Leo) He
# @Email: hqyone@gmail.com
# @Insititute: Hunan Normal Univeristy
# @Date: 2020-11-06 16:08:18
# @LastEditTime: 2020-11-22 02:34:32
# @LastEditors: Quanyuan(Leo) He
# @Description:
# @FilePath: /eccDNA_Explorer/pipeline.sh
# @License: The MIT License (MIT)
###

wdir=$(pwd);
if [[ ! -d "$wdir" ]]; then $wdir="$PWD"; fi
. "$wdir/setup.sh"

echo " #### GenomeFASTA : $GenomeFASTA"
echo " #### GenomeGTF : $GenomeGTF"
echo " #### FastUniq : $FastUniq"
echo " #### Cutadapt : $Cutadapt"
echo " #### SeqPrep : $SeqPrep"
echo " #### BWA : $BWA"
echo " #### Bedtools : $BedTools"
echo " #### SAMTools : $SAMTools"

cd $wdir
echo $wdir
python simulator.py -t g -n test -b $wdir/test_data/hg38_knownGene_seq_s.bed -o $wdir/test_data

fastq1="$wdir/test_data/test_1.fastq"
fastq2="$wdir/test_data/test_2.fastq"

cut_fastq1="$wdir/test_data/cut_test_1.fastq"
cut_fastq2="$wdir/test_data/cut_test_2.fastq"

s1_adapter_5=TTTTTTT
s1_adapter_3=GGGGGGG
s2_adapter_5=`echo $s1_adapter_3 | tr ACGTacgt TGCAtgca |rev`
s2_adapter_3=`echo $s1_adapter_5 | tr ACGTacgt TGCAtgca |rev`

#$Cutadapt -g $adapter_5 $fastq1 > $cut_fastq1 2> report1.txt
$Cutadapt --pair-adapters \
-g $s1_adapter_5 \
-g $s2_adapter_5 \
-G $s2_adapter_5 \
-G $s1_adapter_5 \
-o $cut_fastq1 \
-p $cut_fastq2 \
$fastq1 $fastq2

# Removing the redundant pair reads
echo $cut_fastq1$'\n'$cut_fastq2 > $wdir/test_data/input.file
echo $FastUniq -i $wdir/test_data/input.file -tq -o $wdir/test_data/fu_1.fq -p $wdir/test_data/fu_2.fq
$FastUniq -i $wdir/test_data/input.file -tq -o $wdir/test_data/fu_1.fq -p $wdir/test_data/fu_2.fq

$BWA mem -B 8 -a $GenomeFASTA $cut_fastq1 $cut_fastq2 > aln-pe.sam
# awk '/^*@SQ/ || ($7=="=" && ($9>100 && $9<100000)){print $0;}' aln-pe.sam>temp.sam
# awk '/^*@SQ/ || ($3=="chr1" && $4<169301 && $4>167301){print $0;}' aln-pe.sam>temp.sam
# awk '/^*@SQ/ || ($1~/seq_9\|/ && $3=="chr1" && $4<170301 && $4>166301){print $0;}' aln-pe.sam>temp.sam
# awk '/^*@SQ/ || ($6~/100M/  && $3=="chr1" && $4<170301 && $4>166301){print $0;}' aln-pe.sam>temp.sam
# awk '($6~/[0-9]+[HS][0-9]+M/ || $6~/[0-9]+M[0-9]+[HS]/) && $4~/^16/ {print $0}' aln-pe.sam> temp.sam
$SAMTools view -bS aln-pe.sam > aln-pe.bam
$SAMTools sort aln-pe.bam aln-pe_sort
$SAMTools index aln-pe_sort.bam
#rm aln-pe.bam
#rm aln-pe.sam

python eccSearcher.py -i aln-pe.bam -o ecc_out.bed
# awk '($6~/[0-9]+[HS][0-9]+M/ || $6~/[0-9]+M[0-9]+[HS]/) && $4~/^16/ {print $0}' aln-pe.sam> temp.sam
