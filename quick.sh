#awk '/^*@SQ/ || ($6~/100M/ && $3=="chr1" && $4<170301 && $4>166301){print $0;}' aln-pe.sam>temp.sam###
# @Author: Quanyuan(Leo) He
# @Email: hqyone@gmail.com
# @Insititute: Hunan Normal Univeristy
# @Date: 2020-11-19 15:03:35
# @LastEditTime: 2020-11-22 01:41:55
# @LastEditors: Quanyuan(Leo) He
# @Description:
# @FilePath: /eccDNA_Explorer/quick.sh
# @License: The MIT License (MIT)
###

~/.local/bin/cutadapt --pair-adapters \
-g TTTTTTT \
-g CCCCCCC \
-G CCCCCCC \
-G TTTTTTT \
-o cut_test_1.fastq \
-p cut_test_2.fastq \
/Users/hqyone/PycharmProjects/eccDNA_Explorer/test_data/test_1.fastq \
/Users/hqyone/PycharmProjects/eccDNA_Explorer/test_data/test_2.fastq

#-g TTTTTTT \
#-G AAAAAAA \
#-a GGGGGGG \

#-A TTTTTTT \
#-A GGGGGGG \