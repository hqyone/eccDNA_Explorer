
# coding=utf-8

#  tRNAExplorer v.1.0 - tRNA profiles and regulation networks
#  Copyright (C) 2020  Dr. Quanyuan He
#  School of Medicine, Hunan Normal University
#  Email: hqyone@hotmail.com
#  Freely distributed under the GNU General Public License (GPLv3)
#
# Simulater generate artificial FASTQ file for eccDNA

import numpy as np
import os
import subprocess
import sys, getopt
import time
import pathlib, shutil
from scipy.stats import truncnorm
import random


def get_rc_sequence(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

# Run BLAST for a bam file
def CreateBED(config):
    try:
        wdir = pathlib.Path(__file__).parent.absolute()

        bedtools = config['bedtools']
        name = config['name']
        db_fasta = config['fasta']
        bed = config['bedfile']
        outdir = config['out_dir']

        print("Begin process for genreate BED file ...")
        cmd_bash = wdir+"/"+name+".sh"
        CMD_FILE = open(cmd_bash, "w")
        out_bed_file = outdir + "/"+name+"_seq.bed"
        # bedtools getfasta -fi hg38.fa -bed hg38_knownGene.bed -bedOut>hg38_knownGene_seq.bed
        bedtools_cmd = bedtools+" getfasta -fi "+db_fasta+" -bed "+bed+ " -bedOut > "+out_bed_file
        print(bedtools_cmd)
        CMD_FILE.write(bedtools_cmd)
        CMD_FILE.close()
        #os.popen("bash " + cmd_bash)
        process = subprocess.Popen("bash " + cmd_bash, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print (process.returncode)
        os.remove(cmd_bash)
        print("Finish processing " + name)
        return out_bed_file
    except():
        print("Unexpected error:", sys.exc_info()[0])
        return ""

def eccRegion(seq, mean, sd, low, upp, number=1):
    L = get_truncated_normal(mean=mean, sd=sd, low=low, upp=len(seq))
    region_array=[]
    for i in range(0, number):
        ecc_len = random.choice (list(map(int, list(L.rvs(1)))))
        start = random.choice (range(0, len(seq)-ecc_len))
        end = start+ecc_len
        sub_str = seq[start: end]
        region_array.append({"seq":sub_str, "start":start, "end":end})
    return region_array  

def circ_amplifySeq(seq, length):
    folds = int(length/len(seq)+1)
    ext_seq = seq*folds
    start = random.choice (range(0, len(ext_seq)-length))
    end = start+length
    return (ext_seq[start: end], start, end)


def GeneateFastQ(config):
    bedfile = config['bedfile']
    fastq1 = config['out_dir']+"/"+config['name']+"_1.fastq"
    fastq2 = config['out_dir']+"/"+config['name']+"_2.fastq"
    ecc_tsv = config['out_dir']+"/"+config['name']+".tsv"
    
    
    col_index=12
    #seqs = getSeq(bedfile)
    BED=open(bedfile,'r')
    TSV=open(ecc_tsv,'w')
    FQ1=open(fastq1,'w')
    FQ2=open(fastq2,'w')

    random.seed(30)
    # Genome extraction setting
    ecc_len_mean = config["ecc_len_mean"]
    ecc_len_std = config["ecc_len_std"]
    ecc_len_min = config["ecc_len_min"]
    ecc_len_max = config["ecc_len_max"]  #2M
    ecc_num = config["ecc_num"] 
    
    # Fragmentation settings
    f_mean = config["f_mean"] 
    f_std = config["f_std"]
    f_min = config["f_min"]
    f_max = config["f_max"]
    f_num = random.choice(range(1,100))
    
    # PCR settings
    adapter_5 = config["adapter_5"]
    adapter_3 = config["adapter_3"]
    read_len = config["read_len"]
    
    # Output_settings
    reads_number = config["reads_number"]
    
    seq_index = 0
    for line in BED:
        line = line.strip()
        contents = line.split('\t')
        chro = contents[0]
        start = int(contents[1])
        #end = int(contents[2])
        if len(contents)<col_index+1:
            continue
        seq = contents[col_index]
        head_str="\t".join(contents[:6])
        if len(seq)>ecc_len_min and len(seq)<ecc_len_max:
            ecc_number = random.choice(range(0,ecc_num+1))
            if ecc_number<=0:
                continue
            # Get eccDNA from a read in bed files
            region_array = eccRegion(seq, ecc_len_mean, ecc_len_std, ecc_len_min, ecc_len_max, ecc_number)
            for r in region_array:
                r_start = start+r['start']
                r_end = start+r['end']
                source_infor= chro+":"+str(r_start)+"-"+str(r_end)
                circ_len = random.choice(range(len(r['seq']), len(r['seq'])*10))*10
                (ecc_t_seq, ecc_t_start, ecc_t_end) = circ_amplifySeq(r['seq'], circ_len)
                f_num = random.choice(range(1,100))
                f_seq_array= eccRegion(ecc_t_seq, f_mean, f_std, f_min, f_max, f_num) #Break to fragment
                f_line = head_str+"\t"+str(r_start)+"\t"+str(r_end)+"\t"+str(circ_len)+"\t"+str(len(f_seq_array))
                for f in f_seq_array:
                    f_seq = f['seq']
                    if len(f_seq)>read_len-len(adapter_5):
                        af_seq = adapter_5+f_seq+adapter_3
                        seq1= af_seq[0:read_len+len(adapter_5)]
                        seq2= get_rc_sequence(af_seq[-(read_len+len(adapter_5)):])
                        tsv_line = "@seq_"+str(seq_index)+"\t"+f_line+"\t"+f_seq+"\n"
                        TSV.write(tsv_line)
                        FQ1.write("@seq_"+str(seq_index)+"|"+source_infor+"\n")
                        FQ1.write(seq1+"\n")
                        FQ1.write("+\n")
                        FQ1.write("I"*len(seq1)+"\n")
                        FQ2.write("@seq_"+str(seq_index)+"|"+source_infor+"\n")
                        FQ2.write(seq2+"\n")
                        FQ2.write("+\n")
                        FQ2.write("I"*len(seq1)+"\n")
                        seq_index+=1
                        if seq_index>=reads_number:
                            break
                if seq_index>=reads_number:
                    break
        if seq_index>=reads_number:
            break
    TSV.close()
    FQ1.close()
    FQ2.close()
    BED.close()

def printHelp():
    print('#############################################################')
    print('# Usage: python eccSimulator.py -t b -f fasta -b bedfile, -t bedtools, -o outbed')
    print('# -t : The processing type, b: geneate bed file with seqeunces')
    print('# -n : The name of project')
    print('# -f : The path of the genomeic FASTA file')
    print('# -b : The path of the tareget regions (BED file)')
    print('# -e : The path of the bedtools2 (avaiable at https://github.com/arq5x/bedtools2), default: >')
    print('# -o : The path of the output dir ')

    print('#############################################################')
    print('# Usage: python eccSimulator.py -t g -n test -b bedfile, -i 12 .... ')
    print('# -t : The processing type, g: geneate FASTQs')
    print('# -n : The name of project')
    print('# -b : Absolute path of tareget regions (BED file) with sequence information')
    print('# -i : The index of sequences in the bed file, default 12')
    print('# -em : The mean of the length of ecc regions, default:400 ')
    print('# -es : The std of the length of ecc regions, default:800 ')
    print('# -ei : The min of the length of ecc regions, default:200 ')
    print('# -ea : The max of the length of ecc regions, default:200000 ')
    print('# -en : The number of the ecc DNA extracted from a regions> ')
    print('# -fm : The mean of the length of ecc regions>')
    print('# -fs : The mean of the length of ecc regions>')
    print('# -fi : The mean of the length of ecc regions>')
    print('# -fa : The mean of the length of ecc regions>')
    print('# -fn : The mean of the length of ecc regions>')
    print('# -a5 : The adapter sequence of 5 terminal')
    print('# -a3 : The adapter sequence of 3 terminal')
    print('# -el : The min and max of ecc fragments, eg: 180,800')
    print('# -l : The pair-end read length, default: 75')
    print('# -r : The number of reads in output FASTQ file, default : 100000')
    print('# -o : The path of the output dir ')
    print('# Output 1: <out_dir>/<name>_ecc.tsv , Reads numbers, processing time for each sample')
    print('# Output 2: <out_dir>/<name>_1.fastq & fastq_2 , Two pair-end FASTQ files')
    print('# Output 3: <out_dir>/ecc_tsv , Two pair-end FASTQ files')

# Main function
def main(argv):
    #print(eccRegion("ATTTATTAGGGGGAACCCCATTT",4,3,2,10,1))
    #print(circ_amplifySeq("ATTTATTAGGGGGAACCCCATTT",800))
    bed_key_ls= [
        # Input settings
        "name",
        "fastq",
        "bedfile",
        "seq_index",
        "bedtools",
        "out_dir"
    ]
    fastq_key_ls = [
        # Input settings
        "name",
        "bedfile",
        "seq_index",

        # Genome ecc settings
        "ecc_len_mean",
        "ecc_len_std",
        "ecc_len_min",
        "ecc_len_max",
        "ecc_num",
        
        # Fragmentation settings
        "f_mean",
        "f_std",
        "f_min",
        "f_max",
        "f_num",
        
        # Sequencing  settings
        "adapter_5",
        "adapter_3",
        "read_len",

        # Output settings
        "reads_number",
        "out_dir"
    ]
    wdir = str(pathlib.Path(__file__).parent.absolute())
    config = {
        # Default settings
        #########################################################
        # Input settings
        "type":"g",
        "name":"test",
        "bedfile":wdir+'/test_data/hg38_knownGene_seq_s.bed',
        "seq_index":12, # 0 based

        # Genome ecc settings
        "ecc_len_mean":400,
        "ecc_len_std":800,
        "ecc_len_min":200,
        "ecc_len_max":200000,  #200K
        "ecc_num":1,
        
        # Fragmentation settings
        "f_mean":300,
        "f_std":200,
        "f_min":100,
        "f_max":1000,
        "f_num":random.choice(range(1,100)),
        
        # PCR settings
        "adapter_5":"TTTTTTT",
        "adapter_3":"GGGGGGG",
        "read_len":100,
        
        # Output_settings
        "reads_number":1000,
        "out_dir":wdir+"/test_data"
    }

    version=0.1
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:t:n:f:b:e:o:i:", ["em",'es','ei','ea','en','fm','fs',"fi","fa","fn","a5","a3","l","r"])
    except getopt.GetoptError as err:
        print(err)
        printHelp()
        sys.exit(2)
    print('###########################################################################################')
    print('############                      eccSimulator ['+str(version)+"]             #############")
    print('############      A tool to generate ecc FASTQ files                          #############')
    print('############      Author:  Quanyuan He Ph.D  Contact: hqyone@hotmail.com      #############')
    print('############      Institution :  School of Medicine, Hunan Normal University  #############')
    print('############  Freely distributed under the GNU General Public License (GPLv3) #############')
    print('############                             2020/10/15                           #############')
    print('###########################################################################################')

    for opt, arg in opts:
        if opt == '-h':
            printHelp()
            sys.exit(0)
        elif opt == '-t':
            config["type"] = arg
        elif opt == '-n':
            config["name"] = arg
        elif opt == '-f':
            config["fasta"] = arg
        elif opt == '-b':
            config["bedfile"] = arg
        elif opt == '-e':
            config["bedtools"] = arg
        elif opt == '-o':
            config["out_dir"] = arg
        elif opt == '-i':
            config["seq_index"] = int(arg)
        elif opt == '-em':
            config["ecc_len_mean"] = int(arg)
        elif opt == '-es':
            config["ecc_len_std"] = int(arg)
        elif opt == '-ei':
            config["ecc_len_min"] = int(arg)
        elif opt == '-ea':
            config["ecc_len_max"] = int(arg)
        elif opt == '-en':
            config["ecc_num"] = int(arg)
        elif opt == '-fm':
            config["f_mean"] = int(arg)
        elif opt == '-fs':
            config["f_std"] = int(arg)
        elif opt == '-fi':
            config["f_min"] = int(arg)
        elif opt == '-fa':
            config["f_max"] = int(arg)
        elif opt == '-fn':
            config["f_num"] = int(arg)
        elif opt == '-a5':
            config["adapter_5"] = arg
        elif opt == '-a3':
            config["adapter_3"] = arg
        elif opt == '-l':
            config["read_len"] = int(arg)
        elif opt == '-r':
            config["reads_number"] = int(arg)

    print('##  ------------------ The configs are as following ....  -----------------')
    ecc_config = config['out_dir']+"/"+config['name']+".config"
    CONFIG = open(ecc_config,'w')
    if config['type']=="b":
        for key in bed_key_ls:
            item = config[key]
            print("## "+key+"="+str(item)+"\n")
            CONFIG.write(key+"="+str(item)+"\n")
        CreateBED(config)
    elif config['type']=="g":
        for key in fastq_key_ls:
            item = config[key]
            print("## "+key+"="+str(item)+"\n")
            CONFIG.write(key+"="+str(item)+"\n")
        GeneateFastQ(config)
    else:
        print("Unknown processing type. Avaliable type are (b/g)")
        exit(1)
    CONFIG.close()
    exit(0)
    
if __name__ == "__main__":
    main(sys.argv[1:])