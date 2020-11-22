#!/usr/bin/env python
'''
Author: Quanyuan(Leo) He
Email: hqyone@gmail.com
Insititute: Hunan Normal Univeristy
Date: 2020-11-16 00:17:17
LastEditTime: 2020-11-21 23:27:03
LastEditors: Quanyuan(Leo) He
Description: Main function to search the small eccDNAs
FilePath: /eccDNA_Explorer/eccSearcher.py
License: The MIT License (MIT)
'''

import pysam
import re
import pathlib
import random
import sys, getopt
import time
import pathlib, shutil

# https://www.jianshu.com/p/6e99a291f2c8

class ECC_DNA:
    def __init__(self, chro, start, end):
        self.chro = chro
        self.start = start
        self.end = end
        self.evidences = {
          "junctions":[],
          "spans":[]
        }
    
    def getID(self):
        return self.chro+"_"+str(self.start)+"_"+str(self.end)

samfile="temp.sam"
samfile = "aln-pe_sort.bam"

ecc_dic ={}
cur_read_id = ""

def getClosestPair(list1, list2):
    list1.sort()
    list2.sort()
    min_value = float('inf')
    cp1=None
    cp2=None
    if len(list1)==0 or len(list2)==0:
        return (None, None, None)
    if list1[-1]<=list2[0]:
        cp1= list1[-1]
        cp2= list2[0]
        min_value = abs(cp1-cp2)
    elif list1[0]>=list2[-1]:
        cp1= list1[0]
        cp2= list2[-1]
        min_value = abs(cp1-cp2)
    else:
        for v1 in list1:
            temp_dis=float('inf')
            temp_v1=None
            temp_v2=None
            for v2 in list2:
                if abs(v2-v1)<temp_dis:
                    temp_dis = abs(v2-v1)
                    temp_v1= v1
                    temp_v2 = v2
                else:
                    break
            if temp_dis<min_value:
                min_value=temp_dis
                cp1= temp_v1
                cp2 = temp_v2
            else:
                break
    return (min_value, cp1, cp2)

# Search for junction reads
def search_ecc_junctions(hit_dic):
    ecc_dic={}
    read_id = hit_dic['id']
    for chrom in hit_dic['data']:
        read_starts=[]
        for rt in hit_dic['data'][chrom]:
            for cigar_str in hit_dic['data'][chrom][rt]:
                read_starts+=hit_dic['data'][chrom][rt][cigar_str]
                m = re.match(r"^(\d+)[HS](\d+)M$",cigar_str)
                if m:
                    h_len=int(m.groups(0)[0])
                    m_len=int(m.groups(0)[1])
                    re_cigar1 =str(h_len)+"M"+str(m_len)+"H"
                    re_cigar2 =str(h_len)+"M"+str(m_len)+"S"
                    if re_cigar1 in hit_dic['data'][chrom][rt] or re_cigar2 in hit_dic['data'][chrom][rt]:
                        if re_cigar1 in hit_dic['data'][chrom][rt]:
                            re_cigar = re_cigar1
                        if re_cigar2 in hit_dic['data'][chrom][rt]:
                            re_cigar = re_cigar2
                        result = getClosestPair(hit_dic['data'][chrom][rt][cigar_str],list(map(lambda x: x+h_len,hit_dic['data'][chrom][rt][re_cigar])))
                        if result[0]!=None:
                            new_ecc_dna = ECC_DNA(chrom, result[1], result[2])
                            ecc_id = new_ecc_dna.getID()
                            if ecc_id not in ecc_dic:
                                ecc_dic[ecc_id]=new_ecc_dna
        # Get the evidences information for junction reads
        for ecc_id, ecc_dna in ecc_dic.items():
            if ecc_dna.start in read_starts:
                ecc_dna.evidences["junctions"].append(read_id)                     
    return ecc_dic

def combine_ecc(ecc_dic, read_ecc_dic):
    for k, v in read_ecc_dic.items():
        if k in ecc_dic:
            ecc_dic[k].evidences["junctions"]+= v.evidences["junctions"]
        else:
            ecc_dic[k]=v
    return ecc_dic

def isEmpty(dic):
    return not bool(dic)

def generateReports(ecc_dic, outfile):
    OUT=open(outfile, 'w')
    OUT.close()

def filter_span_data(span_data, min_len=100, max_len=100000):
    f_span_data={}
    for chro in span_data:
        if not isEmpty(span_data[chro]['read1']) and not isEmpty(span_data[chro]['read2']):
            read_pair_ls=[]
            while(True):
                read1_starts = list(map(int,span_data[chro]['read1'].keys()))
                read2_ends = list(map(int,span_data[chro]['read2'].keys()))
                (reg_len, read1_start, read2_end) = getClosestPair(read1_starts, read2_ends)
                if (reg_len>min_len and reg_len<max_len):
                    
                    r1_strand = span_data[chro]['read1'][read1_start]
                    r2_strand = span_data[chro]['read2'][read2_end]
                    if (read1_start<read2_end and r1_strand=="-" and r2_strand=="+") or \
                        (read1_start>read2_end and r1_strand=="+" and r2_strand=="-"):
                        read_pair_item={
                            'r1_start':read1_start,
                            'r1_strand':r1_strand,
                            'r2_end':read2_end,
                            'r2_strand':r2_strand
                        }
                        read_pair_ls.append(read_pair_item)
                    else:
                        print('s')
                span_data[chro]['read1'].pop(read1_start)
                span_data[chro]['read2'].pop(read2_end)
                if isEmpty(span_data[chro]['read1']) or isEmpty(span_data[chro]['read2']):
                    break
            if len(read_pair_ls)!=0:
                f_span_data[chro]=read_pair_ls
    return f_span_data

def SearchEccDNA(bam, output, min_len=100, max_len=1000000):
    # Main function to analysis SAM files
    span_hit_dic={}
    hit_dic={}
    SAM = pysam.AlignmentFile(bam, 'rb')
    for read in SAM:
        rid = read.query_name
        print(rid)
        if (cur_read_id==""):
            cur_read_id=rid
            hit_dic={"id":rid, 'full_match':False, "span_data":{},"data":{}}
        elif (cur_read_id!=rid):
            # reset
            read_ecc_dic = search_ecc_junctions(hit_dic)
            if isEmpty(read_ecc_dic):
                if hit_dic['full_match']:
                    f = filter_span_data(hit_dic["span_data"], min_len, max_len)
                    if not isEmpty(f):
                        span_hit_dic[hit_dic['id']]=f
            else:
                ecc_dic = combine_ecc(ecc_dic, read_ecc_dic)
            cur_read_id = rid
            hit_dic={"id":rid, 'full_match':False, "span_data":{},"data":{}}
        
        chrom = read.reference_name
        start = read.reference_start
        end = read.reference_end
        strand = "+"
        if read.is_reverse: strand='-'
        rt = "Unknown"
        if read.is_read1:
            rt = "read1"
        elif read.is_read2:
            rt = "read2"
        cigar_str = read.cigarstring
        m = re.match(r"^(\d+)[HS](\d+)M$",cigar_str)
        n = re.match(r"^(\d+)M(\d+)[HS]$",cigar_str)
        k = re.match(r"^(\d+)M$",cigar_str)
        if k:
            hit_dic['full_match'] = True
            if chrom not in hit_dic['span_data']:
                hit_dic['span_data'][chrom]={'read1':{},'read2':{}}
            loc=start
            if strand=="-":
                loc+=int(k.groups(0)[0])
            if rt not in hit_dic['span_data'][chrom]:
                hit_dic['span_data'][chrom][rt]={}
            hit_dic['span_data'][chrom][rt][loc]=strand
        elif m or n:
            if chrom not in hit_dic['data']:
                hit_dic['data'][chrom]={"read1":{}, "read2":{}}
            
            if read.is_read1:
                if cigar_str not in hit_dic['data'][chrom]["read1"]:
                    hit_dic['data'][chrom]["read1"][cigar_str] = []
                hit_dic['data'][chrom]["read1"][cigar_str].append(start)
            
            if read.is_read2:
                if cigar_str not in hit_dic['data'][chrom]["read2"]:
                    hit_dic['data'][chrom]["read2"][cigar_str] = []
                hit_dic['data'][chrom]["read2"][cigar_str].append(start)        
    # The last reads
    read_ecc_dic = search_ecc_junctions(hit_dic)
    if isEmpty(read_ecc_dic):
        if hit_dic['full_match']:
            f = filter_span_data(hit_dic["span_data"], min_len, max_len)
            if not isEmpty(f):
                span_hit_dic[hit_dic['id']]=f
    else:
        ecc_dic = combine_ecc(ecc_dic, read_ecc_dic)
    SAM.close()

    # Search for span reads
    for ecc_id, ecc in ecc_dic.items():
        chro= ecc.chro
        start= ecc.start
        end= ecc.end
        for span_id, span_item in span_hit_dic.items():
            if chro in span_item:
                for i in span_item[chro]:
                    if (start<min(i["r1_start"], i["r2_end"])) and \
                        (end>max(i["r1_start"], i["r2_end"])):
                        ecc.evidences['spans'].append(span_id)
                        break
    print(ecc_dic)

def printHelp():
    print('#############################################################')
    print('# Usage: python eccSearcher.py -t b -f fasta -b bedfile, -t bedtools, -o outbed')
    print('# -t : The processing type, b: geneate bed file with seqeunces')
    print('# -n : The name of project')
    print('# -f : The path of the genomeic FASTA file')
    print('# -b : The path of the tareget regions (BED file)')
    print('# -e : The path of the bedtools2 (avaiable at https://github.com/arq5x/bedtools2), default: >')
    print('# -o : The path of the output file ')

# Main function
def main(argv):
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
    CONFIG.close()
    exit(0)
    
if __name__ == "__main__":
    main(sys.argv[1:])