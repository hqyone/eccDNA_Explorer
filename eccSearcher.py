#!/usr/bin/env python
'''
Author: Quanyuan(Leo) He
Email: hqyone@gmail.com
Insititute: Hunan Normal Univeristy
Date: 2020-11-16 00:17:17
LastEditTime: 2020-11-22 01:01:33
LastEditors: Quanyuan(Leo) He
Description: Main function to search the small eccDNAs
FilePath: /eccDNA_Explorer/eccSearcher.py
License: The MIT License (MIT)
'''

import pysam
import re
import pathlib
import sys
import getopt

# https://www.jianshu.com/p/6e99a291f2c8


class ECC_DNA:
    def __init__(self, chro, start, end):
        self.chro = chro
        self.start = start
        self.end = end
        self.evidences = {
            "junctions": [],
            "spans": []
        }

    def getID(self):
        return self.chro+"_"+str(self.start)+"_"+str(self.end)

    def toStr(self):
        return self.chro+"\t"+str(self.start)+"\t"+str(self.end)+"\t" + \
            str(len(self.evidences["junctions"]))+"\t" + \
            str(len(self.evidences["spans"]))+"\t" + \
            ";".join(self.evidences["junctions"])+"\t" + \
            ";".join(self.evidences["spans"]) \



samfile = "temp.sam"
samfile = "aln-pe_sort.bam"

ecc_dic = {}
cur_read_id = ""


def getClosestPair(list1, list2):
    list1.sort()
    list2.sort()
    min_value = float('inf')
    cp1 = None
    cp2 = None
    if len(list1) == 0 or len(list2) == 0:
        return (None, None, None)
    if list1[-1] <= list2[0]:
        cp1 = list1[-1]
        cp2 = list2[0]
        min_value = abs(cp1-cp2)
    elif list1[0] >= list2[-1]:
        cp1 = list1[0]
        cp2 = list2[-1]
        min_value = abs(cp1-cp2)
    else:
        for v1 in list1:
            temp_dis = float('inf')
            temp_v1 = None
            temp_v2 = None
            for v2 in list2:
                if abs(v2-v1) < temp_dis:
                    temp_dis = abs(v2-v1)
                    temp_v1 = v1
                    temp_v2 = v2
                else:
                    break
            if temp_dis < min_value:
                min_value = temp_dis
                cp1 = temp_v1
                cp2 = temp_v2
            else:
                break
    return (min_value, cp1, cp2)

# Search for junction reads


def search_ecc_junctions(hit_dic, min_len=100, max_len=100000):
    """ Search eccDNAs supported by function reads

    Args:
        hit_dic (string): The dictionary for map hits of a read
        min_len (int, optional): The minimun of length of eccDNA. Defaults to 100.
        max_len (int, optional): The maxinum of length of eccDNA. Defaults to 100000.

    Returns:                    
        dictionary: eccDNAs supported by function reads
    """
    ecc_dic = {}
    read_id = hit_dic['id']
    for chrom in hit_dic['data']:
        read_starts = []
        for rt in hit_dic['data'][chrom]:
            for cigar_str in hit_dic['data'][chrom][rt]:
                read_starts += hit_dic['data'][chrom][rt][cigar_str]
                m = re.match(r"^(\d+)[HS](\d+)M$", cigar_str)
                if m:
                    h_len = int(m.groups(0)[0])
                    m_len = int(m.groups(0)[1])
                    re_cigar1 = str(h_len)+"M"+str(m_len)+"H"
                    re_cigar2 = str(h_len)+"M"+str(m_len)+"S"
                    if re_cigar1 in hit_dic['data'][chrom][rt] or re_cigar2 in hit_dic['data'][chrom][rt]:
                        if re_cigar1 in hit_dic['data'][chrom][rt]:
                            re_cigar = re_cigar1
                        if re_cigar2 in hit_dic['data'][chrom][rt]:
                            re_cigar = re_cigar2
                        result = getClosestPair(hit_dic['data'][chrom][rt][cigar_str], list(
                            map(lambda x: x+h_len, hit_dic['data'][chrom][rt][re_cigar])))
                        if result[0] and result[0] >= min_len and result[0] <= max_len:
                            new_ecc_dna = ECC_DNA(chrom, result[1], result[2])
                            ecc_id = new_ecc_dna.getID()
                            if ecc_id not in ecc_dic:
                                ecc_dic[ecc_id] = new_ecc_dna
        # Get the evidences information for junction reads
        for ecc_id, ecc_dna in ecc_dic.items():
            if ecc_dna.start in read_starts:
                ecc_dna.evidences["junctions"].append(read_id)
    return ecc_dic


def combine_ecc(ecc_dic, read_ecc_dic):
    for k, v in read_ecc_dic.items():
        if k in ecc_dic:
            ecc_dic[k].evidences["junctions"] += v.evidences["junctions"]
        else:
            ecc_dic[k] = v
    return ecc_dic


def isEmpty(dic):
    return not bool(dic)


def generateReports(ecc_dic, outfile):
    if isEmpty(ecc_dic):
        return
    OUT = open(outfile, 'w')
    for k, v in ecc_dic.items():
        OUT.write(v.toStr()+"\n")
    OUT.close()


def filter_span_data(span_data, min_len=100, max_len=1000000):
    f_span_data = {}
    for chro in span_data:
        if not isEmpty(span_data[chro]['read1']) and not isEmpty(span_data[chro]['read2']):
            read_pair_ls = []
            while(True):
                read1_starts = list(map(int, span_data[chro]['read1'].keys()))
                read2_ends = list(map(int, span_data[chro]['read2'].keys()))
                (reg_len, read1_start, read2_end) = getClosestPair(
                    read1_starts, read2_ends)
                if (reg_len > min_len and reg_len < max_len):
                    r1_strand = span_data[chro]['read1'][read1_start]
                    r2_strand = span_data[chro]['read2'][read2_end]
                    if (read1_start < read2_end and r1_strand == "-" and r2_strand == "+") or \
                            (read1_start > read2_end and r1_strand == "+" and r2_strand == "-"):
                        read_pair_item = {
                            'r1_start': read1_start,
                            'r1_strand': r1_strand,
                            'r2_end': read2_end,
                            'r2_strand': r2_strand
                        }
                        read_pair_ls.append(read_pair_item)
                    # else:
                    #    print('Not span reads')
                span_data[chro]['read1'].pop(read1_start)
                span_data[chro]['read2'].pop(read2_end)
                if isEmpty(span_data[chro]['read1']) or isEmpty(span_data[chro]['read2']):
                    break
            if len(read_pair_ls) != 0:
                f_span_data[chro] = read_pair_ls
    return f_span_data


def SearchEccDNA(bam, outfile, min_len=100, max_len=1000000):
    """ Main function to search eccDNAs

    Args:
        bam (string): [sorted bam file]
        outfile (string): [output tsv file]
        min_len (int, optional): [The minimum of the eccDNA length]. Defaults to 100.
        max_len (int, optional): [The maximum of the eccDNA length]. Defaults to 1000000.
    """
    # Main function to analysis SAM files
    span_hit_dic = {}
    hit_dic = {}
    ecc_dic = {}
    cur_read_id = ""
    SAM = pysam.AlignmentFile(bam, 'rb')
    for read in SAM:
        rid = read.query_name
        # print(rid)
        if (cur_read_id == ""):
            cur_read_id = rid
            hit_dic = {"id": rid, 'full_match': False,
                       "span_data": {}, "data": {}}
        elif (cur_read_id != rid):
            # reset
            read_ecc_dic = search_ecc_junctions(hit_dic, min_len, max_len)
            if isEmpty(read_ecc_dic):
                if hit_dic['full_match']:
                    f = filter_span_data(
                        hit_dic["span_data"], min_len, max_len)
                    if not isEmpty(f):
                        span_hit_dic[hit_dic['id']] = f
            else:
                ecc_dic = combine_ecc(ecc_dic, read_ecc_dic)
            cur_read_id = rid
            hit_dic = {"id": rid, 'full_match': False,
                       "span_data": {}, "data": {}}

        chrom = read.reference_name
        start = read.reference_start
        end = read.reference_end
        strand = "+"
        if read.is_reverse:
            strand = '-'
        rt = "Unknown"
        if read.is_read1:
            rt = "read1"
        elif read.is_read2:
            rt = "read2"
        cigar_str = read.cigarstring
        if cigar_str:
            m = re.match(r"^(\d+)[HS](\d+)M$", cigar_str)
            n = re.match(r"^(\d+)M(\d+)[HS]$", cigar_str)
            k = re.match(r"^(\d+)M$", cigar_str)
            if k:
                hit_dic['full_match'] = True
                if chrom not in hit_dic['span_data']:
                    hit_dic['span_data'][chrom] = {'read1': {}, 'read2': {}}
                loc = start
                if strand == "-":
                    loc += int(k.groups(0)[0])
                if rt not in hit_dic['span_data'][chrom]:
                    hit_dic['span_data'][chrom][rt] = {}
                hit_dic['span_data'][chrom][rt][loc] = strand
            elif m or n:
                if chrom not in hit_dic['data']:
                    hit_dic['data'][chrom] = {"read1": {}, "read2": {}}

                if read.is_read1:
                    if cigar_str not in hit_dic['data'][chrom]["read1"]:
                        hit_dic['data'][chrom]["read1"][cigar_str] = []
                    hit_dic['data'][chrom]["read1"][cigar_str].append(start)

                if read.is_read2:
                    if cigar_str not in hit_dic['data'][chrom]["read2"]:
                        hit_dic['data'][chrom]["read2"][cigar_str] = []
                    hit_dic['data'][chrom]["read2"][cigar_str].append(start)
        else:
            continue
    # The last reads
    read_ecc_dic = search_ecc_junctions(hit_dic)
    if isEmpty(read_ecc_dic):
        if hit_dic['full_match']:
            f = filter_span_data(hit_dic["span_data"], min_len, max_len)
            if not isEmpty(f):
                span_hit_dic[hit_dic['id']] = f
    else:
        ecc_dic = combine_ecc(ecc_dic, read_ecc_dic)
    SAM.close()

    # Search for span reads
    for ecc_id, ecc in ecc_dic.items():
        chro = ecc.chro
        start = ecc.start
        end = ecc.end
        for span_id, span_item in span_hit_dic.items():
            if chro in span_item:
                for i in span_item[chro]:
                    if (start < min(i["r1_start"], i["r2_end"])) and \
                            (end > max(i["r1_start"], i["r2_end"])):
                        ecc.evidences['spans'].append(span_id)
                        break
    generateReports(ecc_dic, outfile)


# SearchEccDNA("aln-pe.sam", "/home/hqyone/eccDNA_Explorer/ecc_out.bed")

def printHelp():
    print('#############################################################')
    print('# Usage: python eccSearcher.py -i bamfile, --min 100, --max 1000000 -o out.bed')
    print('# -i : The path of the BAM file generate by bwa with -a mode')
    print('# -min : The minium of the length of eccDNAs >')
    print('# -max : The maxium of the length of eccDNAs >')
    print('# -o : The path of the output bed file ')

# Main function


def main(argv):
    wdir = str(pathlib.Path(__file__).parent.absolute())
    config = {
        # Default settings
        #########################################################
        "bamfile": "aln-pe.sam",
        "max_len": 1000000,
        "min_len": 100,
        "out_bed": wdir+"/ecc_out.bed"
    }

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:", ["max=", 'min='])
    except getopt.GetoptError as err:
        print(err)
        printHelp()
        sys.exit(2)
    print('###########################################################################################')
    print('############                      eccSearcher [0.1]                           #############')
    print('############      A tool to identify the small eccDNAs                        #############')
    print('############      Author:  Quanyuan He Ph.D.  Contact: hqyone@hotmail.com     #############')
    print('############      Institution :  School of Medicine, Hunan Normal University  #############')
    print('############      Freely distributed under The MIT License (MIT)              #############')
    print('############                             2020/10/15                           #############')
    print('###########################################################################################')

    for opt, arg in opts:
        if opt == '-h':
            printHelp()
            sys.exit(0)
        elif opt == '-i':
            config["bamfile"] = arg
        elif opt == '--min':
            print("min"+arg)
            config["min_len"] = int(arg)
        elif opt == '--max':
            print("max"+arg)
            config["max_len"] = int(arg)
        elif opt == '-o':
            config["out_bed"] = arg

    print('##  ------------------ The configs are as following ....  -----------------')
    for k, v in config.items():
        print(":".join([k, str(v)])+"\n")
    try:
        SearchEccDNA(config["bamfile"], config["out_bed"],
                     config["min_len"], config["max_len"])
        print("Pipeline finished successfully.")
    except:
        print("Unexpected error:", sys.exc_info()[0])
        exit(1)
    exit(0)


if __name__ == "__main__":
    main(sys.argv[1:])