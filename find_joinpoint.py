import pysam
import re

class ECC_DNA:
    def __init__(self, chro, start, end):
        self.chro = chro
        self.start = start
        self.end = end
        self.evidences = {
          "forward_junc":[],
          "reverse_junc":[],
          "bridges":[]
        }
    
    def getID(self):
        return self.chro+"_"+str(self.start)+"_"+str(self.end)


samfile="temp.sam"
samfile = "aln-pe.sam"
SAM = pysam.AlignmentFile(samfile, "rb")

ecc_dic ={}
cur_read_id = ""
'''
hit_dic = {
    "id":id
    "data":{
        "chr1":{
            read1:{"50H50M":[168301]}
            read2:{"50H50M":[168301]}
        }
    }
}
'''
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
                    break;
            if temp_dis<min_value:
                min_value=temp_dis
                cp1= temp_v1
                cp2 = temp_v2
            else:
                break;
    return (min_value, cp1, cp2)

print(getClosestPair([1,2,3.60,46],[]))

def search_ecc(hit_dic):
    ecc_dic={}
    read_id = hit_dic['id']
    for chrom in hit_dic['data']:
        read1_starts=[]
        read2_starts=[]
        for rt in hit_dic['data'][chrom]:
            for cigar_str in hit_dic['data'][chrom][rt]:
                if rt=="read1":
                    read1_starts+=hit_dic['data'][chrom][rt][cigar_str]
                if rt=="read2":
                    read2_starts+=hit_dic['data'][chrom][rt][cigar_str]
                m = re.match("^(\d+)[HS](\d+)M$",cigar_str)
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
        # Get the evidences information 
        for ecc_id in ecc_dic:
            ecc_dna = ecc_dic[ecc_id]
            if ecc_dna.start in read1_starts:
                ecc_dic[ecc_id].evidences["forward_junc"].append(read_id)
            if ecc_dna.start in read2_starts:
                ecc_dic[ecc_id].evidences["reverse_junc"].append(read_id)                     
    return ecc_dic



def combine_ecc(ecc_dic, read_ecc_dic):
    for k, v in read_ecc_dic.items():
        if k in ecc_dic:
            ecc_dic[k].evidences["forward_junc"]+= v.evidences["forward_junc"]
            ecc_dic[k].evidences["reverse_junc"]+= v.evidences["reverse_junc"]
        else:
            ecc_dic[k]=v
    return ecc_dic

hit_dic={}
for read in SAM:
    rid = read.query_name
    print(rid)
    if (cur_read_id==""):
        cur_read_id=rid
        hit_dic={"id":rid, "data":{}}
    elif (cur_read_id!=rid):
        # reset
        read_ecc_dic = search_ecc(hit_dic)
        if (len(list(read_ecc_dic.keys()))>0):
            print("ok")
        ecc_dic = combine_ecc(ecc_dic, read_ecc_dic)
        cur_read_id = rid
        hit_dic={"id":rid, "data":{}}
    
    chrom = read.reference_name
    start = read.reference_start
    end = read.reference_end
    cigar_str = read.cigarstring
    m = re.match("^(\d+)[HS](\d+)M$",cigar_str)
    n = re.match("^(\d+)M(\d+)[HS]",cigar_str)
    if m or n:
        if chrom not in hit_dic['data']:
            hit_dic['data'][chrom]={"read1":{}, "read2":{}}
        
        if read.is_read1:
            if cigar_str not in hit_dic['data'][chrom]["read1"]:
                hit_dic['data'][chrom]["read1"][cigar_str] = []
            hit_dic['data'][chrom]["read1"][cigar_str].append(start)
        
        if read.is_read2:
            if cigar_str not in hit_dic['data'][chrom]["read2"]:
                hit_dic['data'][chrom]["read2"][cigar_str] = []
            if n:
                hit_dic['data'][chrom]["read2"][cigar_str].append(start)
read_ecc_dic = search_ecc(hit_dic)
if (len(list(read_ecc_dic.keys()))>0):
    print("ok")
ecc_dic = combine_ecc(ecc_dic, read_ecc_dic)
SAM.close()
print(ecc_dic)
    