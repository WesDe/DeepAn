import sys
import csv
import statistics
import re
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
from Bio.Blast import NCBIWWW
import os
csv.field_size_limit(100000000)
dictionary = defaultdict(list)
#Script looks for homology at 10 bp at breakpoint junction with based 0 annotation

def is_valid(inser):
    allowed = "ATCGatcg"
    if all(c in allowed for c in inser):
        return True
    else:
        return False
def converted_nuc(nuc) :
    return nuc.upper()
def get_pos_seq(vcf,dictionary):
    for element in vcf :
        if "#" not in element[0] and "@" not in element[0] :
            if not is_valid(element[4]) :
                parser = element[7].split(';')
                svseq = parser[6].split("=")[1].upper()
            if is_valid(element[4]) :
                svseq=element[4]
                svseq_start = element[4][0:10].upper()
                svseq_end = element[4][-10:].upper()
                chrom = re.sub("chr", "", element[0])
                position = int(element[1])
                couple = (chrom, position, len(svseq))
                dictionary[couple].append(svseq_start)
                dictionary[couple].append(svseq_end)
                dictionary[couple].append(element[3][-1])
            elif svseq != "0":
                svseq_start=parser[6].split("=")[1][1:11].upper()
                svseq_end = parser[6].split("=")[1][-10:].upper()
                chrom = re.sub("chr", "", element[0])
                position=int(element[1])
                couple=(chrom,position,len(svseq))
                dictionary[couple].append(svseq_start)
                dictionary[couple].append(svseq_end)
                dictionary[couple].append(element[3][-1])
    return dictionary


def get_ref_seq(ref,dictionary) :
    reader = SeqIO.parse(ref, "fasta")
    for element in reader :
        chrom = re.sub("chr", "", element.description).split(" ")[0]
        sequence=str(element.seq)
        for key in dictionary :
            if key[0] ==chrom :
                if sequence[key[1]-1] == dictionary[key][2] :
                    start=sequence[key[1]-11:key[1]-1].upper()
                    end = sequence[key[1]:key[1]+10].upper()
                    dictionary[key].append(start)
                    dictionary[key].append(end)
    return dictionary
def get_based(minus,zero,plus,size):
    if minus==size :
        return -1
    elif zero==size :
        return 0
    elif plus==size :
        return 1
def verif_based (ref,dictionary) :
    reader = SeqIO.parse(ref, "fasta")
    count_base_0=0
    count_base_1=0
    count_base_less = 0
    for element in reader:
        chrom = re.sub("chr", "", element.description).split(" ")[0]
        sequence = str(element.seq)
        for key in dictionary:
            if key[0] == chrom:
                based_0 = sequence[key[1]].upper()
                based_1 = sequence[key[1]+1].upper()
                based_less = sequence[key[1]-1].upper()
                if based_0 ==dictionary[key][2].upper() :
                    count_base_0+=1
                if based_1 == dictionary[key][2].upper():
                    count_base_1 += 1
                if based_less == dictionary[key][2].upper:
                    count_base_less += 1
                #else :
                    #print("error")
                    #return 0
    print(count_base_less, count_base_0, count_base_1, len(dictionary))
    return get_based(count_base_less,count_base_0, count_base_1,len(dictionary))

def get_max_start(a,b) :
    for i in range (0,10) :
        if a[0:10-i]==b[0:10-i]:
            return 10-int(i)
    return int(0)


def get_max_end(a, b):
    for i in range(0, 10):
        if a[i:10] == b[i:10]:
            return 10-int(i)
    return int(0)
def align_blastn(dictionary) :
    otp=csv.writer(open("Microhomology.csv","w"),delimiter="\t")
    count=0
    for element in dictionary :

        header =element[0]+"_"+str(element[1])+"_"+str(element[2])
        if len(dictionary[element])==5 :
            inser_start = dictionary[element][0]
            inser_end = dictionary[element][1]
            left = dictionary[element][3]
            right = dictionary[element][4]
            max_align = get_max_start(inser_start, right)+ get_max_end(inser_end, left)
            count += 1
        else : 
            max_align="NaN"
        otp.writerow([header,max_align])
    #print(count,len(dictionary))


ref = open(sys.argv[1], 'r')
vcf_file = csv.reader(open(sys.argv[2], 'r'), delimiter='\t')
dictionary = get_pos_seq(vcf_file, dictionary)
#a=verif_based(ref, dictionary)
#print(a)
dictionary = get_ref_seq(ref, dictionary)
align_blastn(dictionary)
