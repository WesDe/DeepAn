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
vcf_file=csv.reader(open(sys.argv[1],'r'),delimiter='\t')
otp=csv.writer(open(sys.argv[2],"w"),delimiter="\n")
dictionary = defaultdict(list)


def get_fasta(vcf, otp):
    ide=0
    for element in vcf :
        svseq="0"
        if "#" not in element[0] and "@" not in element[0] :
            chrom = re.sub("chr", "", element[0])
            position = int(element[1])
            if not is_valid(element[4]) :
                parser = element[7].split(';')
                svseq = parser[6].split("=")[1].upper()
            if is_valid(element[4]) :
                svseq=element[4]
            if svseq != "0":
                idef=">"+chrom+"_"+"id"+str(ide)+"_"+position+"_"+str(len(svseq))
                liste=(idef,svseq)
                otp.writerow(liste)
                ide+=1
get_fasta(vcf_file,otp)