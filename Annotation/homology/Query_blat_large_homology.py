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
ref=open(sys.argv[1],'r')
vcf_file=csv.reader(open(sys.argv[2],'r'),delimiter='\t')
path_blat=sys.argv[3]
dictionary = defaultdict(list)


def get_pos_seq(vcf, dictionary):
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
                dictionary[couple].append(svseq)
    return dictionary

def get_ref_seq(ref,dictionary) :
    reader = SeqIO.parse(ref, "fasta")
    for element in reader :
        chrom = re.sub("chr", "", element.description).split(" ")[0]
        sequence=str(element.seq)
        for key in dictionary :
            if key[0] ==chrom :
               #choice = max([250, len(dictionary[key][0])])
                #print(250, len(dictionary[key][0]),choice)
                start=sequence[key[1]-len(dictionary[key][0]):key[1]].upper()
                end = sequence[key[1]:key[1]+len(dictionary[key][0])].upper()
                dictionary[key].append(start)
                dictionary[key].append(end)
                #print (dictionary)
                #break
            #break
    return dictionary

def write_sequence(element,a) :
    outp = csv.writer(open("tmp_seq_1.fa", "w"), delimiter="\n")
    outp_start = csv.writer(open("tmp_seq_start.fa", "w"), delimiter="\n")
    outp_end = csv.writer(open("tmp_seq_end.fa", "w"), delimiter="\n")
    header = ">"+element[0]+"_"+str(element[1])
    inser = a[0]
    start = a[1]
    end = a[2]
    couple_1 = [header+"_inser_length_"+str(len(inser)), inser]
    couple_start = [header+"_ref_start", start]
    couple_end = [header+"_ref_end", end]
    outp.writerow(couple_1)
    outp_start.writerow(couple_start)
    outp_end.writerow(couple_end)

def align_blastn(dictionary,path_blat) :
    #os.system("echo qseqid  sseqid  pident  length  mismatch    gapopen qstart  qend    sstart  send    evalue  bitscore > concatenate_result.tsv")
    for element in dictionary :
        write_sequence(element,dictionary[element])
        cmd2=path_blat+" tmp_seq_1.fa tmp_seq_start.fa -q=dna -t=dna tmp.psl -noHead"
        os.system(cmd2)
        os.system("cat tmp.psl >> concatenate_result_v2_blast.tsv")
        cmd3 =path_blat+ " tmp_seq_1.fa tmp_seq_end.fa -q=dna -t=dna tmp.psl -noHead"
        os.system(cmd3)
        os.system("cat tmp.psl >> concatenate_result_v2_blast.tsv")
        #break



dictionary=get_pos_seq(vcf_file,dictionary)
dictionary=get_ref_seq(ref,dictionary)
align_blastn(dictionary,path_blat)
