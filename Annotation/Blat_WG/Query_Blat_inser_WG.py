#!/usr/bin/env python3

"""*******************************************************************************
	Name: DeepAn 
	Description: DeepAn aims to provide an automatic annotation of insertion type in vcf file.
	Author: Wesley Delage
	Contact: wesley.delage@irisa.fr, IRISA/Univ Rennes/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
	
	Copyright (C) 2020 Inria
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU Affero General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Affero General Public License for more details.
	You should have received a copy of the GNU Affero General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************"""

import sys
import csv
import statistics
import re
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
from Bio.Blast import NCBIWWW
import os
import os.path
csv.field_size_limit(100000000)
ref=open(sys.argv[1],'r')
vcf_file=csv.reader(open(sys.argv[2],'r'),delimiter='\t')
dictionary = defaultdict(list)
path_blat=sys.argv[3]
path_ref=sys.argv[4]
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
        chrom = re.sub("chr", "", element.description)
        sequence=str(element.seq)
        for key in dictionary :
            if key[0] ==chrom :
                start=sequence[key[1]-300:key[1]+300].upper()
                dictionary[key].append(dictionary[key][0])

    return dictionary


def write_sequence(element,a) :
    outp = csv.writer(open("tmp_seq.fa", "a"), delimiter="\n")
    header = ">"+element[0]+"_"+"id"+"_"+str(element[1])
    inser = a[1]
    couple_1 = [header+"_"+str(len(inser)), inser]
    outp.writerow(couple_1)

def align_blastn(dictionary,path_blat,path_ref_genome) :
    if os.path.isfile('tmp_seq.fa') :
        cmd_a="rm tmp_seq.fa"
        os.sys(cmd_a)
    if os.path.isfile('concatenate_result.tsv') :
        cmd_b="rm concatenate_result.tsv"
        os.sys(cmd_b)
    for element in dictionary :
        write_sequence(element,dictionary[element])
    cmd = path_blat+" "+ path_ref_genome+" tmp_seq.fa -q=dna -t=dna tmp.psl -noHead"
    os.system(cmd)
    os.system("cat tmp.psl >> concatenate_result.tsv")

dictionary=get_pos_seq(vcf_file,dictionary)
dictionary=get_ref_seq(ref,dictionary)
align_blastn(dictionary,path_blat,path_ref)
