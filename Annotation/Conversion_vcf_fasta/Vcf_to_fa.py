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