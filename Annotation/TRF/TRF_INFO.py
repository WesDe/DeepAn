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

vcf_file_insertion = csv.reader(open(sys.argv[1], 'r'), delimiter='\t')
path_trf=sys.argv[2]
def write_header (outp):
    with open(outp, 'w') as txt:
        txt.write('Contig StartPos EndPos PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy(0-2) Motif Sequence\n')
def convert_readable(txtfile,datfile) :
    with open(txtfile, 'a') as txt:
        chrom = ""
        with open(datfile, 'r') as dat:
            for line in dat:
                splitline = line.split()
                if line.startswith("Sequence:"):
                    chrom = line.split()[1]
                else:
                    # Catch index errors when line is blank
                    try:
                        # Check if in header sequence (all non-header lines start with an int: start pos)
                        try:
                            int(splitline[0])
                        except ValueError:
                            continue
                        txt.write(chrom + ' ' + line)
                    except IndexError:
                        pass

def is_valid(inser):
	allowed = "ATCGatcgNn"
	if all(c in allowed for c in inser):
		return True
	else:
		return False

def find_seq(liste) :
    parser = re.findall(r"[\w']+", liste)
    #print(liste)
    #print(parser[0])
    if is_valid(parser[0]) :
        return parser[0].upper()
    else :
        for elt in parser :
            if elt.upper() == "SEQ" or elt.upper() == "SVSEQ" and is_valid(parser[parser.index(elt)+1]):
                return parser[parser.index(elt)+1].upper()
    raise ValueError("Sequence not found in : ",liste)


def temp_file (couple) :
    outp = csv.writer(open("tmp_seq_1.fa", "w"), delimiter="\n")
    outp.writerow(couple)


def find_trf(vcf_file,trf):
    txt_file="readable_format_TRF.txt"
    for element in vcf_file:
        if "#" not in element[0] and "@" not in element[0]:
            chrom = re.sub("chr", "", element[0])
            position = int(element[1])
            couple = (chrom, position)
            liste=element[4]+";"+element[7]
            svseq = find_seq(liste)
            size=len(svseq)
            #print(len(svseq))
            if len(svseq)>=50 :
                header = ">"+chrom+"_"+str(element[1])+"_"+str(size)
                couple = [header, svseq]
                temp_file(couple)
                cmd2 = trf + " tmp_seq_1.fa 2 5 7 80 10 50 2000 -l 6 -d -h"
                os.system(cmd2)
                convert_readable(txt_file,"tmp_seq_1.fa.2.5.7.80.10.50.2000.dat")
                cmd1 = "rm tmp_seq_1.*"
                os.system(cmd1)
                #break


find_trf(vcf_file_insertion, path_trf)
