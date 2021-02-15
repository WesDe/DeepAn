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
dictionary = defaultdict(list)

inpt=csv.reader(open(sys.argv[1],"r"),delimiter="\t")
otp= csv.writer(open(sys.argv[2]+"_illumina.vcf", "w"), delimiter="\t")
otp_two= csv.writer(open(sys.argv[2]+"_other_tech.vcf", "w"), delimiter="\t")
found=[]
for elt in inpt :
    liste=[]
    tagged=(elt[0],element[1])

    if "HG2_Ill" in elt[2] :
        liste.append(elt[2])
    sp = re.split('; |=|:|;', elt[7])
    for e in sp : 
        if "Ill" in e :
            if e not in liste and "refine" not in e and "HG2" in e and e != "Illcalls" and e != "Illexactcalls":
                 liste.append(e)
            else :
                if tagged not in found :
                    otp_two.writerow(elt)
                    found.append(tagged)
        else :
            if tagged not in found :
                otp_two.writerow(elt)
                found.append(tagged)
    if len(liste)>=1 :
        otp.writerow(elt)
    else :
        if tagged not in found :
            otp_two.writerow(elt)
            found.append(tagged)
    
