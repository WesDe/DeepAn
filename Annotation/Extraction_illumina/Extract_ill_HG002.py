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
    
