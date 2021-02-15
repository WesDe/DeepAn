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
import csv
import numpy as np 
import matplotlib.mlab as mlab 
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
import re
import pandas as pd
import collections
from collections import defaultdict
import os
import statistics
from matplotlib import rc
import matplotlib.pyplot as plt
import Function
########################## Localisation #############################
def get_liste_insertion(inpt):
    lec=csv.reader(open(inpt,"r"),delimiter="\t")
    count=0
    liste=[]
    for elt in lec:
        if "#" not in elt[0] and "@" not in elt[0]:
            chrom_inser=re.sub("chr", "",elt[0])
            if Function.is_valid(elt[4]) :
                svseq=elt[4]
            else :
                parser = elt[7].split(';')
                svseq=parser[6].split("=")[1].upper()
            pos=int(elt[1])
            head=chrom_inser+"_"+str(pos)#+"_"+str(len(svseq))
            liste.append(head)
    return liste

def get_info(liste,inpt,dic) :
    liste_info=set() 
    lec=csv.reader(open(inpt,"r"),delimiter=",")
    for elt in lec :
        heads=re.sub(r'[ \(\) \' ]', "",elt[0]).split(",")
        #print(elt,heads)
        head=heads[0]+"_"+heads[1]
        liste_info=set()
        for i in range (1,len(elt)):
            liste_info.add(elt[i])
        if liste_info and head in liste :
            dic[head]=list(liste_info)
    return dic
def get_other(liste,dic,types) :
    count=0
    for elt in liste :
        if elt not in dic :
            #print(elt)
            dic[elt]=[types]
            count+=1
    print(count)
    return dic

def concatenate_dic (dic_1,dic_4):
    concat_dic=defaultdict(list)
    for elt in dic_4 :
        #print(elt)
        head=elt.split("_")[0]+"_"+elt.split("_")[1]
        size=elt.split("_")[2]
        #print(elt,dic_1[head])
        total=[elt.split("_")[0],elt.split("_")[1],size]+dic_1[head]+dic_4[elt]
        #total.append(techno)
        #print(total)
        
        concat_dic[elt]=total
    return concat_dic