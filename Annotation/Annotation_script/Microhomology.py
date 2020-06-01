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

########################### Microhomology ##########################################
def check_validity_hom(qpos,tpos,asize,size,types):
    if types =="start" :
        if qpos+asize>=size-10 and tpos+asize>=size-10 :
            return True
        else :
            return False
    if types=="end":
        if qpos<=10 and tpos<=10 :
            return True
        else :
            return False
def get_mh_blat(align_file):
    inpts=csv.reader(open(align_file,"r"), delimiter='\t')
    dic=defaultdict(list)
    i=0
    for elt in inpts :
        query=elt[9]
        target=elt[13]
        first_match_size=int(elt[-3].split(",")[0])
        type_q=elt[9].split("_")[-1]
        size_t=int(elt[13].split("_")[-1])
        if len(elt[-3].split(","))>2:
            first_match_size=int(elt[-3].split(",")[0])
            last_match_size=int(elt[-3].split(",")[-2])
            qfirst_match_pos=int(elt[-2].split(",")[0])
            tfirst_match_pos=int(elt[-1].split(",")[0])
            qlast_match_pos=int(elt[-2].split(",")[-2])
            tlast_match_pos=int(elt[-1].split(",")[-2])
            #print(elt)
            #print(first_match_size,last_match_size)
            #print(qfirst_match_pos,qlast_match_pos)
            #print(tfirst_match_pos,tlast_match_pos)
            first=check_validity_hom(qfirst_match_pos,tfirst_match_pos,first_match_size,size_t,type_q)
            last=check_validity_hom(qlast_match_pos,tlast_match_pos,last_match_size,size_t,type_q)
            #print(first,last)
            if type_q=="start" and last:
                dic[target].append((qlast_match_pos,tlast_match_pos,last_match_size,type_q))
            if type_q=="end" and first :
                dic[target].append((qfirst_match_pos,tfirst_match_pos,first_match_size,type_q))
        else :
            first_match_size=int(elt[-3].split(",")[0])
            qfirst_match_pos=int(elt[-2].split(",")[0])
            tfirst_match_pos=int(elt[-1].split(",")[0])
            first=check_validity_hom(qfirst_match_pos,tfirst_match_pos,first_match_size,size_t,type_q)
            if first:
                dic[target].append((qfirst_match_pos,tfirst_match_pos,first_match_size,type_q))
            #print(elt)
            #print(first_match_size)
            #print(qfirst_match_pos)
            #print(tfirst_match_pos)
            #print(first)
    return dic

    
def check_pos_mh(liste,types) :
    first_liste=[]
    if types=="start" :
        size=len(liste)
        i=0
        for elt in liste :
            if i==0 :
                first_match=int(elt[0])+int(elt[2])
                first_liste=elt
                i+=1
                continue
            else :
                other_match=int(elt[0])+int(elt[2])
                if other_match> first_match :
                    first_match=other_match
                    first_liste=elt
    if types=="end" :
        i=0
        for elt in liste :
            if i==0 :
                first_match=int(elt[0])
                first_liste=elt
                i+=1
                continue
            else :
                other_match=int(elt[0])
                if other_match< first_match :
                    first_match=other_match
                    first_liste=elt
    return first_liste
def clean_dic(dic):
    dic_cleanse=defaultdict(list)
    for elt in dic :
        liste_start=[]
        liste_end=[]
        for a in dic[elt] : 
            if "start" in a :
                liste_start.append(a)
            if "end" in a :
                liste_end.append(a)
        if len(liste_start)==1:
            liste_def_start=liste_start[0]
            dic_cleanse[elt].append(liste_def_start)
        if len(liste_end)==1 :
            liste_def_end=liste_end[0]
            dic_cleanse[elt].append(liste_def_end)
        if len(liste_start)>1 :
            liste_def_start=check_pos_mh(liste_start,"start")
            dic_cleanse[elt].append(liste_def_start)
        if len(liste_end)>1 :
            liste_def_end=check_pos_mh(liste_end,"end")
            dic_cleanse[elt].append(liste_def_end)
        #break
    return dic_cleanse
def concatenate_mh(dic_clean) :
    dic_size={}
    for elt in dic_clean :
        head=elt.split("_")[0]+"_"+elt.split("_")[1]+"_"+elt.split("_")[-1]
        if len(dic_clean[elt])==2 :
            if dic_clean[elt][-1]=="start":
                liste_start=dic_clean[elt][0]
                liste_end=dic_clean[elt][1]
            else :
                liste_end=dic_clean[elt][0]
                liste_start=dic_clean[elt][1]
            if int(liste_start[0])+int(liste_start[2])>=int(elt.split("_")[-1])-10 and int(liste_end[0])<=10 :
                if int(liste_start[0]) >= int(liste_end[0]) and int(liste_end[0])+int(liste_end[2])<=int(liste_start[0]) :
                    dic_size[head]=int(liste_start[2])+int(liste_end[2])
            else :
                if int(liste_start[2])> int(liste_end[2]) :
                    dic_size[head]=int(liste_start[2])
                else :
                    dic_size[head]=int(liste_end[2])
        elif len(dic_clean[elt])==1 :
            dic_size[head]=int(dic_clean[elt][0][2])
    return dic_size

def get_small_mh(csv_mh) :
    lec=csv.reader(open(csv_mh,'r'),delimiter='\t')
    dic=defaultdict()
    for elt in lec :
        dic[elt[0]]=elt[1]
    return dic
def fill_type(dic_small,dic_large,dic_all,dic_out) :
    dic_convert=get_convert(dic_all)
    for elt in dic_large :
        #dic_out[elt].append(dic_large[elt][0])
        dic_out[elt].append(dic_large[elt])

    for a in dic_small :
        #print(a)
        if a not in dic_all :
            head=get_converted(dic_convert,a)
        else :
            head=a
        if head not in dic_out :
            dic_out[head].append(dic_small[a])
    for v in dic_all :
        if v not in dic_out :
            dic_out[v].append(0)
    return dic_out


def Assign_type(dic,liste_type,types) :
    tmp=[]
    for elt in liste_type :
        if elt in dic :
            dic[elt].append(types)
            tmp.append(dic[elt][0])
    #print(tmp)
    return dic
def get_convert(dic) :
    dic_b=defaultdict()
    for elt in dic :
        head_dic=elt.split("_")[0]+"_"+elt.split("_")[1]
        dic_b[head_dic]=elt
    return dic_b
def get_converted(dic,a):
    head=a.split("_")[0]+"_"+a.split("_")[1]
    if head in dic :
        return dic[head]
