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


def is_valid(inser):
    allowed = "ATCGatcg"
    if all(c in allowed for c in inser):
        return True
    else:
        return False
def check_size_align(a,b) :
    size=0
    if a > b :
        size =a-b
    else :
        size =b-a
    return size
def check_presence(elt,dic):
    if elt in dic :
        return 1
    else :
        return 0
def get_occu(liste,threshold):
    count=0
    for i in liste :
        if i>= threshold :
            count+=1
    return count
def get_dic_all(inpt) :
    dic_all=[]
    inpt_seq=csv.reader(open(inpt,"r"), delimiter="\t")
    for element in inpt_seq :
        svseq=""
        if "#" not in element[0] and "@" not in element[0] :
            if is_valid(element[4]) :
                svseq=element[4]
            else :
                parser = element[7].split(';')
                for elt in parser :
                    #print(elt)
                    if "SEQ" in elt.split("=")[0] and is_valid(elt.split("=")[1].upper()) :
                        svseq=elt.split("=")[1].upper()
            if svseq!="" :
                chrom = re.sub("chr", "", element[0])
                position=int(element[1])
                headers=chrom+"_"+str(position)+"_"+str(len(svseq))
                if headers not in dic_all and is_valid(svseq) and len(svseq)>=50 :
                    dic_all.append(headers)
    return dic_all
def get_seq(element):
    svseq=""
    if is_valid(element[4]) :
        svseq=element[4]
    else :
        parser = element[7].split(';')
        for elt in parser :
            #print(elt)
            if "SEQ" in elt.split("=")[0] and is_valid(elt.split("=")[1].upper()) :
                svseq=elt.split("=")[1].upper()
    return svseq

def get_list_occu(info_file,threshold,dic_all): 
    inpt=csv.reader(open(info_file,"r"),delimiter="\t")
    count=0
    liste_all=defaultdict()
    dic=defaultdict(list)
    count_td=0
    next(inpt)
    next(inpt)
    next(inpt)
    next(inpt)
    next(inpt)
    for elt in inpt :
        if len(elt[9].split("_"))<=4 :
            match=int(elt[0])
            size_align=check_size_align(int(elt[12]),int(elt[11]))
            size_align_ref=check_size_align(int(elt[16]),int(elt[15]))
            size_inser=int(elt[9].split("_")[-1])
            chrom_inser=re.sub("chr", "",elt[9].split("_")[0])
            inser_pos=int(elt[9].split("_")[2])
            chrom_ref=re.sub("chr", "",elt[13])
            ref_pos=int(elt[15])
            heads=chrom_inser+"_"+str(inser_pos)+"_"+str(size_inser)
            #print(size_align,size_inser,chrom_inser,inser_pos,chrom_ref,ref_pos)
            percentage=float(size_align/size_inser)
            percentage_match=float(match/size_inser)
            dic[heads].append(percentage)
    for head in dic :
        if head in dic_all :
            counts=get_occu(dic[head],threshold)
            liste_all[head]=counts
        else :
            liste_all[head]=0
    return liste_all

def check_duplication(info_file,threshold,dic_occurrence): 
    inpt=csv.reader(open(info_file,"r"),delimiter="\t")
    next(inpt)
    next(inpt)
    next(inpt)
    next(inpt)
    next(inpt)
    count=0
    seen_validate=set()
    dic=defaultdict(list)
    dic_last=defaultdict(float)
    black_list=[]
    dic_segm=defaultdict(float)
    dic_td=defaultdict(float)
    dic_not_valid=defaultdict(float)
    count_td=0
    liste_segm_dup=[]
    count_segm=0
    dic_tes=defaultdict(list)
    for elt in inpt :
        if len(elt[9].split("_"))<=4 :
            match=int(elt[0])
            size_align=check_size_align(int(elt[12]),int(elt[11]))
            size_align_ref=check_size_align(int(elt[16]),int(elt[15]))
            size_inser=int(elt[9].split("_")[-1])
            chrom_inser=re.sub("chr", "",elt[9].split("_")[0])
            inser_pos=int(elt[9].split("_")[2])
            chrom_ref=re.sub("chr", "",elt[13])
            ref_pos=int(elt[15])
            head=chrom_inser+"_"+str(inser_pos)+"_"+str(size_inser)
            #print(head)
            #print(size_align,size_istsnser,chrom_inser,inser_pos,chrom_ref,ref_pos)
            #percentage=float(size_align/size_inser)
            percentage=float(match/size_inser)
            #print(percentage)
            if chrom_ref==chrom_inser and ref_pos>=inser_pos-size_inser and ref_pos<=inser_pos+size_inser :
                #print(elt)
                black_list.append(head)
            dic[head].append(percentage)
            #if head in dic_occurrence and size_inser>=1000 and head not in liste_segm_dup and dic_occurrence[head]>=1 and dic_occurrence[head]<=50 :
            #    liste_segm_dup.append(head)
            #    count_segm+=1
            iden=elt[9]
            size_align=elt[18].split(",")
            size_inse=elt[19].split(",")
            for i in range (0,len(size_align)-1):
                couple=(int(size_inse[i]),int(size_inse[i])+int(size_align[i]))
                dic_tes[iden].append(couple)
    for elt in dic :
        #print(elt)
        best=max(dic[elt])
        if elt not in black_list and elt not in liste_segm_dup and dic_occurrence[elt]<=50 and best >=threshold :
                dic_last[elt]=best
                seen_validate.add(elt)
                count+=1
        if elt in black_list and best >=threshold  and elt not in liste_segm_dup and best >=threshold :
            dic_td[elt]=best
            #seen_validate.add(elt)
            count_td+=1
        if elt in liste_segm_dup and best >=threshold :
            dic_segm[elt]=best
            
        if elt not in liste_segm_dup  and best >=threshold :
            dic_not_valid[elt]=best
                #print(best,elt)
            #print(best)
        #break
        #break
    #for e in black_list :
    #    if e in seen_validate :
    #        seen_validate.remove(e)
    #        count-=1
    #print(count,count_td,count_segm,len(dic_not_valid))
    return(dic_last,dic_td,dic_segm,dic_not_valid,dic_tes)


def get_dic_composite(dic_tes):
    dic_size_align=defaultdict()
    for elt in dic_tes :
        #print("before",dic_tes[elt])
        work_list=sorted(dic_tes[elt], key=lambda tup:tup[0]) 
        #print("after",work_list)
        size=0
        st=0
        for elts in work_list :
            if st==0 :
                begin=elts[0]
                end=elts[1]
                st+=1
                continue
            else :
                if elts[1]>end and elts[0]<=end :
                    end=elts[1]
                elif elts[1]<=end or (elts[0]==begin and elts[1]==end):
                    continue
                elif elts[0]>end :
                    size+=end-begin
                    begin=elts[0]
                    end=elts[1]
                else :
                    print("probleme",elts,begin,end)
        size+=end-begin
        dic_size_align[elt]=size
    return dic_size_align

def make_db(info_file,threshold):
    TRF=csv.reader(open(info_file,"r"),delimiter=" ")
    seen_validate=set()
    dic=defaultdict(list)
    #next(TRF)
    count=0
    dic_last=defaultdict(float)
    dic_seed=defaultdict(list)
    dic_seed_unique=defaultdict(list)
    for elt in TRF :
        #print(elt)
        if len(elt[0].split("_"))<=4 :
            #print(elt)
            size_trf =check_size_align(int(elt[2]),int(elt[1]))+1
            parser=elt[0].split("_")
            chromosome=re.sub("chr", "",parser[0])
            position=int(parser[-2])
            size=int(parser[-1])
            percentage=float(size_trf)/float(size)
            head=chromosome+"_"+str(position)+"_"+str(size)
            #print (percentage)
            if percentage >=threshold :
                dic[head].append(percentage)
                compound=(float(elt[3]),float(elt[4]))
                dic_seed[head].append(compound)
            
    for elt in dic :
        #print(dic[elt],elt)
        best=max(dic[elt])
        if best >=threshold :
            dic_last[elt]=best
            seen_validate.add(elt)
            count+=1
            #print(best)
    for elt in dic_seed :
        if len(dic_seed[elt])!=1 :
            lower_seed=dic_seed[elt][0][0]
            max_rep=dic_seed[elt][0][1]
            for e in dic_seed[elt]:
                if e[0]< lower_seed and e[1]>max_rep :
                    lower_seed=e[0]
                    max_rep=e[1]
            dic_seed_unique[elt]=(lower_seed,max_rep)
        else :
            dic_seed_unique[elt]=(dic_seed[elt][0][0],dic_seed[elt][0][1])
    #print(count)
    return(dic_last,dic_seed_unique)

    
def make_db_me(info_file,threshold,delimiters):
    seen_validate=set()
    dic=defaultdict(list)
    seen=[]
    #next(info_file)
    ME=csv.reader(open(info_file,"r"),delimiter=delimiters)
    count=0
    liste_found=[]
    dic_last=defaultdict(float)
    for elt in ME :
        if "#" not in elt[0] :
            clean_list = [x for x in elt if x != '']
            parser=clean_list[2].split("_")
            if len(parser)<=4 :
                #print(clean_list)
                size_aligned =check_size_align(int(clean_list[12]),int(clean_list[11]))+1
                parser=clean_list[2].split("_")
                chromosome=re.sub("chr", "",parser[0])
                position=int(parser[-2])
                size=int(parser[-1])
                #print(clean_list)
                #print(size_aligned,size)
                percentage=float(size_aligned)/float(size)
                #print(size,size_aligned,percentage)
                #print(clean_list)
                head=chromosome+"_"+str(position)+"_"+str(size)
                if percentage >=threshold :
                    dic[head].append(percentage)

    for elt in dic :
        #print(dic[elt],elt)
        best=max(dic[elt])
        if best >=threshold :
            dic_last[elt]=best
            seen_validate.add(elt)
            count+=1
            #print(dic[elt],elt,best)
    #print(count)
    return(dic_last)

def get_list(dic_all,dic_td,dic_tr,dic_me,dic_dup,segm_dup,not_valid,threshold,compo):
    count_novo=0
    count_novo_2=0
    count_me=0
    count_tr=0
    count_td=0
    count_dup=0
    count_segm=0
    count_order_me=0
    count_order_tr=0
    count_order_td=0
    count_order_dup=0
    count_order_segm=0
    liste_novo=[]
    liste_td=[]
    liste_me=[]
    liste_tr=[]
    liste_dup=[]
    liste_segm=[]
    liste_td_o=[]
    liste_me_o=[]
    liste_tr_o=[]
    liste_dup_o=[]
    liste_segm_o=[]
    #print(len(segm_dup))
    for elt in dic_all :
        if elt not in dic_td and elt not in segm_dup and elt not in dic_tr and elt not in dic_me and elt not in dic_dup and elt not in segm_dup and elt not in not_valid :
            #print(elt)
            count_novo+=1
            liste_novo.append(elt)
            #if elt not in compo :
            #    count_novo+=1
            #    liste_novo.append(elt)
            #elif compo[elt]/int(elt.split("_")[-1])<threshold :
            #    count_novo_2+=1
            #    liste_novo.append(elt)
            #else :
            #    print(elt,compo[elt])
        #if threshold!=1 :
        #    if best_td[1]<1-threshold and best_tr[1]<1-threshold and best_me[1]<1-threshold and best_dup[1]<1-threshold :
        #        count_novo+=1
        if elt in dic_me :
            count_order_me+=1
            liste_me_o.append(elt)
            continue
        if elt in dic_tr :
            count_order_tr+=1
            liste_tr_o.append(elt)
            continue
        if elt in segm_dup :
            count_order_segm+=1
            liste_segm_o.append(elt)
            continue
        if elt in dic_td :
            count_order_td+=1
            liste_td_o.append(elt)
            continue
        if elt in dic_dup :
            count_order_dup+=1
            liste_dup_o.append(elt)
            continue
        #break
    total_order=count_novo+count_order_td+count_order_tr+count_order_me+count_order_dup
    print ("Threshold set at :", threshold,"novo : ", count_novo, "TD : ", count_order_td, "TR : ", count_order_tr, "ME : ", count_order_me, "DUP : ", count_order_dup)#, "SEGM : ",count_order_segm)
    #print(" TOTAL order : ",total_order)
    return liste_novo,liste_td_o,liste_me_o,liste_tr_o,liste_dup_o,liste_segm_o
def get_complex(dic,td,tr,me,dup,novo,segm,threshold):
    count=0
    liste_complex=[]
    for head in dic :
        if head not in novo and head not in td and head not in me and head not in tr and head not in dup and head not in segm:
            count+=1
            liste_complex.append(head)
    print("Unassigned at : ", threshold, " :", count)
    return liste_complex

def count_unassigned(complex_8_o,liste,types) :
    count=0
    for elt in complex_8_o :
        if elt in liste :
            count+=1
    print("Unassigned : potential "+types+" :", count)


#################################################
def write_new_vcf(vcf,dic,otp) :
    new=csv.writer(open(otp,"w"),delimiter="\t")
    lec=csv.reader(open(vcf,"r"),delimiter="\t")
    increm=0
    for elt in lec :
        temp_in=elt
        if "#" in elt[0] or "@" in elt[0]:
            new.writerow(elt)
        else :
            if increm==0 :
                header_one=["##INFO=<ID=Repeat_loc,Number=1,Type=String,Description=\"Location of insertion in potential repeat\">"]
                header_two=["##INFO=<ID=MH_size,Number=1,Type=String,Description=\"Size of junctional homology\">"]
                header_three=["##INFO=<ID=SVTYPE_Refined,Number=1,Type=String,Description=\"Refined insertion type\">"]
                new.writerow(header_one)
                new.writerow(header_two)
                new.writerow(header_three)

                increm+=1
            else :
                head=re.sub("chr", "",elt[0])+"_"+elt[1]+"_"+str(len(get_seq(elt)))
                if head in dic :
                    temp_in[7]=temp_in[7]+"Repeat_loc="+str(dic[head][3])+";MH_size="+str(dic[head][4])+";SVTYPE_Refined="+str(dic[head][5])
                    #temp_in[7]=temp_in[7]+";MH_size="+dic[head][3]+";SVTYPE_Refined="+dic[head][4]
                    new.writerow(temp_in)
                else :
                    temp_in[7]=temp_in[7]+"Repeat_loc=NA;MH_size=NA;SVTYPE_Refined=NA"
                    #temp_in[7]=temp_in[7]+";MH_size="+dic[head][3]+";SVTYPE_Refined="+dic[head][4]
                    new.writerow(temp_in)