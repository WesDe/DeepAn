from collections import defaultdict
import sys
import csv
import re
import pandas as pd
import pyupset as pyu
csv.field_size_limit(100000000)
bed_1 = sys.argv[1]


def make_dictionary(bed, dic):
    bed_file=csv.reader(open(bed, 'r'), delimiter='\t')
    for elt in bed_file :
        if '#' not in elt[0] and '@' not in elt[0] :
            dic.setdefault(re.sub('chr','',elt[0]),[]).append((int(elt[1]),int(elt[2])))
    print(bed)
    return dic


def parsing(bed, dic,types):
    bed_file = csv.reader(open(bed, 'r'), delimiter='\t')
    for elt in bed_file:
        if '#' not in elt[0] and '@' not in elt[0]:
            chrom = re.sub('chr', '', elt[1])
            start = int(elt[2])
            end = int(elt[3])
            dic=check_in(dic,chrom,start,end,types)
    return dic

def check_elt_in_dic(dic, chrom, pos,types):
    if chrom in dic:
        for liste in dic[chrom]:
            if pos >= int(liste[0]) and pos <= int(liste[1]) :
                #outp_file.writerow(element)
                return 1
        return 0

def incremente_dic(couple,dic):
    if couple not in dic_resume :
        dic[couple]=1
    else  :
        dic[couple]+=1
    return dic


def check_in(dic,chrom,pos_start,pos_end,types) :
    for elt in dic :
        if chrom==elt[0] and elt[1]>=pos_start and elt[1]<=pos_end :
            dic[elt].append(types)
            #print(types,chrom,pos_start,pos_end,elt,dic[elt])
            #print(types)
        elif chrom==elt[1] and elt[1]>=pos_end :
            #print(types, chrom, pos_start, pos_end, elt, dic[elt])
            break
    
    return dic




type_1 = bed_1.split("/")[-1].split("_")[0]
dic_resume={}


found_1=0

tot=0
truth_file = csv.reader(open(sys.argv[2], 'r'), delimiter='\t')
outp_file = csv.writer(open(sys.argv[3], 'w'), delimiter=',')
#header = [type_1, type_2,type_3,type_4,type_5,type_6,type_7,type_8,type_9,type_10,"Other"]
#outp_file.writerow(header)
dic=defaultdict(list)
for element in truth_file:
    other=0
    sum_found=0
    #print(element)
    if '#' not in element[0] and '@' not in element[0]:
        tot += 1
        pos = int(element[1])
        chrom = re.sub('chr', '', element[0])
        dic[(chrom,pos)]

print("parsing 1")
dic = parsing(bed_1, dic, type_1)

for elt in dic :
    liste=[elt]
    for el in dic[elt] :
        liste.append(el)
    outp_file.writerow(liste)
#print(dic)
