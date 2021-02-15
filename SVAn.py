#!/usr/bin/env python3
import tempfile
from itertools import islice
import getopt
import sys
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
import multiprocessing as mp
import csv
import numpy as np
#import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt
#import pandas as pd
#from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
#import re
import pandas as pd
#import collections
from collections import defaultdict
import os
from operator import itemgetter
import statistics
#from matplotlib import rc
#import matplotlib.pyplot as plt
csv.field_size_limit(100000000)
import random
import time

def main():
    #print(sys.argv[1:])
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:b:t:m:h:s:r:c:", ["vcf=", "blat=", "trf=", "mobile=", "mobile_hmm=", "threshold=", "ref_genome=", "cpu="])
    except getopt.GetoptError:
    # print help information and exit:
    #print ('error')  # will print something like "option -a not recognized"
        sys.exit(2)

# Default parameters
#print(opts)
    global vcf
    global blat_path
    global trf
    global mobile
    global threshold
    global hmm
    global ref_genome
    global cpu
    global low
    global liste_info
    global output_name


    for opt, arg in opts:
        #print(opt, arg)
        if opt in ('-v', "--vcf"):
            vcf = arg
        #print(i)
        elif opt in ('-b', "--blat" ):
            blat_path = arg
        elif opt in ('-m', "--mobile"):
            mobile = arg
        #print(r)
        elif opt in ('-h', "--mobile_hmm"):
            hmm = arg
        #print(r)
        elif opt in ('-t', "--trf"):
            trf = arg
        elif opt in ('-s', "--threshold"):
            threshold = float(arg)
            low = 1-threshold
        elif opt in ('-r', "--ref_genome"):
            ref_genome = arg
        elif opt in ('-c', "--cpu"):
            cpu = int(arg)
        #print(r)
        else:
            print("error")
            assert False, "unhandled option"
    if vcf == "" or blat_path == "" or mobile == "" or trf == "" or hmm == "" or ref_genome == "":
        print("Mandatory file missing")
    else:
    #print("oui")
    #job(blat_path,vcf,trf,mobile,threshold,hmm,ref_genome)
        unassigned=[]
        output_name = vcf.split(".")[0].split("/")[-1]+"_annotated.vcf"
        output_name_final = vcf.split(".")[0].split("/")[-1]+"_final_annotated.vcf"
        list_ids = get_list(vcf, output_name)
        list_list_ids = splitter(list_ids, cpu)
        t1 = time.time()
            # splitting and controlling process count (Faster)
        get_info(list_list_ids, cpu)
        #print(time.time() - t1)
        #print(len(liste_info))
        #print("DISP check")
        unassigned=create_fasta_disp( output_name)
        align_disp_test(blat_path, ref_genome)
        unassigned_update=parse_multi_blat(threshold, "Tmp_Disp.psl", unassigned)
        #print(unassigned_update)
        write_final_vcf(output_name,unassigned_update,output_name_final)
def splitter(list_ids, NUM_WORKERS):
    split = int((len(list_ids)/NUM_WORKERS))+1
    list_list_ids = []


# 	or i in range(0, NUM_WORKERS):
# 		liste_f.append(split)
# #list_list_ids=splitter(list_ids,NUM_WORKERS)
# 	Output = [list(islice(iter(list_ids), elem))for elem in liste_f]
# 	print (Output)
# 	return Outputf
    for i in np.array_split(list_ids, NUM_WORKERS):  # numpy required
        list_list_ids.append(list(i))
        #print(list_list_ids)
        return list_list_ids


def get_list(vcf_file,output_name):
    list_ids = []
    csv_file = csv.reader(open(vcf_file, "r"), delimiter="\t")
    output = csv.writer(open(output_name, "w"), delimiter="\t")
    i=0
    for elt in csv_file:
        if "#" not in elt[0] and "@" not in elt[0]:
            list_ids.append(elt)
        else :
            output.writerow(elt)
    if i == 0:
        output.writerow(["##INFO=<ID=Annotation_Type,Number=1,Type=String,Description=\"Refined insertion type\">"])
        i += 1
    else :
        output.writerow(elt)
    return list_ids

def get_pgm():
    return blat_path, vcf, trf, mobile, threshold, hmm, ref_genome, output

def mp_process(sub_list):
#blat_path, vcf, trf, mobile, threshold, hmm, ref_genome=get_pgm()
    return job(sub_list, blat_path, vcf, trf, mobile, threshold, hmm, ref_genome, output_name)


def get_info(ids, poolsize):
    ''' A function to initialize the pool of workers and call the wrapper function mp_process. 
    The input 'ids' is a list of  sublists'''
    pool = mp.Pool(poolsize)  # multiprocessing package required
    s = pool.map(mp_process, ids)
    #print('Active children count: %d' % len(mp.active_children()))
    pool.close()
    pool.join()
    return 'OK'


def job(sub_list, blat_path, vcf, trf, mobile, threshold, hmm, ref_genome, output_name):
#vcf_reader = csv.reader(open(vcf, "r"), delimiter="\t")
#print(sub_list)
    output = csv.writer(open(output_name, "a"), delimiter="\t")
    low = 1-threshold
    a=random.randint(0,5000)
    #print(os.system("pwd"))
    cmd_dir = "mkdir tmp_" + str(a) #+ " && cd tmp_"+str(a)
    tmp_dir = "tmp_" + str(a)+"/"
    #print(tmp_dir)
    #print(cmd_dir)
    os.system(cmd_dir)
    nam_tmp_fa = "Tmp_SV.fa"
    i = 0
    for elt in sub_list:
        temp_in = list(elt)
    #print(temp_in)
        best = False
        ref, alt = get_seq(elt)
        #print(ref,alt)
        if len(ref) != len(alt):
            types = identification_type(ref, alt,elt)
            
        if types == "Gain":
            sequence = alt
        else:
            sequence = ref
        write_tmp_seq(sequence,tmp_dir)
        run_command(nam_tmp_fa, mobile, hmm, tmp_dir)
        best,percentage_ME = parse_me(sequence, threshold, tmp_dir)
        if best:
            temp_in[7] = temp_in[7]+";Annotation_Type="+types + "_"+"Mobile_element"
        #writer new annot
        else :
            align_trf(nam_tmp_fa, trf, tmp_dir)
            best,percentage_TRF = parse_trf(sequence, threshold, tmp_dir)
            if best:
                temp_in[7] = temp_in[7]+";Annotation_Type="+types + "_"+"Tandem_Repeat"
            else :
                align_td(sequence, elt[0], int(elt[1]), ref_genome, blat_path, tmp_dir)
                best,percentage_TD = parse_blat(threshold,  "Tmp_TD.psl", sequence, tmp_dir)
                if best:
                    temp_in[7] = temp_in[7]+";Annotation_Type=" + types + "_"+"Tandem_Duplication"
                else:
                    if percentage_ME>=1-threshold or percentage_TD>=1-threshold or percentage_TRF>=1-threshold :
                        temp_in[7] = temp_in[7]+";Annotation_Type="+types+"_NoNew"
                    else :
                        temp_in[7] = temp_in[7]+";Annotation_Type="+types+"_TMP"
        output.writerow(temp_in)
    

    #os.system("rm -rf  "+cmd_dir)

def create_fasta_disp(output_name) :
    reader=csv.reader(open(output_name,"r"),delimiter="\t")
    temp_fasta=csv.writer(open("Tmp_SV.fa",'w'),delimiter="\n")
    unassigned=defaultdict(list)
    for elt in reader :
        if "#" not in elt[0]:
            if  elt[7].split(";")[-1].split("_")[-1]=="TMP" or elt[7].split(";")[-1].split("_")[-1]=="NoNew":
                print ("Looking",elt[7].split(";")[-1] )
                ref, alt = get_seq(elt)
                if len(ref) != len(alt):
                    types = identification_type(ref, alt,elt)
                #print(ref,alt,types)
                if types == "Gain":
                    sequence = alt
                else:
                    sequence = ref
                #print(elt)
                liste=[">"+str(elt[0])+"_"+str(elt[1]),sequence]
                temp_fasta.writerow(liste)
                key=str(elt[0])+"_"+str(elt[1])
                unassigned[key].append(elt)
    return unassigned
def align_disp_test(blat_path, ref_genome):
    try:
        tmp_fa = csv.reader(open("Tmp_SV.fa", "r"), delimiter=" ") 
        cmd = blat_path+" "+ref_genome +" -q=dna -t=dna "+"Tmp_SV.fa "+"Tmp_Disp.psl -noHead" #-ooc="+tmp_dir+"Tmp_SV.fa.11.ooc"  # -ooc=11.ooc not working
        os.system(cmd)
    except IOError:
        print("file does not exist")
        return False
    
def parse_multi_blat(threshold, files, unassigned):
    unassigned_update=defaultdict()
    try :
        inpt = csv.reader(open(files, "r"), delimiter="\t")
        dictionnary_match=defaultdict(float)
        for elt in inpt:
            id_seq=elt[9]
            if id_seq not in dictionnary_match :
                match = int(elt[0])
                
                #print(elt)
            #size_align_ref = check_size_align(int(elt[16]), int(elt[15]))
            #size_inser = int(elt[9].split("_")[-1])
                percentage = float(match/int(elt[10]))
                dictionnary_match[id_seq]=percentage
                #if percentage >= threshold:
                #    return True,percentage
            #return False,percentage
        #print(dictionnary_match)
        #print(unassigned)
        for elt in unassigned :
            key=elt
            value=unassigned[elt][0]
            #print("TEST",unassigned[elt][0][7].split(';'))
            if elt not in dictionnary_match :
                if unassigned[elt][0][7].split(';')[-1].split("_")[-1]!="NoNew" :
                    types=value[7].split(';')[-1].split("=")[-1].split("_")[0]
                    previous=value[7].split(';')[-1]
                    res = [item.replace(previous, "Annotation_Type="+types+"_"+"Novel_sequence") for item in value]
                    unassigned_update[key]=res
                else :
                    types=value[7].split(';')[-1].split("=")[-1].split("_")[0]
                    previous=value[7].split(';')[-1]
                    res = [item.replace(previous, "Annotation_Type="+types+"_"+"Unassigned") for item in value]
                    unassigned_update[key]=res
            else :
                if dictionnary_match[elt]>=threshold :
                    types=value[7].split(';')[-1].split("=")[-1].split("_")[0]
                    previous=value[7].split(';')[-1]
                    res = [item.replace(previous, "Annotation_Type="+types+"_"+"Dispersed_duplication") for item in value]
                    unassigned_update[key]=res
                elif dictionnary_match[elt]<=1-threshold :
                    types=value[7].split(';')[-1].split("=")[-1].split("_")[0]
                    previous=value[7].split(';')[-1]
                    res = [item.replace(previous, "Annotation_Type="+types+"_"+"Novel_sequence") for item in value]
                    unassigned_update[key]=res
                else :
                    types=value[7].split(';')[-1].split("=")[-1].split("_")[0]
                    previous=value[7].split(';')[-1]
                    res = [item.replace(previous, "Annotation_Type="+types+"_"+"Unassigned") for item in value]
                    unassigned_update[key]=res
        return unassigned_update
        #print(unassigned)
    except IOError:
        return False
def write_final_vcf(vcf_modified,unassigned_update,output_name_final) :
    read_vcf=csv.reader(open(vcf_modified,"r"),delimiter="\t")
    write_vcf=csv.writer(open(output_name_final,"w"),delimiter="\t")
    for elt in read_vcf :
        if "#" in elt[0] :
            write_vcf.writerow(elt)
        else :
            key=str(elt[0])+"_"+str(elt[1])
            if key not in unassigned_update :
                write_vcf.writerow(elt)
            else :
                write_vcf.writerow(unassigned_update[key])
def check_size_align(a, b):
    size = 0
    if a > b:
        size = a-b
    else:
        size = b-a
    return size


def is_valid(inser):
    allowed = "ATCGatcg"
    if all(c in allowed for c in inser):
        return True
    else:
        return False


def get_seq(elts):
    ref = elts[3]
    alt = ""
    if is_valid(elts[4]):
        alt = elts[4]
    else:
        parser = elts[7].split(';')
    for elt in parser:
        if "SEQ" in elt.split("=")[0] and is_valid(elt.split("=")[1].upper()):
            alt = elt.split("=")[1].upper()
    if is_valid(alt):
        return (ref, alt)
    if not is_valid(alt):
        return "NA", "NA"


def identification_type(ref, alt,elt):
    if "DEL" in elt[4] :
        return "Loss"
    elif "INS" in elt[4] :
        return "Gain"
    if len(ref) < len(alt):
        return "Gain"
    else:
        return "Loss"

# Mobile Element


def write_tmp_seq(seq, tmp_dir):
    write_tmp = csv.writer(open(tmp_dir+"Tmp_SV.fa", "w"), delimiter="\n")
    write_tmp.writerow([">Tmp_SV_a", seq])


def run_command(tmp_file, mobile, hmm, tmp_dir):
    command_me = "perl " + mobile + "  --fastafile " + tmp_dir+tmp_file + " --hmmfile " + hmm + " --dfam_outfile " + tmp_dir+"Temp_ME"
    #print(command_me)
    os.system(command_me)


def parse_me(sequence, threshold, tmp_dir):
    try:
        ME = csv.reader(open(tmp_dir+"Temp_ME", "r"), delimiter=" ")
        for elt in ME:
            if "#" not in elt[0]:
                clean_list = [x for x in elt if x != '']
                size_aligned = check_size_align(int(clean_list[12]), int(clean_list[11]))+1
                size = len(sequence)
                percentage = float(size_aligned)/float(size)
                if percentage >= threshold:
                    return True,percentage
                return False,percentage    
        return False,0
    except IOError:
        return False,0

    
    # Tandem Repeat

def convert_readable(txtfile, datfile, tmp_dir):
	with open(tmp_dir+txtfile, 'a') as txt:
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
	cmd_rm="rm "+datfile
	os.system(cmd_rm)

def align_trf(fa_file, trf, tmp_dir):
    txt_file = "readable_format_TRF.txt"
    cmd2 = trf + " " + tmp_dir+fa_file+" 2 5 7 80 10 50 2000 -l 6 -d -h "
    os.system(cmd2)
    convert_readable(txt_file, "Tmp_SV.fa.2.5.7.80.10.50.2000.dat", tmp_dir)
#cmd1 = "rm tmp_seq_1.*"
#os.system(cmd1)


def parse_trf(sequence, threshold, tmp_dir):
    try :
        TRF = csv.reader(open(tmp_dir+"readable_format_TRF.txt", "r"), delimiter=" ")
        for elt in TRF:
            #print(elt)
            size_trf = check_size_align(int(elt[2]), int(elt[1]))+1
            size = len(sequence)
            percentage = float(size_trf)/float(size)
            #print (percentage)
            if percentage >= threshold:
                return True,percentage
            return False,percentage
        return False,0
    except IOError:
        return False,0


# Tandem Duplication
def extract_bkpt(ref_genome, pos, chrom, sequence, tmp_dir):
    reader = SeqIO.parse(ref_genome, "fasta")
    writer_bkpt = csv.writer(open(tmp_dir+"Tmp_bkpt.fa", "w"), delimiter="\n")
    start = ""
    end = ""
    for element in reader:
        chromosome = element.description
        ref_seq = str(element.seq)
        if chromosome == chrom:
            start = ref_seq[pos-len(sequence):pos].upper()
            end = ref_seq[pos:pos+len(sequence)].upper()
            writer_bkpt.writerow([">Left", start, ">Right", end])
#if start == "" or end == "":
#	print("Error reference genome")


def align_td(sequence, chrom, pos, ref_genome, blat_path, tmp_dir):
    extract_bkpt(ref_genome, pos, chrom, sequence, tmp_dir)
    #path_global = str(os.system("pwd"))+"/"
    cmd2 = blat_path+" "+tmp_dir+"Tmp_bkpt.fa "+tmp_dir +"Tmp_SV.fa "+"-q=dna -t=dna "+tmp_dir+"Tmp_TD.psl -noHead" #-ooc="+tmp_dir+"Tmp_SV.fa.11.ooc"  # not working
    os.system(cmd2)


def parse_blat(threshold, files, sequence,tmp_dir):
    try :
        inpt = csv.reader(open(tmp_dir+files, "r"), delimiter="\t")
        for elt in inpt:
            print("RESULT",elt)
            match = int(elt[0])
        #size_align_ref = check_size_align(int(elt[16]), int(elt[15]))
        #size_inser = int(elt[9].split("_")[-1])
            percentage = float(match/len(sequence))
            if percentage >= threshold:
                return True,percentage
            return False,percentage
        return False,0
    except IOError:
        return False,0
# Dispersed duplicationls


def align_disp(Tmp_fa, ref_genome, blat_path, threshold, sequence, tmp_dir):
    #path_global = os.system("pwd")+"/"
    cmd = blat_path+" "+tmp_dir+"Tmp_SV.fa " + ref_genome + \
    " -q=dna -t=dna "+tmp_dir+"Tmp_Disp.psl -noHead >&" # -ooc="+tmp_dir+"Tmp_SV.fa.11.ooc"  # -ooc=11.ooc not working
    os.system(cmd)
    best = parse_blat(threshold, "Tmp_Disp.psl", sequence, tmp_dir)
    os.system("cat "+tmp_dir+"Tmp_Disp.psl >> "+tmp_dir+"Tmp_all_Disp.psl")
    if best:
        return True
    return False


def check_novo(tmp_me, tmp_trf, tmp_td, tmp_disp, low, sequence, tmp_dir):
    best_me = False
    best_trf = False
    best_td = False
    best_disp = False
    best_me = parse_me(sequence, low, tmp_dir)
    best_trf = parse_trf(sequence, low, tmp_dir)
    best_td = parse_blat(low, "Tmp_TD.psl", sequence, tmp_dir)
    best_disp = parse_blat(low, "Tmp_Disp.psl", sequence, tmp_dir)
    if not best_me and not best_trf and not best_td and not best_disp:
        return True
    return False


if __name__ == "__main__":
    vcf = ""
    blat_path = ""
    trf = ""
    mobile = ""
    threshold = 0.8
    hmm = ""
    ref_genome = ""
    cpu = 1
    low=.2
    liste_info=[]
    output_name=""
    main()
