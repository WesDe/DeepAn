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
from operator import itemgetter
import statistics
from matplotlib import rc
import matplotlib.pyplot as plt
csv.field_size_limit(100000000)
import Function, Microhomology, Location
import sys
import getopt
def main():
	print(sys.argv[1:])
	try:
		opts, args = getopt.getopt(sys.argv[1:], "v:a:t:m:u:l:o:h:", ["vcf=","alignment=","trf=","mobile=","micro=","large=","output=","threshold="])
	except getopt.GetoptError:
	# print help information and exit:
	#print ('error')  # will print something like "option -a not recognized"
		sys.exit(2)

    # Default parameters
	#print(opts)
	vcf=""
	alignment_wg=""
	trf=""
	mobile=""
	micro=""
	large=""
	output=""
	threshold=0.8
	low=0.2
	for opt, arg in opts:
		print(opt,arg)
		if opt in ('-v',"--vcf"):
			vcf= arg
			#print(i)
		elif opt in ('-a',"--alignment"):
			alignment_wg= arg
		elif opt in ('-m',"--mobile"):
			mobile=arg
			#print(r)
		elif opt in ('-t',"--trf"):
			trf=arg
		elif opt in ('-u',"--micro"):
			micro= arg
		elif opt in ('-l',"--large"):
			large= arg
		elif opt in ('-o',"--output"):
			output= arg
		elif opt in ('-h',"--threshold"):
			threshold= float(arg)
			low=1-threshold

			#print(r)
		else:
	            assert False, "unhandled option"
	if vcf=="" or alignment_wg=="" or mobile=="" or trf=="" or micro=="" or large=="" or output=="":
		print("Mandatory file missing")
	else :
		dic_all=Function.get_dic_all(vcf)
		dic_occu=Function.get_list_occu(alignment_wg,threshold,dic_all)
		dic_occu_less=Function.get_list_occu(alignment_wg,low,dic_all)
		dup,TD,Segm_Dup,not_valid,dic_ali_comp=Function.check_duplication(alignment_wg,threshold,dic_occu)
		dup_b,TD_b,Segm_Dup_b,not_valid_b,dic_ali_comp_b=Function.check_duplication(alignment_wg,low,dic_occu_less)
		compo=Function.get_dic_composite(dic_ali_comp)
		compo_b=Function.get_dic_composite(dic_ali_comp_b)

		TRF,TRF_seed=Function.make_db(trf,threshold)
		TRF_b,TRF_b_seed=Function.make_db(trf,low)

		ME=Function.make_db_me(mobile,0.80," ")
		ME_b=Function.make_db_me(mobile,0.20," ")


		novo,td_o,me_o,tr_o,dup_o,Segm_o=Function.get_list(dic_all,TD,TRF,ME,dup,Segm_Dup,not_valid,threshold,compo)
		novo_b,td_o_b,me_o_b,tr_o_b,dup_o_b,Segm_o_b=Function.get_list(dic_all,TD_b,TRF_b,ME_b,dup_b,Segm_Dup_b,not_valid_b,low,compo_b)

		complex_o=Function.get_complex(dic_all,td_o,tr_o,me_o,dup_o,novo_b,Segm_o,threshold)
		complex_o_b=Function.get_complex(dic_all,td_o_b,tr_o_b,me_o_b,dup_o_b,novo,Segm_o_b,low)


		d_test=defaultdict(list)   
		d_test=Microhomology.get_mh_blat(large)
		clean_dictionary=defaultdict(list)
		clean_dictionary=Microhomology.clean_dic(d_test)

		dic_size_mh={}
		dic_size_mh=Microhomology.concatenate_mh(clean_dictionary)
		dic_small_mh=Microhomology.get_small_mh(micro)


		dic_all_info=defaultdict(list)

		dic_all_info=Microhomology.fill_type(dic_small_mh,dic_size_mh,dic_all,dic_all_info)


		dic_all_info=Microhomology.Assign_type(dic_all_info,me_o,"MobileElement")
		dic_all_info=Microhomology.Assign_type(dic_all_info,tr_o,"TandemRepeat")
		dic_all_info=Microhomology.Assign_type(dic_all_info,td_o,"TandemDuplication")
		dic_all_info=Microhomology.Assign_type(dic_all_info,dup_o,"DispersedDuplication")
		dic_all_info=Microhomology.Assign_type(dic_all_info,novo_b,"NovelSequence")
		dic_all_info=Microhomology.Assign_type(dic_all_info,Segm_o,"SegmentalDuplication")
		dic_all_info=Microhomology.Assign_type(dic_all_info,complex_o,"Unassigned")
		to_remove=[]
		for elt in dic_all_info :
			if len(dic_all_info[elt])<2 :
				to_remove.append(elt)
		for a in to_remove :
			del dic_all_info[a]



		#liste_all=get_liste_insertion(inpt_all)

		liste_all=Location.get_liste_insertion(vcf)

		dic_other_RM=defaultdict(list)
		dic_other_RM=Location.get_other(liste_all,dic_other_RM,"NA")
		dicon_other=Location.concatenate_dic(dic_other_RM,dic_all_info)
		Function.write_new_vcf(vcf,dicon_other,output)


if __name__ == "__main__":

    main()


    

