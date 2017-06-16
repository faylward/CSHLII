#!/usr/bin/python
import sys
from collections import defaultdict
import re
from Bio import SeqIO
import subprocess
import numpy as np

abundance = open("scaffold_combined_abundance_for_plotting.txt", "r")
ratio = open("scaffold_ratio_for_plotting.txt", "r")
index = open("plotting_index2.txt", "r")
cast_index = open("cast_index.txt", "r")

group1 = ['scaffold1|size70376', 'scaffold6|size43974', 'scaffold12|size37975', 'scaffold153|size11910']

def scaffold2group(scaffold):
	if scaffold in group1:
		group = "VG1"
	else:
		list1 = scaffold.split("|")
		scaf = list1[0]
		num = re.sub("scaffold", "", scaf)
		group = "VS"+num
	return(group)

# get viral group indices
group2index = {}
for i in index.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	group = tabs[0]
	index = tabs[5]
	group2index[group] = index

# get cast indices
cast2index = {}
for i in cast_index.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	cast = tabs[0]
	index = tabs[1]
	cast2index[cast] = index

# iterate through abundance matrix
abundance_dict = defaultdict(list)
for i in abundance.readlines():
	line = i.rstrip()
	if line.startswith("S06"):
		casts = line.split("\t")
	else:
		abund = line.split("\t")
		scaffold = abund[0]
		values = abund[1:len(abund)]
		group = scaffold2group(scaffold)
		
		for index, abundance in enumerate(values):
			cast = casts[index]
			group_cast = group +"__"+ cast
			abundance_dict[group_cast].append(float(abundance))

# iterate through ratio matrix
ratio_dict = defaultdict(list)
for i in ratio.readlines():
	line = i.rstrip()
	if line.startswith("S06"):
		casts = line.split("\t")
	else:
		abund = line.split("\t")
		scaffold = abund[0]
		values = abund[1:len(abund)]
		group = scaffold2group(scaffold)
		
		for index, ratio in enumerate(values):
			cast = casts[index]
			group_cast = group +"__"+ cast
			if ratio == "NA":
				ratio_dict[group_cast].append(ratio)
			else:
				ratio_dict[group_cast].append(float(ratio))


# get aggregate group values
final_abund = {}
final_ratio = {}
for j in abundance_dict:
	if len(abundance_dict[j]) > 1:
		mean_abund = np.mean(abundance_dict[j])
		mean_ratio = np.mean(ratio_dict[j])
	else:
		abund = abundance_dict[j]
		ratio = ratio_dict[j]
		#print(abund, ratio, type(abund), type(ratio))
		mean_abund = abund[0]
		mean_ratio = ratio[0]

	final_abund[j] = mean_abund
	final_ratio[j] = mean_ratio
	#print(mean_abund, type(mean_abund))

o = open("longform_data_newdata.txt", "w")
o.write("name\tgroup\tcast\tgroup_index\tcast_index\tabundance\tvc\n")
### output long form data tables
for j in final_abund:
	names = j.split("__")
	group = names[0]
	cast = names[1]
	#print(group, cast, j)
	if group in group2index:
		#o.write(j +"\t"+ group +"\t"+ cast +"\t"+ str(group2index[group]) +"\t"+ str(cast2index[cast]) +"\t"+ str(mean_abund[j]) +"\t"+ str(mean_ratio[j]) +"\n")
		o.write(j +"\t"+ group +"\t"+ cast +"\t"+ str(group2index[group]) +"\t"+ str(cast2index[cast]) +"\t"+ str(final_abund[j]) +"\t"+ str(final_ratio[j]) +"\n")






