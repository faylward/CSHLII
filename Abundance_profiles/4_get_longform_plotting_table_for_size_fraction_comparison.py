#!/usr/bin/python
import sys
from collections import defaultdict
import re
from Bio import SeqIO
import subprocess
import numpy as np

metagenome = open("scaffold_full_abundance_cellular_fraction.txt", "r")
virome = open("scaffold_full_abundance_viral_fraction.txt", "r")

group1 = ['scaffold1|size70376', 'scaffold6|size43974', 'scaffold12|size37975', 'scaffold153|size11910']

def scaffold2group(scaffold):
	if scaffold in group1:
		group = "VG1"
		size = 164000
	else:
		list1 = scaffold.split("|")
		scaf = list1[0]
		num = re.sub("scaffold", "", scaf)
		group = "VS"+num
		size = list1[1]
		size = int(re.sub("size", "", size))
	return(group, size)

size_dict = {}
# iterate through abundance matrix
virome_dict = defaultdict(list)
for i in virome.readlines():
	line = i.rstrip()
	if line.startswith("S06"):
		casts = line.split("\t")
	else:
		abund = line.split("\t")
		scaffold = abund[0]
		values = abund[1:len(abund)]
		group, size = scaffold2group(scaffold)

		size_dict[group] = size
		
		for index, abundance in enumerate(values):
			virome_dict[group].append(float(abundance))


# iterate through ratio matrix
mg_dict = defaultdict(list)
for i in metagenome.readlines():
	line = i.rstrip()
	if line.startswith("S06"):
		casts = line.split("\t")
	else:
		abund = line.split("\t")
		scaffold = abund[0]
		values = abund[1:len(abund)]
		group, size = scaffold2group(scaffold)
		
		for index, abundance in enumerate(values):
			mg_dict[group].append(float(abundance))

#size_dict["VG1"] = 164000
o = open("longform_data_newdata_for_comparison.txt", "w")
o.write("group\tsize\tvirome\tmetagenome\n")
# get aggregate group values
final_virome = {}
final_mg = {}
for j in virome_dict:
	print(j)
	#print(mg_dict[j])
	mean_virome = np.sum(virome_dict[j])
	mean_mg = np.sum(mg_dict[j])
	size = size_dict[j]

	o.write(j +"\t"+ str(size) +"\t"+ str(mean_virome) +"\t"+ str(mean_mg)  +"\n")






