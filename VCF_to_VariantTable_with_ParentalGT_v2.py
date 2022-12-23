#!/usr/bin/python

###################################
# Purpose: Filter annovar-annotated vcf file, extract mutation carriers, and get rid of 0/0 and ./.
# Author: Sheng Chih Jin <shengchih.jin@yale.edu>
# Version: Python3
# Last Modified Date: 03-15-2020
###################################

##Setup the important stuff
#Import necessary packages
import sys, os, numbers, linecache

if len(sys.argv) != 3:
	sys.stderr.write("Usage: python VCF_to_VariantTable_with_ParentalGT_v2.py <vcf file> <ped file>\n")
	sys.stderr.write("Usage: python VCF_to_VariantTable_with_ParentalGT_v2.py exome_calls_pass_step2_normalized_anno.hg19_multianno_0.001.vcf 29Trios.ped \n")
	sys.exit(-1)

cwd = os.getcwd()

ped_list = []

# Convert DP or GQ into 0 if '.' is found
def convertQual(DP):
	if DP == '.':
		return 0
	else:
		return float(DP)

trio_list = {}

pedfile = sys.argv[2]
with open(pedfile,'r') as file:
	for line in file:
		data = line.strip().split('\t')
		proband = data[1]
		father = data[2]
		mother = data[3]
		if father != '0' and mother != '0':
			ped_list.append([proband,father,mother])
		trio_list[proband] = []
		trio_list[father] = []
		trio_list[mother] = []

# count how many unique key-value pairs
#n_keys = 0
#for key,value in trio_list.items():
#	print(key,value)
#	n_keys = n_keys + 1

#print("n_keys = %d " % n_keys)

# get the full sample list from the VCF
vcf_file = sys.argv[1]
with open(vcf_file, 'r') as file:
	for line in file:
		if line[0] == '#':
			if '#CHROM' in line:
				data = line.strip().split('\t')
				list = data[9:]	

# get the singleton list by excluding samples from ped (trios)
singleton_list = {}
for ID in list:
	if ID not in trio_list:
		singleton_list[ID] = []


out = open(vcf_file[:-4]+'_variantTable.txt','w')
header = []
header_flag = 'F'
order = []

with open(vcf_file, 'r') as file:
	for line in file:
		data = line.strip().split('\t')
		if line.startswith('#CHROM'):
			store = data	# store saves each header
		elif line.startswith('##'):
			continue
		else:
			chr = data[0]
			pos = data[1] 
			info = data[7].split(';')  # info stores the info of each line
			Info = {}	# Info is a dictionary; key is the metric; value is the variant info
			for item in info:
				if len(item.split('=')) == 2:
					Info[item.split('=')[0]] = item.split('=')[1]
					order.append(item.split('=')[0])

			# print out header
			if header_flag == 'F':
				for key in order:
					header.append(key)
				out.write('ProbandID\tProband_GT\tProband_AD\tProband_VAF\tProband_DP\tProband_GQ\tFatherID\tFather_GT\tFather_AD\tFather_VAF\tFather_DP\tFather_GQ\tMotherID\tMother_GT\tMother_AD\tMother_VAF\tMother_DP\tMother_GQ\t'+'\t'.join(store[0:7])+'\t'+'\t'.join(header)+'\n')
				header_flag = 'T'
			
			# exonic
			Func = Info['Func.refGene']
			ExonicFunc = Info['ExonicFunc.refGene']
			
			# fill in the information for trios
			format = data[store.index('FORMAT')].split(':')
			for element in ped_list:
				proband = element[0]
				father = element[1]
				mother = element[2]	
				if proband not in store:
					continue
				proband_info = data[store.index(proband)]
				father_info = data[store.index(father)]
				mother_info = data[store.index(mother)]

				proband_GT = proband_info.split(':')[0]
				father_GT = father_info.split(':')[0]
				mother_GT = mother_info.split(':')[0]

				if (proband_GT == '0/0' or proband_GT == './.') and (father_GT == '0/0' or father_GT == './.') and (mother_GT == '0/0' or mother_GT == './.'):
					continue

				# Deal with Proband
				if len(proband_info.split(':')) < len(format):
					proband_AD = 'NA'
					proband_DP = 'NA'
					proband_GQ = 'NA'
				else:
					proband_AD = proband_info.split(':')[format.index('AD')]
					proband_DP = proband_info.split(':')[format.index('DP')]
					proband_GQ = proband_info.split(':')[format.index('GQ')]
				# Deal with father
				if len(father_info.split(':')) < len(format): 
					father_AD = 'NA'
					father_DP = 'NA'
					father_GQ = 'NA'
				else:
					father_AD = father_info.split(':')[format.index('AD')]
					father_DP = father_info.split(':')[format.index('DP')]
					father_GQ = father_info.split(':')[format.index('GQ')]
				# Deal with mother
				if len(mother_info.split(':')) < len(format):
					mother_AD = 'NA'
					mother_DP = 'NA'
					mother_GQ = 'NA'
				else:
					mother_AD = mother_info.split(':')[format.index('AD')]	
					mother_DP = mother_info.split(':')[format.index('DP')]
					mother_GQ = mother_info.split(':')[format.index('GQ')]

				if proband_DP == '0' or proband_DP == '.' or proband_DP == 'NA':
					proband_VAF = 'NA'
				else:
					proband_VAF = str(float(proband_AD.split(',')[1])/float(proband_DP))

				if father_DP == '0' or father_DP == '.' or father_DP == 'NA':
					father_VAF = 'NA'
				else:
					father_VAF = str(float(father_AD.split(',')[1])/float(father_DP))

				if mother_DP == '0' or mother_DP == '.' or mother_DP == 'NA':
					mother_VAF = 'NA'
				else:
					mother_VAF = str(float(mother_AD.split(',')[1])/float(mother_DP))

				out.write('\t'.join([proband,proband_GT,proband_AD,proband_VAF,proband_DP,proband_GQ,father,father_GT,father_AD,father_VAF,father_DP,father_GQ,mother,mother_GT,mother_AD,mother_VAF,mother_DP,mother_GQ])+'\t')

				#print out line contents
				out.write('\t'.join(data[0:7]))
			

				for i in header:
					if i not in Info:
						out.write('\t'+'NA')
					else:
						out.write('\t'+Info[i])
				out.write('\n')
		
			# for each singleton proband
			for proband in singleton_list:
				proband_info = data[store.index(proband)]

				proband_GT = proband_info.split(':')[0]
				if (proband_GT == '0/0' or proband_GT == './.'):
					continue
				proband_AD = proband_info.split(':')[format.index('AD')]
				proband_DP = proband_info.split(':')[format.index('DP')]
				if proband_DP == '.' or proband_DP == '0':
					continue
				proband_GQ = proband_info.split(':')[format.index('GQ')]
				if proband_GQ == '.' or proband_GQ == '0':
					continue
				proband_VAF = str(float(proband_AD.split(',')[1])/float(proband_DP))	

				out.write('\t'.join([proband,proband_GT,proband_AD,proband_VAF,proband_DP,proband_GQ,'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA'])+'\t')
				out.write('\t'.join(data[0:7]))
				for i in header:
					if i not in Info:
						out.write('\t'+'NA')
					else:
						out.write('\t'+Info[i])
				out.write('\n')


out.close()
