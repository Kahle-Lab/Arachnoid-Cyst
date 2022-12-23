#!/bin/python

#Usage: python ga352_denovo_variants.py --vcf filtered_snps_cohort.vcf --fabam Samp_11114.fa_Sorted_MD.bam --mobam Samp_11114.mo_Sorted_MD.bam
#For the record, I had to do this in a miniconda python virtual environment to get pysam to work. I tried all the modules with pysam using 'module spider pysam' but kept getting errors when using pysam to read in the bam files. 
#Note: in order to avoid a weird hexadecimal error, I had to unzip the input vcf before running

# import libraries
import sys
import pysam
import argparse
import pandas as pd
import numpy as np
import io
import os

#Open empty object
parser=argparse.ArgumentParser()

#Add argument for the input vcf
parser.add_argument("--vcf", type=str)

#Add argument for the input parental bam files
parser.add_argument("--fabam", type=str)
parser.add_argument("--mobam", type=str)

#Set-up input args
args=parser.parse_args()
IN_VCF = args.vcf
IN_FAbam = args.fabam
IN_MObam = args.mobam

#Read in the parental bam files
faBAM=pysam.Samfile(IN_FAbam,'rb')
moBAM=pysam.Samfile(IN_MObam,'rb')

############################
# Set up read in functions #
############################

#Define function to read vcf file into pandas 
def ReadInVCF(Input):
    #Open the vcf file 
    with open(Input, 'r') as file:
    	#Remove header (need to do so for filtering in pandas. Will reinstate header at end with )
        Rows=[line for line in file if not line.startswith('##')]
    #Define column value types
    return pd.read_csv(
        io.StringIO(''.join(Rows)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})    
    
#Define function to store the header information 
def StoreVCFHeader(Input):
    
    #Open VCF file and return lines staring w/ ## (header)
    with open(Input, 'r') as f:
        lines=[l for l in f if l.startswith('##')]
    return lines  
    
############################
# Apply read in functions  #
############################

#Read in the jointcalling vcf file and store its header
VCFHeader=StoreVCFHeader(IN_VCF)
JointcallingVCF=ReadInVCF(IN_VCF)

############################
#Define all other functions#
############################

#Define heterozygous calling to accept any format of 0/1
def Heterozygosity(df):
    Heterozygous=((df[0] == '0|1') | (df[0] == '0/1'))
    return Heterozygous

#Define sample not alt to accept any format of 0/0, ./.
def sample_not_alt(df):
    not_alt=((df[0] == '0|0') | (df[0] == '0/0') | (df[0] == '.|.') | (df[0] == './.'))
    return not_alt

#Define heterozygous verification to ensure the allele frequency is between 0.25 and 0.75 and that the depth of reads is 9 or greater (Warns against false-calls on heterozygotes)
def HetVerify(df):
    AlleleDepth=df[1]
    AlleleDepth=AlleleDepth.str.split(',', expand=True).iloc[:,0:2]
    AFrq=AlleleDepth[1].astype(int)/np.sum(AlleleDepth.astype(int), axis=1)
    Verify=((AFrq > 0.25) & (AFrq < 0.75) & (AlleleDepth[1].astype(int) > 8))
    return Verify

#Define filter for DNM!
def DNMFilter(Proband, Father, Mother):
	#Isolate hets (DNM typically hets)
    ProbandHets=Heterozygosity(Proband)
    #Subset proband to remove any non-heterozygous variants
    ProbandSub=Proband[ProbandHets]
    #Verify het calls in proband
    ProbandSub=ProbandSub[HetVerify(ProbandSub)]
    #Index father and mother
    FatherSub=Father.loc[ProbandSub.index,:]
    MotherSub=Mother.loc[ProbandSub.index,:]
    #Set up filters for Father/Mother based on if the allele is not altered in a spot where is is altered in the proband
    FatherFilter=sample_not_alt(FatherSub)
    MotherFilter=sample_not_alt(MotherSub)
    #Subset ProbandSub again removing anything that does not pass Mother/FatherFilters (e.g. remove all transmitted variants)
    ##Note: This will NOT catch DNM where the an unrelated variant in the same location. 
    ##E.g. If for a given location, the mother is 0/0, father is 0/1 with a common, benign variant, the proband is 0/1 or 1/1 with one varian being a rare damaging de novo, I would miss that call.
    ##I don't think that situation and others like them are high enough yield to code around for this problem set as the chance it happens in only two samples is low; but one should be aware this isn't getting everything on all samples
    ProbandSub=ProbandSub[((FatherFilter) & (MotherFilter))]
    return(ProbandSub)

#Define DepthFilter to subset the filtered de novo variants that have more than 9 reads in each parent
def DepthFilter(VCFdf):
    Check=pd.Series(index=VCFdf.index)
    for x in VCFdf.index:
        chromosome=VCFdf['CHROM'].loc[x]
        start=VCFdf['POS'].astype(int).loc[x]
        stop=start+1
        FatherDepth=faBAM.count(chromosome, start, stop)
        MotherDepth=moBAM.count(chromosome, start, stop)
        if FatherDepth > 9 and MotherDepth > 9:
            Check.loc[x]=True
        elif FatherDepth < 10 or MotherDepth < 10:
            Check.loc[x]=False
    DepthFilteredVCFdf=VCFdf[Check]
    return(DepthFilteredVCFdf)
    
##############################
# Set up files for filtering #
##############################

JointcallingVCF=JointcallingVCF[(([len(i) < 7 for i in JointcallingVCF['CHROM']]) & (JointcallingVCF['FILTER'] == 'PASS'))]

#Separate the joint calling VCF data fram into 4 unique dataframes
P1df=JointcallingVCF['11114.p1'].str.split(":", expand=True)
MOdf=JointcallingVCF['11114.mo'].str.split(":", expand=True)
FAdf=JointcallingVCF['11114.fa'].str.split(":", expand=True)
S1df=JointcallingVCF['11114.s1'].str.split(":", expand=True)

##############################
#Filter files using functions#
##############################

#Filter proband and sibling using DNMFilter against mother and father
P1DNM=DNMFilter(P1df, FAdf, MOdf)
S1DNM=DNMFilter(S1df, FAdf, MOdf)
P1DNM_VCF=JointcallingVCF.loc[P1DNM.index,:]
S1DNM_VCF=JointcallingVCF.loc[S1DNM.index,:]

#Further filter proband and sibling using DepthFilter to weed any calls w/ 9 or fewer reads in either parent
P1DNM_DepthFiltered_VCF=DepthFilter(P1DNM_VCF)
S1DNM_DepthFiltered_VCF=DepthFilter(S1DNM_VCF)

##############################
#  Write final output files  #
##############################

#Create final output file for proband
OutP1=open('Samp_11114.p1.DNM.vcf', 'w')
#Write in stored header
for line in VCFHeader:
    OutP1.write(line)
#Close file
OutP1.close()
#Append file with DepthFiltered DNM calls
P1DNM_DepthFiltered_VCF.to_csv('Samp_11114.p1.DNM.vcf', mode='a', sep='\t')

OutS1=open('Samp_11114.s1.DNM.vcf', 'w')
for line in VCFHeader:
    OutS1.write(line)
OutS1.close()
S1DNM_DepthFiltered_VCF.to_csv('Samp_11114.s1.DNM.vcf', mode='a', sep='\t')