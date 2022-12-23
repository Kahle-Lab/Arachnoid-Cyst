library(devtools)
library(denovolyzeR)
library(reshape)
library(dplyr)
options(stringsAsFactors = F)
setwd("~/Desktop/Scripts/DenovolyzeR")

## Read in data to be analyzed
raw_case = read.table("Input_CHv3.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)
raw_case$Cases10GeneName <- toupper(raw_case$Cases10GeneName)
cases = raw_case
case_num = 2416 + 281 #GDX + Yale


raw_control = read.table("1798-Controls-DeNovo-Input.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)
raw_control$Ctrls10GeneName <- toupper(raw_control$Ctrls10GeneName)
ctrls = raw_control
ctrls_num = 1798


cases10 = read.delim("1210_hs37d5_coding_idt_med_v2_spikein_padded_Mar2018_19347-unique-none-0-gene_mpc-modified.txt",sep="\t",stringsAsFactors = F)
controls10 = read.delim("Padded_control_adj_19347-unique-gene_mpc-modified.txt",sep="\t",stringsAsFactors = F)

reformat_pDNM = function(x, mis_filter="Mis_mpc2_or_MetaSVM"){  # MetaSVM only
  names(x)[names(x)==mis_filter] <- "misD"
  x$prot <- x[,names(x) %in% c("lof","mis")] %>% apply(MARGIN=1, sum)
  x$protD <- x[,names(x) %in% c("lof","misD")] %>% apply(MARGIN=1,sum)
  x <- melt(x)                                          #use anything that is 'chr' as id variables
  names(x)[names(x)=="variable"] <- "class"             #change whichever column with the name 'variable' to 'class'
  return(x)
}
unique(reformat_pDNM(cases10)[,"class"])
pDNM_cases10 <- reformat_pDNM(cases10)

unique(reformat_pDNM(controls10)[,"class"])
pDNM_controls10<-reformat_pDNM(controls10)

HBE_genes <- read.table(file='HBE_Genes.txt',sep="\t",header=TRUE,stringsAsFactors = FALSE)
HBE_genes  <- toupper(HBE_genes[[2]])
index <- which(HBE_genes=="N/A")
HBE_genes <- HBE_genes[-index]
HBE_genes <- unique(HBE_genes)

Int_genes <- read.table(file="LoF-Intolerant-Genes-gnomAD2.1.1.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
Int_genes <- toupper(Int_genes[[2]])
index <- which(Int_genes=="N/A")
#Int_genes <- Int_genes[-index]
Int_genes <- unique(Int_genes)

OMIM_genes <- read.table(file="OMIM.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
OMIM_genes <- toupper(OMIM_genes[[2]])
index <- which(OMIM_genes=="N/A")
#OMIM_genes <- OMIM_genes[-index]
OMIM_genes <- unique(OMIM_genes)

# Intolerant HBE genes
HBE_Int_genes <- intersect(HBE_genes, Int_genes)

# Intolerant OMIM genes
Int_OMIM_genes <- intersect(OMIM_genes, Int_genes)

Int_OMIM_HBE_genes <- intersect(Int_OMIM_genes, HBE_genes)

#####################
### GO Gene Lists ###
#####################

#GO:0034647
GO47_genes <- read.table(file="34647.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
GO47_genes <- toupper(GO47_genes[[2]])
index <- which(GO47_genes=="N/A")
#Int_genes <- Int_genes[-index]
GO47_genes <- unique(GO47_genes)

#nBAF
nbaf_genes <- read.table(file="nbaf.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
nbaf_genes <- toupper(nbaf_genes[[2]])
index <- which(nbaf_genes=="N/A")
#Int_genes <- Int_genes[-index]
nbaf_genes <- unique(nbaf_genes)

#histone
histone_genes <- read.table(file="histone.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
histone_genes <- toupper(histone_genes[[2]])
index <- which(histone_genes=="N/A")
#Int_genes <- Int_genes[-index]
histone_genes <- unique(histone_genes)

#Chromatine Redmodelling
chrom_genes <- read.table(file="ChromatinRem.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
chrom_genes <- toupper(chrom_genes[[2]])
index <- which(chrom_genes=="N/A")
#Int_genes <- Int_genes[-index]
chrom_genes <- unique(chrom_genes)

#PosReg
PR_genes <- read.table(file="PosReg.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
PR_genes <- toupper(PR_genes[[2]])
index <- which(PR_genes=="N/A")
#Int_genes <- Int_genes[-index]
PR_genes <- unique(PR_genes)

#TOR
TOR_genes <- read.table(file="TOR.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
TOR_genes <- toupper(TOR_genes[[2]])
index <- which(TOR_genes=="N/A")
TOR_genes <- unique(TOR_genes)


##Chromatin LOF
Chrom_Int_genes <- intersect(chrom_genes, Int_genes)

##Chrom_OMIM_LoF
COL_genes <- intersect(Chrom_Int_genes, OMIM_genes)

setwd("~/Desktop/")
casesByGene = denovolyzeByGene(cases$Cases10GeneName,cases$denovo_call_metaSVM_MPC,case_num, geneId="gene",probTable=pDNM_cases10, includeGenes="all", signifP=3, roundExpected =15, includeClasses=c("protD", "lof", "misD","prot"))
casesByGene = cbind.data.frame(rownames(casesByGene),casesByGene)
colnames(casesByGene)[1] = 'GeneOrder'
write.table(file="Gene_Level_Significance_MetaSVM_MPC-DMis_2697_CH_Cases.txt",casesByGene,col.names=T,row.names=F,sep="\t",append=F,quote=F)
##### Generate enrichment analysis for each functional class
# All genes in 617 cases using the column denovo_call_metaSVM_MPC
DenovolyzeClass <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGene <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGene, file = "DenovolyzeGenes.txt")
write.table(DenovolyzeClass, file = "DenovolyzeClass.txt")


#### HBE
DenovolyzeClassHBE <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=HBE_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGeneHBE <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=HBE_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGeneHBE, file = "DenovolyzeGenesHBE.txt")
write.table(DenovolyzeClassHBE, file = "DenovolyzeClassHBE.txt")

#### INT
DenovolyzeClassInt <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=Int_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGeneInt <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=Int_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGeneInt, file = "DenovolyzeGenesInt.txt")
write.table(DenovolyzeClassInt, file = "DenovolyzeClassInt.txt")

#### OMIM
DenovolyzeClassOMIM <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=OMIM_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGeneOMIM <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=OMIM_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGeneOMIM, file = "DenovolyzeGenesOMIM.txt")
write.table(DenovolyzeClassOMIM, file = "DenovolyzeClassOMIM.txt")
#### HBE+INT
DenovolyzeClassHBEInt <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=HBE_Int_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGeneHBEInt <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=HBE_Int_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGeneHBEInt, file = "DenovolyzeGenesHBEInt.txt")
write.table(DenovolyzeClassHBEInt, file = "DenovolyzeClassHBEInt.txt")

#### OMIM+INT
DenovolyzeClassOMIMInt <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=Int_OMIM_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGeneOMIMInt <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=Int_OMIM_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGeneOMIMInt, file = "DenovolyzeGenesOMIMInt.txt")
write.table(DenovolyzeClassOMIMInt, file = "DenovolyzeClassOMIMInt.txt")

#### OMIM+INT+HBE
DenovolyzeClassOMIMIntHBE <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=Int_OMIM_HBE_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGeneOMIMIntHBE <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=Int_OMIM_HBE_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGeneOMIMIntHBE, file = "DenovolyzeGenesOMIMIntHBE.txt")
write.table(DenovolyzeClassOMIMIntHBE, file = "DenovolyzeClassOMIMIntHBE.txt")


##### Generate enrichment analysis for each functional class

# All genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_Class <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_Class, file = 'DenovolyzeClass.txt')

# HBE genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_Class_HBE <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=HBE_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_Class_HBE, file = 'DenovolyzeClassCTRLHBE.txt')

# LoF_Intolerant genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_Class_INT <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=Int_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_Class_INT, file = 'DenovolyzeClassCTRINT.txt')

# OMIM genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_Class_OMIM <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=OMIM_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_Class_OMIM, file = 'DenovolyzeClassCTRLOMIM.txt')

# HBE_LoF_Intolerant genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_Class_HBE_INT <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=HBE_Int_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_Class_HBE_INT, file = 'DenovolyzeClassCTRLHBEINT.txt')

# OMIM_LoF_Intolerant genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_Class_OMIM_INT <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=Int_OMIM_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_Class_OMIM_INT, file = 'DenovolyzeClassCTRLOMIMINT.txt')

#OMOM_HBE_highPLI
CTRL_Class_OMIM_INT_HBE <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=Int_OMIM_HBE_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_Class_OMIM_INT_HBE, file = 'DenovolyzeClassCTRLOMIMINTHBE.txt')



#########################################
#########################################
#### TOR
DenovolyzeClassTOR <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=TOR_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))


write.table(DenovolyzeClassTOR, file = "DenovolyzeClassTOR.txt")

#### GO:0067647
DenovolyzeClass47 <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=GO47_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGene47 <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=GO47_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGene47, file = "DenovolyzeGenes47.txt")
write.table(DenovolyzeClass47, file = "DenovolyzeClass47.txt")

#### Chromatin Remodeler
DenovolyzeClassCR <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=chrom_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGeneCR <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=chrom_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGeneCR, file = "DenovolyzeGeneszCR.txt")
write.table(DenovolyzeClassCR, file = "DenovolyzeClassCR.txt")

#### histone
DenovolyzeClassHis <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=histone_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGeneHis <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=histone_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGeneHis, file = "DenovolyzeGenesHis.txt")
write.table(DenovolyzeClassHis, file = "DenovolyzeClassHis.txt")

#### nBAF
DenovolyzeClassnBAF <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=nbaf_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGenenBAF <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=nbaf_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGenenBAF, file = "DenovolyzeGenesBAF.txt")
write.table(DenovolyzeClassnBAF, file = "DenovolyzeClassBAF.txt")

### PR 
DenovolyzeClassPR <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=PR_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGenePR <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=PR_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGenePR, file = "DenovolyzeGenesPR.txt")
write.table(DenovolyzeClassPR, file = "DenovolyzeClassPR.txt")

############
### Chrom_Int
DenovolyzeClassCI <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=Chrom_Int_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGeneCI <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=Chrom_Int_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGeneCI, file = "DenovolyzeGenesCI.txt")
write.table(DenovolyzeClassCI, file = "DenovolyzeClassCI.txt")


#### COL
DenovolyzeClassCOL <- denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=COL_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
DenovolyzeGeneCOL <- denovolyzeByGene(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", includeGenes=COL_genes, nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

write.table(DenovolyzeGeneCOL, file = "DenovolyzeGenesCOL.txt")
write.table(DenovolyzeClassCOL, file = "DenovolyzeClassCOL.txt")

#####################

# 47 genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_47 <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=GO47_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_47, file = 'CNTRL_47.txt')

# PR genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_PR <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=PR_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_PR, file = 'CNTRL_PR.txt')

# nbaf genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_nbaf <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=nbaf_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_nbaf, file = 'CNTRL_nbaf.txt')

# histone genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_histone <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=histone_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_histone, file = 'CNTRL_histone.txt')

# Chromatin genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_CR <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=chrom_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_CR, file = 'CNTRL_CR.txt')

# Chromatin + pLI genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_CI <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=Chrom_Int_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_CI, file = 'CNTRL_CI.txt')

# Chromatin + pLI genes in 1798 controls using the column denovo_call_metaSVM_CADD
CTRL_COL <- denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=COL_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(CTRL_COL, file = 'CNTRL_COL.txt')





