ssh -Y ga352@ruddle.hpc.yale.edu

srun -p interactive --x11=batch --pty --mem=120000 bash

TMMnorm<-read.csv("/Users/garrettallington/Desktop/WGCNA/Younger_Than_Three_mRNA-seq_PsychENCODE.csv",head=T,row.name=1,stringsAsFactors=FALSE)

TMMnorm<-read.csv("Younger_Than_Three_mRNA-seq_PsychENCODE.csv",head=T,row.name=1,stringsAsFactors=FALSE)

setwd('/Users/garrettallington/Desktop/WGCNA')

# create target.csv files
conds.all<-read.csv("traits_CS_U3.csv",head=T,row.name=1,stringsAsFactors=FALSE)
rownames(conds.all)<-conds.all$sample
trait<-conds.all[colnames(TMMnorm),];trait
conds<-trait[,1:2];conds


#target<-read.csv("target.csv",head=T,row.name=1,stringsAsFactors=FALSE)
#rownames(target)<-target$sample;target
#target<-target[conds$sample,];target
#TMMnorm<-TMMnorm[,rownames(target)];head(TMMnorm)
############################
cv <- as.matrix(conds)
dim(cv) <- c(1,prod(dim(cv)))
rownames(table(cv))
names <- rownames(table(cv))
col_default <- c("red ","blue","green","orange","bisque4 ","black ","brown ","cyan ","darkgreen ","darkgrey ","darkmagenta ","darkolivegreen ","darkorange ","darkred ","darkslateblue ","darkturquoise ","floralwhite ","greenyellow ","grey ","lightcyan ","lightcyan1 ","lightgreen ","lightsteelblue1 ","lightyellow ","magenta ","mediumpurple3 ","midnightblue ","paleturquoise ","pink ","plum1 ","plum2 ","royalblue ","saddlebrown ","salmon ","sienna3 ","skyblue ","skyblue3 ","steelblue ","tan ","thistle1 ","thistle2 ","turquoise ","violet ","white ","yellowgreen","grey60 ","orangered4 ","brown4 ","darkorange2 ","ivory ")
col <- rainbow(length(names))
names(col)<- names
clab <- matrix(col[as.matrix(conds)],nrow=nrow(conds),ncol=ncol(conds))
colnames(clab) <- colnames(conds) 
rownames(clab)<-rownames(conds)
a=1
for (i in unique(clab[,2])){
  clab[which(clab[,2]==i),2]=col_default[a]
  a=a+1
}
clab

library(WGCNA)
library(flashClust)
sdout <- 3
## Remove outliers
##Calculate signed, weighted biweight midcorrelation
normExpr=TMMnorm
normadj <- (0.5+0.5*bicor(normExpr)^2)

## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- ku-(mean(ku))/sqrt(var(ku))
## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep=""))
print(colnames(normExpr)[outliers])
print(table(outliers))


conds<-conds[rownames(conds)%in%colnames(normExpr),];conds
conds.all<-conds.all[rownames(conds.all)%in%colnames(normExpr),];conds.all

pdf("boxplot_nromExpr.pdf",height=8,width=8)
temp<-normExpr
colnames(temp)[colnames(temp)%in%colnames(normExpr)[outliers]]<-paste("**",colnames(temp)[colnames(temp)%in%colnames(normExpr)[outliers]],sep="")
par(cex.axis=0.5) #
boxplot(temp,las=2,cex.x=0.4)
dev.off()

#normExpr<-normExpr[,!colnames(normExpr)%in%colnames(normExpr)[outliers]];dim(normExpr)
datExpr=as.data.frame(t(normExpr))

conds.all<-conds.all[colnames(normExpr),];dim(conds.all)
conds<-conds[colnames(normExpr),];dim(conds)

#target<-target[colnames(normExpr),]
trait<-trait[colnames(normExpr),]
# read trait file
dim(trait)

pdf("1.1_power.pdf", height=8, width=10)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType="signed",corFnc="bicor")
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence Control All regions",sep=""));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
abline(h=0.80,col="blue")
abline(h=0.70,col="orange")
abline(h=0.60,col="green")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity Control All regions",sep=""))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


dev.off()

## Co-expression similarity and adjacency

softPower = 10;
#dir.create(paste("Power_",softPower,sep=""))
#setwd(paste("Power_",softPower,sep=""))
for (i in ncol(datExpr)){
  datExpr[,i]<-as.numeric(datExpr[,i])
}
adjacency = adjacency(datExpr, power = softPower, type = "signed");

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = flashClust(as.dist(dissTOM), method = "average");

#save(TOM,dissTOM,geneTree,adjacency,softPower,file=paste("TOM.rda",sep=""))


###################################

### Relating dendrogram with traits
Condition=as.numeric(factor(conds.all$group))
#PC1_Sequencing=target$Seq.PC1
#PC2_Sequencing=target$Seq.PC2


geneSigs=matrix(NA,nrow=3,ncol=ncol(datExpr)) # create a vector to hold the data

for(i in 1:ncol(geneSigs)) {
  
  exprvec=as.numeric(datExpr[,i]) # get the expression vector for ith gene
  conditionr=bicor(Condition,exprvec,use="pairwise.complete.obs")
  #pc1r=bicor(PC1_Sequencing,exprvec,use="pairwise.complete.obs")
  #pc2r=bicor(PC2_Sequencing,exprvec,use="pairwise.complete.obs")
  
  
  geneSigs[,i]=c(conditionr) #,pc1r,pc2r)
  
  cat('Done for gene...',i,'\n')
}

clab<-clab[colnames(normExpr),]

geneSigs[1,] =numbers2colors(as.numeric(geneSigs[1,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[2,] =numbers2colors(as.numeric(geneSigs[2,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[3,] =numbers2colors(as.numeric(geneSigs[3,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))

rownames(geneSigs)=c("Condition","PC1_Sequencing","PC2_Sequencing")


mColorh <- mLabelh <- colorLabels <- NULL
for (minModSize in c(40,100,160)) {
  for (dthresh in c(0.1,0.2,0.25)) {
    for (ds in c(2,4)) {
      print("Trying parameters:")
      print(c(minModSize,dthresh,ds))
      tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE, minClusterSize = minModSize, cutHeight = 0.9999, deepSplit = ds, distM = as.matrix(dissTOM))
      merged <- mergeCloseModules(exprData = datExpr,colors = tree$labels,cutHeight = dthresh)
      mColorh <- cbind(mColorh,labels2colors(merged$colors))
      mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
    }
  }
}


mColorh1=cbind(mColorh,geneSigs[1,])
mLabelh1=c(mLabelh,rownames(geneSigs))


pdf("1.2_SignedDendro.pdf",height=25,width=20)
plotDendroAndColors(geneTree,mColorh1,groupLabels=mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Different Module Cutting Parameters")
dev.off()

mms=40
ds =4
dthresh=0.1

tree = cutreeHybrid(dendro = geneTree, pamStage=F, minClusterSize =mms, cutHeight = 0.9999, deepSplit = ds, distM = as.matrix(dissTOM))
merge <- mergeCloseModules(exprData = datExpr,colors = tree$labels, cutHeight = dthresh)
mColorh <- cbind(labels2colors(merge$colors),geneSigs[1,])
mLabelh <- c("Merged Colors",rownames(geneSigs))

pdf("1.3_FinalDendro_40_4_0.1.pdf",height=10,width=16)
plotDendroAndColors(geneTree, mColorh, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power = ",softPower,"mms=",mms,"ds=",ds,"dthresh=",dthresh));
dev.off()

mergedColors = labels2colors(merge$colors);

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
MEList=moduleEigengenes(datExpr, colors = mergedColors,softPower= softPower, nPC=1)
MEs=MEList$eigengenes
MEs=orderMEs(MEs)
moduleColors = mergedColors

write.table(table(moduleColors), file="final_modules_40_4_0.1.csv", sep=",")
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");

# Plot the result

pdf("1.4_Clustering_of_module_eigengene_40_4_0.1.pdf", height=10, width=15)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()
##############################
#Alldegrees1=intramodularConnectivity(adjacency, moduleColors)
head(Alldegrees1)
rownames(MEs)<-rownames(datExpr)
write.table(MEs, file="1.4_ME_output_40_4_0.1.csv", sep=",")

meInfo<-data.frame(rownames(datExpr), MEs)
colnames(meInfo)[1]= "SampleID"
KMEs<-signedKME(datExpr, MEs,outputColumnName = "kME",corFnc = "bicor")

refSeq <-read.csv("GeneList.csv", header=TRUE, row.names=1)
refSeq <- as.data.frame(refSeq)
colnames(refSeq)[1]="gene"
rownames(refSeq)<-refSeq$gene

table(colnames(datExpr)%in%refSeq$gene)
refSeq=refSeq[na.omit(match(colnames(datExpr),refSeq$gene)),]
refSeq <- as.data.frame(refSeq)
colnames(refSeq)[1]="gene"
###Changed from refSeq$gene to refSeq[1] to avoid atomic vector error
geneInfo=as.data.frame(cbind(refSeq$gene,refSeq$gene,moduleColors, KMEs,Alldegrees1))

# merged gene symbol column
colnames(geneInfo)[1]= "refSeq.Gene.ID"
colnames(geneInfo)[2]= "GeneSymbol"
colnames(geneInfo)[3]= "Initially.Assigned.Module.Color"

write.csv(geneInfo,file=paste('1.5_geneInfoSigned_40_4_0.1.csv',sep=''))


PCvalues=MEs[,-ncol(MEs)]

pdf('1.6_ME_trajectory_Plot_40_4_0.1.pdf',width=14,height=12)
##
toplot=t(PCvalues)
cols=substring(colnames(MEs),3,20)
par(mfrow=c(4,3))
par(mar=c(5,6,4,2))
for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(as.vector(as.factor(conds.all$group)),unique(conds.all$group)),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
}  #verboseScatterplot(x=as.numeric(target$Seq.PC1),y=toplot[i,],xlab="PC1_Sequencing",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
  #verboseScatterplot(x=as.numeric(factor(target$Seq.PC2)),y=toplot[i,],xlab="PC2_Sequencing",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
#}

dev.off()


###############################################
###output with scaled##########################
modu<- as.character(unique(geneInfo$Initially.Assigned.Module.Color))
out1<- matrix(nrow=nrow(normExpr),ncol=13)
for (j in 1:length(modu)){
  out1temp<- geneInfo [geneInfo$Initially.Assigned.Module.Color==modu[j],]
  xtra<- out1temp$kWithin/max(out1temp$kWithin)
  out1<-as.data.frame(cbind(out1temp[,1:8],scaled_kWithin=xtra, out1temp[,9:10]))
  if (j==1)
    out<-out1
  else
    out<-rbind(out,out1)
  print(j)
}


datExpr2<-t(normExpr)
nGenes = ncol(datExpr2);
nSamples = nrow(datExpr2);

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

outmm<-cbind(geneInfo[,1:5], geneModuleMembership, MMPvalue)
dim(outmm)
write.table(outmm, "1.7_module_memebership_40_4_0.1.csv", sep=",",row.names=F)

##output KE and module info

ke<-merge(out,outmm, by.x="refSeq.Gene.ID",by.y="refSeq.Gene.ID")
colnames(ke)
ke<-ke[,c(1:11,16:dim(ke)[2])] ##take extra columns out
colnames(ke)
#######################################




a<-tapply(as.character(ke$Probe),ke$moduleColor.x,c)

datTraits<-trait[,-c(1:2)]
#datTraits<-datTraits[,-2]
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
#datTraits1<- datTraits[,1:5]

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

colnames(moduleTraitPvalue) = paste("p.value.", colnames(moduleTraitCor), sep="");

out3<-cbind(Module=rownames(moduleTraitCor ), moduleTraitCor, moduleTraitPvalue)
dim(out3)
write.table(out3, "1.8_moduleTraitCor_40_4_0.1.csv", sep=",",row.names=F)

##plot traits correlation 

textMatrix = paste( signif(moduleTraitCor, 2), '\n(', 
                    signif(moduleTraitPvalue, 1), ')', sep = ''
);
dim(textMatrix) = dim(moduleTraitCor);
pdf("1.9_heatmap_moduleTraitCor_40_4_0.1.pdf", width=15,height=15)

par( mar = c(8, 12, 3, 3) );
labeledHeatmap( Matrix = moduleTraitCor,
                xLabels = names(datTraits),
                yLabels = names(MEs),
                ySymbols = names(MEs),
                colorLabels = F,
                colors = greenWhiteRed(50),
                textMatrix = textMatrix,
                setStdMargins = F,
                cex.text = 0.4,
                zlim = c(-1, 1),
                main = paste("module-trait relationships")
); 
sig = moduleTraitPvalue < (.05 / (ncol(moduleTraitPvalue) * nrow(moduleTraitPvalue)));
dev.off()

## Generating files for GSEA
#I removed rownames from tempMM to get it to run
tempMM<-read.csv("1.7_module_memebership_40_4_0.1.csv",head=T,row.names=1,stringsAsFactors = F);head(tempMM)
geneModuleMembership<-tempMM[,-c(1:4)]
rownames(geneModuleMembership)<-toupper(rownames(geneModuleMembership))
for (i in colnames(geneModuleMembership)){
  temp<-geneModuleMembership[order(-geneModuleMembership[,i]),]
  temp$gene<-rownames(temp)
  temp<-temp[,c("gene",i)]
  write.table(temp,file=paste("Dir_MM_",i,"_40_4_0.1.rnk",sep=""),col.names=FALSE,sep="\t",row.name=FALSE,quote=FALSE)
}

##############################
### Co-expression PPI Plot ###
##############################
library(igraph)

uniquemodcolors = unique(moduleColors);
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']

#load("~/Dropbox (UCLA Mednet)/R codes/BGandIWcombinedPPI_5-8-2014.Rdata")

pdf("2.0_PPI_ModuleCirclePlots.pdf",width=12,height=12)
par(mar=c(3,3,3,3))

for (mod in uniquemodcolors)  {
  
  thismod=geneInfo[geneInfo$Initially.Assigned.Module.Color==mod,]
  modgenes=thismod$refSeq.Gene.ID
  geneNames=thismod$GeneSymbol
  gnS=intersect(modgenes,colnames(datExpr))
  
  thisExpr=datExpr[,match(gnS,colnames(datExpr))]
  thismod=thismod[match(gnS,modgenes),]
  colnames(thisExpr) = toupper(thismod$GeneSymbol)
  
  
  ## Get PPI data ppiMat
  
  keepgenes <- intersect(rownames(ppiMat),colnames(thisExpr))
  
  coexpInd <- na.omit(match(keepgenes,colnames(thisExpr)))
  ppiInd <- na.omit(match(keepgenes,rownames(ppiMat)))
  
  seedExpr <- thisExpr[,coexpInd]
  #modColors <- colorMat[,coexpInd]
  seedPPI <- ppiMat[ppiInd,ppiInd]
  colnames(seedExpr) <- rownames(seedPPI)
  
  ## Make an adjacency matrix
  adjMat <- bicor(seedExpr) ## adjacency matrix for this module
  combMat=adjMat*seedPPI
  
  maxsize <- min(50,nrow(combMat))
  
  kME <- apply(combMat,2,sum) ## within-module kME
  cutoff <- sort(kME,decreasing=TRUE)[maxsize] ## restrict to max module size
  tophubs <- kME>=cutoff
  
  combMat <- combMat[tophubs,tophubs]
  numcors <- min(300,(maxsize^2-maxsize)/2)
  topcors <- sort(as.numeric(combMat),decreasing=TRUE)[numcors]
  combMat[combMat<=topcors] <- 0
  
  ## Different network layouts
  gA <- graph.adjacency(as.matrix(combMat[1:10,1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
  gB <- graph.adjacency(as.matrix(combMat[11:maxsize,11:maxsize]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutCircle <- rbind(layout.circle(gA)/2,layout.circle(gB))
  
  g1 <- graph.adjacency(as.matrix(combMat),mode="undirected",weighted=TRUE,diag=FALSE)
  
  
  plot.igraph(g1,vertex.label=as.character(rownames(combMat)), ##From library R.utils
              vertex.label.dist=0.8,
              #vertex.shape="pie",
              #vertex.pie=valueList,
              #vertex.pie.color=colorList,
              vertex.size=3,#kME[tophubs]*10,
              vertex.label.color="black",
              vertex.color=mod,
              #vertex.frame.color=module,
              vertex.label.cex=1.3,
              vertex.frame.color="black",
              layout=jitter(layoutCircle),
              #layout=layoutFR,
              edge.color="green",
              main=paste("Co-expression PPI plot of",mod))
  
  cat('Done ....',mod,'\n')
  
}

dev.off()



#########Co-expression only Plot
library(WGCNA)
library(igraph)
library(RColorBrewer);

TOM.matrix=as.matrix(TOM)
uniquemodcolors = unique(moduleColors);
pdf("2.1_Coexpression_ModuleCirclePlots.pdf",width=12,height=12,useDingbats=F)
par(mar=c(3,3,3,3))

for (mod in uniquemodcolors)  {
  # mod="greenyellow"
  numgenesingraph = 50;
  numconnections2keep = 500;
  cat('module:',mod,'\n');
  geneInfo=geneInfo[geneInfo$GeneSymbol!="NA",]
  colind = which(colnames(geneInfo)==paste('kME',mod,sep=''));
  rowind = which(geneInfo[,3]==mod);
  cat(' ',length(rowind),'probes in module\n');
  submatrix = geneInfo[rowind,];
  orderind = order(submatrix[,colind],decreasing=TRUE);
  if (length(rowind) < numgenesingraph) {
    numgenesingraph = length(rowind);
    numconnections2keep = numgenesingraph * (numgenesingraph - 1);
  }
  cat('Making network graphs, using top',numgenesingraph,'probes and',numconnections2keep,'connections of TOM\n');
  submatrix = submatrix[orderind[1:numgenesingraph],];
  #Identify the columns in the TOM that correspond to these hub probes
  matchind = match(submatrix$refSeq.Gene.ID,colnames(datExpr));
  reducedTOM = TOM.matrix[matchind,matchind];
  
  orderind = order(reducedTOM,decreasing=TRUE);
  connections2keep = orderind[1:numconnections2keep];
  reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
  reducedTOM[connections2keep] = 1;
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[1:10,1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMata <- layout.circle(g0)
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[11:50,11:50]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMatb <- layout.circle(g0)
  
  g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMat <- rbind(layoutMata*0.5,layoutMatb)
  
  plot(g1,edge.color="grey",vertex.color=mod,vertex.label=as.character(submatrix$GeneSymbol),vertex.label.cex=1.3,vertex.label.dist=0.9,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat,vertex.size=submatrix[,colind]^2*8,main=paste(mod,"module"))
  cat('Done ....',mod,'\n')
  
}


dev.off()


############################
############################

#rm(list=ls())
#load('WGCNA_DRG.rda')
library(enrichR)
dbs<-c('GO_Biological_Process_2018','GO_Cellular_Component_2018','GO_Molecular_Function_2018',"Human_Gene_Atlas","WikiPathways_2019_Human",'KEGG_2019_Human','TRANSFAC_and_JASPAR_PWMs','LINCS_L1000_Chem_Pert_up','LINCS_L1000_Chem_Pert_down')


uniquemodcolors=as.character(unique(geneInfo$Initially.Assigned.Module.Color))
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']


for(i in 1:length(uniquemodcolors)){
  thismod= uniquemodcolors[i]
  
  pdf(paste(thismod,"1_EnrichR_output.pdf",sep=''),height=8,width=10)
  par(mar=c(4,22,4,4))
  thisInfo=geneInfo[geneInfo$Initially.Assigned.Module.Color==thismod,]
  geneNames=as.character(thisInfo$GeneSymbol)
  enriched <- enrichr(geneNames,dbs)
  printEnrich(enriched, file=paste(thismod,"EnrichR_output.csv",sep=','), sep = ",", columns = c(1,3,4,8))
  

  tmp.GO1<-enriched[[1]][1:10,c(1,8)]
  tmp.GO1=tmp.GO1[order(tmp.GO1$Combined.Score,decreasing=F),]
  tmp.GO2<-enriched[[2]][1:10,c(1,8)]
  tmp.GO2=tmp.GO2[order(tmp.GO2$Combined.Score,decreasing=F),]
  
  tmp.GO3<-enriched[[3]][1:10,c(1,8)]
  tmp.GO3=tmp.GO3[order(tmp.GO3$Combined.Score,decreasing=F),]
  
  tmp.HGA<-enriched[[4]][1:10,c(1,8)]
  tmp.HGA=tmp.HGA[order(tmp.HGA$Combined.Score,decreasing=F),]
  
  tmp.WIKI<-enriched[[5]][1:10,c(1,8)]
  tmp.WIKI=tmp.WIKI[order(tmp.WIKI$Combined.Score,decreasing=F),]
  
  tmp.TF<-enriched[[7]][1:10,c(1,8)]
  tmp.TF=tmp.TF[order(tmp.TF$Combined.Score,decreasing=F),]
  
  tmp.KEGG<-enriched[[6]][1:10,c(1,8)]
  tmp.KEGG=tmp.KEGG[order(tmp.KEGG$Combined.Score,decreasing=F),]
  
  tmp.LincsUP<-enriched[[8]][1:10,c(1,8)]
  tmp.LincsUP=tmp.LincsUP[order(tmp.LincsUP$Combined.Score,decreasing=F),]
  
  tmp.LincsDOWN<-enriched[[9]][1:10,c(1,8)]
  tmp.LincsDOWN=tmp.LincsDOWN[order(tmp.LincsDOWN$Combined.Score,decreasing=F),]
  
  wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
  }
  
  barplot(tmp.GO1$Combined.Score,horiz=T,col="blue",names.arg= wrapText(tmp.GO1$Term,45),cex.names=1.2,las=1,main=paste("GO_Biological_Process of",thismod,"Module"),xlab="Z-Score")
  abline(v=2,col="red")
  
  
  barplot(tmp.GO2$Combined.Score,horiz=T,col="blue",names.arg= wrapText(tmp.GO2$Term,45),cex.names=1.2,las=1,main=paste("GO_Cellular_Component of",thismod,"Module"),xlab="Z-Score")
  abline(v=2,col="red")
  
  barplot(tmp.GO3$Combined.Score,horiz=T,col="blue",names.arg= wrapText(tmp.GO3$Term,45),cex.names=1.2,las=1,main=paste("GO_Molecular_Function of",thismod,"Module"),xlab="Z-Score")
  abline(v=2,col="red")
  
  barplot(tmp.WIKI$Combined.Score,horiz=T,col="blue",names.arg= wrapText(tmp.WIKI$Term,45),cex.names=1.2,las=1,main=paste("Wiki Pathways 2019 Human of",thismod,"Module"),xlab="Z-Score")
  abline(v=2,col="red")
  
  barplot(tmp.HGA$Combined.Score,horiz=T,col="blue",names.arg= wrapText(tmp.HGA$Term,45),cex.names=1.2,las=1,main=paste("Human Gene Atlas of",thismod,"Module"),xlab="Z-Score")
  abline(v=2,col="red")
  barplot(tmp.TF$Combined.Score,horiz=T,col="blue",names.arg= wrapText(tmp.TF$Term,50),cex.names=1.2,las=1,main=paste("TRANSFAC_and_JASPAR of",thismod,"Module"),xlab="Z-Score")
  abline(v=2,col="red")
  
  barplot(tmp.KEGG$Combined.Score,horiz=T,col="blue",names.arg= wrapText(tmp.KEGG$Term,50),cex.names=1.2,las=1,main=paste("KEGG_2019 of",thismod,"Module"),xlab="Z-Score")
  abline(v=2,col="red")
  
  barplot(tmp.LincsUP$Combined.Score,horiz=T,col="blue",names.arg= wrapText(tmp.LincsUP$Term,50),cex.names=1.2,las=1,main=paste("LINCS_L1000_UP of",thismod,"Module"),xlab="Z-Score")
  abline(v=2,col="red")
  
  barplot(tmp.LincsDOWN$Combined.Score,horiz=T,col="blue",names.arg= wrapText(tmp.LincsDOWN$Term,50),cex.names=1.2,las=1,main=paste("LINCS_L1000_DOWN of",thismod,"Module"),xlab="Z-Score")
  abline(v=2,col="red")
  
  dev.off()
  
}

########## intramodule correlation ##############
eigmat<- MEs
ind=c(4,3,18,5,20,8,6)
eigmat<- eigmat[,ind]

mds = cmdscale(dist(t(eigmat)), eig = T);
pc1 = mds$eig[1]^2 / sum(mds$eig^2);
pc2 = mds$eig[2]^2 / sum(mds$eig^2)
c =colnames(eigmat)
cols=substring(c,3,30)

adj = bicor(eigmat)
library(igraph)
g1 <- graph.adjacency(as.matrix(adj),mode="undirected",weighted=T,diag=FALSE)
layoutFR <- mds$points

pdf("2.2_ModuleEigengene-MDS.pdf",height=8,width=8,useDingbats=FALSE)
edgecolors = numbers2colors(E(g1)$weight, colors = redWhiteGreen(10, gamma=4), signed=T, centered=T, lim=c(-1,1))
plot.igraph(g1, vertex.label = c,
            vertex.label.dist=1.2,
            vertex.size=6,
            vertex.label.color="black",
            vertex.label.family = "sans",
            vertex.color = cols,
            vertex.label.cex=0.6,
            layout=layoutFR,
            edge.color=edgecolors,
            edge.width=4,asp=1, main="Module Eigengene MDS")
labs = seq(1,-1,by=-.25)
str = paste(labs,sep="\n")
#text(-1.25,-1, labels = paste(labs,collapse='\n'),pos = 4,cex = .5)
p = matrix(NA,nrow=9,ncol=4)
p[,1]=-1.35; p[,2]=-1.3
p[,3]=p[,4] = -.6-.6*seq(0,1,by=.12)
for(i in 1:9) {
  lines(x=p[i,1:2],y=p[i,3:4], lwd = 2, col=numbers2colors(as.numeric(labs[i]), colors = redWhiteGreen(10, gamma=1), signed=T, centered=T, lim=c(-1,1)))
  text(x=-1.3,y=p[i,3],labs[i],cex=.4,pos=4)
}
dev.off()
################ ################ ################ ################ ################
library("anRichment")
GOcollection = buildGOcollection(organism = "mouse")
########## anErichment ################
########## GO Enrichment ##############
moduleC=outmm$Initially.Assigned.Module.Color
symbol=as.character(outmm$GeneSymbol);symbol
entrez = convert2entrez(organism = "mouse", symbol = symbol);# How many conversions were successful
table(is.finite(entrez))
#entrez<-entrez[!is.na(entrez)]
GOenrichment = enrichmentAnalysis(classLabels = moduleC, identifiers = entrez,refCollection = GOcollection,useBackground = "given",threshold = 1e-4,thresholdType = "Bonferroni",getOverlapEntrez = TRUE,getOverlapSymbols = TRUE,ignoreLabels = "grey");

collectGarbage()
names(GOenrichment)
names(GOenrichment$enrichmentTable)
table.display = GOenrichment$enrichmentTable;table.display$overlapGenes = shortenStrings(table.display$overlapGenes, maxLength = 70,split = "|");head(table.display);
write.csv(GOenrichment$enrichmentTable, file = "GOenrichment-enrichmentTable.csv",row.names = FALSE)
########## anErichment ################
########## GSEA Enrichment ############
internalColl = internalCollection(organism = "mouse")
knownGroups(internalColl, sortBy = "size")
dataSetNames(internalColl, groups = "StemCellLists.Lee")
ids = dataSetIDs(internalColl, groups = "BrainLists.MO")
dataSetNames(internalColl, dataSets = ids)
biosysCollection = BioSystemsCollection("mouse")
genomicPosCollection = genomicPositionCollection(organism = "mouse",spacings = 5e6,overlapFactor = 2)
HDSigCollection = HDSigDBCollection(organism = "mouse")
msdbColl = MSigDBCollection(file = "~/Dropbox (UCLA Mednet)/R codes/msigdb_v6.2.xml", organism = "mouse")


combinedCollection = mergeCollections(GOcollection,internalColl,biosysCollection,HDSigCollection,genomicPosCollection,msdbColl)
knownGroups(combinedCollection, sortBy = "size")##  [1] "GO"##  [2] "GO.BP"##  [3] "GO.MF"##  [4] "GO.CC"##  [5] "REACTOME"##  [6] "HDSigDB"##  [7] "KEGG"##  [8] "High-throughput association analysis"##  [9] "BIOCYC"## [10] "Differential expression induced by Htt CAG length expansion"## [11] "JAM"## [12] "Brain region marker enriched gene sets"## [13] "JA Miller"## [14] "Pathway Interaction Database"## [15] "Brain region markers"11


combinedEnrichment =  enrichmentAnalysis(classLabels = moduleC, identifiers = entrez,refCollection = combinedCollection,useBackground = "given",threshold = 1e-4,thresholdType = "Bonferroni")
head(combinedEnrichment$enrichmentTable[, -ncol(combinedEnrichment$enrichmentTable)])

write.csv(combinedEnrichment$enrichmentTable,"EnrichmentTable.csv")

#Creating enrichment labels for classes
eLabels = enrichmentLabels(combinedEnrichment$enrichmentTable,focusOnGroups = c("all", "GO", "Cell type markers", "Brain region markers", "HDSigDB"),groupShortNames = c("all", "GO", "CT", "BR", "HD"),minSize = 0.05,numericClassLabels = FALSE);
write.csv(eLabels,"Enrichment_labels_classes.csv")

################ ################ ################ ################ ################

