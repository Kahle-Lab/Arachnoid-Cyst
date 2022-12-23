library(lattice)
library(readxl)
setwd("~/Desktop/")
QQ <- read.table(file = 'QQ_Input_CHv3.txt',sep="\t",header=TRUE,stringsAsFactors = FALSE)

AC_Case_Control_Final <- read_excel("AC_Case-Control_Final.xlsx")

#QQ <- AC_Case_Control_Final
#OMIMINT <- read_excel("OMIMINT.xlsx")
#OMIM <- read_excel("OMIM.xlsx")
#INT <- read_excel("INT.xlsx")
#HBEINT <- read_excel("HBEINT.xlsx")
#HBE <- read_excel("HBE.xlsx")

qqunif.plot<-function(pvalues, 
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="fill", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ...);
           panel.abline(0,1);
           panel.abline(h=6.06, col = "red")
         }, par.settings=par.settings, ...
  )
}

pdf("CH_QQ_Syn.pdf")
qqunif.plot(QQ$syn_pValue, main = "Synonymous De Novo Variant Burden")
dev.off()

pdf("CH_QQ_MisD.pdf")
qqunif.plot(QQ$misD_pValue, main = "Deleterious Missense De Novo Variant Burden")
dev.off()

pdf("CH_QQ_lof.pdf")
qqunif.plot(QQ$lof_pValue, main = "Predictive Loss-of-Function De Novo Variant Burden")
dev.off()

pdf("CH_QQ_ProtA.pdf")
qqunif.plot(QQ$prot_pValue, main = "Protein-Altering De Novo Variant Burden")
dev.off()

pdf("CH_QQ_ProtD.pdf")
qqunif.plot(QQ$protD_pValue, main = "Protein-Damaging De Novo Variant Burden")
dev.off()

pdf("CH_QQ_All.pdf")
qqunif.plot(QQ$all_pValue, main = "De Novo Variant Burden (All Variants)")
dev.off()

##########
##########


pdf("DNR_QQ_Syn.pdf")
qqunif.plot(QQ$Syn, main = "DenovolyzeR Synonymous De Novo Variant Burden")
dev.off()

pdf("DNR_QQ_MisT.pdf")
qqunif.plot(QQ$MisT, main = "DenovolyzeR Tolerated-Missense De Novo Variant Burden")
dev.off()

pdf("DNR_QQ_MisD.pdf")
qqunif.plot(QQ$MisD, main = "DenovolyzeR Damaging-Missense De Novo Variant Burden")
dev.off()

pdf("DNR_QQ_LoF.pdf")
qqunif.plot(QQ$LoF, main = "DenovolyzeR Loss-of-Function De Novo Variant Burden")
dev.off()

pdf("DNR_QQ_prot.pdf")
qqunif.plot(QQ$Prot, main = "DenovolyzeR Protein-Altering De Novo Variant Burden")
dev.off()

pdf("DNR_QQ_protD.pdf")
qqunif.plot(QQ$ProtD, main = "DenovolyzeR Protein-Damaging De Novo Variant Burden")
dev.off()

pdf("DNR_QQ_All.pdf")
qqunif.plot(QQ$All, main = "DenovolyzeR De Novo Variant Burden")
dev.off()

pdf("DNW_QQ_Syn.pdf")
qqunif.plot(QQ$Syn, main = "DeNovoWEST Synonymous De Novo Variant Burden")
dev.off()

pdf("DNW_QQ_MisT.pdf")
qqunif.plot(QQ$MisT, main = "DeNovoWEST Tolerated-Missense De Novo Variant Burden")
dev.off()

pdf("DNW_QQ_MisD.pdf")
qqunif.plot(QQ$MisD, main = "DeNovoWEST Damaging-Missense De Novo Variant Burden")
dev.off()

pdf("DNW_QQ_LoF.pdf")
qqunif.plot(QQ$LoF, main = "DeNovoWEST Loss-of-Function De Novo Variant Burden")
dev.off()

pdf("DNW_QQ_prot.pdf")
qqunif.plot(QQ$Prot, main = "DeNovoWEST Protein-Altering De Novo Variant Burden")
dev.off()

pdf("DNW_QQ_protD.pdf")
qqunif.plot(QQ$ProtD, main = "DeNovoWEST Protein-Damaging De Novo Variant Burden")
dev.off()

pdf("DNW_QQ_All.pdf")
qqunif.plot(QQ$All, main = "DeNovoWEST De Novo Variant Burden")
dev.off()

pdf("QQ_HBE.pdf")
qqunif.plot(HBE$HBE)
dev.off()

pdf("QQ_INT.pdf")
qqunif.plot(INT$INT)
dev.off()

pdf("QQ_HBEINT.pdf")
qqunif.plot(HBEINT$HBEINT)
dev.off()

pdf("QQ_OMIMINT.pdf")
qqunif.plot(OMIMINT$OMIMINT)
dev.off()

pdf("QQ_OMIM.pdf")
qqunif.plot(OMIM$OMIM)
dev.off()





