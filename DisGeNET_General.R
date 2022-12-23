library(devtools)
#install_bitbucket("ibi_group/disgenet2r")
library(disgenet2r)
disgenet_api_key <- get_disgenet_api_key(
  email = "garrett.allington@yale.edu", 
  password = "YaleM2P2" )
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)
myListOfGenes <- c("CASK", "ARID1B", "KMT5B", "KANSL1", "SON", "KIFC2", "MMP7", 
                    "MRPL50", "PMFBP1", "PRKG1", "SLC22A17", "TMX2", "KCNQ3", 
                    "POLR2A", "CTNNA1", "CHD7", "RDX", "COL11A1", "COX7A2L", 
                    "ERN2", "GPAA1", "KCNH8", "MIB1", "PIK3R4", "SLC4A1", 
                    "SLC6A2", "USP9X", "SORCS3", "ESX1", "MMP25", "RARA", 
                    "CDK13", "CHD3", "MCPH1", "MED13L", "NFIX", "EEF1AKMT4", 
                    "KDM5B", "PSMD12", "PTPN11", "SLC25A25", "ABCA1", "ABCA7", 
                    "ACTA1", "ADAM28", "ADCY8", "AMT", "ARHGEF1", "ARPC3", 
                    "AZU1", "CACNA1E", "CHRNA4", "CLDN5", "CNTNAP1", "COL2A1", 
                    "COL3A1", "COL4A1", "CPXM2", "CTRB2", "DDR2", "DHX36", 
                    "ECHS1", "FBN2", "FLNB", "FOXD3", "FSD2", "GDAP1L1", "GFPT1", 
                    "GNAI1", "GPR135", "GSS", "HCN4", "HDAC8", "HTR3A", "IGFALS", 
                    "IGFN1", "IKZF2", "KAT8", "KCNA7", "KCND3", "KCNJ8", "KCNQ1", 
                    "KCNS2", "KDM4A", "KIF3B", "LONP1", "LRP1", "LRP2", "LRP3", 
                    "MAP3K13", "MAST1", "MFN2", "MGAT5", "NAA15", "NDUFV1", "NF1", 
                    "NKX2-6", "OTOF", "PDHA1", "PRDM15", "PSMC5", "RBPJL", "RPN1", 
                    "RPS10", "RPS10-NUDT3", "SCN1A", "SCN4A", "SEC23IP", "SLC15A4", 
                    "SLC20A2", "SLC4A10", "SLIT1", "SMAD4", "SMARCA2", "SMIM4", 
                    "SNF8", "SOX11", "SOX4", "TAAR5", "TBL1X", "TENM3", "TMEM165", 
                    "TNIK", "TRIO", "TRRAP", "TTC28", "VARS", "VDR", "ADNP", "BPTF", 
                    "TLK2", "MUC16", "SETBP1", "AHDC1", "ARID1A", "AUTS2", "CPED1", 
                    "CTNND1", "DLGAP1", "DST", "DUS1L", "EMSY", "IGF2BP1", "INO80D", 
                    "KAT6B", "KMT2C", "PHF23", "POLG2", "SECISBP2", "SMARCD1", 
                    "TANC2", "TLL2", "VEZF1", "ZNF254")
data2 <- gene2disease(
  gene     = myListOfGenes,
  score =c(0.4, 1),
  verbose  = TRUE
)
data2

plot( data2,
      class = "Network",
      prop = 10)

setwd('/Users/garrettallington/Desktop/')
pdf("DiseaseClass_Gene_Heatmap.pdf", height=10, width=20)
plot( data2,
      class="DiseaseClass", nchars=60)
dev.off()

  plot( data2,
        class  ="Heatmap",
        limit  = 100, nchars = 50 )







diseasesOfInterest <- c("C0039292",	"C0220685",	"C2745959",	"C0268338",	"C0220668")
  
data10 <- disease2disease_by_gene(
  disease = diseasesOfInterest,
  database = "CURATED",
  ndiseases =  5)
plot( data10,
      class = "Venn",
      ndiseases =  5)
plot( data10,
      class = "Diseasome",
      prop  = 0.1)


data11 <- disease2disease_by_gene(
  disease =  c("C1535926", "C1510586", "C0152021"),
  database = "CURATED",
  ndiseases = 3)

pdf("Disease_Venn.pdf", height=20, width=20)
plot( data11,
      class="Venn")
dev.off()

pdf("Disease_Network.pdf", height=20, width=20)
plot( data11,
      class = "Network",
      prop  = 0.1)
dev.off()

pdf("Diseasome.pdf", height=20, width=20)
plot( data10,
      class = "Diseasome",
      prop  = 0.1)
dev.off()

diseasesOfInterest <- c("C1510586", "C0152021", "C1535926")
data11 <- disease2disease_by_gene(
  disease =  diseasesOfInterest,
  database = "ALL",
  ndiseases = 75)
plot( data11,
      class="Venn"
)


