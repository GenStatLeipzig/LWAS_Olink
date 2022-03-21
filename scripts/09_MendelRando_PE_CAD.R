#' ---
#' title: "Mendelian Randomization: PE --> CAD"
#' subtitle: "LWAS Olink"
#' author: "Janne Pott"
#' date: "February 2022"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Intro ####
#' ***
#' I want to test for causality between protein levels and CAD for all proteins with genome-wide significant instruments.
#' 
#' Step 1: get best pQTL per protein as instrument
#' 
#' Step 2: run MR
#' 
#' Step 3: some plots
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../mySourceFile.R")
.libPaths()
setwd(projectpath_main)

#' # Load pQTL data ####
#' ***
load("../temp/08_OlinkGenes.RData")
load("../results/02_pQTLs_all_cis.RData")
load("../results/02_pQTLs_males_cis.RData")
load("../results/02_pQTLs_females_cis.RData")

pQTLs = rbind(pQTLs_a,pQTLs_m,pQTLs_f)
pQTLs[,GeneSetting:=paste(gene,setting,sep=":")]
pQTLs[,SNPGeneSetting:=paste(variant_id_hg38,gene,setting,sep=":")]
pQTLs[,SNPGene:=paste(variant_id_hg38,gene,sep=":")]
length(unique(pQTLs$GeneSetting))

myGenes = pQTLs[hierarch_fdr5proz ==T, unique(GeneSetting)]
pQTLs<-pQTLs[GeneSetting %in% myGenes,]

#' # Load CAD data ####
#' ***
#' Data is given in hg19
myTab_CAD<-fread(path_CAD)
head(myTab_CAD)

myTab_CAD[,rsID:=oldID]
myTab_CAD[,chr:=CHR]
myTab_CAD[,pos:=BP]
myTab_CAD[,effect_allele:=toupper(Allele1)]
myTab_CAD[,other_allele:=toupper(Allele2)]
myTab_CAD[,eaf:=Freq1]
myTab_CAD[,beta:=Effect]
myTab_CAD[,se:=StdErr]
myTab_CAD[,pval:=`P-value`]
myTab_CAD[,n_samples:=34541 + 261984 + 88192 + 162544]
myTab_CAD[,chr_pos_oa_ea:=paste(chr,pos,other_allele,effect_allele,sep=":")]
myTab_CAD[,chr_pos_ea_oa:=paste(chr,pos,effect_allele,other_allele,sep=":")]
myTab_CAD[,chr_pos:=paste(chr,pos,sep=":")]

myNames <- c("rsID","chr_pos_oa_ea","chr_pos_ea_oa","chr_pos","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","n_samples")
colsOut<-setdiff(colnames(myTab_CAD),myNames)
myTab_CAD[,get("colsOut"):=NULL]
setcolorder(myTab_CAD,myNames)
head(myTab_CAD)

#' # Match and Merge CAD and Olink data ####
#' ***
dummy = unlist(strsplit(pQTLs$variant_id_hg19,":"))
pos_hg19 = dummy[seq(2,length(dummy),4)]
pQTLs[,chrPos19:=paste(chr,pos_hg19,sep=":")]
pQTLs[,chrPos19:=gsub("chr","",chrPos19)]
table(is.element(myTab_CAD$chr_pos,pQTLs$chrPos19))
table(is.element(pQTLs$chrPos19,myTab_CAD$chr_pos))

myTab_CAD = myTab_CAD[chr_pos %in% pQTLs$chrPos19,]
table(is.element(myTab_CAD$chr_pos,pQTLs$chrPos19))
table(is.element(pQTLs$chrPos19,myTab_CAD$chr_pos))

data4 = copy(pQTLs)
data4 = data4[chrPos19 %in% myTab_CAD$chr_pos,]

matched = match(data4$chrPos19,myTab_CAD$chr_pos)
myTab_CAD = myTab_CAD[matched,]
table(myTab_CAD$chr_pos == data4$chrPos19)

table(myTab_CAD$effect_allele == data4$effect_allele,
      myTab_CAD$other_allele == data4$other_allele)
table(myTab_CAD$other_allele == data4$effect_allele,
      myTab_CAD$effect_allele == data4$other_allele)

#' I only want to keep those SNPs that have matching effect and other allele (227,217 + 180,672)
#' 
filt1 = myTab_CAD$effect_allele == data4$effect_allele & myTab_CAD$other_allele == data4$other_allele
filt2 = myTab_CAD$other_allele == data4$effect_allele & myTab_CAD$effect_allele == data4$other_allele
table(filt1)
table(filt2)
myTab_CAD = myTab_CAD[filt1 | filt2,]
data4 = data4[filt1 | filt2,]

table(myTab_CAD$effect_allele == data4$effect_allele,
      myTab_CAD$other_allele == data4$other_allele)
table(myTab_CAD$other_allele == data4$effect_allele,
      myTab_CAD$effect_allele == data4$other_allele)
filt3 = myTab_CAD$other_allele == data4$effect_allele

data4[,beta_CAD := myTab_CAD[,beta]]
data4[filt3,beta_CAD := beta_CAD*(-1)]
data4[,se_CAD:=myTab_CAD[,se]]
data4[,pval_CAD:=myTab_CAD[,pval]]
data4[,n_CAD:=myTab_CAD[,n_samples]]
data4[,eaf_CAD:=myTab_CAD[,eaf]]
data4[filt3,eaf_CAD:=1-eaf_CAD]

plot(data4$eaf,data4$eaf_CAD)

myMRDat_PE_CAD = copy(data4)
save(myMRDat_PE_CAD,file="../results/09_a_MR_matchedData.RData")

leadpQTLs<-myMRDat_PE_CAD[,.SD[pval==min(pval)],by=.(GeneSetting)]
leadpQTLs = leadpQTLs[pval<=5e-8,]
table(duplicated(leadpQTLs$GeneSetting))
leadpQTLs = leadpQTLs[!duplicated(GeneSetting),]

#' # MR Function ####
#' ***
myMRfunction_ratio<-function(betaX,seX,betaY,seY){
  betaIV<-betaY/betaX
  se.ratio.st<- seY/sqrt(betaX^2)
  se.ratio.nd<- sqrt(seY^2/betaX^2 + betaY^2*seX^2/betaX^4)
  p1<-2*pnorm(-abs(betaIV/se.ratio.st))
  p2<-2*pnorm(-abs(betaIV/se.ratio.nd))
  res<-data.table(
    beta_IV=betaIV,
    se_IV1=se.ratio.st,
    se_IV2=se.ratio.nd,
    p_IV1=p1,
    p_IV2=p2)
  return(res)
}

test = leadpQTLs[,myMRfunction_ratio(betaX = beta,seX = se,betaY = beta_CAD,seY = se_CAD)]

myDat = data.table(GeneSetting = leadpQTLs$GeneSetting,
                   gene = leadpQTLs$gene,
                   setting = leadpQTLs$setting,
                   exposure = paste0("PE_",leadpQTLs$gene),
                   outcome = "CAD",
                   lead_pQTL = leadpQTLs$variant_id_hg38, 
                   beta_MR=test$beta_IV,
                   se_MR=test$se_IV2, 
                   pval_MR=test$p_IV2)
myDat

#' # FDR ####
#' ***
myTab_MR_a = copy(myDat)
myTab_MR_a = myTab_MR_a[setting == "combined",]
myTab_MR_m = copy(myDat)
myTab_MR_m = myTab_MR_m[setting == "males",]
myTab_MR_f = copy(myDat)
myTab_MR_f = myTab_MR_f[setting == "females",]

myTab_MR_a[,FDR := p.adjust(pval_MR,method = "bonferroni")]
myTab_MR_m[,FDR := p.adjust(pval_MR,method = "bonferroni")]
myTab_MR_f[,FDR := p.adjust(pval_MR,method = "bonferroni")]

table(myTab_MR_a$FDR<=0.05)
table(myTab_MR_m$FDR<=0.05)
table(myTab_MR_f$FDR<=0.05)

myTab_MR<-rbind(myTab_MR_a,myTab_MR_m,myTab_MR_f)
myTab_MR[,table(setting, FDR<0.05)]
myTab_MR[FDR<0.05,]

#' change from long format to wide format
names(myTab_MR)
myTab_MR3<-dcast(myTab_MR, formula = gene ~ setting,
                 value.var = c("lead_pQTL","beta_MR","se_MR","pval_MR","FDR"),
                 sep = "_")
head(myTab_MR3)
myTab_MR3[,exposure := paste0("PE_",gene)]
myTab_MR3[,outcome := "CAD"]
myTab_MR3[,table(lead_pQTL_combined == lead_pQTL_females)]
myTab_MR3[,table(lead_pQTL_combined == lead_pQTL_males)]

myTab_MR3[,table(FDR_males<=0.05,FDR_females<=0.05,FDR_combined<=0.05)]
myTab_MR3[,table(pval_MR_males<=0.05,pval_MR_females<=0.05,pval_MR_combined<=0.05)]

#' ## Plots ####
#' I want a positive SNP effect on protein levels 
myPlotData = copy(leadpQTLs)
filt = myPlotData$beta < 0
table(filt)
myPlotData[filt,beta_CAD := beta_CAD *(-1)]
myPlotData[filt,beta := beta *(-1)]
myPlotData[filt,eaf := 1-eaf]
myPlotData[filt,eaf_CAD := 1-eaf_CAD]

table(myPlotData$GeneSetting == myTab_MR$GeneSetting)
myPlotData[,FDR_MR := myTab_MR$FDR]
myPlotData[,beta_MR := myTab_MR$beta_MR]
myPlotData[,se_MR := myTab_MR$se_MR]
myPlotData[,pval_MR := myTab_MR$pval_MR]

myPlotData[pval_MR>=0.05,sig:="no"]
myPlotData[pval_MR<0.05,sig:="yes"]
myPlotData[FDR_MR<0.05,sig:="yes (FDR 5%)"]
myPlotData[,gene2:=""]
myPlotData[FDR_MR<0.05,gene2:=gene]
myPlotData[,gene3:=""]
myPlotData[pval_MR<0.05,gene3:=gene]

myPlot = ggplot(myPlotData[sig!="no",], aes(x=beta, y=beta_CAD, color=gene3,label=gene2)) +
  facet_wrap(facets = vars(setting)) + 
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1.25)+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_point(data=myPlotData[sig=="no",], aes(x=beta, y=beta_CAD),col="black",size=1.5,alpha=0.5,shape = 16)+
  geom_point(size=3.5)+ 
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"))+
  labs(x="Effect on protein level", y = "Effect on CAD",
       color = "Genes with \n causal effect")+#, shape = "Type of \nsex interaction")+
  guides(size="none",fill="none")+
  scale_shape_manual(values=c(15, 17, 18,16))+
  geom_label_repel(data = subset(myPlotData, sig!="no"),# & type=="beta estimate"),
                   aes(fill=gene3),
                   box.padding = 1.15,
                   max.overlaps = 25,
                   show.legend=F,color="black")
myPlot

tiff(filename = "../results/09_ScatterPlot_MR_GE_CAD_allSettings.tif", 
     width = 1400, height = 800, res=125, compression = 'lzw')
myPlot
dev.off()

myPlot2 = ggplot(myPlotData[setting == "combined" & sig!="no",], aes(x=beta, y=beta_CAD, color=gene3,label=gene2)) +
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1.25)+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_point(data=myPlotData[setting == "combined" & sig=="no",], aes(x=beta, y=beta_CAD),
             col="black",size=1.5,alpha=0.5,shape = 16)+
  geom_point(size=3.5)+ 
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"))+
  labs(x="Effect on protein level", y = "Effect on CAD",
       color = "Genes with \n causal effect")+#, shape = "Type of \nsex interaction")+
  guides(size="none",fill="none")+
  scale_shape_manual(values=c(15, 17, 18,16))+
  geom_label_repel(data = subset(myPlotData, sig!="no" & setting == "combined"),# & type=="beta estimate"),
                   aes(fill=gene3),
                   box.padding = 1.15,
                   max.overlaps = 25,
                   show.legend=F,color="black")
myPlot2

tiff(filename = "../results/09_ScatterPlot_MR_GE_CAD_combined.tif", 
     width = 1200, height = 900, res=125, compression = 'lzw')
myPlot2
dev.off()

#' # Save ####
#' ***
myTab[,MR_PE_CAD:=NA]
myTab[Gene %in% myTab_MR_a$gene, MR_PE_CAD:=F]
myTab[Gene %in% myTab_MR_a[FDR<=0.05,gene], MR_PE_CAD:=T]
table(myTab$MR_PE_CAD)

write.table(myTab_MR,file = "../results/09_MR_PE_CAD_long.txt",
            col.names = T,row.names = F,quote = F,dec = ",",sep="\t")
write.table(myTab_MR3,file = "../results/09_MR_PE_CAD_combined.txt",
            col.names = T,row.names = F,quote = F,dec = ",",sep="\t")
save(myTab_MR,file="../results/09_MR_PE_CAD_long.RData")
save(myTab,file="../temp/09_OlinkGenes.RData")

#' # sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,formateTimediff(Sys.time()-time0))
