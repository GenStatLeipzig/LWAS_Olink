#' ---
#' title: "Mendelian Randomization: GE --> PE"
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
#' I want to test for causality between gene expression and its protein levels in all available tissues.
#' 
#' Step 1: get best eQTL per gene as instrument
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

#' # Load data ####
#' ***
load("../temp/07_OlinkGenes.RData")
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

#' # Set parameter ####
#' ***
datalist_trait_fn2<-dir("../temp/GTEx_v8_filtered/")
traitlist_trait2<-gsub(".RData","",datalist_trait_fn2)
traitlist_trait2<-gsub("_"," ",traitlist_trait2)

#' # Loop 1: get eQTL data ####
#' ***
dumTab<-foreach(i=1:length(traitlist_trait2))%do%{
  # dumTab<-foreach(i=1:5)%do%{
  #i=1
  trait2<-traitlist_trait2[i]
  message("\nWorking on trait ",trait2,"\n")
  
  # Load eQTL data
  loaded1<-load(paste0("../temp/GTEx_v8_filtered/",datalist_trait_fn2[i]))
  data2<- get(loaded1)
  
  # Filt eQTL data in relevant columns and rows
  filt<-!is.na(data2$beta) & !is.na(data2$se) & !is.na(data2$n_samples)
  data2<-data2[filt,]
  
  # Filt eQTL data for genomewide significant hits
  data2 = data2[pval <= 5e-8,]
  
  # Match eQTL data to pQTL data
  data4 = copy(pQTLs)
  data4 = data4[variant_id_hg38 %in% data2$variant_id,]
  table(is.element(data4$variant_id_hg38, data2$variant_id))
  table(is.element(data2$variant_id,data4$variant_id_hg38))
  
  matched = match(data4$variant_id_hg38, data2$variant_id)
  table(is.na(matched))
  data5 = data2[matched,]
  table(data5$variant_id==data4$variant_id_hg38)
  stopifnot(data5$effect_allele==data4$effect_allele)
  
  data4[,beta_GTEx := data5[,beta]]
  data4[,se_GTEx:=data5[,se]]
  data4[,pval_GTEx:=data5[,pval]]
  data4[,n_GTEx:=data5[,n_samples]]
  data4[,maf_GTEx:=data5[,maf]]
  data4[,tissue_GTEx:=data5[,tissue]]
  
  data4[,GeneSettingTissue:=paste(gene,setting,tissue_GTEx,sep=":")]
  data4
}

myMRDat_GE_PE = rbindlist(dumTab)
save(myMRDat_GE_PE,file="../results/08_a_MR_matchedData.RData")

leadeQTLs<-myMRDat_GE_PE[,.SD[pval_GTEx==min(pval_GTEx)],by=.(GeneSettingTissue)]
table(duplicated(leadeQTLs$GeneSettingTissue))
leadeQTLs = leadeQTLs[!duplicated(GeneSettingTissue),]

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

test = leadeQTLs[,myMRfunction_ratio(betaX = beta_GTEx,seX = se_GTEx,betaY = beta,seY = se)]

myDat = data.table(GeneSettingTissue = leadeQTLs$GeneSettingTissue,
                   gene = leadeQTLs$gene,
                   setting = leadeQTLs$setting,
                   tissue = leadeQTLs$tissue,
                   exposure = paste0("GE_",leadeQTLs$gene),
                   outcome = paste0("PE_",leadeQTLs$gene),
                   lead_eQTL = leadeQTLs$variant_id_hg38, 
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

HierFDR_Holger = function(data){
  # data = pQTLs_a
  
  data2 = copy(data)
  
  myFDR<- addHierarchFDR(pvalues = data2[,pval_MR], categs = data2[,gene],quiet = T)
  data2[,hierarch_fdr5proz:=myFDR$hierarch_fdr5proz]
  
  x1 = myFDR[,.(min_level1 = min(fdr_level1)), by = list(category,fdr_level2)]
  x1[,hierarch_fdr5proz:=fdr_level2<0.05]
  x0 = myFDR[hierarch_fdr5proz==T,.N,.(category,fdr_level2)]
  matched<-match(x1$category,x0$category)
  x1[,N:=x0$N[matched]]
  x1[is.na(N),N:=0]
  
  return(x1)
} 

myFDR_a<- addHierarchFDR(pvalues = myTab_MR_a[,pval_MR], categs = myTab_MR_a[,gene],quiet = F)
myFDR_m<- addHierarchFDR(pvalues = myTab_MR_m[,pval_MR], categs = myTab_MR_m[,gene],quiet = F)
myFDR_f<- addHierarchFDR(pvalues = myTab_MR_f[,pval_MR], categs = myTab_MR_f[,gene],quiet = F)

myTab_MR_a[,hierFDR:=myFDR_a$hierarch_fdr5proz]
myTab_MR_m[,hierFDR:=myFDR_m$hierarch_fdr5proz]
myTab_MR_f[,hierFDR:=myFDR_f$hierarch_fdr5proz]

x1 = HierFDR_Holger(data = myTab_MR_a)
x2 = HierFDR_Holger(data = myTab_MR_m)
x3 = HierFDR_Holger(data = myTab_MR_f)

dim(x1[hierarch_fdr5proz==T])
dim(x2[hierarch_fdr5proz==T])
dim(x3[hierarch_fdr5proz==T])


myGenes = unique(c(x1$category,x2$category,x3$category))
myFDR_short = data.table(category = myGenes)
setorder(myFDR_short,"category")

matched1 = match(myFDR_short$category,x1$category)
matched2 = match(myFDR_short$category,x2$category)
matched3 = match(myFDR_short$category,x3$category)

myFDR_short = cbind(myFDR_short,x1[matched1,-1],x2[matched2,-1],x3[matched3,-1])
names(myFDR_short)[2:13]<-c("FDRl2_a","SimesP_a","Sig_a","NSig_a",
                            "FDRl2_m","SimesP_m","Sig_m","NSig_m",
                            "FDRl2_f","SimesP_f","Sig_f","NSig_f")
head(myFDR_short)
dim(myFDR_short)

dummy<-myTab[Gene %in% myFDR_short$category,c(1,2,5)]  
stopifnot(dummy$Gene == myFDR_short$category)
myFDR2<-cbind(dummy,myFDR_short[,-1])
head(myFDR2)

myFDR2[,table(Sig_m,Sig_f,Sig_a)]
myFDR2[Sig_f==F & Sig_a==F & Sig_m==T,] 
myFDR2[Sig_f==F & Sig_a==T & Sig_m==T,] 
myFDR2[Sig_f==T & Sig_a==T & Sig_m==F,] 

myTab_MR<-rbind(myTab_MR_a,myTab_MR_m,myTab_MR_f)
myTab_MR[,table(setting, hierFDR)]

#' change from long format to wide format
names(myTab_MR)
myTab_MR[,dumID := paste(gene,tissue,sep="::")]
myTab_MR3<-dcast(myTab_MR, formula = dumID ~ setting,
                 value.var = c("lead_eQTL","beta_MR","se_MR","pval_MR","hierFDR"),
                 sep = "_")
head(myTab_MR3)
myTab_MR3[,gene := gsub("::.*","",dumID)]
myTab_MR3[,tissue := gsub(".*::","",dumID)]
myTab_MR3[,exposure := paste0("GE_",gene)]
myTab_MR3[,outcome := paste0("PE_","",gene)]
myTab_MR3[,table(lead_eQTL_combined == lead_eQTL_females)]
myTab_MR3[,table(lead_eQTL_combined == lead_eQTL_males)]
myTab_MR3[,lead_eQTL := lead_eQTL_combined]
myTab_MR3[is.na(lead_eQTL_combined),lead_eQTL := lead_eQTL_males]
table(is.na(myTab_MR3$lead_eQTL))
myTab_MR3[,lead_eQTL_combined := NULL]
myTab_MR3[,lead_eQTL_females := NULL]
myTab_MR3[,lead_eQTL_males := NULL]

myTab_MR3[,table(hierFDR_males,hierFDR_females,hierFDR_combined)]

#' ## Plots ####
myPlotMR<-ggplot(myTab_MR3[hierFDR_females==T | hierFDR_males==T,], aes(x=beta_MR_males, y=beta_MR_females,col=gene)) +
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1.25)+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_point(size=2.5)+ 
  theme_bw(base_size = 10)+
  ggtitle(label = "Scatter plot of Mendelian Randomization (hier. FDR ==T for at least one strata)")+
  theme(plot.title = element_text(hjust = 0.5))+ 
  xlab("male")+
  ylab("female")+
  guides(size="none",fill="none",col="none")
myPlotMR


#' # Save ####
#' ***
myTab[,MR_GE_PE:=NA]
myTab[Gene %in% myTab_MR_a$gene, MR_GE_PE:=F]
myTab[Gene %in% myTab_MR_a[hierFDR==T,gene], MR_GE_PE:=T]
table(myTab$MR_GE_PE)

write.table(myTab_MR,file = "../results/08_MR_GE_PE_long.txt",
            col.names = T,row.names = F,quote = F,dec = ",",sep="\t")
write.table(myTab_MR3,file = "../results/08_MR_GE_PE_combined.txt",
            col.names = T,row.names = F,quote = F,dec = ",",sep="\t")
save(myTab_MR,file="../results/08_MR_GE_PE_long.RData")
save(myTab,file="../temp/08_OlinkGenes.RData")

#' # sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,formateTimediff(Sys.time()-time0))
