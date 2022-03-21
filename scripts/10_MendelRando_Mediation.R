#' ---
#' title: "Mendelian Randomization: GE --> PE --> CAD"
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
#' I want to test for mediation between gene expression levels and CAD by protein levels for four proteins with significant causal effect on CAD.
#' 
#' Step 1: get best eQTL per protein as instrument
#' 
#' Step 2: run MR & mediation
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
load("../temp/09_OlinkGenes.RData")
load("../results/08_a_MR_matchedData.RData")
load("../results/09_a_MR_matchedData.RData")

#' # Match CAD data ####
#' ***
table(is.element(myMRDat_PE_CAD$variant_id_hg38,myMRDat_GE_PE$variant_id_hg38))
table(is.element(myMRDat_GE_PE$variant_id_hg38,myMRDat_PE_CAD$variant_id_hg38))
data4 = copy(myMRDat_GE_PE)
data4 = data4[variant_id_hg38 %in% myMRDat_PE_CAD$variant_id_hg38,]

matched = match(data4$variant_id_hg38,myMRDat_PE_CAD$variant_id_hg38)
myMRDat_PE_CAD = myMRDat_PE_CAD[matched,]
table(myMRDat_PE_CAD$variant_id_hg38 == data4$variant_id_hg38)

table(myMRDat_PE_CAD$effect_allele == data4$effect_allele,
      myMRDat_PE_CAD$other_allele == data4$other_allele)
table(myMRDat_PE_CAD$other_allele == data4$effect_allele,
      myMRDat_PE_CAD$effect_allele == data4$other_allele)

data4[,beta_CAD := myMRDat_PE_CAD[,beta_CAD]]
data4[,se_CAD:=myMRDat_PE_CAD[,se_CAD]]
data4[,pval_CAD:=myMRDat_PE_CAD[,pval_CAD]]
data4[,n_CAD:=myMRDat_PE_CAD[,n_CAD]]
data4[,eaf_CAD:=myMRDat_PE_CAD[,eaf_CAD]]

plot(data4$eaf,data4$eaf_CAD)
names(data4)

myMRDat_GE_CAD = copy(data4)
save(myMRDat_GE_CAD,file="../results/10_a_MR_matchedData.RData")

myCandidates = c("PCSK9:combined","AXL:combined","IL6R:combined","TFPI:combined",
                 "PCSK9:males","IL6R:males","IL6R:females","TFPI:females")

leadeQTLs<-myMRDat_GE_CAD[,.SD[pval_GTEx==min(pval_GTEx)],by=.(GeneSettingTissue)]
leadeQTLs = leadeQTLs[is.element(GeneSetting,myCandidates),]
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

test = leadeQTLs[,myMRfunction_ratio(betaX = beta_GTEx,seX = se_GTEx,betaY = beta_CAD,seY = se_CAD)]

myDat = data.table(GeneSettingTissue = leadeQTLs$GeneSettingTissue,
                   gene = leadeQTLs$gene,
                   setting = leadeQTLs$setting,
                   tissue = leadeQTLs$tissue,
                   exposure = paste0("GE_",leadeQTLs$gene),
                   outcome = "CAD",
                   lead_eQTL = leadeQTLs$variant_id_hg38, 
                   beta_MR=test$beta_IV,
                   se_MR=test$se_IV2, 
                   pval_MR=test$p_IV2)
head(myDat)

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

myTab_MR<-rbind(myTab_MR_a,myTab_MR_m,myTab_MR_f)
myTab_MR[, table(hierFDR,gene,setting)]

#' # Check Mediation ####
#' *** 
#' * myTab_MR_GE_CAD: MR results for GE --> CAD (tissue dependent)
#' * myTab_MR_PE_CAD: MR results for PE --> CAD (tissue independent)
#' * myTab_MR_GE_PE: MR results for GE --> PE (tissue dependent)

myTab_MR_GE_CAD = copy(myTab_MR)
load("../results/08_MR_GE_PE_long.RData")
myTab_MR_GE_PE = copy(myTab_MR)
load("../results/09_MR_PE_CAD_long.RData")
myTab_MR_PE_CAD = copy(myTab_MR)

#' Filter for relevant settings
#' 
myTab_MR_GE_PE = myTab_MR_GE_PE[is.element(GeneSettingTissue,myTab_MR_GE_CAD$GeneSettingTissue)]
myTab_MR_PE_CAD = myTab_MR_PE_CAD[is.element(GeneSetting,myCandidates)]

myTab_MR_mediation = copy(myTab_MR_GE_CAD)
myTab_MR_mediation = myTab_MR_mediation[,c(1:4)]

#' Merge
#'
matched = match(myTab_MR_mediation$GeneSettingTissue,myTab_MR_GE_PE$GeneSettingTissue)
table(matched == 1:52)
myTab_MR_mediation[,lead_GE_PE := myTab_MR_GE_PE$lead_eQTL]
myTab_MR_mediation[,beta_GE_PE := myTab_MR_GE_PE$beta_MR]
myTab_MR_mediation[,se_GE_PE := myTab_MR_GE_PE$se_MR]
myTab_MR_mediation[,pval_GE_PE := myTab_MR_GE_PE$pval_MR]
myTab_MR_mediation[,hierFDR_GE_PE := myTab_MR_GE_PE$hierFDR]

myTab_MR_mediation[,lead_GE_CAD := myTab_MR_GE_CAD$lead_eQTL]
myTab_MR_mediation[,beta_GE_CAD := myTab_MR_GE_CAD$beta_MR]
myTab_MR_mediation[,se_GE_CAD := myTab_MR_GE_CAD$se_MR]
myTab_MR_mediation[,pval_GE_CAD := myTab_MR_GE_CAD$pval_MR]
myTab_MR_mediation[,hierFDR_GE_CAD := myTab_MR_GE_CAD$hierFDR]

dumID = paste(myTab_MR_mediation$gene,myTab_MR_mediation$setting,sep=":")
matched = match(dumID,myTab_MR_PE_CAD$GeneSetting)
myTab_MR_mediation[,lead_PE_CAD := myTab_MR_PE_CAD$lead_pQTL[matched]]
myTab_MR_mediation[,beta_PE_CAD := myTab_MR_PE_CAD$beta_MR[matched]]
myTab_MR_mediation[,se_PE_CAD := myTab_MR_PE_CAD$se_MR[matched]]
myTab_MR_mediation[,pval_PE_CAD := myTab_MR_PE_CAD$pval_MR[matched]]
myTab_MR_mediation[,FDR_PE_CAD := myTab_MR_PE_CAD$FDR[matched]]

head(myTab_MR_mediation)

#' Check 1: sig causal effect of GE on PE?
table(myTab_MR_mediation$hierFDR_GE_PE)
myTab_MR_mediation2 = copy(myTab_MR_mediation)
myTab_MR_mediation2 = myTab_MR_mediation2[hierFDR_GE_PE==T,]

#' Check 2: sig causal effect of PE on CAD? 
table(myTab_MR_mediation2$FDR_PE_CAD<0.05)

#' Check 3: sig causal effect of GE on CAD?
table(myTab_MR_mediation2$hierFDR_GE_CAD)

#' Check 4: sig (in)direct effect?
interactionTest  = function(mean1, se1, mean2, se2) { ## interaction test see http://www.bmj.com/content/326/7382/219.long
  meandiff_se = sqrt(se1^2 + se2^2)
  meandiff = mean2 - mean1
  meandiff_cilow = meandiff - 1.96*meandiff_se
  meandiff_cihigh = meandiff + 1.96*meandiff_se
  meandiff_z = meandiff/meandiff_se
  meandiff_p = pnorm(abs(meandiff_z), lower.tail = F)*2
  if(meandiff_p>1) meandiff_p = 1
  data.table::data.table(mean1, se1, mean2, se2, meandiff, meandiff_se, meandiff_cilow, meandiff_cihigh, meandiff_z, meandiff_p)
}

alpha<-myTab_MR_mediation2[,beta_GE_PE]
SE_alpha<-myTab_MR_mediation2[,se_GE_PE]
beta<-myTab_MR_mediation2[,beta_PE_CAD]
SE_beta<-myTab_MR_mediation2[,se_PE_CAD]

indir<-alpha*beta
SE_indir <- sqrt(alpha^2 * SE_beta^2 + beta^2 * SE_alpha^2 )
p_indir <- 2*pnorm(q = abs(indir/SE_indir), mean = 0, sd = 1, lower.tail = F)

myTab_MR_mediation2[,beta_GE_CAD_indirect:=indir]
myTab_MR_mediation2[,se_GE_CAD_indirect:=SE_indir]
myTab_MR_mediation2[,pval_GE_CAD_indirect:=p_indir]

test2 = myTab_MR_mediation2[,interactionTest(mean1 = beta_GE_CAD_indirect,
                                             se1 = se_GE_CAD_indirect,
                                             mean2 = beta_GE_CAD, 
                                             se2 = se_GE_CAD)]

myTab_MR_mediation2[,beta_GE_CAD_direct:=test2$meandiff]
myTab_MR_mediation2[,se_GE_CAD_direct:=test2$meandiff_se]
myTab_MR_mediation2[,pval_GE_CAD_direct:=test2$meandiff_p]

table(myTab_MR_mediation2$pval_GE_CAD_indirect<0.05)
table(myTab_MR_mediation2$pval_GE_CAD_direct<0.05,myTab_MR_mediation2$hierFDR_GE_CAD==T)

myTab_MR_a = copy(myTab_MR_mediation2)
myTab_MR_a = myTab_MR_a[setting == "combined",]
myTab_MR_m = copy(myTab_MR_mediation2)
myTab_MR_m = myTab_MR_m[setting == "males",]
myTab_MR_f = copy(myTab_MR_mediation2)
myTab_MR_f = myTab_MR_f[setting == "females",]

myFDR_a1<- addHierarchFDR(pvalues = myTab_MR_a[,pval_GE_CAD_direct], 
                          categs = myTab_MR_a[,gene],quiet = F)
myFDR_m1<- addHierarchFDR(pvalues = myTab_MR_m[,pval_GE_CAD_direct], 
                          categs = myTab_MR_m[,gene],quiet = F)
myFDR_f1<- addHierarchFDR(pvalues = myTab_MR_f[,pval_GE_CAD_direct], 
                          categs = myTab_MR_f[,gene],quiet = F)
myFDR_a2<- addHierarchFDR(pvalues = myTab_MR_a[,pval_GE_CAD_indirect], 
                          categs = myTab_MR_a[,gene],quiet = F)
myFDR_m2<- addHierarchFDR(pvalues = myTab_MR_m[,pval_GE_CAD_indirect], 
                          categs = myTab_MR_m[,gene],quiet = F)
myFDR_f2<- addHierarchFDR(pvalues = myTab_MR_f[,pval_GE_CAD_indirect], 
                          categs = myTab_MR_f[,gene],quiet = F)

myTab_MR_a[,hierFDR_GE_CAD_direct:=myFDR_a1$hierarch_fdr5proz]
myTab_MR_m[,hierFDR_GE_CAD_direct:=myFDR_m1$hierarch_fdr5proz]
myTab_MR_f[,hierFDR_GE_CAD_direct:=myFDR_f1$hierarch_fdr5proz]
myTab_MR_a[,hierFDR_GE_CAD_indirect:=myFDR_a2$hierarch_fdr5proz]
myTab_MR_m[,hierFDR_GE_CAD_indirect:=myFDR_m2$hierarch_fdr5proz]
myTab_MR_f[,hierFDR_GE_CAD_indirect:=myFDR_f2$hierarch_fdr5proz]

myTab_MR_mediation3 = rbind(myTab_MR_a,myTab_MR_m,myTab_MR_f)
table(myTab_MR_mediation3$hierFDR_GE_CAD_direct,
      myTab_MR_mediation3$hierFDR_GE_CAD_indirect)


#' # Save ####
#' ***
write.table(myTab_MR_mediation3,file = "../results/10_MR_Mediation.txt",
            col.names = T,row.names = F,quote = F,dec = ",",sep="\t")
save(myTab_MR_mediation3,file="../results/10_MR_Mediation.RData")


#' # sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,formateTimediff(Sys.time()-time0))








