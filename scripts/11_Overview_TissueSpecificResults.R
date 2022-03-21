#' ---
#' title: "Combine tissue-specific data"
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
#' I want to combine all tissue-specific data (eQTLs, Coloc, MetaXcan, MR).
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
load("../results/04_c_eQTL_lookup_Summary_LD.RData")
load("../results/06_Coloc_Genewise.RData")
load("../results/07_MetaXcan_long.RData")
load("../results/08_MR_GE_PE_long.RData")

load("../temp/09_OlinkGenes.RData")

#' # Filter ####
#' ***
#' I only want the combined setting, only GTEx tissues, only proteins with significant effect in LWAS
#' 
myTab = myTab[GWAS_sig ==T,]

leadpQTLs4 = leadpQTLs4[setting == "combined",]
leadpQTLs4 = leadpQTLs4[tissue_GTEx != "Whole_Blood_LIFE",]
leadpQTLs4 = leadpQTLs4[gene %in% myTab$Gene,]

Coloc = Coloc_genewise[setting == "combined",]
Coloc<-Coloc[comment=="good coloc",]
Coloc = Coloc[gene %in% myTab$Gene,]

MetaXcan2 = MetaXcan2[setting == "combined",]
MetaXcan2 = MetaXcan2[gene %in% myTab$Gene,]

myTab_MR = myTab_MR[setting == "combined",]
myTab_MR = myTab_MR[gene %in% myTab$Gene,]

#' # Get unique ID ####
#' ***
#' I use gene::tissue as ID
leadpQTLs4[,gene_tissue:=paste(gene,tissue_GTEx,sep="::")]
table(duplicated(leadpQTLs4$gene_tissue))

Coloc[,tissue:=gsub("GE in ","",trait2)]
Coloc[,tissue:=gsub(" ","_",tissue)]
Coloc[,gene_tissue:=paste(gene,tissue,sep="::")]
table(duplicated(Coloc$gene_tissue))

MetaXcan2[,gene_tissue:=paste(gene,tissue,sep="::")]
table(duplicated(MetaXcan2$gene_tissue))

myTab_MR[,gene_tissue:=paste(gene,tissue,sep="::")]
table(duplicated(myTab_MR$gene_tissue))

myID<-unique(c(leadpQTLs4$gene_tissue,Coloc$gene_tissue, MetaXcan2$gene_tissue, myTab_MR$gene_tissue))
dummy<-unlist(strsplit(myID,"::"))
myData<-data.table(myID=myID,
                   gene=dummy[seq(1,length(myID)*2,2)],
                   tissue=dummy[seq(2,length(myID)*2,2)])
setorder(myData,myID)
myData

#' # Merge info ####
#' ***
#' * Protein, Cyto, Gene, Tissue
matched1<-match(myData$gene,myTab$Gene)
table(is.na(matched1))
myData[,protein:=myTab[matched1,Parameter]]
myData[,cyto:=myTab[matched1,Cytoband]]

#' # Merge best eQTL ####
#' ***
#' * lead pQTL, sameDirpQTLs, best eQTL, effect, pval, tissue, LD r2
leadpQTLs4[,same_direction2:=NA]
leadpQTLs4[sign(beta) ==sign(beta_GTEx),same_direction2:=T]
leadpQTLs4[sign(beta) !=sign(beta_GTEx),same_direction2:=F]
leadpQTLs4[pval_GTEx>0.05,same_direction2:=NA]
table(leadpQTLs4$same_direction2, leadpQTLs4$pval_GTEx<0.05)
table(leadpQTLs4$same_direction, leadpQTLs4$pval_GTEx<0.05)

matched2<-match(myData$myID,leadpQTLs4$gene_tissue)
table(is.na(matched2))
myData[,eQTLs_leadpQTL_hg19:=leadpQTLs4[matched2,variant_id_hg19]]
myData[,eQTLs_leadpQTL_hg38:=leadpQTLs4[matched2,variant_id_hg38]]
myData[,eQTLs_sameDir:=leadpQTLs4[matched2,same_direction2]]
myData[,eQTLs_besteQTL_hg19:=leadpQTLs4[matched2,SNPID_GTEx_best_hg19]]
myData[,eQTLs_besteQTL_hg38:=leadpQTLs4[matched2,SNPID_GTEx_best]]
myData[,eQTLs_besteQTL_LDr2:=leadpQTLs4[matched2,LD_r2]]

#' # Merge Coloc ####
#' ***
#' * Coloc (nSNPs, H0_H4)
matched3<-match(myData$myID,Coloc$gene_tissue)
table(is.na(matched3))
myData[,Coloc_nsnps:=Coloc[matched3,nsnps]]
myData[,Coloc_PPH3:=Coloc[matched3,PP.H3.abf]]
myData[,Coloc_PPH4:=Coloc[matched3,PP.H4.abf]]

#' # Merge MetaXcan ####
#' ***
#' * MetaXcan (effect size, pval, hierFDR, nSNPs)
matched4<-match(myData$myID,MetaXcan2$gene_tissue)
table(is.na(matched4))
myData[,MetaX_nsnps:=MetaXcan2[matched4,n_snps_used]]
myData[,MetaX_beta:=MetaXcan2[matched4,effect_size]]
myData[,MetaX_pval:=MetaXcan2[matched4,pval]]
myData[,MetaX_hierFDR:=MetaXcan2[matched4,hierFDR]]

#' # Merge MR ####
#' ***
#' * MR (beta, se, pval, hierFDR)
matched5<-match(myData$myID,myTab_MR$gene_tissue)
table(is.na(matched5))
myData[,MRGEP_beta:=myTab_MR[matched5,beta_MR]]
myData[,MRGEP_se:=myTab_MR[matched5,se_MR]]
myData[,MRGEP_pval:=myTab_MR[matched5,pval_MR]]
myData[,MRGEP_hierFDR:=myTab_MR[matched5,hierFDR]]

#' # Take a look ####
#' ***
myData[MetaX_hierFDR==T,table(eQTLs_sameDir,MetaX_beta>0)]
myData[MRGEP_hierFDR==T,table(eQTLs_sameDir,MRGEP_beta>0)]

myData[MetaX_hierFDR==T,table(Coloc_PPH4>0.75,MetaX_beta>0)]
myData[MetaX_hierFDR==T,table(Coloc_PPH3>0.75,MetaX_beta>0)]

myData[MRGEP_hierFDR==T,table(Coloc_PPH4>0.75,MRGEP_beta>0)]
myData[MRGEP_hierFDR==T,table(Coloc_PPH3>0.75,MRGEP_beta>0)]

myData[,table(Coloc_PPH4>0.75,eQTLs_besteQTL_LDr2>0.1)]
myData[,table(Coloc_PPH3>0.75,eQTLs_besteQTL_LDr2>0.1)]

#' # Save ####
#' ***
save(myData,file="../results/11_allResults_tissueSpec.RData")
write.table(myData,file="../results/11_allResults_tissueSpec.txt",col.names = T, row.names = F,quote = F,sep="\t", dec=",")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in hours): " ,round(difftime(Sys.time(), time0, tz,units = "hours"),2))
