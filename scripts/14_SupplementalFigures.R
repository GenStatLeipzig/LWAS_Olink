#' ---
#' title: "Supplemental Figures"
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
#' I want to create the supplemental figures S5 and S7. 
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../mySourceFile.R")
.libPaths()
setwd(projectpath_main)

#' # Load ####
#' ***
load("../results/11_allResults_tissueSpec.RData")
head(myData)

#' # Prep ####
#' ***
genes_lowLD<-myData[!is.na(eQTLs_besteQTL_LDr2) & eQTLs_besteQTL_LDr2<0.1,myID]
genes_Coloc3<-myData[!is.na(Coloc_PPH3) & Coloc_PPH3>=0.75,myID]
genes_Coloc4<-myData[!is.na(Coloc_PPH4) & Coloc_PPH4>=0.75,myID]
genes_MetaX_neg<-myData[!is.na(MetaX_hierFDR) & MetaX_hierFDR==T & MetaX_beta<0,myID]
genes_MetaX_pos<-myData[!is.na(MetaX_hierFDR) & MetaX_hierFDR==T & MetaX_beta>0,myID]

#' # Plot 1: Venn Diagramm (SupFig S5) ####
#' ***
#' 
dummy1<-venn4(genes_lowLD, genes_Coloc3, genes_Coloc4, genes_MetaX_pos,
              mytitle = "Venn Diagram Olink (pos MX)",
              mylabels = c("low LD","indep.","shared","pos. MX"))
dummy2<-venn4(genes_lowLD, genes_Coloc3, genes_Coloc4, genes_MetaX_neg,
              mytitle = "Venn Diagram Olink (neg MX)",
              mylabels = c("low LD","indep.","shared","neg. MX"))

#' # Plot 2: Euler Diagramm (SupFig S7) ####
#' ***
myDat5<-myData[!is.na(MetaX_hierFDR) & MetaX_hierFDR==T,.(myID,MetaX_beta)]
myDat6<-myData[!is.na(MRGEP_hierFDR) & MRGEP_hierFDR==T,.(myID,MRGEP_beta)]

names(myDat5)[2]<-"beta"
names(myDat6)[2]<-"beta"

makeEuler(mytab1 = myDat5 , mytab2 = myDat6 , 
          object_colname = "myID",
          direction_colname = "beta",
          titletext = "Effect MetaXcan - MR",
          legendtext = c("MetaXcan","MR"))

myDat7 = myDat5[myID %in% myDat6$myID]
myDat8 = myDat6[myID %in% myDat7$myID]

stopifnot(myDat7$myID == myDat8$myID)
table(myDat7$beta>0, myDat8$beta>0)
table(myDat7$beta>0)
table(myDat8$beta>0)

#' There are 7 GE-PE pairs with positive causality but negative association.
#' 
#' There are 14 GE-PE pairs with positive association but negative causality.
#' 
#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))

