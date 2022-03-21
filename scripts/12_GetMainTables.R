#' ---
#' title: "Create Main Tables"
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
#' I want to create the main tables. 
#' 
#' Of course, I cannot simply but the .txt files into the Word manuscript (a bit manual editing will be necessary). But I can extract the main findings. 
#' 
#' Again, I only do this for the combined setting. 
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
#' # Filter ####
#' ***
#' I only want the combined setting, only GTEx tissues, only proteins with significant effect in LWAS
#' 
myTab = myTab[GWAS_sig ==T,]

leadpQTLs4 = leadpQTLs4[setting == "combined",]
leadpQTLs3 = leadpQTLs4[tissue_GTEx == "Whole_Blood_LIFE",]
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

#' # Get Tables ####
#' ***
#' ## eQTLs ####
z1<-leadpQTLs4[pval_GTEx<=0.05 & same_direction==F & !is.na(same_direction),.N,by=gene]
z2<-leadpQTLs4[pval_GTEx<=0.05 & same_direction==T & !is.na(same_direction),.N,by=gene]
table(is.element(z1$gene,z2$gene))
z3<-unique(c(z1$gene,z2$gene))
z3<-data.table(gene=z3)
z3[,N_same:=z2[match(z3$gene,gene),N]]
z3[,N_diff:=z1[match(z3$gene,gene),N]]
z3[is.na(N_same),N_same:=0]
z3[is.na(N_diff),N_diff:=0]
z3[,N:=N_same + N_diff]
z3[,ratio:=N_diff/N]
print(z3[ratio>=0.75,])
message("     Found ",dim(z3[ratio>=0.75,])[1]," proteins with different effect direction compared to eQTLs ...")
z4<-z3[ratio>=0.75,]

#' ## MetaXcan ####
y1<-MetaXcan2[hierFDR==T & effect_size<0,.N,by=gene]
y2<-MetaXcan2[hierFDR==T & effect_size>0,.N,by=gene]
table(is.element(y1$gene,y2$gene))
y3<-unique(c(y1$gene,y2$gene))
y3<-data.table(gene=y3)
y3[,N_same:=y2[match(y3$gene,gene),N]]
y3[,N_diff:=y1[match(y3$gene,gene),N]]
y3[is.na(N_same),N_same:=0]
y3[is.na(N_diff),N_diff:=0]
y3[,N:=N_same + N_diff]
y3[,ratio:=N_diff/N]
print(y3[ratio>=0.75,])
message("     Found ",dim(y3[ratio>=0.75,])[1]," proteins with mainly negative associations in MetaXcan ...")
y4<-y3[ratio>=0.75,]

#' ## MendelRando ####
x1<-myTab_MR[hierFDR==T & beta_MR<0,.N,by=gene]
x2<-myTab_MR[hierFDR==T & beta_MR>0,.N,by=gene]
table(is.element(x1$gene,x2$gene))
x3<-unique(c(x1$gene,x2$gene))
x3<-data.table(gene=x3)
x3[,N_same:=x2[match(x3$gene,gene),N]]
x3[,N_diff:=x1[match(x3$gene,gene),N]]
x3[is.na(N_same),N_same:=0]
x3[is.na(N_diff),N_diff:=0]
x3[,N:=N_same + N_diff]
x3[,ratio:=N_diff/N]
print(x3[ratio>=0.75,])
message("     Found ",dim(x3[ratio>=0.75,])[1]," proteins with mainly negative causality in Mendelian Randomization ...")
x4<-x3[ratio>=0.75,]

#' # Venn (Figure 3) ####
#' Will be a template for Figure 3. That Figure is created in PowerPoint, but I add the following genes stored in variable *dummy*. The function *venn3* originates from tha package *toolboxH*. 

dummy<-venn3(z4$gene,y4$gene,x4$gene, 
               mytitle = paste0("Venn Diagram Olink - pred. neg. link - combined"),
               mylabels = c("QTLs","MetaXcan","MR"))

#' Negative links in all three approaches:
dummy$q1
#' Different effect direction & negative causality:
dummy$q2
#' Negative association & negative causality:
dummy$q3
#' Different effect direction & negative association:
dummy$q4
#' Only different effect direction:
dummy$q5
#' Only negative association:
dummy$q6
#' Only negative causality:
dummy$q7


#' # Get Main Tables ####
#' ***
#' 
#' ## Table 1: Basic sample description ###
#' Taken from Master Thesis from Tarcyane Garcia (co-author of paper)
#' 
#' Required columns:
#' 
#' * Age (years), BMI (kg/m^2), Current smoker, Hypertension, Type 2 diabetes, Statin therapy
#' * TC (mmol/l), LDL-C (mmol/l), HDL-C (mmol/l)
#' 

tab1 = fread("../exampleData/BasicSampleDescription.txt")
tab1

#' ## Table 2: eQTLs vs pQTLs (different effect direction) ###
#' 
#' Required columns:
#' 
#' * Gene (diff/tot)	
#' * SNP-ID (hg38) 
#' * effect allele
#' * effect allele frequency
#' * beta & p-value as pQTL
#' * beta & p-value as eQTL in whole blood
#' * beta & p-value as eQTL in best tissue
#' * replicated in LIFE? --> TRUE == effect in LIFE eQTL data is also in different direction compared to pQTL
#' 
TissuesAbr = fread("../temp/TissuesGTExv8.txt")
myGenes = z4$gene

leadpQTLs5 = copy(leadpQTLs4)
leadpQTLs5 = leadpQTLs5[is.element(gene,myGenes),]

leadpQTLs6<-leadpQTLs5[same_direction==F & tissue_GTEx != "Whole_Blood",.SD[pval_GTEx==min(pval_GTEx)],by=.(gene)]
leadpQTLs7<-leadpQTLs5[tissue_GTEx == "Whole_Blood",]
leadpQTLs8<-leadpQTLs3[setting == "combined" & is.element(gene,myGenes) & pval_GTEx<0.05,]

table(leadpQTLs6$gene == z4$gene)
matched = match(leadpQTLs6$gene, z4$gene)
z4 = z4[matched,]
table(leadpQTLs6$gene == z4$gene)

tab2 = data.table(gene = leadpQTLs6$gene,
                  ratio = paste0(z4$N_diff,"/",z4$N_diff + z4$N_same),
                  SNP = leadpQTLs6$rs_id,
                  ea = leadpQTLs6$effect_allele,
                  eaf = round(leadpQTLs6$eaf,3),
                  beta_pQTL = round(leadpQTLs6$beta,3),
                  pval_pQTL = signif(leadpQTLs6$pval,3),
                  beta_eQTL = round(leadpQTLs6$beta_GTEx,3),
                  pval_eQTL = signif(leadpQTLs6$pval_GTEx,3),
                  tissue_eQTL = leadpQTLs6$tissue_GTEx)
tab2

matched = match(leadpQTLs7$gene, z4$gene)
z4 = z4[matched,]
table(leadpQTLs7$gene == z4$gene)

tab2_1 = data.table(gene = leadpQTLs7$gene,
                    ratio = paste0(z4$N_diff,"/",z4$N_diff + z4$N_same),
                    SNP = leadpQTLs7$rs_id,
                    ea = leadpQTLs7$effect_allele,
                    eaf = round(leadpQTLs7$eaf,3),
                    beta_pQTL = round(leadpQTLs7$beta,3),
                    pval_pQTL = signif(leadpQTLs7$pval,3),
                    beta_eQTL = round(leadpQTLs7$beta_GTEx,3),
                    pval_eQTL = signif(leadpQTLs7$pval_GTEx,3),
                    tissue_eQTL = leadpQTLs7$tissue_GTEx)
tab2_1

tab2 = rbind(tab2,tab2_1)
setorder(tab2,"pval_pQTL")
tab2

matched = match(tab2$tissue_eQTL,TissuesAbr$Original_Name)
table(is.na(matched))
table(tab2$tissue_eQTL == TissuesAbr[matched,Original_Name])
tab2[,tissue_eQTL:= TissuesAbr[matched,Abbreviation]]

matched = match(tab2$gene,leadpQTLs8$gene)
table(is.na(matched))
table(tab2$gene == leadpQTLs8[matched,gene])
tab2[,LIFE:= !leadpQTLs8[matched,same_direction]]

tab2
setorder(tab2,"pval_pQTL",-"tissue_eQTL")

#' ## Table 3: negative association in MetaXcan ###
#' 
#' Required columns:
#' 
#' * Gene (diff/tot)	
#' * beta & p-value of association in whole blood
#' * beta & p-value of association in best tissue
#' * replicated in LIFE? --> TRUE == correlation of LIFE GE data with LIFE Olink data is also negative
#' 
MetaXcan6 = fread("../results/07_MetaXcan_LIFE.txt",dec = ",",sep="\t")
myGenes = y4$gene

MetaXcan3 = copy(MetaXcan2)
MetaXcan3 = MetaXcan3[setting == "combined",]
MetaXcan3 = MetaXcan3[is.element(gene,myGenes),]
MetaXcan3 = MetaXcan3[hierFDR==T,]

MetaXcan4<-MetaXcan3[effect_size<0 & tissue != "Whole_Blood",.SD[pval==min(pval)],by=.(gene)]
MetaXcan5<-MetaXcan3[tissue == "Whole_Blood",]
MetaXcan7<-MetaXcan6[tissue == "LIFE_Whole.Blood" & is.element(gene,myGenes) & pval_combined<0.05,]

table(MetaXcan4$gene == y4$gene)
matched = match(MetaXcan4$gene, y4$gene)
y4 = y4[matched,]
table(MetaXcan4$gene == y4$gene)

tab3 = data.table(gene = MetaXcan4$gene,
                  ratio = paste0(y4$N_diff,"/",y4$N_diff + y4$N_same),
                  effectSize = round(MetaXcan4$effect_size,3),
                  pval = signif(MetaXcan4$pval,3),
                  n_SNPs = MetaXcan4$n_snps_used,
                  tissue = MetaXcan4$tissue)
tab3

matched = match(MetaXcan5$gene, y4$gene)
y4 = y4[matched,]
table(MetaXcan5$gene_name == y4$gene)

tab3_1 = data.table(gene = MetaXcan5$gene,
                    ratio = paste0(y4$N_diff,"/",y4$N_diff + y4$N_same),
                    effectSize = round(MetaXcan5$effect_size,3),
                    pval = signif(MetaXcan5$pval,3),
                    n_SNPs = MetaXcan5$n_snps_used,
                    tissue = MetaXcan5$tissue)
tab3_1

tab3 = rbind(tab3,tab3_1)
setorder(tab3,"pval")
tab3

matched = match(tab3$tissue,TissuesAbr$Original_Name)
table(is.na(matched))
table(tab3$tissue == TissuesAbr[matched,Original_Name])
tab3[,tissue:= TissuesAbr[matched,Abbreviation]]

matched = match(tab3$gene,MetaXcan7$gene)
table(is.na(matched))
table(tab3$gene == MetaXcan7[matched,gene])
tab3[,LIFE:= MetaXcan7[matched,effect_size_combined<0]]

setorder(tab3,"pval",-"tissue")
tab3

#' ## Table 4: negative causality in MendelRando ###
#' 
#' Required columns:
#' 
#' * Gene (diff/tot)
#' * beta & p-value of causality and association in whole blood
#' * beta & p-value of causality and association in best tissue
#' 
myGenes = x4$gene

myTab_MR2 = copy(myTab_MR)
myTab_MR2 = myTab_MR2[is.element(gene,myGenes),]
myTab_MR2 = myTab_MR2[hierFDR==T,]

myTab_MR3<-myTab_MR2[beta_MR<0 & tissue != "Whole_Blood",.SD[pval_MR==min(pval_MR)],by=.(gene)]
myTab_MR4<-myTab_MR2[tissue == "Whole_Blood",]

table(myTab_MR3$gene == x4$gene)
matched = match(myTab_MR3$gene, x4$gene)
x5 = copy(x4)
x4 = x4[matched,]
table(myTab_MR3$gene == x4$gene)

tab4 = data.table(gene = myTab_MR3$gene,
                  ratio = paste0(x4$N_diff,"/",x4$N_diff + x4$N_same),
                  beta_causal = round(myTab_MR3$beta_MR,3),
                  pval_causal = signif(myTab_MR3$pval_MR,3),
                  tissue = myTab_MR3$tissue)
tab4

matched = match(myTab_MR4$gene, x5$gene)
x5 = x5[matched,]
table(myTab_MR4$gene == x5$gene)

tab4_1 = data.table(gene = myTab_MR4$gene,
                    ratio = paste0(x5$N_diff,"/",x5$N_diff + x5$N_same),
                    beta_causal = round(myTab_MR4$beta_MR,3),
                    pval_causal = signif(myTab_MR4$pval_MR,3),
                    tissue = myTab_MR4$tissue)
tab4_1

tab4 = rbind(tab4,tab4_1)
setorder(tab4,"pval_causal")
tab4

dumID = paste(tab4$gene,tab4$tissue,sep="_")

matched = match(tab4$tissue,TissuesAbr$Original_Name)
table(is.na(matched))
table(tab4$tissue == TissuesAbr[matched,Original_Name])
tab4[,tissue:= TissuesAbr[matched,Abbreviation]]

matched = match(dumID,paste(MetaXcan3$gene,MetaXcan3$tissue,sep="_"))
table(is.na(matched))
table(tab4$gene == MetaXcan3[matched,gene])
tab4[,beta_assoc:= round(MetaXcan3[matched,effect_size],3)]
tab4[,pval_assoc:= signif(MetaXcan3[matched,pval],3)]

tab4

#' ## Save tables ###
save(tab1,tab2,tab3,tab4,file="../results/12_MainTables.RData")
tosave4 = data.table(data = c("tab1","tab2", "tab3", "tab4"), 
                     SheetNames = c("Table1","Table2", "Table3", "Table4"))
excel_fn = "../results/12_MainTables.xlsx"
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=F, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in hours): " ,round(difftime(Sys.time(), time0, tz,units = "hours"),2))

