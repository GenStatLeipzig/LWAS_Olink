#' ---
#' title: "Hierarchical FDR"
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
#' I want to use hierarchical FDR for my Olink results. For this, I use the *addHierarchFDR()* function from the R-package *toolboxH* from Holger Kirsten (co-author of the publication; member of group Genstat at IMISE). You can download it from github via *devtools::install_github('holgerman/toolboxH')*
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../mySourceFile.R")
.libPaths()
setwd(projectpath_main)

require(toolboxH)

#' # Prep ####
#' ***
myTab<-fread("../temp/OlinkGenes.txt")
head(myTab)

#' # Functions for this script ####
#' ***

HierFDR_Holger = function(data){
  # data = pQTLs_a
  
  data2 = copy(data)
  
  myFDR<- addHierarchFDR(pvalues = data2[,pval], categs = data2[,gene],quiet = F)
  data2[,hierarch_fdr5proz:=myFDR$hierarch_fdr5proz]
  
  x1 = myFDR[,.(min_level1 = min(fdr_level1)), by = list(category,fdr_level2)]
  x1[,hierarch_fdr5proz:=fdr_level2<0.05]
  x0 = myFDR[hierarch_fdr5proz==T,.N,.(category,fdr_level2)]
  matched<-match(x1$category,x0$category)
  x1[,N:=x0$N[matched]]
  x1[is.na(N),N:=0]
  
  return(x1)
} 

HierFDR_TreeQTL = function(data,myFDR_short,setting,fn_in,fn_out){
  # data = pQTLs_a
  # myFDR_short = myFDR_short
  # setting = "combined"
  # fn_in = "../results/01_pQTLs_all_cis.txt"
  # fn_out = "../results/01_pQTLs_all_cis_pAssoc_TreeQTL.txt"
  
  data2 = copy(data)
  fdr2 = copy(myFDR_short)
  
  myNames5 = c("FDRl2","SimesP","Sig","NSig")
  if(setting == "combined"){
    myNames5 = paste0(myNames5,"_a")
  }else if(setting == "males"){
    myNames5 = paste0(myNames5,"_m")
  }else if(setting == "females"){
    myNames5 = paste0(myNames5,"_f")
  }
  myNames6 = c("category",myNames5)
  colsOut<-setdiff(colnames(fdr2),myNames6)
  fdr2[,get("colsOut"):=NULL]
  setcolorder(fdr2,myNames6)
  myNames7<-c("category","FDRl2","SimesP","Sig","NSig")
  names(fdr2) = myNames7
  fdr2
  
  # Prep data
  myTest_snps<-data2[,.N,variant_id_hg19]
  names(myTest_snps)<-c("family","n_tests")
  myTest_gene<-data2[,.N,gene]
  names(myTest_gene)<-c("family","n_tests")
  
  #' Get associated genes
  data_gene_BH<- get_eGenes(n_tests_per_gene = myTest_gene, 
                            m_eqtl_out = fn_in,
                            method="BH",
                            level1 = 0.05, 
                            level2 = 0.05, 
                            slice_size = 1e+05,
                            silent = FALSE, 
                            gene_pvals = NA)
  
  data_gene_BH2<-get_eAssociations(eDiscoveries = data_gene_BH, 
                                   n_tests = myTest_gene, 
                                   m_eqtl_out = fn_in, 
                                   out_file = fn_out,
                                   by_snp = FALSE,
                                   slice_size = 1e+05,
                                   silent = FALSE)
  
  head(data_gene_BH2)
  setDT(data_gene_BH)
  head(data_gene_BH)
  
  #' Check 1: are the same genes noted as significant? --> Yes :)
  stopifnot(is.element(data_gene_BH$family,fdr2[Sig==T,category]))
  stopifnot(data_gene_BH$family==fdr2[Sig==T,category])
  message("Check 1 - passed; the same genes are significant")
  
  #' Check 2: are the same number of SNPs significant? --> Yes :)
  stopifnot(data_gene_BH$n_sel==fdr2[Sig==T,NSig])
  message("Check 2 - passed; the same number of SNPs are significant")
  
  #' Check 3: are the Simes p-Values the same? --> Yes :)
  stopifnot(signif(data_gene_BH$fam_p,6)==signif(fdr2[Sig==T,SimesP],6))
  message("Check 3 - passed; Simes p-Values per gene are the same")
  
  #' Check 4: are the SNPs associated?
  data_gene_BH2[,dumID:=paste(gene,SNP,sep="::")]
  data2[,dumID:=paste(gene,variant_id_hg38,sep="::")]
  
  stopifnot(is.element(data2[hierarch_fdr5proz==T,dumID],data_gene_BH2[,dumID]))
  stopifnot(is.element(data_gene_BH2[,dumID],data2[hierarch_fdr5proz==T,dumID]))
  message("Check 4 - passed; the same SNP-gene pairs are associated")
  
  matched<-match(data2[,dumID],data_gene_BH2[,dumID])
  data_gene_BH2<-data_gene_BH2[matched,]
  #stopifnot(data_gene_BH2$dumID == data2$dumID)
  data2[,BBFDR:=data_gene_BH2$BBFDR]
  data2[,dumID:=NULL]
  
  return(data2)
  
}

#' # Load Olink data hg38 ####
#' ***
pQTLs_a = fread("../exampleData/pQTLs_combined_hg38.txt.gz")
pQTLs_f = fread("../exampleData/pQTLs_females_hg38.txt.gz")
pQTLs_m = fread("../exampleData/pQTLs_males_hg38.txt.gz")

pQTLs_a[,length(unique(gene))]
pQTLs_f[,length(unique(gene))]
pQTLs_m[,length(unique(gene))]

#' # Hierarchical FDR ####
#' ***
#' 
#' * Level 1: SNPs per gene
#' * Level 2: lead SNPs per gene 
#' * get k as number of genes associated with at least one SNP with $\alpha = 0.05$
#' * define adjusted $\alpha_2$ as $\alpha_2 = \alpha \cdot k/n$ 
#' * use $\alpha_2$ on Level 1 FDRs 
#' 
myFDR_a<- addHierarchFDR(pvalues = pQTLs_a[,pval], categs = pQTLs_a[,gene],quiet = F)
myFDR_m<- addHierarchFDR(pvalues = pQTLs_m[,pval], categs = pQTLs_m[,gene],quiet = F)
myFDR_f<- addHierarchFDR(pvalues = pQTLs_f[,pval], categs = pQTLs_f[,gene],quiet = F)

pQTLs_a[,hierarch_fdr5proz:=myFDR_a$hierarch_fdr5proz]
pQTLs_m[,hierarch_fdr5proz:=myFDR_m$hierarch_fdr5proz]
pQTLs_f[,hierarch_fdr5proz:=myFDR_f$hierarch_fdr5proz]

x1 = HierFDR_Holger(data = pQTLs_a)
x2 = HierFDR_Holger(data = pQTLs_m)
x3 = HierFDR_Holger(data = pQTLs_f)

dim(x1[hierarch_fdr5proz==T])
dim(x2[hierarch_fdr5proz==T])
dim(x3[hierarch_fdr5proz==T])

#' There are 64 significantly associated genes in the combined setting, 54 in males only, and 48 in females only.
#' 
#' Hence, I can define the BBFDR as FDR * (92/k) in the different settings 
#' 
pQTLs_a[,BBFDR:= FDR * 92/64]
pQTLs_f[,BBFDR:= FDR * 92/48]
pQTLs_m[,BBFDR:= FDR * 92/54]

#' # Summary ####
#' ***
#' In x1, x2, and x3: 
#' 
#' * min_level1 is the Simes p-value. 
#' * fdr_level2 is the FDR corrected Simes p-value -> pheno counted as significant if fdr_level2 < 0.05 x (k/n)
#' * hierarch_fdr5proz is a T/F flag if pheno is significanlty associated with at least one SNP
#' * N is the number of associated SNPs
#'  
myFDR_short = cbind(x1,x2[,-1],x3[,-1])
names(myFDR_short)[2:13]<-c("FDRl2_a","SimesP_a","Sig_a","NSig_a",
                            "FDRl2_m","SimesP_m","Sig_m","NSig_m",
                            "FDRl2_f","SimesP_f","Sig_f","NSig_f")
head(myFDR_short)
dim(myFDR_short)

table(myFDR_short$Sig_f,myFDR_short$Sig_m,myFDR_short$Sig_a)

dummy<-myTab[,c(1,2,5)]  
stopifnot(dummy$Gene == myFDR_short$category)
myFDR2<-cbind(dummy,myFDR_short[,-1])

myTab[,GWAS_sig:=myFDR_short$Sig_a]
myTab[,SimesP:=myFDR_short$SimesP_a]
myTab[,nSNPs:=myFDR_short$NSig_a]

myTab[,GWAS_sig_male:=myFDR_short$Sig_m]
myTab[,SimesP_male:=myFDR_short$SimesP_m]
myTab[,nSNPs_male:=myFDR_short$NSig_m]

myTab[,GWAS_sig_female:=myFDR_short$Sig_f]
myTab[,SimesP_female:=myFDR_short$SimesP_f]
myTab[,nSNPs_female:=myFDR_short$NSig_f]
head(myTab)

#' # Save ####
#' ***
save(myTab,file="../temp/01_OlinkGenes.RData")
write.table(myFDR_short,file="../results/01_pQTLs_sigGenes.txt",
            col.names = T,row.names = F, quote = F, sep="\t",dec = ",")

save(pQTLs_a,file="../results/01_pQTLs_all_cis.RData")
save(pQTLs_m,file="../results/01_pQTLs_males_cis.RData")
save(pQTLs_f,file="../results/01_pQTLs_females_cis.RData")

#' # SessionInfo
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
