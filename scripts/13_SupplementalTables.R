#' ---
#' title: "Supplemental Tables"
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
#' I want to create the supplemental tables in one script instead of doing it by hand. 
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../mySourceFile.R")
.libPaths()
setwd(projectpath_main)

#' # Get Content (tab0) ####
#' ***
tab0 = data.table(Table = paste0("S",c(1:12)),
                  Title = c("Overview (gene based)	all biomarkers and their mean and sd + all relevant results (T/F)",
                            "Hier. FDR by gene	Simes p-value per gene and setting (all, male, female)",
                            "Summary Stats	GWAS Pipeline Output (hier. FDR, tagSNP==T)",
                            "Lead cis pQTLs (1)	One per biomarker, annotated with sex stratified results and PC adjusted",
                            "Lead cis pQTLs (2)	all tissues, all significant biomarker",
                            "Overview (tissue based)	short overview of all results across all analyzed tissues",
                            "Colocalization Results	all tissues, all significant biomarker",
                            "MetaXcan all results	all tissues, all significant biomarker",
                            "MetaXcan top results	best associated tissue, correlation results in LIFE",
                            "MR GE --> P	best eQTLsof GTEx p<5e-8",
                            "MR P --> CAD	best pQTLs p<5e-8",
                            "MR GE --> P --> CAD	for sig. P-->CAD pairs only"),
                  Source = c("../temp/09_OlinkGenes.RData",
                             "../results/02_pQTLs_sigGenes.txt",
                             "../exampleData/topliste_2022-02-02_Olink_cis_neu.txt",
                             "../results/03_leadSNPs_IA.RData",
                             "../results/04_c_eQTL_lookup_Summary_LD.RData",
                             "../results/11_allResults_tissueSpec.RData",
                             "../results/06_Coloc_Genewise.RData",
                             "../results/07_MetaXcan_long.txt",
                             "../results/07_MetaXcan_LIFE.txt",
                             "../results/08_MR_GE_PE_long.RData",
                             "../results/09_MR_PE_CAD_long.RData",
                             "../results/10_MR_Mediation.RData"))

tab0

#' # Get Overview (gene-based) (tab1) ####
#' ***
#' 
#' * Protein Information: Biomarker, Full name, Sample size, Overall, Female, Male, P-value
#' * Gene Information: Gene, Cytoband, ENSG, chr, bp start (hg38), bp stop (hg38)
#' * Genetic Association (SNP level): all, male, female, novel sex IA, sex IA type, diff dir, n tissues with diff dir
#' * Genetic Association (gene level): Coloc - shared, Coloc - indep, MetaXcan, Correlation
#' * Causality: MR - GE, MR - CAD
#' 

load(tab0$Source[1])
tab1 = copy(myTab)
names(tab1)

tab1 = tab1[,c(1,4,12:19, 2,5,3,6,9,10, 20,23,26,36,29,30,32,32, 37,38,49,51,52,50,54,55)]
names(tab1)[23] = "diffDir"
names(tab1)[24] = "n_diffDir"

tab1[n_diffDir>0,diffDir:=T]
tab1[n_diffDir==0 | is.na(n_diffDir),diffDir:=F]
tab1[,diffDir:=as.logical(diffDir)]
tab1[is.na(n_diffDir),diffDir:=NA]
head(tab1)

#' # Get Hier. FDR (tab2) ####
#' ***
#' Supplemental Table S2: Results of the hierarchical FDR over all 92 analyzed biomarkers in three setting (all, male-specific, female-specific). Simes P-value is the minimal FDR corrected p-value per gene (first level). FDR l2 is the second level FDR over all Simes p-values, and proteins are considered significantly associated if FDR l2 < 0.05 (hier. FDR == T). This resulted in k = 64, 54 and 48 proteins with significant associations in all, males, and females, respectively. We then used 0.05 * k /92 as significance level on the FDR corrected p-value per SNP of the first level, and report here how many SNPs are considered associated per gene. 														
#' 
#' Please note: in this example script, the pruning information (how many independent SNPs per protein) is not added, as that required individual level data (not possible to host on github).

tab2 = fread(tab0$Source[2],dec = ",")
head(tab2)
matched = match(tab2$category,myTab$Gene)
table(is.na(matched))
tab2[,protein:= myTab$Parameter[matched]]
tab2[,cytoband:= myTab$Cytoband[matched]]

myNames = names(tab2)
myNames = myNames[c(14,15,1,3,2,4,5,7,6,8,9,11,10,12,13)]
setcolorder(tab2,myNames)
setnames(tab2,"category","gene")
head(tab2)

#' # Get Summary Stats (tab3) ####
#' ***
#' "Supplemental Table S3: Genetic association results, after hierarchical FDR and priority pruning (n=776 unique SNPs, n=792 protein - gene pairs). All other summary statistics can be downloaded at the Leipzig Heath Atlas (Link). Due to overlapping cis-loci and priority pruning, not all lead pQTLs are listed (see Table S4 & S5 for all lead pQTLs per biomarker). As significance level alpha = 0.05 *64/92 was used on the FDR values. EA = effect allele; OA = other allele; EAF = effect allele frequency; info = imputation info score; #SNPs_tagged = number of significant SNPs that are in LD (r2>0.1) with the respective SNP; GWASCatalog = traits from the GWAS catalog including LD to catalog SNP (r2>0.3)"											
tab3 = fread(tab0$Source[3])
tab3 = tab3[setting == "combined",]
dumID_GWAS<-paste(tab3$provided_name,tab3$gene,sep="::")
filt<-duplicated(dumID_GWAS)
tab3<-tab3[!filt,]
dumID_GWAS<-dumID_GWAS[!filt]

names(tab3)[1:60]
myNames = c("markername","variant_id_hg38","gene","BBFDR","chr","effect_allele",
            "other_allele","topeaf", "topinfo","cyto","nearestgenes","CADD_scaled","corinfo_gwas2",
            "cisgene","transgene", "KEGG","reactome","DOSE", "GO")
myNames
colsOut<-setdiff(colnames(tab3),myNames)
tab3[,get("colsOut"):=NULL]
setcolorder(tab3,myNames)
myNames2<-c("markername_hg19","markername_hg38","gene","BBFDR","chr","ea","oa","eaf","info","cytoband","nearestgenes","CADD","GWASCatalog","eQTL_cis","eQTL_trans", "KEGG","reactome","DOSE", "GO")
names(tab3)<-myNames2

load("../results/02_pQTLs_all_cis.RData")
pQTLs = copy(pQTLs_a)
pQTLs[,dumID:=paste(variant_id_hg19,gene,sep="::")]
dumID_GWAS<-paste(tab3$markername_hg19,tab3$gene,sep="::")
matched = match(dumID_GWAS,pQTLs$dumID)
table(is.na(matched))
pQTLs = pQTLs[matched,]
tab3[,biomarker := pQTLs$pheno]
tab3[,beta := pQTLs$beta]
tab3[,se := pQTLs$se]
tab3[,pval := pQTLs$pval]
tab3[,logPval := -log10(pQTLs$pval)]
tab3[,FDR := pQTLs$FDR]
tab3[,ea2 := pQTLs$effect_allele]
tab3[,oa2 := pQTLs$other_allele]
tab3[,eaf2 := pQTLs$eaf]

myNames = names(tab3)[c(1,2,5,26:28,9:10,3,11,20:25,4,12:14,16:19)]
myNames
colsOut<-setdiff(colnames(tab3),myNames)
tab3[,get("colsOut"):=NULL]
setcolorder(tab3,myNames)
head(tab3)

#' # Get lead SNPs (1) (tab4) ####
#' ***
#'
#' * Locus information: Gene, Cyto, SNP ID hg38, SNP ID hg19, EA, info 
#' * PE ~ SNP in combined setting: EAF, beta, se, p-value, FDR_sig
#' * PE ~ SNP in males: EAF, beta, se, pvalue, FDR_sig
#' * PE ~ SNP in females: EAF, beta, se, pvalue, FDR_sig
#' * Sex-Interaction: diff, dff_se, diff_p, diff_p_FDR, type
#' * PC-Adjustment: beta, se, pvalue
#' 
#' Please note: in this example script, the PC adjusted statistics is not added, as it is not part of the example data
 
load(tab0$Source[4])
tab4 = copy(leadSNPs_IA)
names(tab4)
matched = match(tab4$gene, myTab$Gene)
table(tab4$gene == myTab$Gene[matched])
tab4[,protein := myTab$Parameter[matched]]

myNames = names(tab4)[c(36,3,2,4:6,7, 8:11,14, 16:19,22, 24:27,30, 31:35)]
myNames
colsOut<-setdiff(colnames(tab4),myNames)
tab4[,get("colsOut"):=NULL]
setcolorder(tab4,myNames)
head(tab4)

#' # Get lead SNPs (2) (tab5) ####
#' ***
#' 
#' * Locus information: protein, cytoband, gene, setting
#' * SNP information: SNP ID hg38, SNP ID hg19, EA, info
#' * Statistics protein ~ pQTL: eaf	maf beta	se	pval	FDR	
#' * Statistics GTEx eQTLs: beta	se	pval	n	tissue	direction	
#' * Statistics GE ~ best associated eQTL: SNP ID GE	beta	se	pval	n	LD_r2
#' 					

loaded = load(tab0$Source[5])
loaded
tab5 = copy(leadpQTLs4)
names(tab5)
myNames = names(tab5)[c(18,16,17,19, 2,1,5,10, 7:9,11:14, 25:30, 32,34:37,40)]
myNames
colsOut<-setdiff(colnames(tab5),myNames)
tab5[,get("colsOut"):=NULL]
setcolorder(tab5,myNames)
head(tab5)

#' # Get Overview (tissue-based) (tab6) ####
#' ***
#'
#' * Locus information: protein	cytoband	gene	GTEx tissue
#' * eQTL information: lead pQTL	same dir	lead eQTL	beta	pval	r2
#' * Coloc: nsnps	PP.H3	PP.H4	
#' * MetaXcan: nsnps	beta	p-value	hierFDR
#' * MR GE - P: beta	se	p-value	hierFDR	
#' 		

loaded = load(tab0$Source[6])
loaded
tab6 = copy(myData)
names(tab6)
head(tab6)
dumTab = copy(leadpQTLs4)
dumTab[, dumID2 :=paste(gene,tissue_GTEx,sep="::")]
dumTab = dumTab[setting == "combined"]
table(is.element(tab6$myID,dumTab$dumID2))
tab6[!is.element(myID,dumTab$dumID2)]
matched = match(tab6$myID,dumTab$dumID2)
tab6[,eQTLs_besteQTL_beta := dumTab[matched,beta_GTEx_best]]
tab6[,eQTLs_besteQTL_pval := dumTab[matched,pval_GTEx_best]]
myNames = names(tab6)[c(4,5,2,3, 6:10,23,24,11, 12:14, 15:18, 19:22)]
myNames
colsOut<-setdiff(colnames(tab6),myNames)
tab6[,get("colsOut"):=NULL]
setcolorder(tab6,myNames)
head(tab6)

#' # Get Coloc (gene wise) (tab7) ####
#' ***
#' protein	cytoband	gene	setting GTEx tissue	#SNPs	PP.H0.abf	PP.H1.abf	PP.H2.abf	PP.H3.abf	PP.H4.abf
loaded = load(tab0$Source[7])
loaded
tab7 = copy(Coloc_genewise)
names(tab7)
matched = match(tab7$gene, myTab$Gene)
table(tab7$gene == myTab$Gene[matched])
tab7[,protein := myTab$Parameter[matched]]
myNames = names(tab7)[c(13,8,9,12,11,1:6)]
myNames
colsOut<-setdiff(colnames(tab7),myNames)
tab7[,get("colsOut"):=NULL]
setcolorder(tab7,myNames)
head(tab7)

# tab7= tab7[setting =="combined"]
# tab7[,setting:=NULL]
# head(tab7)

#' # Get MetaXcan (tissue wise) (tab8) ####
#' ***
#' protein	cytoband	gene setting	GTEx tissue	Z score	effect size	p-value	hierarch_fdr5proz	var_g	pred_perf_r2	pred_perf_pval	n_snps_used	n_snps_in_cov	n_snps_in_model

tab8 = fread(tab0$Source[8],dec=",")
names(tab8)
matched = match(tab8$gene, myTab$Gene)
table(tab8$gene == myTab$Gene[matched])
tab8[,protein := myTab$Parameter[matched]]
tab8[,cytoband := myTab$Cytoband[matched]]
myNames = names(tab8)[c(19,20,2,15,16,3,4,5,17,6,10,13)]
myNames
colsOut<-setdiff(colnames(tab8),myNames)
tab8[,get("colsOut"):=NULL]
setcolorder(tab8,myNames)
head(tab8)

# tab8= tab8[setting =="combined"]
# tab8[,setting:=NULL]
# head(tab8)

#' # Get MetaXcan (best tissue) (tab9) ####
#' ***
#' protein	cytoband	gene	setting tissue	
#' best: effect size	p-value	hierarch_fdr5proz	
#' whole blood: effect size	p-value	hierarch_fdr5proz	
#' LIFE: part. cor	p-value	N

tab9 = fread(tab0$Source[9],dec=",")
names(tab9)
matched = match(tab9$gene, myTab$Gene)
table(tab9$gene == myTab$Gene[matched])
tab9[,protein := myTab$Parameter[matched]]
tab9[,cytoband := myTab$Cytoband[matched]]

tab9_a = tab9[,c(19,20,17,18,2,5,8,11,14)]
tab9_m = tab9[,c(19,20,17,18,4,7,10,13,16)]
tab9_f = tab9[,c(19,20,17,18,3,6,9,12,15)]
tab9_a[,setting := "combined"]
tab9_m[,setting := "males"]
tab9_f[,setting := "females"]
names(tab9_a) = c("protein","cytoband","gene","tissue","zscore","beta","pval","hierFDR","Nsnp","setting")
names(tab9_m) = c("protein","cytoband","gene","tissue","zscore","beta","pval","hierFDR","Nsnp","setting")
names(tab9_f) = c("protein","cytoband","gene","tissue","zscore","beta","pval","hierFDR","Nsnp","setting")
tab9_2 = rbind(tab9_a,tab9_m,tab9_f)
tab9_2[, dumID:=paste(gene,setting,sep=":")]

tab9_3<-tab9_2[tissue != "Whole_Blood",.SD[pval==min(pval)],by=.(dumID)]
tab9_4<-tab9_2[tissue == "Whole_Blood",]
tab9_5<-tab9_2[tissue == "LIFE_Whole.Blood",]

names(tab9_3)
myNames = names(tab9_3)[c(1:4,11,5,7:9)]
myNames
colsOut<-setdiff(colnames(tab9_3),myNames)
tab9_3[,get("colsOut"):=NULL]
setcolorder(tab9_3,myNames)

matched = match(tab9_3$dumID,tab9_4$dumID)
table(is.na(matched))
tab9_3[,beta_WB:= tab9_4[matched,beta]]
tab9_3[,pval_WB:= tab9_4[matched,pval]]
tab9_3[,hierFDR_WB:= tab9_4[matched,hierFDR]]

matched = match(tab9_3$dumID,tab9_5$dumID)
table(is.na(matched))
tab9_3[,partCor_LIFE:= tab9_5[matched,beta]]
tab9_3[,pval_LIFE:= tab9_5[matched,pval]]
tab9_3[,N_LIFE:= tab9_5[matched,Nsnp]]

head(tab9_3)
tab9 = copy(tab9_3)
tab9[,dumID:=NULL]
head(tab9)

# tab9= tab9[setting =="combined"]
# tab9[,setting:=NULL]
# head(tab9)

#' # Get MR (GE - PE) (tab10) ####
#' ***
#' protein	cytoband	gene	tissue	SNP	chr	pos	ea	info	maf	beta	se	p	eaf	beta	se	p	beta	se	p	hierFDR

loaded = load(tab0$Source[10])
loaded
tab10 = copy(myTab_MR)
names(tab10)
matched = match(tab10$gene, myTab$Gene)
table(tab10$gene == myTab$Gene[matched])
tab10[,protein := myTab$Parameter[matched]]
tab10[,cytoband := myTab$Cytoband[matched]]
myNames = names(tab10)[c(13,14,2,3,4,5:11)]
myNames
colsOut<-setdiff(colnames(tab10),myNames)
tab10[,get("colsOut"):=NULL]
setcolorder(tab10,myNames)
head(tab10)

# tab10= tab10[setting =="combined"]
# tab10[,setting:=NULL]
# head(tab10)

#' # Get MR (PE - CAD) (tab11) ####
#' ***
#' protein	cytoband	gene	tissue	SNP	chr	pos	ea	info	maf	beta	se	p	eaf	beta	se	p	beta	se	p	hierFDR
loaded = load(tab0$Source[11])
loaded
tab11 = copy(myTab_MR)
names(tab11)
matched = match(tab11$gene, myTab$Gene)
table(tab11$gene == myTab$Gene[matched])
tab11[,protein := myTab$Parameter[matched]]
tab11[,cytoband := myTab$Cytoband[matched]]
myNames = names(tab11)[c(11,12,2,3,4:10)]
myNames
colsOut<-setdiff(colnames(tab11),myNames)
tab11[,get("colsOut"):=NULL]
setcolorder(tab11,myNames)
head(tab11)

# tab11= tab11[setting =="combined"]
# tab11[,setting:=NULL]
# head(tab11)

#' # Get MR (GE - PE - CAD) (tab12) ####
#' ***
loaded = load(tab0$Source[12])
loaded
tab12 = copy(myTab_MR_mediation3)
names(tab12)
matched = match(tab12$gene, myTab$Gene)
table(tab12$gene == myTab$Gene[matched])
tab12[,protein := myTab$Parameter[matched]]
tab12[,cytoband := myTab$Cytoband[matched]]
myNames = names(tab12)[c(28,29,2:22,27,23:26)]
myNames
colsOut<-setdiff(colnames(tab12),myNames)
tab12[,get("colsOut"):=NULL]
setcolorder(tab12,myNames)
head(tab12)

# tab12= tab12[setting =="combined"]
# tab12[,setting:=NULL]
# head(tab12)

#' # Check colnames ####
#' ***
#' I want the first colnames to be the same 
#' protein - cytoband - gene - setting - tissue
names(tab1)
names(tab1)[1] = "Protein"
names(tab2)
names(tab3)
myNames = names(tab3)[c(11,8,9,1:7,12:17,10,18:24)]
setcolorder(tab3,myNames)
names(tab3)[1]= "protein"
names(tab3)[7]= "EA"
names(tab3)[8]= "OA"
names(tab3)[9]= "EAF"
names(tab3)
names(tab4)
names(tab4)[2]= "cytoband"
names(tab4) = gsub("_a","_combined",names(tab4))
names(tab4) = gsub("_f","_females",names(tab4))
names(tab4) = gsub("_m","_males",names(tab4))
names(tab4)
names(tab5)
names(tab5)[1:2]= c("protein","cytoband")
names(tab6)
names(tab6)[2]= "cytoband"
names(tab7)
names(tab7)[2]= "cytoband"
names(tab7)[5]= "tissue"
names(tab8)
names(tab9)
names(tab10)
names(tab11)
names(tab12)

#' # Save tables ###
#' ***
tosave4 = data.table(data = c("tab0","tab1","tab2", "tab3", "tab4","tab5","tab6",
                              "tab7","tab8","tab9","tab10","tab11","tab12"), 
                     SheetNames = c("Content","TableS1","TableS2", "TableS3", 
                                    "TableS4","TableS5",
                                    "TableS6","TableS7","TableS8","TableS9","TableS10",
                                    "TableS11","TableS12"))
excel_fn = "../results/13_SupplementalTables.xlsx"
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))
