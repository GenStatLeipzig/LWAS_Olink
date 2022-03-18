#' ---
#' title: "MetaXcan"
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
#' I want to test the association between genetically regulated gene expression and its protein levels in all settings for all proteins with significant associations. 
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
load("../temp/06_OlinkGenes.RData")
myTab2 = copy(myTab)
myTab2 = myTab2[GWAS_sig==T | GWAS_sig_female==T | GWAS_sig_male==T,]

GECorrData_LIFE = fread("../exampleData/GECorrData_LIFE.txt")

todofile = data.table(gene = rep(myTab2$Gene,3),
                      protein = rep(myTab2$Parameter,3),
                      setting = c(rep("combined",66),rep("females",66),rep("males",66)),
                      GWAS_sig = c(myTab2$GWAS_sig,myTab2$GWAS_sig_female,myTab2$GWAS_sig_male))
todofile[,dumID:=paste(protein,setting,sep="_")]

myHarmonizedObjects = dir(path = "../temp/harmonize//",pattern = "hg38_*")
myHarmonizedObjects_pheno = myHarmonizedObjects
myHarmonizedObjects_pheno = gsub("gwas_hg38_","",myHarmonizedObjects_pheno)
myHarmonizedObjects_pheno = gsub(".txt","",myHarmonizedObjects_pheno)
table(is.element(myHarmonizedObjects_pheno,todofile$dumID))
matched = match(todofile$dumID,myHarmonizedObjects_pheno)
stopifnot(todofile$dumID == myHarmonizedObjects_pheno[matched])

todofile[,harmonizedfn := myHarmonizedObjects[matched]]
todofile[,harmonizedpath := paste0("../temp/harmonize/",harmonizedfn)]

#' # Check Python environment ####
#' ***
use_python(path_phyton)
x<-conda_list(conda=path_conda)
stopifnot(is.element("metaxcan",x$name))
y<-x[x$name=="metaxcan",]
y2<-gsub("/bin/python","",y$python)
use_condaenv(condaenv = y2,conda = path_conda,required = TRUE)
py_config()
message("using conda enviroment ",y$name," at ",y2)

#' # Loop ####
#' ***
path_tissuewise = "../temp/tissuewise/"
myeQTLs<-dir(path = paste0(path_metaxcan_data, "/models/eqtl/mashr/"),pattern = "txt.gz")
myeQTLs<-gsub("mashr_","",myeQTLs)
myeQTLs<-gsub(".txt.gz","",myeQTLs)

registerDoMC(cores=10)

dumTab = foreach(i = 1:dim(todofile)[1])%do%{
  #i=1
  myRow = todofile[i,]
  message("\nWorking on protein ",myRow$protein," coded by gene ",myRow$gene, " in setting ",myRow$setting, " (",i," of ",dim(todofile)[1],") ...")
  
  time2<-Sys.time()
  
  dumTab2 = foreach(j = 1:length(myeQTLs))%dopar%{
    #j=16
    myTissue = myeQTLs[j]
    message("     ... working on tissue ",myTissue,", ",j, " of ",length(myeQTLs))
    fn = paste(myRow$gene,myRow$setting,sep="_")
    
    mycall = paste0(path_metaxcan_tools,"/SPrediXcan.py --gwas_file ",myRow$harmonizedpath," ",
                    "--snp_column panel_variant_id ",
                    "--effect_allele_column effect_allele ",
                    "--non_effect_allele_column non_effect_allele ",
                    "--beta_column effect_size ",
                    "--se_column standard_error ",
                    "--model_db_path ",path_metaxcan_data,
                    "/models/eqtl/mashr/mashr_",myTissue,".db ",
                    "--covariance ",path_metaxcan_data,
                    "/models/eqtl/mashr/mashr_",myTissue,".txt.gz ",
                    "--keep_non_rsid ",
                    "--additional_output ",
                    "--model_db_snp_key varID ",
                    "--throw ",
                    "--output_file ",path_tissuewise,fn,"__",myTissue,".csv")
    
    print(mycall)
    #system(mycall)
    
    check = file.exists(paste0(path_tissuewise,fn,"__",myTissue,".csv"))
    if(check == T){
      MetaXcan = read.csv(paste0(path_tissuewise,fn,"__",myTissue,".csv"))
      setDT(MetaXcan)
    }else{
      MetaXcan = data.table(gene_name = myRow$gene,
                            comment = "no SNPs available in this tissue")
    }
    
    
    MetaXcan[,setting:=myRow$setting]
    MetaXcan = MetaXcan[gene_name == myRow$gene,]
    MetaXcan[,tissue:=myTissue]
    MetaXcan
    }
  dumTab2 = rbindlist(dumTab2,fill=T)
  
  # Step 3: return run time
  message("     TOTAL TIME of gene ",myRow$gene, " (in mins): " ,round(difftime(Sys.time(), time2, tz,units = "mins"),2))
  
  dumTab2
}

MetaXcan<-rbindlist(dumTab,fill=T)
MetaXcan = MetaXcan[is.na(comment),]
MetaXcan[,comment := NULL]
MetaXcan

setnames(MetaXcan,"gene","ENSG")
setnames(MetaXcan,"gene_name","gene")
setnames(MetaXcan,"pvalue","pval")

MetaXcan_a<-MetaXcan[setting=="combined",]
MetaXcan_m<-MetaXcan[setting=="males",]
MetaXcan_f<-MetaXcan[setting=="females",]

#' # FDR pro Setting ####
#' ***
#' Similar to script 02
HierFDR_Holger = function(data){
  # data = pQTLs_a
  
  data2 = copy(data)
  
  myFDR<- addHierarchFDR(pvalues = data2[,pval], categs = data2[,gene],quiet = T)
  data2[,hierarch_fdr5proz:=myFDR$hierarch_fdr5proz]
  
  x1 = myFDR[,.(min_level1 = min(fdr_level1)), by = list(category,fdr_level2)]
  x1[,hierarch_fdr5proz:=fdr_level2<0.05]
  x0 = myFDR[hierarch_fdr5proz==T,.N,.(category,fdr_level2)]
  matched<-match(x1$category,x0$category)
  x1[,N:=x0$N[matched]]
  x1[is.na(N),N:=0]
  
  return(x1)
} 

myFDR_a<- addHierarchFDR(pvalues = MetaXcan_a[,pval], categs = MetaXcan_a[,gene],quiet = F)
myFDR_m<- addHierarchFDR(pvalues = MetaXcan_m[,pval], categs = MetaXcan_m[,gene],quiet = F)
myFDR_f<- addHierarchFDR(pvalues = MetaXcan_f[,pval], categs = MetaXcan_f[,gene],quiet = F)

MetaXcan_a[,hierFDR:=myFDR_a$hierarch_fdr5proz]
MetaXcan_m[,hierFDR:=myFDR_m$hierarch_fdr5proz]
MetaXcan_f[,hierFDR:=myFDR_f$hierarch_fdr5proz]

x1 = HierFDR_Holger(data = MetaXcan_a)
x2 = HierFDR_Holger(data = MetaXcan_m)
x3 = HierFDR_Holger(data = MetaXcan_f)

dim(x1[hierarch_fdr5proz==T])
dim(x2[hierarch_fdr5proz==T])
dim(x3[hierarch_fdr5proz==T])

#' ## Summary #### 
#' ***
#' 
#' There are 58 significant gGE associations in the combined setting, 50 in the male-stratified setting and 48 in the female-stratifed setting. 
#' 
#' In x1, x2, and x3: 
#' 
#' * min_level1 is the Simes p-value. 
#' * fdr_level2 is the FDR corrected Simes p-value -> pheno counted as significant if fdr_level2 < 0.05 x (k/n)
#' * hierarch_fdr5proz is a T/F flag if pheno is significanty associated with at least one SNP
#' * N is the number of associated SNPs
#'  
myFDR_short = cbind(x1,x2[,-1],x3[,-1])
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

#' * 8 genes are not associated in any setting
#' * 42 genes are associated in all settings  
#' * 6 genes are associated in the combined and females only settings (--> FABP4, PLAUR, SCGB3A2, TNFRSF11B, TNFRSF14, TNFSF13B)
#' * 8 gene are associated in the combined and males only settings (--> ALCAM, CNTN1, EGFR, FAS, MMP9, MPO, SELP, SPP1)
#' * 2 gene are only associated in the combined setting (--> RETN, TFF3)

myFDR2[Sig_f==F & Sig_a==F & Sig_m==F,] 
myFDR2[Sig_f==T & Sig_a==T & Sig_m==F,] 
myFDR2[Sig_f==F & Sig_a==T & Sig_m==T,] 
myFDR2[Sig_f==F & Sig_a==T & Sig_m==F,] 

MetaXcan2<-rbind(MetaXcan_a,MetaXcan_m,MetaXcan_f)
MetaXcan2[,table(setting, hierFDR)]

#' change from long format to wide format
names(MetaXcan2)
MetaXcan2[,dumID := paste(gene,tissue,sep="::")]
MetaXcan3<-dcast(MetaXcan2, formula = dumID ~ setting,
                 value.var = c("zscore","effect_size","pval","hierFDR","n_snps_used"),
                 sep = "_")
head(MetaXcan3)
MetaXcan3[,gene := gsub("::.*","",dumID)]
MetaXcan3[,tissue := gsub(".*::","",dumID)]
MetaXcan3[,table(hierFDR_males,hierFDR_females,hierFDR_combined)]

#' ## Plot #### 
#' ***

myPlotMetaXcan<-ggplot(MetaXcan3[hierFDR_females==T | hierFDR_males==T,], aes(x=zscore_males, y=zscore_females,col=gene)) +
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1.25)+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_point(size=2.5)+ 
  theme_bw(base_size = 10)+
  ggtitle(label = "Scatter plot of MetaXcan associations (hier. FDR ==T for at least one strata)")+
  theme(plot.title = element_text(hjust = 0.5))+ 
  xlab("male")+
  ylab("female")+
  guides(size="none",fill="none",col="none")
myPlotMetaXcan

#' ## Merge LIFE data to MetaXcan ####
head(MetaXcan3)

MetaXcan4<-rbind(MetaXcan3,GECorrData_LIFE,fill=T)
MetaXcan4

#' # Save ####
#' ***
dummy1<-unique(MetaXcan4[hierFDR_combined==T,gene])
dummy4<-MetaXcan4[grepl("LIFE",tissue) & (pval_combined<=0.05 | pval_males<=0.05 | pval_females<=0.05),gene]
dummy5 = MetaXcan2[,.SD[pval==min(pval)],by=.(gene)]
bestTissue = dummy5[,paste(tissue,collapse = "|"),by=gene]
bestPvalue = dummy5[!duplicated(gene),pval,by=gene]
bestSetting = dummy5[!duplicated(gene),setting,by=gene]

myTab[GWAS_sig==F,MetaXcan:=NA]
myTab[GWAS_sig==T,MetaXcan:=F]
myTab[Gene %in% dummy1,MetaXcan:=T]
table(myTab$MetaXcan)

myTab[GE_ok=="no" | GWAS_sig==F,CorrLIFE:=NA]
myTab[GE_ok=="yes" & GWAS_sig==T,CorrLIFE:=F]
myTab[Gene %in% dummy4,CorrLIFE:=T]
table(myTab$CorrLIFE)

myTab[,table(CorrLIFE,MetaXcan)]
myTab[CorrLIFE==T & MetaXcan==F,]

matched = match(myTab$Gene,bestTissue$gene)
myTab[,MultiXcan_bestTissue:=bestTissue[matched,V1]]
matched = match(myTab$Gene,bestPvalue$gene)
myTab[,MultiXcan_bestPvalue:=bestPvalue[matched,pval]]
matched = match(myTab$Gene,bestSetting$gene)
myTab[,MultiXcan_bestSetting:=bestSetting[matched,setting]]
head(myTab)

write.table(MetaXcan2,file = "../results/07_MetaXcan_long.txt",
            col.names = T,row.names = F,quote = F,dec = ",",sep="\t")
write.table(MetaXcan4,file = "../results/07_MetaXcan_LIFE.txt",
            col.names = T,row.names = F,quote = F,dec = ",",sep="\t")
save(MetaXcan2,file="../results/07_MetaXcan_long.RData")
save(myTab,file="../temp/07_OlinkGenes.RData")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in hours): " ,round(difftime(Sys.time(), time0, tz,units = "hours"),2))
