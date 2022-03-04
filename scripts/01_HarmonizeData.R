#' ---
#' title: "Harmonize Olink data to hg38"
#' subtitle: "LWAS Olink"
#' author: "Janne Pott"
#' date: "November 2021"
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
#' The goal of this script is to uplift the genetic association data from hg19 to hg38. 
#' 
#' To do this, I use the harmonization function from the github repository [summary-gwas-imputation](https://github.com/hakyimlab/summary-gwas-imputation). 
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../mySourceFile.R")
.libPaths()
setwd(projectpath_main)

#' # Prep ####
#' ***
myTab<-data.table(read_excel("../temp/OlinkGenes.xlsx"))
head(myTab)

todofiles = data.table(setting = c("combined","females","males"),
                       fn = paste0("../exampleData/",
                                   c("pQTLs_combined_hg19.txt.gz",
                                     "pQTLs_females_hg19.txt.gz",
                                     "pQTLs_males_hg19.txt.gz")))
todofiles[,outfn := gsub("hg19","hg38",fn)]
todofiles[,outfn := gsub(".gz","",outfn)]
head(todofiles)

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

#' # Python calls to harmonize ####
#' ***
myPhenos = unique(myTab$Parameter)

registerDoMC(cores=10)

dumTab = foreach(i = 1:dim(todofiles)[1])%do%{
  # i=1
  time1<-Sys.time()
  myRow = todofiles[i,]
  myData = fread(myRow[,fn])
  
  dumTab2 = foreach(j = 1:length(myPhenos))%dopar%{
    # dumTab2 = foreach(j = 1:5)%dopar%{
    # j=1
    myPheno = myPhenos[j]
    myData2 = copy(myData)
    myData2 = myData2[pheno == myPheno]
    
    temp_fn = paste0("../temp/harmonize/gwas_hg19_",myPheno,"_",myRow[,setting],".txt")
    fwrite(myData2, file = temp_fn,quote = F,sep="\t")
    temp_fn_out = gsub("hg19","hg38",temp_fn)
    
    mycall1 = paste0(path_metaxcan_gwastools,"/gwas_parsing.py -gwas_file ",temp_fn," ",
                     "-liftover ",path_metaxcan_data,"/liftover/hg19ToHg38.over.chain.gz ",
                     "-snp_reference_metadata ",path_metaxcan_data,
                     "/reference_panel_1000G/variant_metadata.txt.gz METADATA ",
                     "-output_column_map rs_id variant_id ",
                     "-output_column_map other_allele non_effect_allele ",
                     "-output_column_map effect_allele effect_allele ",
                     "-output_column_map beta effect_size ",
                     "-output_column_map se standard_error ",
                     "-output_column_map pval pvalue ",
                     "-output_column_map chr chromosome ",
                     "--chromosome_format ",
                     "-output_column_map pos_hg19 position ",
                     "-output_column_map eaf frequency ",
                     "-output_column_map n_samples sample_size ",
                     "-output_column_map info info ",
                     "-output_column_map pheno pheno ",
                     "-output_column_map gene gene ",
                     "-output_column_map cyto cyto ",
                     "-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size info pheno cyto gene ",
                     "-output ",temp_fn_out)
    
    print(mycall1)
    system(mycall1)
    
    myData3 = fread(temp_fn_out)
    
    # Filter for gene region
    myRow2 = copy(myTab)
    myRow2 = myRow2[Parameter == myPheno,]
    myData3 = myData3[position <= myRow2$stop_hg38 +500000]
    myData3 = myData3[position >= myRow2$start_hg38 -500000]
    
    # Add FDR 
    myData3[,FDR:=p.adjust(pvalue,method="fdr")]
    myData3[,maf:=frequency]
    myData3[frequency>0.5,maf:=1-frequency]
    myData3[,setting :=myRow$setting]
    dim(myData3)
    
    # Check Columns
    myNames1<-c("variant_id","panel_variant_id","chromosome","position","effect_allele",
                "non_effect_allele","frequency","maf","sample_size","info",
                "effect_size","standard_error","pvalue","FDR","zscore",
                "cyto","gene","pheno","setting")
    setcolorder(myData3,myNames1)
    myNames2<-c("variant_id_hg19","variant_id_hg38","chr","pos_hg38","effect_allele",
                "other_allele","eaf","maf","n_samples","info",
                "beta","se","pval","FDR","zscore",
                "cyto","gene","pheno","setting")
    names(myData3)<-myNames2
    
    # Return
    myData3
  }
  
  dumTab2 = rbindlist(dumTab2)
  
  myData[,dumID := paste(rs_id,pheno,sep="_")]
  dumTab2[,dumID := paste(variant_id_hg19,pheno,sep="_")]
  
  message("Checking how many SNPs are still there after harmonization")
  print(table(is.element(myData$dumID,dumTab2$dumID)))
  
  myData[,dumID := NULL]
  dumTab2[,dumID := NULL]
  
  fwrite(dumTab2, file = myRow[,outfn],quote = F,sep="\t")
  R.utils::gzip(paste0(myRow[,outfn]), destname = paste0(myRow[,outfn],".gz"))
  
  dif = round(difftime(Sys.time(), time1, tz,units = "mins"),2)
  message("Finished working on setting ",myRow[,setting]," (",dif," minutes)")
  
  myRow[,time := dif]
  myRow
}

dumTab = rbindlist(dumTab)
dumTab

#' # SessionInfo
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
