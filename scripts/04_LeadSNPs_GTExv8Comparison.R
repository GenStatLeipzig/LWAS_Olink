#' ---
#' title: "Lead SNPs: GTEx v8 look-up"
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
#' I want first to define the lead SNP per gene and setting (see also script 03) and then to compare the effect estimates to those in GTEx v8 (same SNP across all tissues). 
#' 
#' In addition, I want to find the best eQTL per tissue (defined by lowest p-value) and test the LD between the best eQTL and the best pQTL. 
#' 
#' In the first step, I load and filter the GTEx v8 data by tissue. I also save the matched data set, as I need this one for the co-localization analyses later. 
#'  
#' I assume the GTEx data is named *GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_TISSUE.allpairs.txt.gz*. If the files are named differently, please change code lines x & y accordingly. 
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
load("../temp/03_OlinkGenes.RData")

load("../results/02_pQTLs_all_cis.RData")
load("../results/02_pQTLs_males_cis.RData")
load("../results/02_pQTLs_females_cis.RData")
pQTLs = rbind(pQTLs_a,pQTLs_m,pQTLs_f)

load("../results/03_leadSNPs_IA.RData")

#' # Get lead SNPs ####
#' ***
names(pQTLs)
pQTLs[,chrPosGeneSetting:=paste(chr,pos_hg38,gene,setting,sep=":")]
pQTLs[,chrPosGene:=paste(chr,pos_hg38,gene,sep=":")]
pQTLs[,chrPos:=paste(chr,pos_hg38,sep="::")]

dummy = unlist(strsplit(pQTLs$variant_id_hg19,split = ":"))
pos_hg19 = dummy[seq(2,length(dummy),4)]
pQTLs[,pos_hg19:=as.numeric(pos_hg19)]
pQTLs[,rs_id:=dummy[seq(1,length(dummy),4)]]
pQTLs[rs_id %in% c(1:22),rs_id:=paste(rs_id,pos_hg19,sep=":")]
pQTLs[,chrPos_hg19:=paste(chr,pos_hg19,sep="::")]

leadpQTLs<-pQTLs[hierarch_fdr5proz==T,.SD[pval==min(pval)],by=.(gene,setting)]
leadpQTLs[,table(gene,setting)]
dummy<-pQTLs[variant_id_hg19 == "rs78811301:43740523:G:A" & setting=="males",]
leadpQTLs = leadpQTLs[!(gene =="TFF3" & setting =="males"),]
leadpQTLs = rbind(leadpQTLs,dummy)

leadpQTLs[,table(setting)]

setcolorder(leadpQTLs,names(pQTLs_a))
head(leadpQTLs)

#' # Get GTEx eQTL data ####
#' ***
myeQTLs<-dir(path = path_GTExv8,pattern = "GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8")

dumTab1<-foreach (i=c(1:length(myeQTLs)))%do%{
  #i=1
  myTissue<-gsub(".allpairs.txt.gz","",myeQTLs[i])
  myTissue<-gsub("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_","",myTissue)
  time1 = Sys.time()
  
  # Loading data
  message("Loading eQTL data ",myTissue," number ",i," of ",length(myeQTLs))
  data1<-fread(paste0(basicpath, "/01_daten/2007_GTEx_v8/",myeQTLs[i]))
  time2 = Sys.time()
  message("          Finished loading eQTL data in ",round(difftime(time2, time1, tz,units = "min"),2)," minutes")
  
  # Change GeneID
  data1[,ENSG:=gsub(gene_id, pattern = "\\..*", replacement = "")]
  # data2<-copy(data1)
  # data1 = copy(data2)
  
  # Filtering data
  message("          Filtering eQTL data ",myTissue," number ",i," of ",length(myeQTLs))
  filt<-is.element(data1$ENSG,myTab$ENSG)
  data1<-data1[filt,]
  time3 = Sys.time()
  message("          Finished filtering eQTL data in ",round(difftime(time3, time2, tz,units = "min"),2)," minutes")
  
  # Harmonisation of column names
  message("          Harmonizing column names for ",myTissue," number ",i," of ",length(myeQTLs))
  # table(is.element(data1$variant_id,pQTLs_a$variant_id_hg38))
  # table(is.element(pQTLs_a$variant_id_hg38,data1$variant_id))
  dummy<-unlist(strsplit(data1$variant_id,"_"))
  chr<-dummy[seq(1,length(dummy),by=5)]
  chr<-gsub("chr","",chr)
  chr<-as.numeric(chr)
  pos<-as.numeric(dummy[seq(2,length(dummy),by=5)])
  other_allele<-dummy[seq(3,length(dummy),by=5)]
  effect_allele<-dummy[seq(4,length(dummy),by=5)]
  
  data1[,chr:=chr,]
  data1[,pos:=pos,]
  data1[,other_allele:=other_allele,]
  data1[,effect_allele:=effect_allele,]
  data1[,n_samples:=round(ma_count/maf/2,1)]
  data1[,chrPos:=paste(chr,pos,sep="::")]
  setnames(data1,"slope_se","se")
  setnames(data1,"pval_nominal","pval")
  setnames(data1,"slope","beta")
  
  # Add gene name and cytoband
  matched<-match(data1$ENSG,myTab$ENSG)
  data1[,gene:=myTab[matched,Gene]]
  data1[,cyto:=myTab[matched,Cytoband]]
  data1[,tissue:=myTissue]
  data1[,chrPosGene:=paste(chr,pos,gene,sep=":")]
  
  # data set 1: no duplicates, MAF>=0.01
  # data set 2: SNPs in pQTLs
  table(duplicated(data1$chrPosGene),is.element(data1$variant_id,pQTLs$variant_id_hg38))
  dups = data1[duplicated(data1$chrPosGene),]
  dups2 = data1[is.element(chrPosGene,dups$chrPosGene),]
  dups3 = dups2[is.element(variant_id,pQTLs$variant_id_hg38),]
  data1 = rbind(data1[!is.element(chrPosGene,dups$chrPosGene),],dups3)
  data1 = data1[maf>=0.01,]
  
  filt<-is.element(data1$variant_id,pQTLs$variant_id_hg38)
  table(filt)
  data3 = data1[filt,]
  
  table(is.element(pQTLs$variant_id_hg38,data3$variant_id))
  table(is.element(data3$variant_id,pQTLs$variant_id_hg38))
  
  x1 = dim(data1)[1]
  x2 = dim(data3)[1]
  message("          Saving ",x1, " unfiltered SNPs, and ",x2," filtered SNPs")
  
  # Save eqtl data
  outfn1<-paste0("../temp/GTEx_v8_unfiltered/",myTissue,".RData")
  outfn2<-paste0("../temp/GTEx_v8_filtered/",myTissue,".RData")
  save(data1,file=outfn1)
  save(data3,file=outfn2)
  
  # Check if lead SNP is in data
  check1<-table(is.element(myTab$ENSG, data1$ENSG))
  check2<-table(is.element(myTab$ENSG, data3$ENSG))
  message("          Found ", check1[2], " of 92 candidate genes in eQTL data (all)")
  message("          Found ", check2[2], " of 92 candidate genes in eQTL data (matched)")
  
  # Get lead pQTLs
  filt<-is.element(data3$chrPosGene,leadpQTLs$chrPosGene)
  data6<-data3[filt,]
  data6[,comment := "matched_pQTL"]
  
  # Get lead eQTLs
  data7<-data3[,.SD[pval==min(pval)],by=.(gene)]
  data7[,comment := "best_eQTL"]
  
  # Merge and save
  data8<-rbind(data6,data7)
  outfn3<-paste0("../temp/GTEx_v8_matched_best/",myTissue,".RData")
  save(data8,file=outfn3)
  message("\n          Finished tissue ",myTissue," in ",round(difftime(Sys.time(), time1, tz,units = "min"),2)," minutes\n")
  
  # remove temp and data1
  tmp_dir <- tempdir()
  files <- list.files(tmp_dir, full.names = T)
  file.remove(files)
  file.exists(files)
  rm(data1)
  
  #return something
  data8[,x1 :=x1]
  data8[,x2 :=x2]
  data8
}
eQTL_GTEx<-rbindlist(dumTab1)
head(eQTL_GTEx)
names(eQTL_GTEx)

x1 = eQTL_GTEx[,unique(x1),by=tissue]
x1
x2 = eQTL_GTEx[,unique(x2),by=tissue]
x2
save(x1,x2,file = "../temp/control_freadReadings.RData")

eQTL_GTEx[,x1:=NULL]
eQTL_GTEx[,x2:=NULL]
save(eQTL_GTEx,file="../results/04_a_eQTL_lookup_allTissues.RData")

#' # sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in hours): " ,round(difftime(Sys.time(), time0, tz,units = "hours"),2))
