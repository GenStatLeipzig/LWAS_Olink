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
pQTLs[,SNPGene:=paste(variant_id_hg38,gene,sep=":")]

eQTL_LIFE = fread("../exampleData/eQTLs_LIFE_hg19.txt")

leadpQTLs<-pQTLs[hierarch_fdr5proz==T,.SD[pval==min(pval)],by=.(gene,setting)]
leadpQTLs[,table(gene,setting)]
dummy<-pQTLs[variant_id_hg19 == "rs78811301:43740523:G:A" & setting=="males",]
leadpQTLs = leadpQTLs[!(gene =="TFF3" & setting =="males"),]
leadpQTLs = rbind(leadpQTLs,dummy)
leadpQTLs[,table(setting)]
setcolorder(leadpQTLs,names(pQTLs_a))
leadpQTLs[,SNPGeneSetting:=paste(variant_id_hg38,gene,setting,sep=":")]
leadpQTLs_a = leadpQTLs[setting == "combined"]
leadpQTLs_m = leadpQTLs[setting == "males"]
leadpQTLs_f = leadpQTLs[setting == "females"]

#' # Get GTEx eQTL data ####
#' ***
#' Please note: I commented out some lines of code as the interim results were saved in the *../temp/GTEx_v8_filtered directory.
#' 
#' In addition, I assume the GTEx data is named *GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_TISSUE.allpairs.txt.gz*. If the files are named differently, please change code lines 60 & 61 accordingly. 
#' 
myeQTLs<-dir(path = path_GTExv8,pattern = "GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8")

dumTab1<-foreach (i=c(1:length(myeQTLs)))%do%{
  #i=1
  myTissue<-gsub(".allpairs.txt.gz","",myeQTLs[i])
  myTissue<-gsub("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_","",myTissue)
  time1 = Sys.time()

  # # Loading data
  # message("Loading eQTL data ",myTissue," number ",i," of ",length(myeQTLs))
  # data0<-fread(paste0(basicpath, "/01_daten/2007_GTEx_v8/",myeQTLs[i]))
  # time2 = Sys.time()
  # message("          Finished loading eQTL data in ",round(difftime(time2, time1, tz,units = "min"),2)," minutes")
  # 
  # # Change GeneID
  # data0[,ENSG:=gsub(gene_id, pattern = "\\..*", replacement = "")]
  # # data2<-copy(data1)
  # # data1 = copy(data2)
  # 
  # # Filtering data
  # message("          Filtering eQTL data ",myTissue," number ",i," of ",length(myeQTLs))
  # filt<-is.element(data0$ENSG,myTab$ENSG)
  # data0<-data0[filt,]
  # time3 = Sys.time()
  # message("          Finished filtering eQTL data in ",round(difftime(time3, time2, tz,units = "min"),2)," minutes")
  # 
  # # Harmonisation of column names
  # message("          Harmonizing column names for ",myTissue," number ",i," of ",length(myeQTLs))
  # dummy<-unlist(strsplit(data0$variant_id,"_"))
  # chr<-dummy[seq(1,length(dummy),by=5)]
  # chr<-gsub("chr","",chr)
  # chr<-as.numeric(chr)
  # pos<-as.numeric(dummy[seq(2,length(dummy),by=5)])
  # other_allele<-dummy[seq(3,length(dummy),by=5)]
  # effect_allele<-dummy[seq(4,length(dummy),by=5)]
  # 
  # data0[,chr:=chr,]
  # data0[,pos:=pos,]
  # data0[,other_allele:=other_allele,]
  # data0[,effect_allele:=effect_allele,]
  # data0[,n_samples:=round(ma_count/maf/2,1)]
  # setnames(data0,"slope_se","se")
  # setnames(data0,"pval_nominal","pval")
  # setnames(data0,"slope","beta")
  # 
  # # Add gene name and cytoband
  # matched<-match(data0$ENSG,myTab$ENSG)
  # data0[,gene:=myTab[matched,Gene]]
  # data0[,cyto:=myTab[matched,Cytoband]]
  # data0[,tissue:=myTissue]
  # data0[,SNPGene:=paste(variant_id,gene,sep=":")]
  # 
  # # Filtering: no duplicates, MAF>=0.01, SNPs in pQTLs
  # table(duplicated(data0$SNPGene),is.element(data0$variant_id,pQTLs$variant_id_hg38))
  # dups = data0[duplicated(data0$SNPGene),]
  # dups2 = data0[is.element(SNPGene,dups$SNPGene),]
  # dups3 = dups2[is.element(variant_id,pQTLs$variant_id_hg38),]
  # 
  # data1 = rbind(data0[!is.element(SNPGene,dups$SNPGene),],dups3)
  # data1 = data1[maf>=0.01,]
  # data1 = data1[is.element(variant_id,pQTLs$variant_id_hg38),]
  # 
  # table(is.element(pQTLs$variant_id_hg38,data1$variant_id))
  # table(is.element(data1$variant_id,pQTLs$variant_id_hg38))
  # 
  # x1 = dim(data1)[1]
  # message("          Saving ",x1," filtered SNPs")
  # 
  # # Save eqtl data
  outfn1<-paste0("../temp/GTEx_v8_filtered/",myTissue,".RData")
  # save(data1,file=outfn1)
  load(outfn1)
  
  # Check if lead SNP is in data
  check1<-table(is.element(myTab$ENSG, data1$ENSG))
  message("          Found ", check1[2], " of 92 candidate genes in eQTL data (matched)")
  
  # Get lead pQTLs
  filt<-is.element(data1$SNPGene,leadpQTLs$SNPGene)
  data2a<-data1[filt,]
  data2a[,comment := "matched_pQTL"]
  
  # Get lead eQTLs
  data2b<-data1[,.SD[pval==min(pval)],by=.(gene)]
  data2b[,comment := "best_eQTL"]
  
  # Merge and save
  data2<-rbind(data2a,data2b)
  outfn2<-paste0("../temp/GTEx_v8_matched_best/",myTissue,".RData")
  save(data2,file=outfn2)
  message("\n          Finished tissue ",myTissue," in ",round(difftime(Sys.time(), time1, tz,units = "min"),2)," minutes\n")
  
  # remove temp and data1
  tmp_dir <- tempdir()
  files <- list.files(tmp_dir, full.names = T)
  file.remove(files)
  file.exists(files)
  rm(data1)
  
  #return something
  data2
}
eQTL_GTEx<-rbindlist(dumTab1)
head(eQTL_GTEx)
names(eQTL_GTEx)
save(eQTL_GTEx,file="../results/04_a_eQTL_lookup_allTissues.RData")

matched = match(eQTL_GTEx$SNPGene, pQTLs$SNPGene)
table(is.na(matched))
table(eQTL_GTEx$effect_allele ==pQTLs$effect_allele[matched])
eQTL_GTEx[,variant_id_hg19:=pQTLs[matched,variant_id_hg19]]
eQTL_GTEx[,eaf:=maf]
plot(eQTL_GTEx$eaf,pQTLs$eaf[matched])
filt = pQTLs_a$eaf[matched]> 0.5
table(filt)
eQTL_GTEx[filt,eaf:=1-maf]
eQTL_GTEx[is.na(variant_id_hg19),eaf:=NA]
plot(eQTL_GTEx$eaf,pQTLs_a$eaf[matched])

#' # Merge eQTL data ####
#' ***
names(eQTL_LIFE)
names(eQTL_GTEx)

names(eQTL_LIFE)[!is.element(names(eQTL_LIFE),names(eQTL_GTEx))]
setnames(eQTL_GTEx,"pos","pos_hg38")
setnames(eQTL_GTEx,"variant_id","variant_id_hg38")
eQTL_LIFE[,SNPGene:=paste(variant_id_hg38,gene,sep=":")]
names(eQTL_LIFE)[!is.element(names(eQTL_LIFE),names(eQTL_GTEx))]
names(eQTL_GTEx)[!is.element(names(eQTL_GTEx),names(eQTL_LIFE))]

colsOut<-setdiff(colnames(eQTL_GTEx),names(eQTL_LIFE))
eQTL_GTEx[,get("colsOut"):=NULL]
colsOut<-setdiff(colnames(eQTL_LIFE),names(eQTL_GTEx))
eQTL_LIFE[,get("colsOut"):=NULL]

setcolorder(eQTL_GTEx,names(eQTL_LIFE))
table(names(eQTL_LIFE) == names(eQTL_GTEx))
eQTLs<-rbind(eQTL_GTEx,eQTL_LIFE)

length(unique(eQTLs$gene))
myGenes = myTab[GWAS_sig==T | GWAS_sig_male==T | GWAS_sig_female==T,Gene]
filt<-is.element(eQTLs$gene,myGenes)
table(filt)
eQTLs<-eQTLs[filt,]
save(eQTLs,file="../results/04_b_eQTL_GTEx_LIFE.RData")

#' # Split best and matched eQTL ####
#' ***
names(eQTLs)

table(is.element(eQTLs$SNPGene,leadpQTLs$SNPGene))
table(is.element(leadpQTLs$SNPGene,eQTLs$SNPGene))

data_matched<-eQTLs[comment=="matched_pQTL"]
data_besteQTL<-eQTLs[comment=="best_eQTL"]

#' # Map eQTLs to lead pQTLs ###
#' ***
myMappingFunction = function(data1,data2,data3){
  # data1 = leadpQTLs_a
  # data2 = data_matched
  # data3 = data_besteQTL
  
  # step 0: eQTL data filter
  data2 = data2[SNPGene %in% data1$SNPGene,]
  data2[,dumID:=paste(gene,tissue,sep="::")]
  data3[,dumID:=paste(gene,tissue,sep="::")]
  data3 = data3[dumID %in% data2$dumID,]
  data3 = data3[!duplicated(dumID),]
  
  # step 1: Matched eQTL
  matched1<-match(data2$SNPGene,data1$SNPGene)
  table(is.na(matched1))
  leadpQTLs1<-data1[matched1,]
  
  table(leadpQTLs1$SNPGene==data2$SNPGene)
  table(leadpQTLs1$effect_allele==data2$effect_allele, 
        leadpQTLs1$other_allele==data2$other_allele)
  filt<-leadpQTLs1$effect_allele!=data2$effect_allele
  table(filt)
  
  leadpQTLs1[,SNPID_GTEx:=data2[,variant_id_hg38]]
  leadpQTLs1[,beta_GTEx:=data2[,beta]]
  leadpQTLs1[filt,beta_GTEx:=(-1)*beta_GTEx]
  leadpQTLs1[,se_GTEx:=data2[,se]]
  leadpQTLs1[,pval_GTEx:=data2[,pval]]
  leadpQTLs1[,n_GTEx:=data2[,n_samples]]
  leadpQTLs1[,tissue_GTEx:=data2[,tissue]]
  leadpQTLs1[,same_direction:=NA]
  leadpQTLs1[beta>0 & beta_GTEx>0,same_direction:=T]
  leadpQTLs1[beta<0 & beta_GTEx<0,same_direction:=T]
  leadpQTLs1[beta>0 & beta_GTEx<0,same_direction:=F]
  leadpQTLs1[beta<0 & beta_GTEx>0,same_direction:=F]
  head(leadpQTLs1)
  table(leadpQTLs1$same_direction)
  leadpQTLs1[,dumID:=paste(gene,tissue_GTEx,setting,sep=":")]
  table(duplicated(leadpQTLs1$dumID))
  
  # step 2: Best eQTL
  matched2<-match(data3$gene,data1$gene)
  table(is.na(matched2))
  leadpQTLs2<-data1[matched2,]
  leadpQTLs2[,SNPID_GTEx_best:=data3[,variant_id_hg38]]
  leadpQTLs2[,SNPID_GTEx_best_hg19:=data3[,variant_id_hg19]]
  leadpQTLs2[,beta_GTEx_best:=data3[,beta]]
  leadpQTLs2[,se_GTEx_best:=data3[,se]]
  leadpQTLs2[,pval_GTEx_best:=data3[,pval]]
  leadpQTLs2[,n_GTEx_best:=data3[,n_samples]]
  leadpQTLs2[,tissue_GTEx_best:=data3[,tissue]]
  leadpQTLs2[,effect_allele_GTEx_best:=data3[,effect_allele]]
  leadpQTLs2[,other_allele_GTEx_best:=data3[,other_allele]]
  head(leadpQTLs2)
  leadpQTLs2[,dumID:=paste(gene,tissue_GTEx_best,setting,sep=":")]
  table(duplicated(leadpQTLs2$dumID))
  
  # step 3: merge matched and best eQTL by dumID
  matched3<-match(leadpQTLs1$dumID,leadpQTLs2$dumID)
  table(is.na(matched3))
  leadpQTLs2<-leadpQTLs2[matched3,]
  table(leadpQTLs2$dumID==leadpQTLs1$dumID)
  
  leadpQTLs3 = copy(leadpQTLs1)
  leadpQTLs3[,SNPID_GTEx_best:=leadpQTLs2[,SNPID_GTEx_best]]
  leadpQTLs3[,SNPID_GTEx_best_hg19:=leadpQTLs2[,SNPID_GTEx_best_hg19]]
  leadpQTLs3[,beta_GTEx_best:=leadpQTLs2[,beta_GTEx_best]]
  leadpQTLs3[,se_GTEx_best:=leadpQTLs2[,se_GTEx_best]]
  leadpQTLs3[,pval_GTEx_best:=leadpQTLs2[,pval_GTEx_best]]
  leadpQTLs3[,n_GTEx_best:=leadpQTLs2[,n_GTEx_best]]
  leadpQTLs3[,tissue_GTEx_best:=leadpQTLs2[,tissue_GTEx_best]]
  head(leadpQTLs3)
  table(leadpQTLs3$SNPID_GTEx==leadpQTLs3$SNPID_GTEx_best)
  leadpQTLs3[,LeadIsBest:=leadpQTLs3$SNPID_GTEx==leadpQTLs3$SNPID_GTEx_best]
  
  # step 4: return merged object
  return(leadpQTLs3)
  
}

leadpQTLs_a2 = myMappingFunction(data1 = leadpQTLs_a,data2 = data_matched,data3 = data_besteQTL)
leadpQTLs_m2 = myMappingFunction(data1 = leadpQTLs_m,data2 = data_matched,data3 = data_besteQTL)
leadpQTLs_f2 = myMappingFunction(data1 = leadpQTLs_f,data2 = data_matched,data3 = data_besteQTL)

leadpQTLs4 = rbind(leadpQTLs_a2,leadpQTLs_m2,leadpQTLs_f2)
table(duplicated(leadpQTLs4$dumID))
table(is.element(leadpQTLs$SNPGeneSetting,leadpQTLs4$SNPGeneSetting))
leadpQTLs[!is.element(SNPGeneSetting,leadpQTLs4$SNPGeneSetting)]

save(leadpQTLs4,file="../results/04_c_eQTL_lookup_Summary.RData")

#' # LD estimation ####
#' 
leadpQTLs4[,table(is.na(SNPID_GTEx_best_hg19))]

dumTab5<-data.table(pQTL=leadpQTLs4$variant_id_hg19,
                    eQTL=leadpQTLs4$SNPID_GTEx_best_hg19)
dumTab5 = dumTab5[!is.na(eQTL),]
dumTab5 = dumTab5[eQTL != pQTL,]
dumTab5[,out:=paste0("../temp/PLINK_LD/",pQTL,"__",eQTL)]
dumTab5[,out:=gsub(":","_",out)]
table(duplicated(dumTab5$out))
dumTab5 = dumTab5[!duplicated(out),]

registerDoMC(cores=10)

dumTab6<-foreach(i = 1:dim(dumTab5)[1])%dopar%{
  #dumTab6<-foreach (i=c(1:3))%do%{
  #i=1
  message("Working on i=",i,", of ",dim(dumTab5)[1])
  myRow<-dumTab5[i,]
  eQTL<-myRow[,eQTL]
  pQTL<-myRow[,pQTL]
  outfile<-myRow[,out]
    
  if(!file.exists(paste0(outfile,".log"))){
    # plink call for LD
    myCall = paste(path_plink2, 
                   " --bfile ", path_1000Genomes,
                   " --ld ", pQTL, " ",eQTL,
                   " --out" , outfile)
    myCall
    system(myCall)
  }
    
  # load PLINK output
  line <- readLines(paste0(outfile,".log"))
  res<-line[grep("D' = ",line)]
  res2<-unlist(strsplit(x = res,split = "D' = "))
  res2<-gsub(" ","",res2)
  res2<-gsub(".*=","",res2)
  res2<-as.numeric(res2)
  
  myRow[,LD_r2:=res2[1]]
  myRow[,LD_Dprime:=res2[2]]
  myRow
}
dumTab6<-rbindlist(dumTab6)
table(dumTab6$pQTL==dumTab5$pQTL)
table(dumTab6$eQTL==dumTab5$eQTL)

dumTab6[is.na(LD_r2),]
dumTab6[,matchID := paste(pQTL,eQTL,sep="__")]
leadpQTLs4[,matchID := paste(variant_id_hg19,SNPID_GTEx_best_hg19,sep="__")]

matched4<-match(leadpQTLs4$matchID,dumTab6$matchID)
table(leadpQTLs4$matchID==dumTab6$matchID[matched4])
leadpQTLs4[,LD_r2:=dumTab6[matched4,LD_r2]]
leadpQTLs4[,LD_Dprime:=dumTab6[matched4,LD_Dprime]]
leadpQTLs4[LeadIsBest==T,LD_r2:=1]
leadpQTLs4[LeadIsBest==T,LD_Dprime:=1]
head(leadpQTLs4)
leadpQTLs4[,matchID := NULL]

save(leadpQTLs4,file="../results/04_c_eQTL_lookup_Summary_LD.RData")

#' # Check effect directions ####
#' ***
leadpQTLs5 = copy(leadpQTLs4)
leadpQTLs5 = leadpQTLs5[setting == "combined",]
leadpQTLs5 = leadpQTLs5[tissue_GTEx != "Whole_Blood_LIFE",]

z1<-leadpQTLs5[pval_GTEx<=0.05 & same_direction==F & !is.na(same_direction),.N,by=gene]
z2<-leadpQTLs5[pval_GTEx<=0.05 & same_direction==T & !is.na(same_direction),.N,by=gene]
table(is.element(z1$gene,z2$gene))
z3<-unique(c(z1$gene,z2$gene))
z3<-data.table(gene=z3)
z3[,N_same:=z2[match(z3$gene,gene),N]]
z3[,N_diff:=z1[match(z3$gene,gene),N]]
z3[is.na(N_same),N_same:=0]
z3[is.na(N_diff),N_diff:=0]
z3[,N:=N_same + N_diff]
z3[,ratio:=N_diff/N]
z3[ratio>=0.75,]

#' There are 13 proteins with different effect direction compared to eQTLs (in the combined setting only). 
#'  

#' # Plotting ####
#' ***
#' Plotting Figure 2 ...
myPlotData<-copy(leadpQTLs4)
myPlotData<-myPlotData[setting =="combined",]
myPlotData<-myPlotData[tissue_GTEx != "Whole_Blood_LIFE",]
myPlotData<-myPlotData[pval_GTEx<=0.05,]
length(unique(myPlotData$gene))
table(is.element(myPlotData$gene,z3[ratio>=0.75,gene]))
table(is.element(z3[ratio>=0.75,gene],myPlotData$gene))

ylim1<-min(myPlotData$beta,na.rm = T)-0.1
ylim2<-max(myPlotData$beta,na.rm = T)+0.1
xlim1<-min(myPlotData$beta_GTEx,na.rm = T)-0.1
xlim2<-max(myPlotData$beta_GTEx,na.rm = T)+0.1

myPlotData[,tissue_gene:= paste(tissue_GTEx,gene,sep=":")]
x5<-myPlotData[same_direction==F,.SD[pval_GTEx==min(pval_GTEx)],by=.(gene)]
x5<-x5[is.element(gene,z3[ratio>=0.75,gene])]

matched1<-match(myPlotData$gene,x5$gene)
table(is.na(matched1))
matched2<-match(myPlotData$tissue_gene,x5$tissue_gene)
table(is.na(matched2))

myPlotData[,gene_color:=""]
myPlotData[,flag2:=F]
myPlotData[!is.na(matched1),gene_color:=gene]
myPlotData[!is.na(matched1),flag2:=T]
myPlotData[,gene_label:=""]
myPlotData[!is.na(matched2),gene_label:=gene]
table(myPlotData$same_direction,myPlotData$flag2)

myPlot2<-ggplot(myPlotData[flag2==T,], aes(x=beta_GTEx, y=beta, color=gene_color,label=gene_label)) +
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1.25)+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_point(data=myPlotData[flag2==F,], aes(x=beta_GTEx, y=beta),col="black",size=1.5,alpha=0.5)+
  geom_point(size=2.5)+ 
  theme_bw(base_size = 10)+
  ggtitle(label = "Scatter plot of effect estimates")+
  theme(plot.title = element_text(hjust = 0.5))+ 
  xlab("Effect on mRNA levels")+
  ylab("Effect on protein levels")+
  guides(size="none",fill="none",color="none")+
  geom_label_repel(data = subset(myPlotData, flag2==T),
                   aes(fill=gene_label),
                   box.padding = 1.15,
                   max.overlaps = 80,
                   show.legend=F,color="black")
myPlot2

tiff(filename = "../results/04_ScatterPlot_GEvsPE.tif", 
     width = 1200, height = 900, res=125, compression = 'lzw')
myPlot2
dev.off()

#' # Update myTab ###
#' ****
#' Merge meta data of GTEx onto myTab
leadpQTLs_GTEx = copy(leadpQTLs4)
leadpQTLs_GTEx <- leadpQTLs_GTEx[!grepl("LIFE",tissue_GTEx) & pval_GTEx<0.05,]

x1<-leadpQTLs_GTEx[setting =="combined",.N,by=c("gene")]
x2<-leadpQTLs_GTEx[setting =="combined" & same_direction==T,.N,by=c("gene")]
x3<-leadpQTLs_GTEx[setting =="combined" & same_direction==F,.N,by=c("gene")]

head(myTab)
filt1<-is.element(myTab$Gene,x1$gene)
matched2<-match(myTab$Gene,x2$gene)
matched3<-match(myTab$Gene,x3$gene)
myTab[,n_sig_eQTLs_sameDir:=x2[matched2,N]]
myTab[,n_sig_eQTLs_diffDir:=x3[matched3,N]]
myTab[filt1 & is.na(n_sig_eQTLs_sameDir),n_sig_eQTLs_sameDir:=0]
myTab[filt1 & is.na(n_sig_eQTLs_diffDir),n_sig_eQTLs_diffDir:=0]

dummy<-myTab[,n_sig_eQTLs_diffDir/(n_sig_eQTLs_sameDir + n_sig_eQTLs_diffDir)]
hist(dummy, breaks=20)
filt2<-dummy>0.8 & !is.na(dummy)
myTab[filt2,]

#' Merge meta data of LIFE onto myTab
leadpQTLs_LIFE = copy(leadpQTLs4)
leadpQTLs_LIFE <- leadpQTLs_LIFE[grepl("LIFE",tissue_GTEx) & pval_GTEx<0.05,]
x4<-leadpQTLs_LIFE[setting =="combined" & same_direction==F,.N,by=c("gene")]
myTab[,n_sig_eQTLs_diffDir_LIFE:=NA]
myTab[GE_ok=="yes" & GWAS_sig==T,n_sig_eQTLs_diffDir_LIFE:=F]
myTab[Gene %in% x4$gene,n_sig_eQTLs_diffDir_LIFE:=T]
myTab[filt2 & n_sig_eQTLs_diffDir_LIFE==T,]

#' Save myTab
save(myTab,file="../temp/04_OlinkGenes.RData")

#' # sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
