#' ---
#' title: "Co-localization"
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
#' I want to perform a co-localization test between pQTLs and eQTLs for all significant biomarkers in all settings. 
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
load("../temp/05_OlinkGenes.RData")
load("../results/02_pQTLs_all_cis.RData")
load("../results/02_pQTLs_males_cis.RData")
load("../results/02_pQTLs_females_cis.RData")

pQTLs = rbind(pQTLs_a,pQTLs_m,pQTLs_f)
pQTLs[,SNPGene:=paste(variant_id_hg38,gene,sep=":")]
pQTLs[,SNPGeneSetting:=paste(variant_id_hg38,gene,setting,sep=":")]
pQTLs[,GeneSetting:=paste(gene,setting,sep=":")]

myGenes_combined = myTab[GWAS_sig == T,paste(Gene,"combined",sep=":")]
myGenes_males = myTab[GWAS_sig_male == T,paste(Gene,"males",sep=":")]
myGenes_females = myTab[GWAS_sig_female == T,paste(Gene,"females",sep=":")]
myGeneSettings = c(myGenes_combined,myGenes_females,myGenes_males)
myGeneSettings = myGeneSettings[order(myGeneSettings)]
pQTLs<-pQTLs[GeneSetting %in% myGeneSettings,]

length(unique(pQTLs$gene))
length(unique(pQTLs$GeneSetting))

#' # Define function ####
#' ***
myColocFunction<-function(tab1,tab2,trait1,trait2,locus,locus_name,
                          plot=F,SD=F,sd1=0,sd2=0){
  # x1<-myColocFunction(myTab_Olink2, data2,trait1,trait2,myLocus,myGene)
  # tab1<-myTab_Olink_filtered
  # tab2<-data3
  # trait1<-trait1
  # trait2<-trait2
  # locus<-myLocus
  # locus_name<-myGene
  # plot=T
  # SD=T
  # sd1 = sd_m
  # sd2 = sd_f
  
  table(is.element(tab1$variant_id_hg38,tab2$variant_id_hg38))
  table(is.element(tab2$variant_id_hg38,tab1$variant_id_hg38))
  
  dumTab1<-tab1[is.element(tab1$variant_id_hg38,tab2$variant_id_hg38),]
  dumTab2<-tab2[is.element(tab2$variant_id_hg38,tab1$variant_id_hg38),]
  
  matched<-match(dumTab1$variant_id_hg38,dumTab2$variant_id_hg38)
  dumTab2<-dumTab2[matched,]
  stopifnot(dumTab1$variant_id_hg38==dumTab2$variant_id_hg38)
  
  if(SD==T){
    my_res<- coloc.abf(dataset1=list(beta=dumTab1$beta,sdY = sd1,
                                     varbeta=(dumTab1$se)^2, MAF=dumTab1$maf,
                                     N=dumTab1$n_samples,snp=dumTab1$variant_id_hg38,
                                     type="quant"),
                       dataset2=list(beta=dumTab2$beta, sdY = sd2,
                                     varbeta=(dumTab2$se)^2, MAF=dumTab2$maf,
                                     N=dumTab2$n_samples,snp=dumTab2$variant_id_hg38,
                                     type="quant")) 
  }else{
    my_res<- coloc.abf(dataset1=list(beta=dumTab1$beta,
                                     varbeta=(dumTab1$se)^2, 
                                     N=dumTab1$n_samples,snp=dumTab1$variant_id_hg38,
                                     type="quant"),
                       dataset2=list(beta=dumTab2$beta, 
                                     varbeta=(dumTab2$se)^2, 
                                     N=dumTab2$n_samples,snp=dumTab2$variant_id_hg38,
                                     type="quant"),
                       MAF=dumTab1$maf)
  }
  
  
  my_res2<-my_res[[1]]
  
  #plotting
  if(plot==T){
    myXlab<-locus
    
    dumTab1[,plot(pos_hg38, -log10(pval),col = rgb(0,0,1,0.3),pch=19,
                  xlab = myXlab, 
                  ylab = "")]
    axis(side = 2, col = "blue", col.ticks = "blue",col.axis="blue")
    mtext(side = 2, line = 2.1, 
          bquote(.(trait1) ~ -log[10](italic(p-value))), col = "blue", cex = 0.7)
    
    par(new = T)
    dumTab2[,plot(pos_hg38, -log10(pval),axes=F, xlab=NA, ylab=NA,  
                  col = rgb(1,0,0,0.2), pch = 17)]
    axis(side = 4, col = "red", col.ticks = "red",col.axis="red")
    mtext(side = 4, line = 2.1, 
          bquote(.(trait2) ~ -log[10](italic(p-value))), col = "red", cex = 0.7)
    
    #legendtext = paste0("APOB gene\nmore GWAS-risk,\nmore expression") 
    #legend("topright",legend=legendtext,  text.col=c("black"), bty="n")
    
    title(paste0(trait1," vs. ",trait2," for ", locus_name), cex.main = 1)
    
    mtext(paste0("H0 | H1 | ... | H4:  ", paste(formatC(round(my_res2[2:6],3), digits=3, format="f" ), collapse = " | ")),cex = 0.7)
  }
  
  return(my_res2)
}

#' # Loop 1: pQTLs vs eQTLs ####
#' ***
par(mfrow=c(1,1))
par(mar=c(5.1 ,4.1 ,5.1 ,3.5))

datalist_trait_fn2<-dir("../temp/GTEx_v8_filtered/")
traitlist_trait2<-gsub(".RData","",datalist_trait_fn2)
traitlist_trait2<-gsub("_"," ",traitlist_trait2)

genelist<-unique(pQTLs$GeneSetting)

registerDoMC(cores=10)

dumTab<-foreach(i=1:length(traitlist_trait2))%dopar%{
  #dumTab<-foreach(i=1:5)%do%{
  #i=1
  trait2<-traitlist_trait2[i]
  message("\nWorking on trait ",trait2,"\n")
  
  # Load eQTL data
  loaded1<-load(paste0("../temp/GTEx_v8_filtered/",datalist_trait_fn2[i]))
  data2<- get(loaded1)
  
  # Filt eQTL data in relevant columns and rows
  filt<-!is.na(data2$beta) & !is.na(data2$se) & !is.na(data2$n_samples)
  data2<-data2[filt,]
  setnames(data2,"variant_id","variant_id_hg38")
  setnames(data2,"pos","pos_hg38")
  
  dumTab2<-foreach(j=1:length(genelist))%do%{
    #j=1
    myToDo = genelist[j]
    myGene<-unlist(strsplit(myToDo,":"))[1]    
    mySetting<-unlist(strsplit(myToDo,":"))[2]    
    myLocus<-myTab[Gene==myGene,Cytoband]
    trait1<-myTab[Gene==myGene,Parameter]
    message("Working on locus ",myLocus," of phenotype ",trait1, " with candidate gene(s) ",myGene," in setting ",mySetting," \n(",j," of ",length(genelist),")")
    
    myTab_Olink_filtered<-copy(pQTLs)
    myTab_Olink_filtered<-myTab_Olink_filtered[pheno==trait1,]
    myTab_Olink_filtered<-myTab_Olink_filtered[setting==mySetting,]
    myTab_Olink_filtered<-myTab_Olink_filtered[!duplicated(variant_id_hg19),]
    stopifnot(duplicated(myTab_Olink_filtered$variant_id_hg19)==F)
    stopifnot(length(unique(myTab_Olink_filtered$cyto))==1)
    
    data3<-copy(data2)
    data3<-data3[is.element(gene,myGene)]
    data3[,chr := paste0("chr",chr)]
    
    if(dim(data3)[1]==0){
      x4<-data.table(nsnps=0,PP.H0.abf=NA,PP.H1.abf=NA,
                     PP.H2.abf=NA,PP.H3.abf=NA,PP.H4.abf=NA,
                     comment="no eQTL data available")
    }else if(unique(data3$chr)!=unique(myTab_Olink_filtered$chr)){
      x4<-data.table(nsnps=0,PP.H0.abf=NA,PP.H1.abf=NA,
                     PP.H2.abf=NA,PP.H3.abf=NA,PP.H4.abf=NA,
                     comment="candidate is a trans eQTL - no coloc possible")
    }else {
      x1<-myColocFunction(tab1 = myTab_Olink_filtered, tab2 = data3,trait1 = "PE levels",trait2 = "GE levels",
                          locus = myLocus,locus_name = myGene,plot=F,SD = F)
      x2<-as.data.table(x1)
      x3<-t(x2)
      x4<-as.data.table(x3)
      names(x4)<-names(x1)
      x4[,comment:="good coloc"]
    }
    x4[,locus:=myLocus]
    x4[,gene:=myGene]
    x4[,trait1:=trait1]
    x4[,trait2:=paste0("GE in ",trait2)]
    x4[,setting:=mySetting]
    x4
  }
  Coloc<-rbindlist(dumTab2)
  Coloc
}

Coloc_genewise<-rbindlist(dumTab)
Coloc_genewise[PP.H4.abf>=0.75,table(gene,setting)]
Coloc_genewise[PP.H3.abf>=0.75,table(gene,setting)]
Coloc_genewise[,table(comment)]

#' # Loop 2: sex comparison ####
#' ***
par(mfrow=c(1,1))
par(mar=c(5.1 ,4.1 ,5.1 ,3.5))

genelist<-unique(pQTLs$gene)

dumTab<-foreach(i=1:length(genelist))%do%{
  # dumTab<-foreach(i=1:5)%do%{
  #i=1
  myGene<-genelist[i]
  myLocus = myTab[Gene == myGene,Cytoband]
  myPheno = myTab[Gene == myGene,Parameter]
  message("\nWorking on gene ",myGene," at ",myLocus," \n")
  
  # Get sdY
  sd_m = myTab[Gene == myGene,sd_males]
  sd_f = myTab[Gene == myGene,sd_females]
  
  # Get pQTL data
  data_males = copy(pQTLs_m)
  data_males = data_males[gene == myGene,]
  data_females = copy(pQTLs_f)
  data_females = data_females[gene == myGene,]
  
  # Coloc
  x1<-myColocFunction(tab1 = data_males, tab2 = data_females,trait1 = "males",trait2 = "females",
                      locus = myLocus,locus_name = myGene,plot=T,SD = T,sd1 = sd_m,sd2 = sd_f)
  x2<-as.data.table(x1)
  x3<-t(x2)
  x4<-as.data.table(x3)
  names(x4)<-names(x1)
  x4[,comment:="good coloc"]
  x4[,locus:=myLocus]
  x4[,gene:=myGene]
  x4[,trait1:="males"]
  x4[,trait2:="females"]
  x4
}

Coloc_sex<-rbindlist(dumTab)
Coloc_sex[PP.H4.abf>=0.75,]
Coloc_sex[PP.H3.abf>=0.75,]
Coloc_sex[,table(comment)]

#' # Save ####
#' ***
Coloc_a = copy(Coloc_genewise)
Coloc_a = Coloc_a[setting == "combined"]
dummy1<-Coloc_a[PP.H4.abf>=0.75 ,.N,by=.(gene)]
dummy2<-Coloc_a[PP.H3.abf>=0.75 ,.N,by=.(gene)]
matched1<-match(myTab$Gene,dummy1$gene)
matched2<-match(myTab$Gene,dummy2$gene)
myTab[,Coloc_H4:=NA]
myTab[Gene %in% Coloc_a$gene, Coloc_H4:=F]
myTab[Gene %in% dummy1$gene, Coloc_H4:=T]
myTab[,Coloc_H3:=NA]
myTab[Gene %in% Coloc_a$gene, Coloc_H3:=F]
myTab[Gene %in% dummy2$gene, Coloc_H3:=T]
myTab[,Coloc_H4_NR_tissues:=dummy1[matched1,N]]
myTab[,Coloc_H3_NR_tissues:=dummy2[matched2,N]]

Coloc_m = copy(Coloc_genewise)
Coloc_m = Coloc_m[setting == "males"]
dummy1<-Coloc_m[PP.H4.abf>=0.75 ,.N,by=.(gene)]
dummy2<-Coloc_m[PP.H3.abf>=0.75 ,.N,by=.(gene)]
matched1<-match(myTab$Gene,dummy1$gene)
matched2<-match(myTab$Gene,dummy2$gene)
myTab[,Coloc_H4_males:=NA]
myTab[Gene %in% Coloc_m$gene, Coloc_H4_males:=F]
myTab[Gene %in% dummy1$gene, Coloc_H4_males:=T]
myTab[,Coloc_H3_males:=NA]
myTab[Gene %in% Coloc_m$gene, Coloc_H3_males:=F]
myTab[Gene %in% dummy2$gene, Coloc_H3_males:=T]
myTab[,Coloc_H4_NR_tissues_males:=dummy1[matched1,N]]
myTab[,Coloc_H3_NR_tissues_males:=dummy2[matched2,N]]

Coloc_f = copy(Coloc_genewise)
Coloc_f = Coloc_f[setting == "females"]
dummy1<-Coloc_f[PP.H4.abf>=0.75 ,.N,by=.(gene)]
dummy2<-Coloc_f[PP.H3.abf>=0.75 ,.N,by=.(gene)]
matched1<-match(myTab$Gene,dummy1$gene)
matched2<-match(myTab$Gene,dummy2$gene)
myTab[,Coloc_H4_females:=NA]
myTab[Gene %in% Coloc_f$gene, Coloc_H4_females:=F]
myTab[Gene %in% dummy1$gene, Coloc_H4_females:=T]
myTab[,Coloc_H3_females:=NA]
myTab[Gene %in% Coloc_f$gene, Coloc_H3_females:=F]
myTab[Gene %in% dummy2$gene, Coloc_H3_females:=T]
myTab[,Coloc_H4_NR_tissues_females:=dummy1[matched1,N]]
myTab[,Coloc_H3_NR_tissues_females:=dummy2[matched2,N]]

save(Coloc_genewise,file="../results/06_Coloc_Genewise.RData")
write.table(Coloc_genewise,file="../results/06_Coloc_Genewise.txt",col.names = T,row.names = F, quote = F, sep="\t",dec = ",")

save(Coloc_sex,file="../results/06_Coloc_sexTest.RData")
write.table(Coloc_sex,file="../results/06_Coloc_sexTest.txt",col.names = T,row.names = F, quote = F, sep="\t",dec = ",")

save(myTab,file="../temp/06_OlinkGenes.RData")

#' # Plotting ####
#' ***
Coloc_a = Coloc_a[gene %in% myTab[GWAS_sig == T,Gene]]
Coloc_m = Coloc_m[gene %in% myTab[GWAS_sig_male == T,Gene]]
Coloc_f = Coloc_f[gene %in% myTab[GWAS_sig_female == T,Gene]]

#' ## Build Matrix ####
#' ***
#' I want one row per tissue and one column per gene

myBuildMatrixFunction = function(data){
  #data = Coloc_a
  Tab_PP4<-dcast(data, formula = trait2 ~ gene, value.var = c("PP.H4.abf"), sep = "_")
  M4<-as.matrix(Tab_PP4[,-1])
  
  Tab_PP3<-dcast(data, formula = trait2 ~ gene, value.var = c("PP.H3.abf"), sep = "_")
  M3<-as.matrix(Tab_PP3[,-1])
  
  x1 = dim(M4)[1]
  x2 = dim(M4)[2]
  
  M<-matrix(0,x1,x2)
  for (i in 1:x1){
    for (j in 1:x2){
      m4<-M4[i,j]
      m3<-M3[i,j]
      
      if(is.na(m3)==T){
        M[i,j] = 0
      }else if(m4>m3){
        M[i,j]<-m4
      }else{
        M[i,j]<- -m3
      }
    }
  }
  rownames(M)<-Tab_PP4$trait2
  colnames(M)<-names(Tab_PP4)[-1]
  
  return(M)
}

MA = myBuildMatrixFunction(data = Coloc_a)
MM = myBuildMatrixFunction(data = Coloc_m)
MF = myBuildMatrixFunction(data = Coloc_f)

#' ## Plots ####
#' ***
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(MA, method="circle",tl.col="black",#tl.srt=45, 
         cl.ratio = 0.7, cl.align = "l",tl.cex=0.25,
         col = brewer.pal(n = 8, name = "RdBu"))
corrplot(MM, method="circle",tl.col="black",#tl.srt=45, 
         cl.ratio = 0.7, cl.align = "l",tl.cex=0.25,
         col = brewer.pal(n = 8, name = "RdBu"))
corrplot(MF, method="circle",tl.col="black",#tl.srt=45, 
         cl.ratio = 0.7, cl.align = "l",tl.cex=0.25,
         col = brewer.pal(n = 8, name = "RdBu"))

tiff(filename = "../results/06_ColocPlot_combined.tif",
     width = 1280, height = 620, res = 200, compression = 'lzw')
corrplot(MA, method="circle",tl.col="black",#tl.srt=45, 
         cl.ratio = 0.7, cl.align = "l",tl.cex=0.25,
         col = brewer.pal(n = 8, name = "RdBu"))
dev.off()

tiff(filename = "../results/06_ColocPlot_males.tif",
     width = 1280, height = 620, res = 200, compression = 'lzw')
corrplot(MM, method="circle",tl.col="black",#tl.srt=45, 
         cl.ratio = 0.7, cl.align = "l",tl.cex=0.25,
         col = brewer.pal(n = 8, name = "RdBu"))
dev.off()

tiff(filename = "../results/06_ColocPlot_females.tif",
     width = 1280, height = 620, res = 200, compression = 'lzw')
corrplot(MF, method="circle",tl.col="black",#tl.srt=45, 
         cl.ratio = 0.7, cl.align = "l",tl.cex=0.25,
         col = brewer.pal(n = 8, name = "RdBu"))
dev.off()

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
