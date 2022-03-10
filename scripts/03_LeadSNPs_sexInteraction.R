#' ---
#' title: "Lead SNPs: sex interaction"
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
#' I want first to define the lead SNP per gene and setting and then to compare the effect estimates in males and females (are they sex-specific or sex-related?).
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../mySourceFile.R")
.libPaths()
setwd(projectpath_main)

#' # Load & Prep Olink data ####
#' ***
#' Created by script 02_Hierarchical_FDR.R
#' 
load("../temp/02_OlinkGenes.RData")
load("../results/02_pQTLs_all_cis.RData")
load("../results/02_pQTLs_males_cis.RData")
load("../results/02_pQTLs_females_cis.RData")

#' # Get lead SNPs ####
#' ***
leadpQTLs1<-pQTLs_a[hierarch_fdr5proz==T,.SD[pval==min(pval)],by=.(gene)]
leadpQTLs2<-pQTLs_m[hierarch_fdr5proz==T,.SD[pval==min(pval)],by=.(gene)]
leadpQTLs3<-pQTLs_f[hierarch_fdr5proz==T,.SD[pval==min(pval)],by=.(gene)]

table(duplicated(leadpQTLs2$gene))
table(leadpQTLs2$gene)
dummy = leadpQTLs2[gene=="TFF3",]
head(dummy)

#' Problem: Two parameters had in the other strata low MAF and were filtered (before harmonization!): 
#' 
#' * TFF3 (females) MAF 0.0087269
#' * AXL (males) MAF 0.00891873
#' 
#' My Plan: In this script, I use  rs140914768:43735844:C:T as TFF3 proxy for the lead SNPs (is in high LD, r2 = 0.908) and rs143715176:41517850:T:A as AXL proxy for the lead SNP (is in low LD, r2 = 0.0652)

dum1 = pQTLs_m[hierarch_fdr5proz==T & gene =="TFF3",]
dum2 = pQTLs_f[ gene =="TFF3",]
table(is.element(dum1$variant_id_hg19,dum2$variant_id_hg19))
dum3 = dum1[is.element(variant_id_hg19,dum2$variant_id_hg19),]
setorder(dum3,pval)
dum3[1,]
leadpQTLs2 = leadpQTLs2[gene!="TFF3",]
dummy<-pQTLs_m[variant_id_hg19 == "rs140914768:43735844:C:T",]
leadpQTLs2 = rbind(leadpQTLs2,dummy)

dum1 = pQTLs_f[hierarch_fdr5proz==T & gene =="AXL",]
dum2 = pQTLs_m[ gene =="AXL",]
table(is.element(dum1$variant_id_hg19,dum2$variant_id_hg19))
dum3 = dum1[is.element(variant_id_hg19,dum2$variant_id_hg19),]
setorder(dum3,pval)
dum3[1,]
leadpQTLs3 = leadpQTLs3[gene!="AXL",]
dummy<-pQTLs_f[variant_id_hg19 == "rs143715176:41517850:T:A",]
leadpQTLs3 = rbind(leadpQTLs3,dummy)

leadpQTLs<-rbind(leadpQTLs1,leadpQTLs2,leadpQTLs3)
leadpQTLs[,dumID:=paste(variant_id_hg38,gene,sep="__")]

save(leadpQTLs, file="../results/03_leadSNPs.RData")

#' # Look-up all lead SNPs in other settings ####
#' ***
mySNPs<-unique(leadpQTLs$variant_id_hg38)

missingSNPs<-unique(c(mySNPs[!is.element(mySNPs,pQTLs_a[,variant_id_hg38])],
                      mySNPs[!is.element(mySNPs,pQTLs_m[,variant_id_hg38])],
                      mySNPs[!is.element(mySNPs,pQTLs_f[,variant_id_hg38])]))

dummy = leadpQTLs[variant_id_hg38 %in% missingSNPs,]
length(unique(dummy$gene))
dummy

leadpQTLs_a<-pQTLs_a[variant_id_hg38 %in% mySNPs,]
leadpQTLs_a[,dumID:=paste(variant_id_hg38,gene,sep="__")]
leadpQTLs_a<-leadpQTLs_a[dumID %in% leadpQTLs[,dumID]]

leadpQTLs_m<-pQTLs_m[variant_id_hg38 %in% mySNPs,]
leadpQTLs_m[,dumID:=paste(variant_id_hg38,gene,sep="__")]
leadpQTLs_m<-leadpQTLs_m[dumID %in% leadpQTLs[,dumID]]

leadpQTLs_f<-pQTLs_f[variant_id_hg38 %in% mySNPs,]
leadpQTLs_f[,dumID:=paste(variant_id_hg38,gene,sep="__")]
leadpQTLs_f<-leadpQTLs_f[dumID %in% leadpQTLs[,dumID]]

table(leadpQTLs_a$variant_id_hg38==leadpQTLs_f$variant_id_hg38)
table(leadpQTLs_a$variant_id_hg38==leadpQTLs_m$variant_id_hg38)

#' # Combine settings ####
#' ***
names(leadpQTLs_a)
myFDR4<-cbind(leadpQTLs_a[,c(22,17,16,2,1,10,5,7,11:13,15,9,20)],
              leadpQTLs_m[,c(5,7,11:13,15,9,20)],
              leadpQTLs_f[,c(5,7,11:13,15,9,20)])
names(myFDR4)
names(myFDR4)<-c("dumID","gene","cyto","variant_id_hg38","variant_id_hg19","info",
                 "ea_a","eaf_a","beta_a","se_a","pval_a","zscore_a","n_a","hfdr_a",
                 "ea_m","eaf_m","beta_m","se_m","pval_m","zscore_m","n_m","hfdr_m",
                 "ea_f","eaf_f","beta_f","se_f","pval_f","zscore_f","n_f","hfdr_f")
myFDR4

plot(myFDR4$eaf_a,myFDR4$eaf_m)
plot(myFDR4$eaf_a,myFDR4$eaf_f)

table(myFDR4$ea_a==myFDR4$ea_m)
table(myFDR4$ea_a==myFDR4$ea_f)

#' # Interaction test ####
#' ***
interactionTest_classic  = function(mean1, se1, mean2, se2) { 
  ## interaction test see http://www.bmj.com/content/326/7382/219.long
  meandiff_se = sqrt(se1^2 + se2^2)
  meandiff = mean2 - mean1
  meandiff_cilow = meandiff - 1.96*meandiff_se
  meandiff_cihigh = meandiff + 1.96*meandiff_se
  meandiff_z = meandiff/meandiff_se
  meandiff_p = pnorm(abs(meandiff_z), lower.tail = F)*2
  if(meandiff_p>1) meandiff_p = 1
  data.table::data.table(mean1, se1, mean2, se2, meandiff, meandiff_se, meandiff_cilow, meandiff_cihigh, meandiff_z, meandiff_p)
}

sexIA<-interactionTest_classic(myFDR4$beta_m,myFDR4$se_m,myFDR4$beta_f,myFDR4$se_f)
table(sexIA$meandiff_p<0.05)
sexIA[,meandiff_p_FDR:=p.adjust(meandiff_p,method = "BH")]
sexIA[,table(meandiff_p<0.05,meandiff_p_FDR<0.05)]

myFDR5<-cbind(myFDR4,sexIA[,c(5,6,10,11)])
myFDR5[meandiff_p_FDR<0.05]
myFDR5[meandiff_p<0.05]
myFDR5[meandiff_p>=0.05,sexIA_type:="sex-unspecific"]
myFDR5[meandiff_p<0.05 & hfdr_m==T & hfdr_f==F,sexIA_type:="male-specific"]
myFDR5[meandiff_p<0.05 & hfdr_m==F & hfdr_f==T,sexIA_type:="female-specific"]
myFDR5[meandiff_p<0.05 & hfdr_m==T & hfdr_f==T,sexIA_type:="sex-related"]

#' # Plotting ####
#' ***
#' I want a beta-beta plot. None-significant interaction should be black, nominal-significant interaction should be colored, and FDR-significant interaction should be labeled.
#' 
 
myPlotData<-data.table(dumID=myFDR5$dumID,
                       rs_id=myFDR5$variant_id_hg38,
                       cyto=myFDR5$cyto,
                       gene=myFDR5$gene,
                       eaf_male = myFDR5$eaf_m,
                       beta_male = myFDR5$beta_m,
                       eaf_female = myFDR5$eaf_f,
                       beta_female = myFDR5$beta_f,
                       meandiff = myFDR5$meandiff,
                       meandiff_p = myFDR5$meandiff_p,
                       meandiff_p_FDR = myFDR5$meandiff_p_FDR,
                       sexIA_type = myFDR5$sexIA_type)
myPlotData[meandiff_p>=0.05,sig:="no"]
myPlotData[meandiff_p<0.05,sig:="yes"]
myPlotData[meandiff_p_FDR<0.05,sig:="yes (FDR 5%)"]
myPlotData[,gene2:=""]
myPlotData[meandiff_p_FDR<0.05,gene2:=gene]
myPlotData[,gene3:=""]
myPlotData[meandiff_p<0.05,gene3:=gene]

myPlot = ggplot(myPlotData[sig!="no",], aes(x=beta_male, y=beta_female, color=gene3,label=gene2,shape = sexIA_type)) +
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1.25)+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_point(data=myPlotData[sig=="no",], aes(x=beta_male, y=beta_female),col="black",size=1.5,alpha=0.5,shape = 16)+
  geom_point(size=3.5)+ 
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"))+
  labs(x="Effect in males", y = "Effect in females",
       color = "Genes with \nsex interaction", shape = "Type of \nsex interaction")+
  guides(size="none",fill="none")+
  scale_shape_manual(values=c(15, 17, 18,16))+
  geom_label_repel(data = subset(myPlotData, sig!="no"),# & type=="beta estimate"),
                   aes(fill=gene3),
                   box.padding = 1.15,
                   max.overlaps = 25,
                   show.legend=F,color="black")
myPlot

tiff(filename = "../results/03_ScatterPlot_sexIA.tif", 
     width = 1200, height = 900, res=125, compression = 'lzw')
myPlot
dev.off()

#' # Save ####
#' ***
leadSNPs_IA<-copy(myFDR5)
save(leadSNPs_IA,file="../results/03_leadSNPs_IA.RData")

sigLoci_IA_sig<-copy(myFDR5)
sigLoci_IA_sig = sigLoci_IA_sig[meandiff_p<0.05,]

myTab[,sexIA:=F]
myTab[Gene %in% sigLoci_IA_sig$gene,sexIA:=T]
head(myTab)
matched = match(myTab$Gene,sigLoci_IA_sig$gene)
table(is.na(matched))
length(unique(sigLoci_IA_sig$gene))
myTab[,sexIA_type:=sigLoci_IA_sig$sexIA_type[matched]]

save(myTab,file="../temp/03_OlinkGenes.RData")

#' # sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
