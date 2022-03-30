#' ---
#' title: "Harmonize Olink data to hg38"
#' subtitle: "LWAS Olink - Heart"
#' author: "Janne Pott"
#' date: "29. March 2022"
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
#' I cannot do this with the same script as in Adult, since the github tool cannot lift chr X. Hence, I use the [UCSC LiftOver unix command line tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver). 
#' 
#' As SNP identifier I want to use chr_pos_effectAllele_otherAllele (this is the GTEx coding - I use the GTEx annotation file to check for same effect allele)
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../mySourceFile.R")
.libPaths()
setwd(projectpath_main_Heart)

#' # Load data ####
#' ***
myTab<-fread("../temp_Heart/OlinkGenes_Heart.txt")
head(myTab)

pQTLs_a = fread("../exampleData_Heart/pQTLs_combined_hg19.txt.gz")
pQTLs_f = fread("../exampleData_Heart/pQTLs_females_hg19.txt.gz")
pQTLs_m = fread("../exampleData_Heart/pQTLs_males_hg19.txt.gz")
pQTLs = rbind(pQTLs_a,pQTLs_f,pQTLs_m)
head(pQTLs)

gtex8annot= fread(path_GTExv8_annotation)
head(gtex8annot)
table(gtex8annot$chr)

chain = import.chain(path_Liftover_Chainfile)
names(chain)

#' # Lifting data ####
#' ***
#' 
#' ## Step 1: Reduce to unique chromosomal position ####
#' 
#' I only need all unique chr_pos. 
#' 
myLiftingTab = copy(pQTLs)
myLiftingTab[,chrPos:= paste(chr,pos_hg19,sep="_")]
myLiftingTab = myLiftingTab[!duplicated(chrPos),]
myLiftingTab = myLiftingTab[,c(1:7)]
head(myLiftingTab)

#' ## Step 2: Check column names ####
#' 
#' I need the columns *snps*, *chr* and *pos*
#' 
setnames(myLiftingTab,"rs_id","snps")
setnames(myLiftingTab,"pos_hg19","pos")
stopifnot(c('snps', 'chr', 'pos') %in% names(myLiftingTab))

#' ## Step 3: Check for NAs ####
#' 
#' The chr should be still numeric, but I need it as character starting with "chr". In addition, chromosome X should be given as X, not 23. I also check that all chr and pos entries are not NA. 
#' 
check1 = myLiftingTab[,stopifnot(all(na.omit(chr) %in% 1:23))]
check1
myLiftingTab[,chr := as.character(chr)]
myLiftingTab[chr =="23", chr := "X"]
myLiftingTab[,chr := paste0("chr",chr)]

check2 = is.na(myLiftingTab$chr) | is.na(myLiftingTab$pos)
stopifnot(sum(check2)==0)

#' ## Step 4: Change to genomic range ####
#' 
#' I use the *GRanges* function from the GenomicRanges package
#' 
tolift_gr = GRanges(seqnames = Rle(myLiftingTab$chr), 
                    ranges = IRanges(start = myLiftingTab$pos, end = myLiftingTab$pos),
                    strand = Rle("*", length(seqnames)))
tolift_gr

values(tolift_gr) <- DataFrame(snps = myLiftingTab$snps)
tolift_gr    

#' ## Step 5: Do lifting ####
#' 
#' I use the chainfile from [UCSC LiftOver tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver) for hg19 to hg38. 
#' 
message("Using for recoding chainfile: ",path_Liftover_Chainfile)
lifted = liftOver(tolift_gr, chain)
class(lifted)
lifted = unlist(lifted)
lifted
genome(lifted) = "hg38"
names(lifted@elementMetadata) = "snps"

#' ## Step 6: Merge with original data ####
#' 
#' I merge the lifted and original data set. Then I do some checks (same chr, how many successfully lifted?)
lifted_dt = data.table(snps = lifted$snps, 
                       chr_hg38 = as.character(seqnames(lifted)), 
                       pos_hg38 = ranges(lifted)@start)

setkey(lifted_dt, snps)
setkey(myLiftingTab, snps)

res = lifted_dt[myLiftingTab]
res[,1:5]

names(res)
setnames(res,"chr","chr_hg19")
setnames(res,"pos","pos_hg19")

table(res$chr_hg19 == res$chr_hg38)

res = res[ ,chr_hg38:= str_replace(chr_hg38, "chr", "")]
res[,chr:=chr_hg38]
res = res[chr_hg38== "X", chr_hg38 := '23']
res = res[ ,chr_hg19:= str_replace(chr_hg19, "chr", "")]
res = res[chr_hg19== "X", chr_hg19 := '23']
res[,liftmethod := "liftover successful"]
res[is.na(pos_hg38),liftmethod := "liftover failed"]
table(res$liftmethod)

message("Succesfully lifted: " ,dim(res[liftmethod == "liftover successful",])[1], " entries...\n")
message("Lifting failed: " ,dim(res[liftmethod == "liftover failed",])[1], " entries...\n")

#' # Adding hg38 info ####
#' ***
#' Now I can rematch the snpwise hg38 position information to my pQTLs data set. 
#' 
pQTLs[,chrPos_hg19 := paste(chr,pos_hg19,sep="_")]
res[,chrPos_hg19 := paste(chr_hg19,pos_hg19,sep="_")]
matched1 = match(pQTLs$chrPos_hg19,res$chrPos_hg19)
table(is.na(matched1))
pQTLs[,pos_hg38 := res[matched1,pos_hg38]]

pQTLs_filtered = copy(pQTLs)
pQTLs_filtered = pQTLs_filtered[is.na(pos_hg38)]
pQTLs = pQTLs[!is.na(pos_hg38)]
pQTLs_filtered[,table(setting,cyto)]
pQTLs[is.element(cyto,pQTLs_filtered$cyto),table(setting,cyto)]

#' Yes, I lose a few variants, but not one locus completely. 
#' 
#' # Comparing to GTEx v8 annotation ####
#' ***
#' I use the GTEx v8 annotation table to get the same effect allele (for those SNPs with matching positions). 
#' 
res[,chrPos_hg38 := paste0("chr",chr,"_",pos_hg38)]
gtex8annot[,chrPos_hg38 := paste0(chr,"_",variant_pos)]
filt = is.element(gtex8annot$chrPos_hg38,res$chrPos_hg38)
table(filt)
gtex8annot = gtex8annot[filt,]
head(gtex8annot)
table(gtex8annot$num_alt_per_site)

#' To make my life easier, I filter all sites with more than one alternative allele, focusing on biallelic sites.  
#' 
gtex8annot = gtex8annot[num_alt_per_site == 1,]

matched2 = match(res$chrPos_hg38,gtex8annot$chrPos_hg38)
table(is.na(matched2))
res[,variant_id_hg38:= gtex8annot[matched2,variant_id]]
res[,variant_id_hg37:= gtex8annot[matched2,variant_id_b37]]
res[,rs_id:= gtex8annot[matched2,rs_id_dbSNP151_GRCh38p7]]
res[,GTEx_effectAllele:= gtex8annot[matched2,ref]]
res[,GTEx_otherAllele:= gtex8annot[matched2,alt]]
head(res)
table(res$liftmethod,is.na(res$variant_id_hg38))
table(res$chr_hg38==23,is.na(res$variant_id_hg38))

matched3 = match(pQTLs$chrPos_hg19,res$chrPos_hg19)
table(is.na(matched3))
table(pQTLs$chrPos_hg19 == res[matched3,chrPos_hg19])
pQTLs[,variant_id_hg38 := res[matched3,variant_id_hg38]]
pQTLs[,variant_id_hg37 := res[matched3,variant_id_hg37]]
pQTLs[,rs_id2 := res[matched3,rs_id]]
pQTLs[,GTEx_effectAllele:= res[matched3,GTEx_effectAllele]]
pQTLs[,GTEx_otherAllele:= res[matched3,GTEx_otherAllele]]

#' Check coding - is my effect allele the same as the GTEx effect allele?
pQTLs[!is.na(GTEx_effectAllele),table(effect_allele == GTEx_effectAllele)]
pQTLs[!is.na(GTEx_otherAllele),table(other_allele == GTEx_otherAllele)]
pQTLs[!is.na(GTEx_effectAllele),table(other_allele == GTEx_effectAllele)]
pQTLs[!is.na(GTEx_otherAllele),table(effect_allele == GTEx_otherAllele)]

pQTLs[,dumID1:=paste0("chr",chr,"_",pos_hg38,"_",effect_allele,"_",other_allele,"_b38")]
pQTLs[,dumID2:=paste0("chr",chr,"_",pos_hg38,"_",other_allele,"_",effect_allele,"_b38")]
pQTLs[,dumID1:=gsub("chr23_","chrX_",dumID1)]
pQTLs[,dumID2:=gsub("chr23_","chrX_",dumID2)]

pQTLs[!is.na(GTEx_effectAllele),table(variant_id_hg38 == dumID1,variant_id_hg38==dumID2)]
pQTLs_filtered2 = copy(pQTLs)
pQTLs_filtered2 = pQTLs_filtered2[variant_id_hg38!=dumID1 & variant_id_hg38!=dumID2]
pQTLs_filtered2[,table(setting,cyto)]
pQTLs[(variant_id_hg38==dumID1 | variant_id_hg38==dumID2) & is.element(cyto,pQTLs_filtered2$cyto),table(setting,cyto)]

filt = pQTLs$dumID2 == pQTLs$variant_id_hg38
table(filt)
pQTLs[filt, effect_allele := GTEx_effectAllele]
pQTLs[filt, other_allele := GTEx_otherAllele]
pQTLs[filt, beta := beta * (-1)]
pQTLs[filt, eaf := 1-eaf]
pQTLs[,dumID3:=paste0("chr",chr,"_",pos_hg38,"_",effect_allele,"_",other_allele,"_b38")]
pQTLs[,dumID3:=gsub("chr23_","chrX_",dumID3)]
table(pQTLs$dumID3 == pQTLs$variant_id_hg38)

#' # Filter for 500 kb Window ####
#' ***
dim(pQTLs)

dumTab = foreach(i = 1:dim(myTab)[1])%do%{
  #i=1
  myRow = myTab[i,]
  message("Working on ",myRow$Parameter, ", number ",i," of 92")
  # Get genetic range
  myDat = copy(pQTLs)
  myDat = myDat[gene == myRow$Gene]
  myDat = myDat[pos_hg38 <= myRow$stop_hg38 + 500000]
  myDat = myDat[pos_hg38 >= myRow$start_hg38 -500000]
  
  # Add FDR & maf & zscore
  myDat[,FDR:=p.adjust(p = pval,method="fdr"),by=setting]
  myDat[,maf:=eaf]
  myDat[eaf>0.5,maf:=1-eaf]
  myDat[,zscore:=beta/se]
  myDat
  
}
pQTLs2 = rbindlist(dumTab)
dim(pQTLs2)

#' # Get necessary columns and save ####
#' ***
#' I want the same columns as in the LIFE Adult analysis:
#' 
#' * variant_id_hg19: ID as used in LIFE 
#' * variant_id_hg38: ID as used by GTEx v8
#' * rs_id: rsID according to dbSNP151 (GRCh38p7)
#' * chr: chromosome as character using chr7 format
#' * pos_hg19: base position accoring to hg19
#' * pos_hg38: base position accoring to hg38
#' * effect_allele: used effect allele (aka risk allele, reference allele, or coded allele)  
#' * other_allele: used other allele (aka alternate allele, or non-effect allele)
#' * eaf: effect allele frequency
#' * maf: minor allele frequency
#' * n_samples: sample size
#' * info: imputation info score
#' * beta: effect estimate
#' * se: standard error of effect estimate
#' * pval: p-value of SNP effect
#' * FDR: adjusted p-value of SNP effect, adjusting for each protein
#' * zscore: Z-scores of SNP effect (beta/se)
#' * cyto: cytoband of considered gene
#' * gene: gene name of considered protein
#' * pheno: protein name as given by Olink
#' * setting: setting of analyses, either combined, males or females
#' 
names(pQTLs2)
myNames1<-c("rs_id","dumID3","rs_id2","chr","pos_hg19","pos_hg38",
            "effect_allele","other_allele","eaf","maf","n_samples","info",
            "beta","se","pval","FDR","zscore",
            "cyto","gene","pheno","setting")
colsOut<-setdiff(colnames(pQTLs2),myNames1)
pQTLs2[,get("colsOut"):=NULL]
setcolorder(pQTLs2,myNames1)
setnames(pQTLs2,"dumID3","variant_id_hg38")
setnames(pQTLs2,"rs_id","variant_id_hg19")
setnames(pQTLs2,"rs_id2","rs_id")
head(pQTLs2)

pQTLs_a<-pQTLs2[setting=="combined",]
pQTLs_m<-pQTLs2[setting=="males",]
pQTLs_f<-pQTLs2[setting=="females",]

fwrite(pQTLs_a, file = "../exampleData_Heart/pQTLs_combined_hg38.txt",quote = F,sep="\t")
R.utils::gzip("../exampleData_Heart/pQTLs_combined_hg38.txt", 
              destname = "../exampleData_Heart/pQTLs_combined_hg38.txt.gz")
fwrite(pQTLs_f, file = "../exampleData_Heart/pQTLs_females_hg38.txt",quote = F,sep="\t")
R.utils::gzip("../exampleData_Heart/pQTLs_females_hg38.txt", 
              destname = "../exampleData_Heart/pQTLs_females_hg38.txt.gz")
fwrite(pQTLs_m, file = "../exampleData_Heart/pQTLs_males_hg38.txt",quote = F,sep="\t")
R.utils::gzip("../exampleData_Heart/pQTLs_males_hg38.txt", 
              destname = "../exampleData_Heart/pQTLs_males_hg38.txt.gz")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
