#' ---
#' title: "Circular Plot of cis-associations"
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
#' I want to create a circular plot for the cis effects. I do this only for the combined setting.
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../mySourceFile.R")
.libPaths()
setwd(projectpath_main)

#' # Load ####
#' ***
load("../results/02_pQTLs_all_cis.RData")
head(pQTLs_a)
load("../results/04_a_eQTL_lookup_allTissues.RData")
head(eQTL_GTEx)
load("../temp/04_OlinkGenes.RData")
head(myTab)

#' Load toplist and GWAS catalog annotation (created by Genstat intern pipeline)
GWAScat<-fread("../exampleData/gwasinfo_2022-02-02_Olink_cis_neu.txt")
head(GWAScat)
TopList = fread("../exampleData/topliste_2022-02-02_Olink_cis_neu.txt")
head(TopList)

#' # Filt for relevant genes ####
#' ***
leadpQTLs<-pQTLs_a[,.SD[pval==min(pval)],by=.(gene)]
leadpQTLs<-leadpQTLs[gene %in% myTab[GWAS_sig==T,Gene]]

data_besteQTL<-eQTL_GTEx[comment=="best_eQTL"]
data_besteQTL_top<-data_besteQTL[,.SD[pval==min(pval)],by=.(gene)]

matched1<-match(myTab$Gene,data_besteQTL_top$gene)
myTab[,top_eQTL_tissue:=data_besteQTL_top$tissue[matched1]]

matched2<-match(myTab$Gene,leadpQTLs$gene)
myTab[,leadSNP:=leadpQTLs$variant_id_hg38[matched2]]

myTab<-myTab[GWAS_sig==T]
head(myTab)
table(is.na(myTab$top_eQTL_tissue))

#' # Search for novel findings ####
#' ***
TopList = TopList[setting == "combined",]
TopList[,dumID:=paste(litsnp,gene,sep=":")]
TopList = TopList[!duplicated(dumID),]
TopList2 = TopList[!is.na(corinfo_gwas2)]
TopList3 = TopList[is.na(corinfo_gwas2)]

genes2 = unique(TopList2$gene)
genes3 = unique(TopList3$gene)
genes_novel = genes3[!is.element(genes3,genes2)]
genes_novel
#' So these 12 genes are not in the GWAS catalog at all. 
#' In addition, I accept as novel also genes that have a GWAS catalog entry with anything but biomarker
#'  
GWAScat = GWAScat[snps %in% TopList2$litsnp]
length(unique(GWAScat$phenotype_gwas))

exclude = c("Blood metabolite levels","Blood protein levels","Blood protein levels in cardiovascular risk",
            "Chronic obstructive pulmonary disease-related biomarkers","Matrix metalloproteinase levels",
            "Matrix metalloproteinase-10 levels", "NT-proBNP levels in acute coronary syndrome",
            "Plasma proprotein convertase subtilisin/kexin type 9 levels in stable coronary artery disease",
            "Protein biomarker", "Protein quantitative trait loci","Soluble interleukin-2 receptor subunit alpha")
GWAScat2 = GWAScat[phenotype_gwas %nin% exclude]
GWAScat3 = GWAScat[phenotype_gwas %in% exclude]

genes_known = unique(TopList2[litsnp %in% GWAScat3$snps,gene])

table(is.element(genes_known,genes_novel))
table(is.element(myTab$Gene,genes_known),is.element(myTab$Gene,genes_novel))

myTab[,novel:="maybe"]
myTab[Gene %in% genes_known,novel:="no"]
myTab[Gene %in% genes_novel,novel:="yes"]
table(myTab$novel)

#' ## Parameter setting ####
#' 
pQTLs<-pQTLs_a[hierarch_fdr5proz==T,]
min(pQTLs$maf)
min(pQTLs$info)
minMaf = 0.01
minInfo = 0.8
max(pQTLs$pval)
max(pQTLs$FDR)
minP = 0.05 * (64/92)
logPCut = -log10(5e-8)
minLogP = -log10(minP)

#' ## pQTLs ####
#' 
setnames(pQTLs,"variant_id_hg38","markername")
table(pQTLs$chr)

#' ## eQTLs ####
#' Here, I can use the data for coloc! 
#' 1) Load tissue
#' 2) Get all SNPs for relevant genes
#' 3) Filter for maf
#' 4) save
#' 
topTissues<-unique(data_besteQTL_top[,tissue])
length(topTissues)

datalist_trait_fn2<-dir("../temp/GTEx_v8_filtered/")
traitlist_trait2<-gsub(".RData","",datalist_trait_fn2)

matched<-match(topTissues,traitlist_trait2)
table(is.na(matched))
datalist_trait_fn2<-datalist_trait_fn2[matched]
traitlist_trait2<-traitlist_trait2[matched]
table(topTissues==traitlist_trait2)
traitlist_trait2<-gsub("_"," ",traitlist_trait2)

dumTab<-foreach(i=1:length(topTissues))%do%{
  #dumTab<-foreach(i=1:5)%do%{
  #i=1
  myTissue<-topTissues[i]
  myPhenos<-myTab[top_eQTL_tissue==myTissue,Parameter]
  myGenes<-myTab[top_eQTL_tissue==myTissue,Gene]
  myCytos<-myTab[top_eQTL_tissue==myTissue,Cytoband]
  
  message("Working on ",myTissue, " (",i," of ",length(topTissues),")")
  
  # Load eQTL data
  loaded1<-load(paste0("../temp/GTEx_v8_filtered/",datalist_trait_fn2[i]))
  data3<- get(loaded1)
  
  # Check & filter overlap
  data3<-data3[gene %in% myGenes,]
  table(data3$gene)
  data3[,tissue := myTissue]
  data3[,type:="eQTL"]
  
  for(i in 1:length(myGenes)){
    data3[gene==myGenes[i],cyto:=myCytos[i]]
    data3[gene==myGenes[i],phenotype:=myPhenos[i]]
  }
  
  # Filter pvalues
  data4<-copy(data3)
  data4<-data4[pval<=minP,]
  
  # Reduce to relevant columns
  myNames<-c("variant_id","chr","pos","maf","beta","se","pval","cyto","gene","phenotype","type","tissue")
  colsOut<-setdiff(colnames(data4),myNames)
  data4[,get("colsOut"):=NULL]
  setcolorder(data4,myNames)
  myNames2<-c("markername","chr","pos","eaf","beta","se","pval","cyto","gene","pheno","type","tissue")
  names(data4)<-myNames2
  head(data4)
  data4
  
}
eQTLs<-rbindlist(dumTab)
eQTLs
table(eQTLs$chr)
eQTLs[,chr:= paste0("chr",chr)]
length(unique(eQTLs$gene))
dumTab2<-eQTLs[,.SD[pval==min(pval)],by=.(pheno)]

#' # Best logP per SNP ####
#' ***
pQTLs[,myID:=paste(chr,pos_hg38,sep=":")]
eQTLs[,myID:=paste(chr,pos,sep=":")]

table(is.element(pQTLs$markername,eQTLs$markername))
table(is.element(eQTLs$markername,pQTLs$markername))
table(is.element(pQTLs$myID,eQTLs$myID))
table(is.element(eQTLs$myID,pQTLs$myID))

pQTLs[,logp:=-log10(pval)]
eQTLs[,logp:=-log10(pval)]

#' matching
plotData = data.table(myID = unique(c(pQTLs[,myID], eQTLs[,myID])))
matched1 = match(plotData[,myID], pQTLs[,myID])
matched2 = match(plotData[,myID], eQTLs[,myID])
matched3 = match(myTab[,leadSNP],leadpQTLs$variant_id_hg38)
myTab[,myID:=leadpQTLs[matched3,paste(chr,pos_hg38,sep=":")]]

#' basic plot data object
plotData[,chr:=pmax(pQTLs[matched1, chr], eQTLs[matched2, chr], na.rm=T)]
plotData[,start:=pmax(pQTLs[matched1, pos_hg38], eQTLs[matched2, pos], na.rm=T)]
plotData[,end:=start+1]
setcolorder(plotData, c("chr","start","end","myID"))
setkeyv(plotData, cols=c("chr","start"))
plotData
plotData[,chr2:=as.numeric(gsub("chr","",chr))]
setkeyv(plotData, cols=c("chr2","start"))
plotData

table(is.na(plotData[,myID]))
table(is.na(plotData[,chr]))
table(is.na(plotData[,start]))

#' get max log p of pQTLs and eQTLs object
p.max = pQTLs[,.SD[pval==min(pval)],by=.(myID)]
e.max = eQTLs[,.SD[pval==min(pval)],by=.(myID)]

matched1 = match(plotData[,myID], p.max$myID)
matched2 = match(plotData[,myID], e.max$myID)
matched3 = match(plotData[,myID], myTab$myID)
table(is.na(matched3))

plotData[,pQTL:=0]
plotData[,pQTL:=p.max$logp[matched1]]
plotData[,eQTL:=0]
plotData[,eQTL:=e.max$logp[matched2]]
plotData[,Protein:=""]
plotData[,Protein:=myTab$Parameter[matched3]]
plotData[,Gene:=""]
plotData[,Gene:=myTab$Gene[matched3]]
plotData[,Cyto:=""]
plotData[,Cyto:=myTab$Cytoband[matched3]]
plotData[,novel:=""]
plotData[,novel:=myTab$novel[matched3]]

table(is.na(plotData[,pQTL]), is.na(plotData[,eQTL]))
table(is.na(plotData[,pQTL]), is.na(plotData[,novel]))

plotData[,maxP:=pmax(pQTL,eQTL, na.rm=T)]

save(plotData, file="../results/05_plotData_forCircPlot.RData")
save(myTab,file="../temp/05_OlinkGenes.RData")

#' # Load & Prepare ####
#' ***
dim(plotData)  

#' rename chromosomes
myChrs<-unique(plotData$chr)

#' seperate gene legend data
legendData = plotData[novel == "yes" | novel =="maybe",]
dim(legendData)
table(is.na(legendData$pQTL))

#' limit y-axis to 15
maxi = ceiling(max(c(plotData[,pQTL],plotData[,eQTL]), na.rm=T))
maxi
xCut = 20
if (maxi>xCut) {
  maxi=xCut+1
  plotData[eQTL>=xCut,eQTL:=xCut]
  plotData[pQTL>=xCut,pQTL:=xCut]
}
gws = -log10(5*10^(-8))

#' get median per chromosome for later chr plotting on inner line
getMiddle = function(x) {
  x = x[!is.na(x)]
  res = round((min(x) + max(x)) / 2)
  return(res)
}
chrData = ddply(plotData, "chr", summarise, myMedian = getMiddle(start))
chrData = as.data.table(chrData)

#' Optional: rearrange quots data for MiamiPlot feeling --> not right now!
#plotData[,eQTL := xCut - eQTL]

#' # Define Function ####
#' ***
myCircosPlot <- function() {
  #' initialize plot
  circos.clear()
  circos.par("start.degree" = 90)
  circos.par("gap.after"=c(rep(2,19),15))
  circos.initializeWithIdeogram(chromosome.index = myChrs,
                                plotType = NULL)
  
  #' Track 1 & 2: novel genes
  myCol = legendData[,novel]
  myCol = gsub("yes","darkblue",myCol)
  myCol = gsub("maybe","black",myCol)
  myCex = rep(x = 0.7,63)
  myCex[myCol=="black"]<-0.6
  
  circos.genomicLabels(legendData, labels.column = "Gene", 
                       side = "outside", col=myCol, 
                       line_col = myCol, padding=0.5, cex=myCex, 
                       line_lwd = 0.5, 
                       connection_height = convert_height(1, "mm"))
  
  #' Track 3: pQTL results in manhattan style
  circos.genomicTrack(data=plotData, 
                      ylim=c(-log10(minP),xCut), 
                      numeric.column = c("pQTL"), 
                      bg.border="grey", 
                      cell.padding = c(0, 1.00, 0.02, 1.00),
                      panel.fun = function(region, value, ...) {
                        circos.genomicPoints(region, value, pch=16, cex=0.5, col="darkgreen",...)
                        xlimit = get.cell.meta.data("xlim")
                        circos.lines(xlimit,c(gws,gws),col="red",lwd=1)
                      },
                      bg.col="darkseagreen1", 
                      track.margin = c(0,0.01))
  
  circos.yaxis(side="left",
               at=c(4,7.3,15),
               labels=c(4,7.3,15), 
               sector.index="chr1", 
               track.index=get.current.track.index(), 
               labels.cex = 0.6, 
               labels.niceFacing = F) 
  
  #' Track 4: eQTL results in manhattan style
  circos.genomicTrack(data=plotData, 
                      ylim=c(-log10(minP),xCut), 
                      numeric.column = c("eQTL"), 
                      bg.border="grey", 
                      cell.padding = c(0, 1.00, 0.02, 1.00),
                      panel.fun = function(region, value, ...) {
                        circos.genomicPoints(region, value, pch=16, cex=0.5, col="blue",...)
                        xlimit = get.cell.meta.data("xlim")
                        circos.lines(xlimit,c(gws,gws),col="red",lwd=1)
                      },
                      bg.col="lightsteelblue1", 
                      track.margin = c(0,0.01))
  
  circos.yaxis(side="left",
               at=c(4,7.3,15),
               labels=c(4,7.3,15), 
               sector.index="chr1", 
               track.index=get.current.track.index(), 
               labels.cex = 0.6, 
               col="black")
  
  #' Track 5: chromosome annotation
  circos.track(ylim = c(0, 3), track.height = uh(3, "mm"), 
               bg.border=NA, cell.padding	= c(0.02, 1.00, 0, 1.00), 
               track.margin = c(0,0.01), 
               panel.fun = function(x, y) {
                 chr = CELL_META$sector.index
                 xlim = CELL_META$xlim
                 ylim = CELL_META$ylim
                 xmin = get.cell.meta.data("xlim")[1]
                 xmax = get.cell.meta.data("xlim")[2]
                 circos.lines(c(xmin, xmax), 
                              c(2.5, 2.5), 
                              col = "black", 
                              lwd = 2, 
                              lty = 1)
                 circos.text(x=mean(xlim), 
                             y=0, 
                             labels=gsub("chr","",chr), 
                             cex = 0.5, col = "black",
                             facing = "inside", niceFacing = TRUE)
               })
  
  #' Track 6 & 7: cytoband of novel genes
  # circos.genomicLabels(legendData, labels.column = "Cyto",
  #                      side = "inside", col=myCol,
  #                      line_col = myCol, padding=1, cex=0.5,
  #                      line_lwd = 1.5,
  #                      connection_height = convert_height(2, "mm"))
  
  
  #finish plot
  circos.clear()
  
} 

#get some infos
circos.info()

#' # Create Legend ####
#' ***
#' I want a legend for my circos plot... according to Katrin, I have to use viewPorts for this. 

point_eQTL = Legend(at = expression(-log[10] (p[eQTL]) (GTEx)), 
                    type = "points", 
                    legend_gp = gpar(col = "blue"), 
                    labels_gp = gpar(fontsize = 9), 
                    title_position = "topleft", 
                    title = NULL, 
                    background = "lightsteelblue1", 
                    grid_height = unit(4, "mm"), 
                    grid_width = unit(4, "mm"))

point_pQTL = Legend(at = expression(-log[10] (p[pQTL])), 
                    type = "points", 
                    legend_gp = gpar(col = "darkgreen"), 
                    labels_gp = gpar(fontsize = 9), 
                    title = NULL, 
                    background = "darkseagreen1", 
                    grid_height = unit(4, "mm"), 
                    grid_width = unit(4, "mm"))

line_blue = Legend(at = c("novel"), 
                   type = "lines", 
                   legend_gp = gpar(col = "darkblue", lwd = 2), 
                   labels_gp = gpar(col = "darkblue", fontsize = 9), 
                   title = "",
                   background = NULL, 
                   grid_height = unit(4, "mm"), 
                   grid_width = unit(4, "mm"))

line_black = Legend(at = c("reported for other traits \n(not blood protein levels)"), 
                    type = "lines", 
                    legend_gp = gpar(col = "black", lwd = 2), 
                    labels_gp = gpar(col = "black", fontsize = 9), 
                    title = "", 
                    background = NULL, 
                    grid_height = unit(8, "mm"), 
                    grid_width = unit(4, "mm"))

lgd_list_vertical = packLegend(line_blue, line_black, 
                                point_pQTL, point_eQTL, 
                                gap = unit(0.2, "mm") )
#' # Test Plotting ####
#' ***
myCircosPlot()
grid.draw(lgd_list_vertical)

#' # Create Plot ####
#' ***
#' Legend placed inside the circles & only yes & maybe genes displayed

# Step 1: open picture
tiff(filename = "../results/05_CircosPlot.tif", width = 1500, height = 1500, res = 200, compression = 'lzw')

# Step 2: plot circos plot
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, 
                      width = 1, height = 1, 
                      just = "center"))
myCircosPlot()
grid.draw(lgd_list_vertical)
upViewport()

# Step 3: close picture
dev.off()


#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
