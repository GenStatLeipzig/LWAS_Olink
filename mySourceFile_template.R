#############################
# this is a template source file
# please change all paths accordingly
#############################

#############################
# Working directory
#############################
projectpath_main = "/PATH/TO/PROJECTS/LWAS_Olink/scripts/"

#############################
# R library and R packages
#############################
.libPaths("/PATH/TO/RLibrary/VERSION_4.x/") 
.libPaths()

suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
library(readxl)
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doMC))
library(reticulate)
library(ggplot2)
library(ggrepel)
library(plyr)
suppressPackageStartupMessages(library(toolboxH))
suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))   
library(gridBase)
suppressPackageStartupMessages(library(corrplot))
library(RColorBrewer)

#############################
# Other Tools 
#############################
path_phyton = "/PATH/TO/TOOLS/anaconda2/bin/python"
path_conda = "/PATH/TO/TOOLS/anaconda2/bin/conda"
path_metaxcan_gwastools = "/PATH/TO/TOOLS/MetaXcan/summary-gwas-imputation-master/src"
path_metaxcan_data = "/PATH/TO/TOOLS/MetaXcan/data"
path_metaxcan_tools = "/PATH/TO/TOOLS/MetaXcan/MetaXcan-master/software/"

path_plink2 = "/PATH/TO/TOOLS/plink2.0/20210203/unix_64/plink2"

#############################
# Downloaded data sets 
#############################
path_GTExv8 = "/PATH/TO/DATA/GTEx_v8/"
path_1000Genomes = "/PATH/TO/DATA/1000genomes_phase3_vs5/plinkformat/EUR"
path_CAD = "/PATH/TO/DATA/SumStatsCAD/CAD_META.gz"
