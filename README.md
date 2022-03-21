# LWAS Olink

Scripts used in my Olink project (submitted to PLOS ONE - currently under revision)

Here, I want to include all relevant R scripts used in the LWAS Olink project. 

Complete data sets including genetic data of LIFE-Adult participants cannot be made publicly available due to ethical and legal restrictions, as they are sufficient to identity study participants. This is not covered by the informed consent. However, access to the LIFE-Adult data is possible via project agreements addressed to:

* LIFE Research Center for Civilization Diseases, Medical Faculty, University of Leipzig, Leipzig, Germany
* E-mail: info@life.uni-leipzig.de
* Homepage: https://www.uniklinikum-leipzig.de/einrichtungen/life/kontakt (currently only available in German)

# Source File

You will need to customize a source file, indicating

* path to R library
* path to Phyton and Conda 
* path to MetaXcan (github repro)
* path to summary-gwas-imputation (github repro)
* path to PLINK
* path to 1000Genomes data
* path to [GTEx v8 Data](https://gtexportal.org/home/protectedDataAccess)
* path to CAD statistics [van der Harst et al.](https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.117.312086?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed), [summary statistics](https://data.mendeley.com/datasets/gbbsrpx6bs/1)

# Included (tbc)

1) Harmonization calls to lift from hg19 to hg38
2) Hierarchical FDR for combined, males and females
3) Lead SNPs: sex-interaction test --> **Figure S4**
4) Lead SNPs: comparison with effect on gene expression (GTEx v8 data) --> **Figure 2**
5) Circular Plot --> **Figure 1**
6) Locus level: Co-localization of pQTLs and eQTLs (GTEx v8 data)  --> **Figure S6**
7) Locus level: Correlation of protein levels and genetically regulated gene expression using MetaXcan
8) Causal analyses: Test causality of gene expression on protein levels
9) Causal analyses: Test causality of protein levels on CAD (van der Haarst data) --> **Figure S8**
10) Causal analyses: Test mediating effect of protein levels on GE-CAD relation

# Not included (reason)

1) Plink Calls to generate association statistics (individual level data)
2) Priority Pruning with LIFE-Adult data (individual level data)
3) Annotation Pipeline (not yet published from GenStatLeipzig Group)
4) Figures:
    * Figure 4: DAGs MR Mediation (generated in PowerPoint, based on Table S12)
    * Figure S1: Sample Overview (individual level data)
    * Figure S2: Flowchart of Analysis Plan (generated in PowerPoint)
    * Figure S3: Genetic PC Plot (individual level data)
    
