# LWAS Olink

Scripts used in my Olink project (submitted to PLOS ONE - currently under revision)

Here, I want to include all relevant R scripts used in the LWAS Olink project. 

Complete data sets including genetic data of LIFE-Adult participants cannot be made publicly available due to ethical and legal restrictions, as they are sufficient to identity study participants. This is not covered by the informed consent. However, access to the LIFE-Adult data is possible via project agreements addressed to:

* LIFE Research Center for Civilization Diseases, Medical Faculty, University of Leipzig, Leipzig, Germany
* E-mail: info@life.uni-leipzig.de
* Homepage: https://www.uniklinikum-leipzig.de/einrichtungen/life/kontakt (currently only available in German)

# Source File

You will need to customize a source file, including

* Path to R library
* Path to Phyton and Conda 
* Path to MetaXcan (github repro)
* Path to summary-gwas-imputation (github repro)
* all relevant R packages

To do: generate template

# Included (tbc)

1) Harmonization calls to lift from hg19 to hg38
2) Hierarchical FDR for combined, males and females
3) Lead SNPs: sex-interaction test
4) Lead SNPs: comparison with effect on gene expression (GTEx v8 data)
5) Circular Plot
6) Locus level: Co-localization of pQTLs and eQTLs (GTEx v8 data)
7) Locus level: Correlation of protein levels and genetically regulated gene expression using MetaXcan
8) Causal analyses: Test causality of gene expression on protein levels
9) Causal analyses: Test causality of protein levels on CAD (van der Haarst data)
10) Causal analyses: Test mediating effect of protein levels on GE-CAD relation

# Not included (reason)

1) Plink Calls to generate association statistics (individual level data)
2) Priority Pruning with LIFE-Adult data (individual level data)
