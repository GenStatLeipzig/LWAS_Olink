# LWAS Olink

We are providing the main scripts used in the locus-wide association study (LWAS) of Olink proteins in LIFE-Adult, to empower other researchers to reproduce our results, starting from the summary statistics (for protein ~ SNP, gene expression ~ SNP and protein ~ gene expression associations in LIFE). These data are intended for research purposes only. The corresponding manuscript "Genetically regulated gene expression and proteins revealed discordant effects" of Janne Pott, Tarcyane Garcia, Stephanie M. Hauck et al. is submitted to PLOS ONE and currently under revision. 

Complete data sets including genetic data of LIFE-Adult participants cannot be made publicly available due to ethical and legal restrictions, as they are sufficient to identity study participants. This is not covered by the informed consent. However, access to the LIFE-Adult data is possible via project agreements addressed to:

-   LIFE Research Center for Civilization Diseases, Medical Faculty, University of Leipzig, Leipzig, Germany
-   E-mail: [info\@life.uni-leipzig.de](mailto:info@life.uni-leipzig.de)
-   Homepage: <https://www.uniklinikum-leipzig.de/einrichtungen/life/kontakt> (currently only available in German)


# Source File

You will need to customize a source file, indicating

-   path to R library (please use [R Version 4.x](https://cran.r-project.org/), all necessary packages are listed in the source file)
-   path to Phyton and Conda (e.g. [Anaconda](https://www.anaconda.com/products/individual))
-   path to [MetaXcan](https://github.com/hakyimlab/MetaXcan)
-   path to [summary-gwas-imputation](https://github.com/hakyimlab/summary-gwas-imputation)
-   path to [PLINK](https://www.cog-genomics.org/plink/2.0/)
-   path to [1000 Genomes Phase 3 EUR data](https://www.internationalgenome.org/data-portal/data-collection/phase-3)
-   path to [GTEx v8 Data](https://gtexportal.org/home/protectedDataAccess)
-   path to CAD statistics [van der Harst et al.](https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.117.312086?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed), [summary statistics](https://data.mendeley.com/datasets/gbbsrpx6bs/1)

# Included (tbc)

1) Harmonization calls to lift from hg19 to hg38 
2) Hierarchical FDR for combined, males and females 
3) Lead SNPs: sex-interaction test --> **Figure S4** 
4) Lead SNPs: comparison with effect on gene expression (GTEx v8 data) --> **Figure 2** 
5) Circular Plot --> **Figure 1** 
6) Locus level: Co-localization of pQTLs and eQTLs (GTEx v8 data) --> **Figure S6** 
7) Locus level: Correlation of protein levels and genetically regulated gene expression using MetaXcan 
8) Causal analyses: Test causality of gene expression on protein levels 
9) Causal analyses: Test causality of protein levels on CAD (van der Haarst data) --> **Figure S8** 
10) Causal analyses: Test mediating effect of protein levels on GE-CAD relation
11) Overview of all tissue-specific results (combined setting only)
12) **Main Tables** & **Figure 3**
13) **Supplemental Tables**
14) **Figures S5 and S7**

# Not included (reason)

1) Plink Calls to generate association statistics (individual level data) 
2) Priority Pruning with LIFE-Adult data (individual level data) 
3) Annotation Pipeline (not yet published from GenStatLeipzig Group) 
4) Figures: 
    * Figure 4: DAGs MR Mediation (generated in PowerPoint, based on Table S12) 
    * Figure S1: Sample Overview (individual level data) 
    * Figure S2: Flowchart of Analysis Plan (generated in PowerPoint) 
    * Figure S3: Genetic PC Plot (individual level data)
