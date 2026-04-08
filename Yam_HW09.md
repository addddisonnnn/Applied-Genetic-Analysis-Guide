## Homework 9: Polygenic (Risk) Scores
### BS859 Applied Genetic Analysis
### Addison Yam
### April 8, 2026

```bash
# load the necessary modules
module load R
module load python2
module load ldsc
# set the proper shortcut paths
export LDSCORES_DIR=/projectnb/bs859/data/ldscore_files
export UKBB_EUR_LD=$LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid
```

1.	We will explore in more detail using PRSice to develop Alzheimer Disease polygenic scores using the TGEN data and the IGAP Alzheimer disease summary statistics.  First change the clumping r2 parameter (see PRSice manual:  https://choishingwan.github.io/PRSice/command_detail/#clumping): Set the r2 threshold for clumping to 0.05, and then 0.15 instead of 0.1.   Then , use the default clump-r2 of 0.1, and change clump-kb 500.  
a.	Fill in the results in the table below.  

|                                                                  |     clump-r2 0.1 (default, in   class)    |     clump-r2 0.05 (new)    |     Clump r2 0.15 (new)    |     Clump r2=0.1, clump-kb =   500     |
|------------------------------------------------------------------|-------------------------------------------|----------------------------|----------------------------|----------------------------------------|
|     # of SNPs after clumping                                     |                                           |                            |                            |                                        |
|     # snps in optimal score                                      |                                           |                            |                            |                                        |
|     proportion of variance in   AD explained by optimal score    |                                           |                            |                            |                                        |
|     optimal p-value threshold                                    |                                           |                            |                            |                                        |
|     p-value of best PRS with   AD phenotype in TGEN              |                                           |                            |                            |                                        |

b.	Which clumping r2 threshold among those with clump-kb=250 explains the largest proportion of variance and has the smallest p-value?  Can we be sure this is the optimal clump-r2 parameter to use?

c.	Do your results for clump-kb 500 and clump-r2 0.1 change from –clump-kb 250 with --clump-r2 0.1?  If so, how? 
d.	Explain the differences in numbers of SNPs after clumping and in the optimal score across the four sets of parameters in the table.


2.	Lower hippocampal volume is associated with Alzheimer Disease.  This manuscript describes a GWAS of hippocampal volume: https://www.nature.com/articles/s41380-018-0262-7.  The GWAS summary statistics from the study were downloaded from here:  ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/vanderMeerD_30279459_GCST006871/Meta_Whole_hippocampus.zip
and are on the scc here: /projectnb/bs859/data/hippo 
zcat /projectnb/bs859/data/hippo/Meta_Whole_hippocampus.txt.gz |head
CHR POS SNP A1 A2 N P Z
1 100000012 rs10875231 T G 20 0.3166 -1.00146924840626
1 100000827 rs6678176 T C 21 0.7427 -0.328279953323096
1 100000843 rs78286437 C T 14 0.05404 1.92651580743038
1 100001201 rs76909621 T G 20 0.5975 -0.527999043244635
1 100002490 rs78642210 T C 14 0.06876 1.81999346559286
1 100002713 rs77140576 T C 20 0.263 -1.11932855055492
1 100002714 rs113470118 G A 14 0.06146 1.87017154359303
1 100002882 rs7545818 G T 21 0.7463 -0.323521932665517
1 100002991 rs75635821 A G 20 0.2449 -1.1628262611733

A1 is the coded allele, Z is the Z-statistic for association, and P is the p-value (determined from the Z statistic).  Ignore the “N” column, as this is NOT the sample size!)

Use the summary statistics to create a PRS using the default PRSice settings EXCEPT change the clumping KB distance to 500 from its default setting.  Include 1000 permutations to assess the significance of the most significant PRS, using seed 1443.

a.	Is the “best” PRS associated with Alzheimer Disease in the TGEN sample?  Explain.

b.	Present plots of the results for the thresholds that were tested.  What is the optimal p-value threshold, and what percentage of the variation in AD does this PRS explain?

c.	How many SNPs are in the best model?

d.	Present and interpret the association statistic for the final model, even if it is not significant.  (What is the direction of effect?)

e.	Summarize your findings:   Are the genetic variants that are associated with hippocampal volume associated with AD in the TGEN sample?

f.	What is the TGEN target sample size?  Do you think it is adequate to allow you to draw conclusions about this research question?  Can you think of a different way to go about testing whether hippocampal volume variants also affect AD that might be more powerful that what we have done here?
