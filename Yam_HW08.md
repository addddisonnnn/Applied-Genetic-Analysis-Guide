## Homework 8: Post-GWAS Analyses - Heritability
### BS859 Applied Genetic Analysis
### Addison Yam
### April 1, 2026

The Global lipids consortium published a meta-analysis of GWAS of several lipid traits; you can find the paper here: https://www.nature.com/articles/ng.2797  (Willer et al, Discovery and Refinement of Loci Associated with Lipid Levels. Nat Genet. 2013;45(11):1274-1283. doi:10.1038/ng.2797
The summary statistics for the blood levels of LDLand HDL cholesterol are in  '/projectnb/bs859/data/lipids'
There is a README file in the directory with information about the columns of the two summary statistics files: 
README 
jointGwasMc_HDL.txt.gz 
jointGwasMc_LDL.txt.gz

```bash
# load the necessary modules
module load R
module load python2
module load ldsc
# set the proper shortcut paths
export LDSCORES_DIR=/projectnb/bs859/data/ldscore_files
export UKBB_EUR_LD=$LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid
```

1)	Use LD score regression and the UK Biobank EUR LD scores to estimate the heritability of LDL cholesterol levels.  Show your code, report your findings and interpret your results.
- Answer: The heritability h² = 0.1918 (SE 0.0301), which means 19.2% of the variance is explained in LDL chloesterol by the common genetic variants in this sample. The LDSC intercept is 0.8893, which is less than 1, which is due to prior genomic control correction. And after correction is performed, the heritability estimate is still reliable. The Lambda GC = 1.0165, which is close to 1, meaning that there is minimal inflation. The mean X² = 1.1938, which is a bit over 1, which is what we expect for an extremely polygenic trait. 

```bash
munge_sumstats.py \
  --sumstats /projectnb/bs859/data/lipids/jointGwasMc_LDL.txt.gz \
  --snp rsid \
  --a1 A1 \
  --a2 A2 \
  --signed-sumstats Beta,0 \
  --p P-value \
  --N-col N \
  --merge-alleles $LDSCORES_DIR/w_hm3.snplist \
  --out LDL
Interpreting column names as follows:
P-value:        p-Value
rsid:   Variant ID (e.g., rs number)
N:      Sample size
A1:     Allele 1, interpreted as ref allele for signed sumstat.
beta:   Directional summary statistic as specified by --signed-sumstats.
A2:     Allele 2, interpreted as non-ref allele for signed sumstat.
Reading list of SNPs for allele merge from /projectnb/bs859/data/ldscore_files/w_hm3.snplis
t
Read 1217311 SNPs for allele merge.
Reading sumstats from /projectnb/bs859/data/lipids/jointGwasMc_LDL.txt.gz into memory 50000
00 SNPs at a time.
.WARNING: 3 SNPs had P outside of (0,1]. The P column may be mislabeled.
 done
Read 2437751 SNPs from --sumstats file.
Removed 1394207 SNPs not in --merge-alleles.
Removed 795 SNPs with missing values.
Removed 0 SNPs with INFO <= 0.9.
Removed 0 SNPs with MAF <= 0.01.
Removed 3 SNPs with out-of-bounds p-values.
Removed 5 variants that were not SNPs or were strand-ambiguous.
1042741 SNPs remain.
Removed 0 SNPs with duplicated rs numbers (1042741 SNPs remain).
Removed 9514 SNPs with N < 59925.3333333 (1033227 SNPs remain).
Median value of SIGNED_SUMSTATS was 0.0049, which seems sensible.
Removed 34 SNPs whose alleles did not match --merge-alleles (1033193 SNPs remain).
Writing summary statistics for 1217311 SNPs (1033193 with nonmissing beta) to LDL.sumstats.
gz.
Metadata:
Mean chi^2 = 1.194
Lambda GC = 1.015
Max chi^2 = 1244.05
1341 Genome-wide significant SNPs (some may have been removed by filtering).
Conversion finished at Tue Mar 31 23:38:01 2026
Total time elapsed: 44.64s

ldsc.py \
  --h2 LDL.sumstats.gz \
  --ref-ld $UKBB_EUR_LD \
  --w-ld $UKBB_EUR_LD \
  --out LDL_h2
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--h2 LDL.sumstats.gz \
--out LDL_h2 \
--w-ld /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--ref-ld /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid 
Beginning analysis at Tue Mar 31 23:38:43 2026
Reading summary statistics from LDL.sumstats.gz ...
Read summary statistics for 1033193 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/
UKBB.EUR.rsid ...
Read reference panel LD Scores for 1094844 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscor
e/UKBB.EUR.rsid ...
Read regression weight LD Scores for 1094844 SNPs.
After merging with reference panel LD, 954004 SNPs remain.
After merging with regression SNP LD, 954004 SNPs remain.
Using two-step estimator with cutoff at 30.
Total Observed scale h2: 0.1918 (0.0301)
Lambda GC: 1.0165
Mean Chi^2: 1.1938
Intercept: 0.8893 (0.0081)
Ratio < 0 (usually indicates GC correction).
Analysis finished at Tue Mar 31 23:39:51 2026
Total time elapsed: 1.0m:8.01s
```

2)	Use LD score regression and the UK Biobank EUR LD scores to estimate the heritability of HDL cholesterol levels. Show your code, report your findings and interpret your results.
- Answer: The heritability h² is = 0.2149 (SE 0.0296), which means 29.6% of the variance is explained in HDL chloesterol by the common genetic variants in this sample. The LDSC intercept is 0.851, which is less than 1, which is due to prior genomic control correction. And after correction is performed, the heritability estimate is still reliable. The Lambda GC = 1.0225, which is close to 1, meaning that there is minimal inflation. The mean X² = 1.2097, which is a bit over 1, which is what we expect for an extremely polygenic trait. 

```bash
munge_sumstats.py \
  --sumstats /projectnb/bs859/data/lipids/jointGwasMc_HDL.txt.gz \
  --snp rsid \
  --a1 A1 \
  --a2 A2 \
  --signed-sumstats Beta,0 \
  --p P-value \
  --N-col N \
  --merge-alleles $LDSCORES_DIR/w_hm3.snplist \
  --out HDL
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./munge_sumstats.py \
--signed-sumstats Beta,0 \
--out HDL \
--merge-alleles /projectnb/bs859/data/ldscore_files/w_hm3.snplist \
--N-col N \
--a1 A1 \
--a2 A2 \
--snp rsid \
--sumstats /projectnb/bs859/data/lipids/jointGwasMc_HDL.txt.gz \
--p P-value 
Interpreting column names as follows:
P-value:        p-Value
rsid:   Variant ID (e.g., rs number)
N:      Sample size
A1:     Allele 1, interpreted as ref allele for signed sumstat.
beta:   Directional summary statistic as specified by --signed-sumstats.
A2:     Allele 2, interpreted as non-ref allele for signed sumstat.
Reading list of SNPs for allele merge from /projectnb/bs859/data/ldscore_files/w_hm3.snplis
t
Read 1217311 SNPs for allele merge.
Reading sumstats from /projectnb/bs859/data/lipids/jointGwasMc_HDL.txt.gz into memory 50000
00 SNPs at a time.
.WARNING: 11 SNPs had P outside of (0,1]. The P column may be mislabeled.
 done
Read 2447441 SNPs from --sumstats file.
Removed 1401331 SNPs not in --merge-alleles.
Removed 803 SNPs with missing values.
Removed 0 SNPs with INFO <= 0.9.
Removed 0 SNPs with MAF <= 0.01.
Removed 11 SNPs with out-of-bounds p-values.
Removed 5 variants that were not SNPs or were strand-ambiguous.
1045291 SNPs remain.
Removed 0 SNPs with duplicated rs numbers (1045291 SNPs remain).
Removed 10925 SNPs with N < 62874.0 (1034366 SNPs remain).
Median value of SIGNED_SUMSTATS was 0.0046, which seems sensible.
Removed 34 SNPs whose alleles did not match --merge-alleles (1034332 SNPs remain).
Writing summary statistics for 1217311 SNPs (1034332 with nonmissing beta) to HDL.sumstats.gz.

Metadata:
Mean chi^2 = 1.21
Lambda GC = 1.02
Max chi^2 = 858.182
1557 Genome-wide significant SNPs (some may have been removed by filtering).
Conversion finished at Tue Mar 31 23:49:59 2026
Total time elapsed: 45.32s

ldsc.py \
  --h2 HDL.sumstats.gz \
  --ref-ld $UKBB_EUR_LD \
  --w-ld $UKBB_EUR_LD \
  --out HDL_h2
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--h2 HDL.sumstats.gz \
--out HDL_h2 \
--w-ld /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--ref-ld /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid 

Beginning analysis at Tue Mar 31 23:50:06 2026
Reading summary statistics from HDL.sumstats.gz ...
Read summary statistics for 1034332 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid ...
Read reference panel LD Scores for 1094844 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid ...
Read regression weight LD Scores for 1094844 SNPs.
After merging with reference panel LD, 954875 SNPs remain.
After merging with regression SNP LD, 954875 SNPs remain.
Using two-step estimator with cutoff at 30.
Total Observed scale h2: 0.2149 (0.0296)
Lambda GC: 1.0225
Mean Chi^2: 1.2097
Intercept: 0.851 (0.0092)
Ratio < 0 (usually indicates GC correction).
Analysis finished at Tue Mar 31 23:51:27 2026
Total time elapsed: 1.0m:20.34s
```

3)	Use LD score regression and the UK Biobank EUR LD scores to estimate the genetic correlations between LDL, HDL, and Alzheimer disease.  Show your code, report your findings, and interpret your result.
- Answer: Between LDL and HDL, the genetic correlation is -0.0835 (SE 0.0486), the Z-score is -1.7157, and the p-value is 0.0862. This isn't statistically significant and the negative correlation means that alleles increasing LDL may cause a decrease for HDL, but this strongly supported and may be potentially de to chance. Between LDL and AD, the genetic correlation is 0.0379 (SE 0.0544), the Z-score is 0.6965, and the p-value is 0.4861. This isn't statistically significant and there isn't evidence of shared genetic correlation between LDL and AD. Between HDL and AD, the genetic correlation is 0.1746 (SE 0.0639), the Z-score is 2.7324, and the p-value is 0.0063. This is statistically signiciant and a positive genetic correlation is shown between HDL and AD, so genetic variants that increase LDL also increase the chance of AD or the other way around.

```bash
munge_sumstats.py \
  --sumstats /projectnb/bs859/data/igap/Kunkle_etal_Stage1_results2019.txt.gz \
  --snp MarkerName \
  --N-cas 21982 \
  --N-con 41944 \
  --a1 Effect_allele \
  --a2 Non_Effect_allele \
  --signed-sumstats Beta,0 \
  --merge-alleles $LDSCORES_DIR/w_hm3.snplist \
  --out AD
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./munge_sumstats.py \
--signed-sumstats Beta,0 \
--out AD \
--merge-alleles /projectnb/bs859/data/ldscore_files/w_hm3.snplist \
--N-con 41944.0 \
--N-cas 21982.0 \
--a1 Effect_allele \
--a2 Non_Effect_allele \
--snp MarkerName \
--sumstats /projectnb/bs859/data/igap/Kunkle_etal_Stage1_results2019.txt.gz 

Interpreting column names as follows:
Effect_allele:  Allele 1, interpreted as ref allele for signed sumstat.
MarkerName:     Variant ID (e.g., rs number)
Beta:   Directional summary statistic as specified by --signed-sumstats.
Pvalue: p-Value
Non_Effect_allele:      Allele 2, interpreted as non-ref allele for signed sumstat.

Reading list of SNPs for allele merge from /projectnb/bs859/data/ldscore_files/w_hm3.snplist
Read 1217311 SNPs for allele merge.
Reading sumstats from /projectnb/bs859/data/igap/Kunkle_etal_Stage1_results2019.txt.gz into memory 5000000 SNPs at a time.
...WARNING: 1 SNPs had P outside of (0,1]. The P column may be mislabeled.
 done
Read 11480632 SNPs from --sumstats file.
Removed 10271272 SNPs not in --merge-alleles.
Removed 2287 SNPs with missing values.
Removed 0 SNPs with INFO <= 0.9.
Removed 0 SNPs with MAF <= 0.01.
Removed 1 SNPs with out-of-bounds p-values.
Removed 0 variants that were not SNPs or were strand-ambiguous.
1207072 SNPs remain.
Removed 0 SNPs with duplicated rs numbers (1207072 SNPs remain).
Median value of SIGNED_SUMSTATS was 0.0, which seems sensible.
Removed 0 SNPs whose alleles did not match --merge-alleles (1207072 SNPs remain).
Writing summary statistics for 1217311 SNPs (1207072 with nonmissing beta) to AD.sumstats.gz.

Metadata:
Mean chi^2 = 1.117
Lambda GC = 1.091
Max chi^2 = 689.858
217 Genome-wide significant SNPs (some may have been removed by filtering).

Conversion finished at Wed Apr  1 00:02:37 2026
Total time elapsed: 1.0m:51.35s

# LDL vs HDL
ldsc.py --rg LDL.sumstats.gz,HDL.sumstats.gz --ref-ld $UKBB_EUR_LD --w-ld $UKBB_EUR_LD --out LDL_HDL_rg
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--out LDL_HDL_rg \
--rg LDL.sumstats.gz,HDL.sumstats.gz \
--w-ld /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--ref-ld /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid 

Beginning analysis at Wed Apr  1 00:09:04 2026
Reading summary statistics from LDL.sumstats.gz ...
Read summary statistics for 1033193 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid ...
Read reference panel LD Scores for 1094844 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid ...
Read regression weight LD Scores for 1094844 SNPs.
After merging with reference panel LD, 954004 SNPs remain.
After merging with regression SNP LD, 954004 SNPs remain.
Computing rg for phenotype 2/2
Reading summary statistics from HDL.sumstats.gz ...
Read summary statistics for 1217311 SNPs.
After merging with summary statistics, 954004 SNPs remain.
953731 SNPs with valid alleles.

Heritability of phenotype 1
---------------------------
Total Observed scale h2: 0.2074 (0.0387)
Lambda GC: 1.0165
Mean Chi^2: 1.193
Intercept: 0.864 (0.0306)
Ratio < 0 (usually indicates GC correction).

Heritability of phenotype 2/2
-----------------------------
Total Observed scale h2: 0.2237 (0.0347)
Lambda GC: 1.0225
Mean Chi^2: 1.21
Intercept: 0.8375 (0.0409)
Ratio < 0 (usually indicates GC correction).

Genetic Covariance
------------------
Total Observed scale gencov: -0.018 (0.0112)
Mean z1*z2: -0.0989
Intercept: -0.0858 (0.008)

Genetic Correlation
-------------------
Genetic Correlation: -0.0835 (0.0486)
Z-score: -1.7157
P: 0.0862


Summary of Genetic Correlation Results
p1               p2      rg      se       z       p  h2_obs  h2_obs_se  h2_int  h2_int_se  gcov_int  gcov_int_se
LDL.sumstats.gz  HDL.sumstats.gz -0.0835  0.0486 -1.7157  0.0862  0.2237     0.0347  0.8375     0.0409   -0.0858        0.008

Analysis finished at Wed Apr  1 00:10:18 2026
Total time elapsed: 1.0m:13.7s

# LDL vs AD
ldsc.py --rg LDL.sumstats.gz,AD.sumstats.gz --ref-ld $UKBB_EUR_LD --w-ld $UKBB_EUR_LD --out LDL_AD_rg
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--out LDL_AD_rg \
--rg LDL.sumstats.gz,AD.sumstats.gz \
--w-ld /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--ref-ld /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid 

Beginning analysis at Wed Apr  1 00:10:41 2026
Reading summary statistics from LDL.sumstats.gz ...
Read summary statistics for 1033193 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid ...
Read reference panel LD Scores for 1094844 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid ...
Read regression weight LD Scores for 1094844 SNPs.
After merging with reference panel LD, 954004 SNPs remain.
After merging with regression SNP LD, 954004 SNPs remain.
Computing rg for phenotype 2/2
Reading summary statistics from AD.sumstats.gz ...
Read summary statistics for 1217311 SNPs.
After merging with summary statistics, 954004 SNPs remain.
953969 SNPs with valid alleles.

Heritability of phenotype 1
---------------------------
Total Observed scale h2: 0.2077 (0.0386)
Lambda GC: 1.0165
Mean Chi^2: 1.1928
Intercept: 0.863 (0.0301)
Ratio < 0 (usually indicates GC correction).

Heritability of phenotype 2/2
-----------------------------
Total Observed scale h2: 0.0458 (0.014)
Lambda GC: 1.0926
Mean Chi^2: 1.1163
Intercept: 1.0545 (0.027)
Ratio: 0.4685 (0.2323)

Genetic Covariance
------------------
Total Observed scale gencov: 0.0037 (0.0054)
Mean z1*z2: 0.0181
Intercept: 0.0136 (0.0102)

Genetic Correlation
-------------------
Genetic Correlation: 0.0379 (0.0544)
Z-score: 0.6965
P: 0.4861


Summary of Genetic Correlation Results
p1              p2      rg      se       z       p  h2_obs  h2_obs_se  h2_int  h2_int_se  gcov_int  gcov_int_se
LDL.sumstats.gz  AD.sumstats.gz  0.0379  0.0544  0.6965  0.4861  0.0458      0.014  1.0545      0.027    0.0136       0.0102

Analysis finished at Wed Apr  1 00:11:54 2026
Total time elapsed: 1.0m:13.03s

# HDL vs AD
ldsc.py --rg HDL.sumstats.gz,AD.sumstats.gz --ref-ld $UKBB_EUR_LD --w-ld $UKBB_EUR_LD --out HDL_AD_rg
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--out HDL_AD_rg \
--rg HDL.sumstats.gz,AD.sumstats.gz \
--w-ld /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--ref-ld /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid 
Beginning analysis at Wed Apr  1 00:11:54 2026
Reading summary statistics from HDL.sumstats.gz ...
Read summary statistics for 1034332 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB
.EUR.rsid ...
Read reference panel LD Scores for 1094844 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid ...
Read regression weight LD Scores for 1094844 SNPs.
After merging with reference panel LD, 954875 SNPs remain.
After merging with regression SNP LD, 954875 SNPs remain.
Computing rg for phenotype 2/2
Reading summary statistics from AD.sumstats.gz ...
Read summary statistics for 1217311 SNPs.
After merging with summary statistics, 954875 SNPs remain.
954841 SNPs with valid alleles.

Heritability of phenotype 1
---------------------------
Total Observed scale h2: 0.2237 (0.0347)
Lambda GC: 1.0225
Mean Chi^2: 1.2095
Intercept: 0.8373 (0.0407)
Ratio < 0 (usually indicates GC correction).

Heritability of phenotype 2/2
-----------------------------
Total Observed scale h2: 0.0454 (0.0141)
Lambda GC: 1.0926
Mean Chi^2: 1.1163
Intercept: 1.0551 (0.0273)
Ratio: 0.4735 (0.2344)

Genetic Covariance
------------------
Total Observed scale gencov: 0.0176 (0.0059)
Mean z1*z2: 0.0066
Intercept: -0.0261 (0.0067)

Genetic Correlation
-------------------
Genetic Correlation: 0.1746 (0.0639)
Z-score: 2.7324
P: 0.0063


Summary of Genetic Correlation Results
p1              p2      rg      se       z       p  h2_obs  h2_obs_se  h2_int  h2_int_se  gcov_int  gcov_int_se
HDL.sumstats.gz  AD.sumstats.gz  0.1746  0.0639  2.7324  0.0063  0.0454     0.0141  1.0551     0.0273   -0.0261       0.0067

Analysis finished at Wed Apr  1 00:13:12 2026
Total time elapsed: 1.0m:18.2s
```

4)	Compare the heritability estimates reported in the genetic correlation analysis to heritability reported when using the –h option for AD, LDL, and HDL. If they are different, identify the differences between the two analyses that might cause different estimates.
- Answer: From the --h2 runs: for LDL, the heritability h² = 0.1918 (0.0301), for HDL, the heritability h² = 0.2149 (0.0296) and for AD, the heritability h² = 0.0458 (0.014). From the --rg outputs, for LDL, the heritability h² = 0.02074 (0.0387), for HDL, the heritability h² = 0.02237 (0.0347) and for AD, the heritability h² = 0.0454 (0.014). The heritability numbers do differ for LDL and HDL. And for AD, there is overlap between the confidence intervals. The differences may be due to using different sets of SNPs as the cross-trait analysis includes SNPs that are present in both summary statistics, which can slightly change which SNPs are used. Differences can also be due to sample size differences, GC correction, and SNP filtering. 

5)	Using LD score regression in class we determined the enrichment for annotations specific to CNS cell types in Alzheimer disease.  Run partitioned ld score regression for the other cell types, and report the enrichment scores, SEs, and p-values for all 10 cell types (including the one we did in class) in a table.  Use population prevalence of 0.10 for AD. Which cell types are most enriched for heritability of AD? Which has the greatest enrichment (look at the enrichment score), and which has the most significant enrichment (look at p-value)?  In your answer, be sure to consider whether a multiple testing correction should be applied.
- Answer: A multiple testing correction should be applied and after applying a multiple testing correction of 0.05/10 = 0.005, the cardiovascular and kidney groups do not reach statistical significance even though the enrichment score of the kidney group is high. The cell types that are the most enriched for heritability of AD are the liver, kidney, adrenal pancreas, and connective bones. The liver has the greatest enrichment at 9.96176003472. The hematopoietic is the most significant enrichment with a p-value of 2.8765959145279837e-07.  

|     Cell type           |     Enrichment score    |     SE    |     p-value (unadjusted)    |
|-------------------------|-------------------------|-----------|-----------------------------|
|     Adrenal_Pancreas    |5.73752919162|1.70698512449|0.0021881685758397842|
|     Cardiovascular      |2.47683904836|1.25954237272|0.26250189310893168|
|     CNS                 |4.45221852312| 0.957346487366 |0.00020583325421419109|
|     Connective_Bone     |5.72367480954 |1.45249272873| 4.8374550100772444e-05|
|     GI                  |4.10062787731| 0.845104692174 |0.00066623764205058823|
|     Hematopoietic       |4.96955736469 |0.853198334309| 2.8765959145279837e-07|
|     Kidney              |6.92913011059 |2.45754590367| 0.0099291283932174731|
|     Liver               |9.96176003472 |2.60554633952| 1.7871904506649469e-05|
|     Other               |4.85626047693 |1.00445465504| 2.5687901929833607e-05|
|     SkeletalMuscle      |5.67939535601 |1.35565187875| 0.00032979073414727102|
```bash
cat $LDSCORES_DIR/1000G_Phase3_cell_type_groups/names
file_num        cell_type
1       Adrenal_Pancreas.bed
2       Cardiovascular.bed
3       CNS.bed
4       Connective_Bone.bed
5       GI.bed
6       Hematopoietic.bed
7       Kidney.bed
8       Liver.bed
9       Other.bed
10      SkeletalMuscle.bed

# I ran this ten times for each group from 1 to 10, this is what it looks like for group 1
ldsc.py \
  --h2 AD.sumstats.gz \
  --ref-ld-chr $LDSCORES_DIR/1000G_Phase3_cell_type_groups/cell_type_group.1.,$LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline. \
  --w-ld-chr $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --overlap-annot \
  --frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
  --pop-prev 0.10 \
  --samp-prev 0.344 \
  --print-coefficients \
  --out AD_celltype_1

# group 1 output
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.1.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline. \
--out AD_celltype_1 \
--overlap-annot  \
--pop-prev 0.10 \
--samp-prev 0.344 \
--frqfile-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients  

Beginning analysis at Wed Apr  1 01:07:07 2026
Reading summary statistics from AD.sumstats.gz ...
Read summary statistics for 1207072 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.1.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Read reference panel LD Scores for 1190321 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.[1-22] ...
Read regression weight LD Scores for 1187349 SNPs.
After merging with reference panel LD, 1186816 SNPs remain.
After merging with regression SNP LD, 1183845 SNPs remain.
Removed 35 SNPs with chi^2 > 80 (1183810 SNPs remain)
Total Liability scale h2: 0.072 (0.0117)
Categories: L2_0 baseL2_1 Coding_UCSC.bedL2_1 Coding_UCSC.extend.500.bedL2_1 Conserved_LindbladToh.bedL2_1 Conserved_LindbladToh.extend.500.bedL2_1 CTCF_Hoffman.bedL2_1 CTCF_Hoffman.extend.500.bedL2_1 DGF_ENCODE.bedL2_1 DGF_ENCODE.extend.500.bedL2_1 DHS_peaks_Trynka.bedL2_1 DHS_Trynka.bedL2_1 DHS_Trynka.extend.500.bedL2_1 Enhancer_Andersson.bedL2_1 Enhancer_Andersson.extend.500.bedL2_1 Enhancer_Hoffman.bedL2_1 Enhancer_Hoffman.extend.500.bedL2_1 FetalDHS_Trynka.bedL2_1 FetalDHS_Trynka.extend.500.bedL2_1 H3K27ac_Hnisz.bedL2_1 H3K27ac_Hnisz.extend.500.bedL2_1 H3K27ac_PGC2.bedL2_1 H3K27ac_PGC2.extend.500.bedL2_1 H3K4me1_peaks_Trynka.bedL2_1 H3K4me1_Trynka.bedL2_1 H3K4me1_Trynka.extend.500.bedL2_1 H3K4me3_peaks_Trynka.bedL2_1 H3K4me3_Trynka.bedL2_1 H3K4me3_Trynka.extend.500.bedL2_1 H3K9ac_peaks_Trynka.bedL2_1 H3K9ac_Trynka.bedL2_1 H3K9ac_Trynka.extend.500.bedL2_1 Intron_UCSC.bedL2_1 Intron_UCSC.extend.500.bedL2_1 PromoterFlanking_Hoffman.bedL2_1 PromoterFlanking_Hoffman.extend.500.bedL2_1 Promoter_UCSC.bedL2_1 Promoter_UCSC.extend.500.bedL2_1 Repressed_Hoffman.bedL2_1 Repressed_Hoffman.extend.500.bedL2_1 SuperEnhancer_Hnisz.bedL2_1 SuperEnhancer_Hnisz.extend.500.bedL2_1 TFBS_ENCODE.bedL2_1 TFBS_ENCODE.extend.500.bedL2_1 Transcribed_Hoffman.bedL2_1 Transcribed_Hoffman.extend.500.bedL2_1 TSS_Hoffman.bedL2_1 TSS_Hoffman.extend.500.bedL2_1 UTR_3_UCSC.bedL2_1 UTR_3_UCSC.extend.500.bedL2_1 UTR_5_UCSC.bedL2_1 UTR_5_UCSC.extend.500.bedL2_1 WeakEnhancer_Hoffman.bedL2_1 WeakEnhancer_Hoffman.extend.500.bedL2_1
Lambda GC: 1.0895
Mean Chi^2: 1.1093
Intercept: 1.0311 (0.0121)
Ratio: 0.2845 (0.1103)
Reading annot matrix from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.1.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Results printed to AD_celltype_1.results
Analysis finished at Wed Apr  1 01:08:46 2026
Total time elapsed: 1.0m:39.55s

# group 2
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.2.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline. \
--out AD_celltype_2 \
--overlap-annot  \
--pop-prev 0.10 \
--samp-prev 0.344 \
--frqfile-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients  

Beginning analysis at Wed Apr  1 01:08:48 2026
Reading summary statistics from AD.sumstats.gz ...
Read summary statistics for 1207072 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.2.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Read reference panel LD Scores for 1190321 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.[1-22] ...
Read regression weight LD Scores for 1187349 SNPs.
After merging with reference panel LD, 1186816 SNPs remain.
After merging with regression SNP LD, 1183845 SNPs remain.
Removed 35 SNPs with chi^2 > 80 (1183810 SNPs remain)
Total Liability scale h2: 0.0724 (0.0117)
Categories: L2_0 baseL2_1 Coding_UCSC.bedL2_1 Coding_UCSC.extend.500.bedL2_1 Conserved_LindbladToh.bedL2_1 Conserved_LindbladToh.extend.500.bedL2_1 CTCF_Hoffman.bedL2_1 CTCF_Hoffman.extend.500.bedL2_1 DGF_ENCODE.bedL2_1 DGF_ENCODE.extend.500.bedL2_1 DHS_peaks_Trynka.bedL2_1 DHS_Trynka.bedL2_1 DHS_Trynka.extend.500.bedL2_1 Enhancer_Andersson.bedL2_1 Enhancer_Andersson.extend.500.bedL2_1 Enhancer_Hoffman.bedL2_1 Enhancer_Hoffman.extend.500.bedL2_1 FetalDHS_Trynka.bedL2_1 FetalDHS_Trynka.extend.500.bedL2_1 H3K27ac_Hnisz.bedL2_1 H3K27ac_Hnisz.extend.500.bedL2_1 H3K27ac_PGC2.bedL2_1 H3K27ac_PGC2.extend.500.bedL2_1 H3K4me1_peaks_Trynka.bedL2_1 H3K4me1_Trynka.bedL2_1 H3K4me1_Trynka.extend.500.bedL2_1 H3K4me3_peaks_Trynka.bedL2_1 H3K4me3_Trynka.bedL2_1 H3K4me3_Trynka.extend.500.bedL2_1 H3K9ac_peaks_Trynka.bedL2_1 H3K9ac_Trynka.bedL2_1 H3K9ac_Trynka.extend.500.bedL2_1 Intron_UCSC.bedL2_1 Intron_UCSC.extend.500.bedL2_1 PromoterFlanking_Hoffman.bedL2_1 PromoterFlanking_Hoffman.extend.500.bedL2_1 Promoter_UCSC.bedL2_1 Promoter_UCSC.extend.500.bedL2_1 Repressed_Hoffman.bedL2_1 Repressed_Hoffman.extend.500.bedL2_1 SuperEnhancer_Hnisz.bedL2_1 SuperEnhancer_Hnisz.extend.500.bedL2_1 TFBS_ENCODE.bedL2_1 TFBS_ENCODE.extend.500.bedL2_1 Transcribed_Hoffman.bedL2_1 Transcribed_Hoffman.extend.500.bedL2_1 TSS_Hoffman.bedL2_1 TSS_Hoffman.extend.500.bedL2_1 UTR_3_UCSC.bedL2_1 UTR_3_UCSC.extend.500.bedL2_1 UTR_5_UCSC.bedL2_1 UTR_5_UCSC.extend.500.bedL2_1 WeakEnhancer_Hoffman.bedL2_1 WeakEnhancer_Hoffman.extend.500.bedL2_1
Lambda GC: 1.0895
Mean Chi^2: 1.1093
Intercept: 1.0308 (0.012)
Ratio: 0.2819 (0.1099)
Reading annot matrix from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.2.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Results printed to AD_celltype_2.results
Analysis finished at Wed Apr  1 01:10:27 2026
Total time elapsed: 1.0m:39.08s

# group 3 
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.
3.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline. \
--out AD_celltype_3 \
--overlap-annot  \
--pop-prev 0.10 \
--samp-prev 0.344 \
--frqfile-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMH
C. \
--print-coefficients  
Beginning analysis at Wed Apr  1 01:10:28 2026
Reading summary statistics from AD.sumstats.gz ...
Read summary statistics for 1207072 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_typ
e_groups/cell_type_group.3.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/basel
ine.[1-22] ...
Read reference panel LD Scores for 1190321 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_weight
s_hm3_no_MHC/weights.hm3_noMHC.[1-22] ...
Read regression weight LD Scores for 1187349 SNPs.
After merging with reference panel LD, 1186816 SNPs remain.
After merging with regression SNP LD, 1183845 SNPs remain.
Removed 35 SNPs with chi^2 > 80 (1183810 SNPs remain)
Total Liability scale h2: 0.0721 (0.0118)
Categories: L2_0 baseL2_1 Coding_UCSC.bedL2_1 Coding_UCSC.extend.500.bedL2_1 Conserved_Lindblad
Toh.bedL2_1 Conserved_LindbladToh.extend.500.bedL2_1 CTCF_Hoffman.bedL2_1 CTCF_Hoffman.extend.5
00.bedL2_1 DGF_ENCODE.bedL2_1 DGF_ENCODE.extend.500.bedL2_1 DHS_peaks_Trynka.bedL2_1 DHS_Trynka
.bedL2_1 DHS_Trynka.extend.500.bedL2_1 Enhancer_Andersson.bedL2_1 Enhancer_Andersson.extend.500
.bedL2_1 Enhancer_Hoffman.bedL2_1 Enhancer_Hoffman.extend.500.bedL2_1 FetalDHS_Trynka.bedL2_1 F
etalDHS_Trynka.extend.500.bedL2_1 H3K27ac_Hnisz.bedL2_1 H3K27ac_Hnisz.extend.500.bedL2_1 H3K27a
c_PGC2.bedL2_1 H3K27ac_PGC2.extend.500.bedL2_1 H3K4me1_peaks_Trynka.bedL2_1 H3K4me1_Trynka.bedL
2_1 H3K4me1_Trynka.extend.500.bedL2_1 H3K4me3_peaks_Trynka.bedL2_1 H3K4me3_Trynka.bedL2_1 H3K4m
e3_Trynka.extend.500.bedL2_1 H3K9ac_peaks_Trynka.bedL2_1 H3K9ac_Trynka.bedL2_1 H3K9ac_Trynka.ex
tend.500.bedL2_1 Intron_UCSC.bedL2_1 Intron_UCSC.extend.500.bedL2_1 PromoterFlanking_Hoffman.be
dL2_1 PromoterFlanking_Hoffman.extend.500.bedL2_1 Promoter_UCSC.bedL2_1 Promoter_UCSC.extend.50
0.bedL2_1 Repressed_Hoffman.bedL2_1 Repressed_Hoffman.extend.500.bedL2_1 SuperEnhancer_Hnisz.be
dL2_1 SuperEnhancer_Hnisz.extend.500.bedL2_1 TFBS_ENCODE.bedL2_1 TFBS_ENCODE.extend.500.bedL2_1
 Transcribed_Hoffman.bedL2_1 Transcribed_Hoffman.extend.500.bedL2_1 TSS_Hoffman.bedL2_1 TSS_Hof
fman.extend.500.bedL2_1 UTR_3_UCSC.bedL2_1 UTR_3_UCSC.extend.500.bedL2_1 UTR_5_UCSC.bedL2_1 UTR
_5_UCSC.extend.500.bedL2_1 WeakEnhancer_Hoffman.bedL2_1 WeakEnhancer_Hoffman.extend.500.bedL2_1
Lambda GC: 1.0895
Mean Chi^2: 1.1093
Intercept: 1.031 (0.0122)
Ratio: 0.284 (0.1119)
Reading annot matrix from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cel
l_type_group.3.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] .
..
Results printed to AD_celltype_3.results
Analysis finished at Wed Apr  1 01:12:08 2026
Total time elapsed: 1.0m:39.97s

# group 4
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.4.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline. \
--out AD_celltype_4 \
--overlap-annot  \
--pop-prev 0.10 \
--samp-prev 0.344 \
--frqfile-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients  

Beginning analysis at Wed Apr  1 01:12:10 2026
Reading summary statistics from AD.sumstats.gz ...
Read summary statistics for 1207072 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.4.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Read reference panel LD Scores for 1190321 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.[1-22] ...
Read regression weight LD Scores for 1187349 SNPs.
After merging with reference panel LD, 1186816 SNPs remain.
After merging with regression SNP LD, 1183845 SNPs remain.
Removed 35 SNPs with chi^2 > 80 (1183810 SNPs remain)
Total Liability scale h2: 0.0712 (0.0119)
Categories: L2_0 baseL2_1 Coding_UCSC.bedL2_1 Coding_UCSC.extend.500.bedL2_1 Conserved_LindbladToh.bedL2_1 Conserved_LindbladToh.extend.500.bedL2_1 CTCF_Hoffman.bedL2_1 CTCF_Hoffman.extend.500.bedL2_1 DGF_ENCODE.bedL2_1 DGF_ENCODE.extend.500.bedL2_1 DHS_peaks_Trynka.bedL2_1 DHS_Trynka.bedL2_1 DHS_Trynka.extend.500.bedL2_1 Enhancer_Andersson.bedL2_1 Enhancer_Andersson.extend.500.bedL2_1 Enhancer_Hoffman.bedL2_1 Enhancer_Hoffman.extend.500.bedL2_1 FetalDHS_Trynka.bedL2_1 FetalDHS_Trynka.extend.500.bedL2_1 H3K27ac_Hnisz.bedL2_1 H3K27ac_Hnisz.extend.500.bedL2_1 H3K27ac_PGC2.bedL2_1 H3K27ac_PGC2.extend.500.bedL2_1 H3K4me1_peaks_Trynka.bedL2_1 H3K4me1_Trynka.bedL2_1 H3K4me1_Trynka.extend.500.bedL2_1 H3K4me3_peaks_Trynka.bedL2_1 H3K4me3_Trynka.bedL2_1 H3K4me3_Trynka.extend.500.bedL2_1 H3K9ac_peaks_Trynka.bedL2_1 H3K9ac_Trynka.bedL2_1 H3K9ac_Trynka.extend.500.bedL2_1 Intron_UCSC.bedL2_1 Intron_UCSC.extend.500.bedL2_1 PromoterFlanking_Hoffman.bedL2_1 PromoterFlanking_Hoffman.extend.500.bedL2_1 Promoter_UCSC.bedL2_1 Promoter_UCSC.extend.500.bedL2_1 Repressed_Hoffman.bedL2_1 Repressed_Hoffman.extend.500.bedL2_1 SuperEnhancer_Hnisz.bedL2_1 SuperEnhancer_Hnisz.extend.500.bedL2_1 TFBS_ENCODE.bedL2_1 TFBS_ENCODE.extend.500.bedL2_1 Transcribed_Hoffman.bedL2_1 Transcribed_Hoffman.extend.500.bedL2_1 TSS_Hoffman.bedL2_1 TSS_Hoffman.extend.500.bedL2_1 UTR_3_UCSC.bedL2_1 UTR_3_UCSC.extend.500.bedL2_1 UTR_5_UCSC.bedL2_1 UTR_5_UCSC.extend.500.bedL2_1 WeakEnhancer_Hoffman.bedL2_1 WeakEnhancer_Hoffman.extend.500.bedL2_1
Lambda GC: 1.0895
Mean Chi^2: 1.1093
Intercept: 1.0316 (0.0121)
Ratio: 0.2889 (0.111)
Reading annot matrix from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.4.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Results printed to AD_celltype_4.results
Analysis finished at Wed Apr  1 01:13:50 2026
Total time elapsed: 1.0m:39.8s

# group 5
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.5.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline. \
--out AD_celltype_5 \
--overlap-annot  \
--pop-prev 0.10 \
--samp-prev 0.344 \
--frqfile-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients  

Beginning analysis at Wed Apr  1 01:13:51 2026
Reading summary statistics from AD.sumstats.gz ...
Read summary statistics for 1207072 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.5.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Read reference panel LD Scores for 1190321 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.[1-22] ...
Read regression weight LD Scores for 1187349 SNPs.
After merging with reference panel LD, 1186816 SNPs remain.
After merging with regression SNP LD, 1183845 SNPs remain.
Removed 35 SNPs with chi^2 > 80 (1183810 SNPs remain)
Total Liability scale h2: 0.072 (0.0117)
Categories: L2_0 baseL2_1 Coding_UCSC.bedL2_1 Coding_UCSC.extend.500.bedL2_1 Conserved_LindbladToh.bedL2_1 Conserved_LindbladToh.extend.500.bedL2_1 CTCF_Hoffman.bedL2_1 CTCF_Hoffman.extend.500.bedL2_1 DGF_ENCODE.bedL2_1 DGF_ENCODE.extend.500.bedL2_1 DHS_peaks_Trynka.bedL2_1 DHS_Trynka.bedL2_1 DHS_Trynka.extend.500.bedL2_1 Enhancer_Andersson.bedL2_1 Enhancer_Andersson.extend.500.bedL2_1 Enhancer_Hoffman.bedL2_1 Enhancer_Hoffman.extend.500.bedL2_1 FetalDHS_Trynka.bedL2_1 FetalDHS_Trynka.extend.500.bedL2_1 H3K27ac_Hnisz.bedL2_1 H3K27ac_Hnisz.extend.500.bedL2_1 H3K27ac_PGC2.bedL2_1 H3K27ac_PGC2.extend.500.bedL2_1 H3K4me1_peaks_Trynka.bedL2_1 H3K4me1_Trynka.bedL2_1 H3K4me1_Trynka.extend.500.bedL2_1 H3K4me3_peaks_Trynka.bedL2_1 H3K4me3_Trynka.bedL2_1 H3K4me3_Trynka.extend.500.bedL2_1 H3K9ac_peaks_Trynka.bedL2_1 H3K9ac_Trynka.bedL2_1 H3K9ac_Trynka.extend.500.bedL2_1 Intron_UCSC.bedL2_1 Intron_UCSC.extend.500.bedL2_1 PromoterFlanking_Hoffman.bedL2_1 PromoterFlanking_Hoffman.extend.500.bedL2_1 Promoter_UCSC.bedL2_1 Promoter_UCSC.extend.500.bedL2_1 Repressed_Hoffman.bedL2_1 Repressed_Hoffman.extend.500.bedL2_1 SuperEnhancer_Hnisz.bedL2_1 SuperEnhancer_Hnisz.extend.500.bedL2_1 TFBS_ENCODE.bedL2_1 TFBS_ENCODE.extend.500.bedL2_1 Transcribed_Hoffman.bedL2_1 Transcribed_Hoffman.extend.500.bedL2_1 TSS_Hoffman.bedL2_1 TSS_Hoffman.extend.500.bedL2_1 UTR_3_UCSC.bedL2_1 UTR_3_UCSC.extend.500.bedL2_1 UTR_5_UCSC.bedL2_1 UTR_5_UCSC.extend.500.bedL2_1 WeakEnhancer_Hoffman.bedL2_1 WeakEnhancer_Hoffman.extend.500.bedL2_1
Lambda GC: 1.0895
Mean Chi^2: 1.1093
Intercept: 1.0311 (0.0121)
Ratio: 0.2847 (0.1107)
Reading annot matrix from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.5.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Results printed to AD_celltype_5.results
Analysis finished at Wed Apr  1 01:15:30 2026
Total time elapsed: 1.0m:39.4s

# group 6
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.6.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline. \
--out AD_celltype_6 \
--overlap-annot  \
--pop-prev 0.10 \
--samp-prev 0.344 \
--frqfile-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients  

Beginning analysis at Wed Apr  1 01:15:32 2026
Reading summary statistics from AD.sumstats.gz ...
Read summary statistics for 1207072 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.6.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Read reference panel LD Scores for 1190321 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.[1-22] ...
Read regression weight LD Scores for 1187349 SNPs.
After merging with reference panel LD, 1186816 SNPs remain.
After merging with regression SNP LD, 1183845 SNPs remain.
Removed 35 SNPs with chi^2 > 80 (1183810 SNPs remain)
Total Liability scale h2: 0.071 (0.0118)
Categories: L2_0 baseL2_1 Coding_UCSC.bedL2_1 Coding_UCSC.extend.500.bedL2_1 Conserved_LindbladToh.bedL2_1 Conserved_LindbladToh.extend.500.bedL2_1 CTCF_Hoffman.bedL2_1 CTCF_Hoffman.extend.500.bedL2_1 DGF_ENCODE.bedL2_1 DGF_ENCODE.extend.500.bedL2_1 DHS_peaks_Trynka.bedL2_1 DHS_Trynka.bedL2_1 DHS_Trynka.extend.500.bedL2_1 Enhancer_Andersson.bedL2_1 Enhancer_Andersson.extend.500.bedL2_1 Enhancer_Hoffman.bedL2_1 Enhancer_Hoffman.extend.500.bedL2_1 FetalDHS_Trynka.bedL2_1 FetalDHS_Trynka.extend.500.bedL2_1 H3K27ac_Hnisz.bedL2_1 H3K27ac_Hnisz.extend.500.bedL2_1 H3K27ac_PGC2.bedL2_1 H3K27ac_PGC2.extend.500.bedL2_1 H3K4me1_peaks_Trynka.bedL2_1 H3K4me1_Trynka.bedL2_1 H3K4me1_Trynka.extend.500.bedL2_1 H3K4me3_peaks_Trynka.bedL2_1 H3K4me3_Trynka.bedL2_1 H3K4me3_Trynka.extend.500.bedL2_1 H3K9ac_peaks_Trynka.bedL2_1 H3K9ac_Trynka.bedL2_1 H3K9ac_Trynka.extend.500.bedL2_1 Intron_UCSC.bedL2_1 Intron_UCSC.extend.500.bedL2_1 PromoterFlanking_Hoffman.bedL2_1 PromoterFlanking_Hoffman.extend.500.bedL2_1 Promoter_UCSC.bedL2_1 Promoter_UCSC.extend.500.bedL2_1 Repressed_Hoffman.bedL2_1 Repressed_Hoffman.extend.500.bedL2_1 SuperEnhancer_Hnisz.bedL2_1 SuperEnhancer_Hnisz.extend.500.bedL2_1 TFBS_ENCODE.bedL2_1 TFBS_ENCODE.extend.500.bedL2_1 Transcribed_Hoffman.bedL2_1 Transcribed_Hoffman.extend.500.bedL2_1 TSS_Hoffman.bedL2_1 TSS_Hoffman.extend.500.bedL2_1 UTR_3_UCSC.bedL2_1 UTR_3_UCSC.extend.500.bedL2_1 UTR_5_UCSC.bedL2_1 UTR_5_UCSC.extend.500.bedL2_1 WeakEnhancer_Hoffman.bedL2_1 WeakEnhancer_Hoffman.extend.500.bedL2_1
Lambda GC: 1.0895
Mean Chi^2: 1.1093
Intercept: 1.0318 (0.0123)
Ratio: 0.291 (0.1122)
Reading annot matrix from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.6.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Results printed to AD_celltype_6.results
Analysis finished at Wed Apr  1 01:17:13 2026
Total time elapsed: 1.0m:41.54s

# group 7
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.
7.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline. \
--out AD_celltype_7 \
--overlap-annot  \
--pop-prev 0.10 \
--samp-prev 0.344 \
--frqfile-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMH
C. \
--print-coefficients  
Beginning analysis at Wed Apr  1 01:21:18 2026
Reading summary statistics from AD.sumstats.gz ...
Read summary statistics for 1207072 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_typ
e_groups/cell_type_group.7.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/basel
ine.[1-22] ...
Read reference panel LD Scores for 1190321 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_weight
s_hm3_no_MHC/weights.hm3_noMHC.[1-22] ...
Read regression weight LD Scores for 1187349 SNPs.
After merging with reference panel LD, 1186816 SNPs remain.
After merging with regression SNP LD, 1183845 SNPs remain.
Removed 35 SNPs with chi^2 > 80 (1183810 SNPs remain)
Total Liability scale h2: 0.0718 (0.0117)
Categories: L2_0 baseL2_1 Coding_UCSC.bedL2_1 Coding_UCSC.extend.500.bedL2_1 Conserved_LindbladToh.bedL2_1 Conserved_LindbladToh.extend.500.bedL2_1 CTCF_Hoffman.bedL2_1 CTCF_Hoffman.extend.500.bedL2_1 DGF_ENCODE.bedL2_1 DGF_ENCODE.extend.500.bedL2_1 DHS_peaks_Trynka.bedL2_1 DHS_Trynka.bedL2_1 DHS_Trynka.extend.500.bedL2_1 Enhancer_Andersson.bedL2_1 Enhancer_Andersson.extend.500.bedL2_1 Enhancer_Hoffman.bedL2_1 Enhancer_Hoffman.extend.500.bedL2_1 FetalDHS_Trynka.bedL2_1 FetalDHS_Trynka.extend.500.bedL2_1 H3K27ac_Hnisz.bedL2_1 H3K27ac_Hnisz.extend.500.bedL2_1 H3K27ac_PGC2.bedL2_1 H3K27ac_PGC2.extend.500.bedL2_1 H3K4me1_peaks_Trynka.bedL2_1 H3K4me1_Trynka.bedL2_1 H3K4me1_Trynka.extend.500.bedL2_1 H3K4me3_peaks_Trynka.bedL2_1 H3K4me3_Trynka.bedL2_1 H3K4me3_Trynka.extend.500.bedL2_1 H3K9ac_peaks_Trynka.bedL2_1 H3K9ac_Trynka.bedL2_1 H3K9ac_Trynka.extend.500.bedL2_1 Intron_UCSC.bedL2_1 Intron_UCSC.extend.500.bedL2_1 PromoterFlanking_Hoffman.bedL2_1 PromoterFlanking_Hoffman.extend.500.bedL2_1 Promoter_UCSC.bedL2_1 Promoter_UCSC.extend.500.bedL2_1 Repressed_Hoffman.bedL2_1 Repressed_Hoffman.extend.500.bedL2_1 SuperEnhancer_Hnisz.bedL2_1 SuperEnhancer_Hnisz.extend.500.bedL2_1 TFBS_ENCODE.bedL2_1 TFBS_ENCODE.extend.500.bedL2_1 Transcribed_Hoffman.bedL2_1 Transcribed_Hoffman.extend.500.bedL2_1 TSS_Hoffman.bedL2_1 TSS_Hoffman.extend.500.bedL2_1 UTR_3_UCSC.bedL2_1 UTR_3_UCSC.extend.500.bedL2_1 UTR_5_UCSC.bedL2_1 UTR_5_UCSC.extend.500.bedL2_1 WeakEnhancer_Hoffman.bedL2_1 WeakEnhancer_Hoffman.extend.500.bedL2_1
Lambda GC: 1.0895
Mean Chi^2: 1.1093
Intercept: 1.0312 (0.0121)
Ratio: 0.2855 (0.1106)
Reading annot matrix from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.7.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Results printed to AD_celltype_7.results
Analysis finished at Wed Apr  1 01:22:58 2026
Total time elapsed: 1.0m:39.74s

# group 8
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.
8.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline. \
--out AD_celltype_8 \
--overlap-annot  \
--pop-prev 0.10 \
--samp-prev 0.344 \
--frqfile-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMH
C. \
--print-coefficients  
Beginning analysis at Wed Apr  1 01:22:59 2026
Reading summary statistics from AD.sumstats.gz ...
Read summary statistics for 1207072 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_typ
e_groups/cell_type_group.8.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/basel
ine.[1-22] ...
Read reference panel LD Scores for 1190321 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_weight
s_hm3_no_MHC/weights.hm3_noMHC.[1-22] ...
Read regression weight LD Scores for 1187349 SNPs.
After merging with reference panel LD, 1186816 SNPs remain.
After merging with regression SNP LD, 1183845 SNPs remain.
Removed 35 SNPs with chi^2 > 80 (1183810 SNPs remain)
Total Liability scale h2: 0.0717 (0.0118)
Categories: L2_0 baseL2_1 Coding_UCSC.bedL2_1 Coding_UCSC.extend.500.bedL2_1 Conserved_LindbladToh.bedL2_1 Conserved_LindbladToh.extend.500.bedL2_1 CTCF_Hoffman.bedL2_1 CTCF_Hoffman.extend.500.bedL2_1 DGF_ENCODE.bedL2_1 DGF_ENCODE.extend.500.bedL2_1 DHS_peaks_Trynka.bedL2_1 DHS_Trynka.bedL2_1 DHS_Trynka.extend.500.bedL2_1 Enhancer_Andersson.bedL2_1 Enhancer_Andersson.extend.500.bedL2_1 Enhancer_Hoffman.bedL2_1 Enhancer_Hoffman.extend.500.bedL2_1 FetalDHS_Trynka.bedL2_1 FetalDHS_Trynka.extend.500.bedL2_1 H3K27ac_Hnisz.bedL2_1 H3K27ac_Hnisz.extend.500.bedL2_1 H3K27ac_PGC2.bedL2_1 H3K27ac_PGC2.extend.500.bedL2_1 H3K4me1_peaks_Trynka.bedL2_1 H3K4me1_Trynka.bedL2_1 H3K4me1_Trynka.extend.500.bedL2_1 H3K4me3_peaks_Trynka.bedL2_1 H3K4me3_Trynka.bedL2_1 H3K4me3_Trynka.extend.500.bedL2_1 H3K9ac_peaks_Trynka.bedL2_1 H3K9ac_Trynka.bedL2_1 H3K9ac_Trynka.extend.500.bedL2_1 Intron_UCSC.bedL2_1 Intron_UCSC.extend.500.bedL2_1 PromoterFlanking_Hoffman.bedL2_1 PromoterFlanking_Hoffman.extend.500.bedL2_1 Promoter_UCSC.bedL2_1 Promoter_UCSC.extend.500.bedL2_1 Repressed_Hoffman.bedL2_1 Repressed_Hoffman.extend.500.bedL2_1 SuperEnhancer_Hnisz.bedL2_1 SuperEnhancer_Hnisz.extend.500.bedL2_1 TFBS_ENCODE.bedL2_1 TFBS_ENCODE.extend.500.bedL2_1 Transcribed_Hoffman.bedL2_1 Transcribed_Hoffman.extend.500.bedL2_1 TSS_Hoffman.bedL2_1 TSS_Hoffman.extend.500.bedL2_1 UTR_3_UCSC.bedL2_1 UTR_3_UCSC.extend.500.bedL2_1 UTR_5_UCSC.bedL2_1 UTR_5_UCSC.extend.500.bedL2_1 WeakEnhancer_Hoffman.bedL2_1 WeakEnhancer_Hoffman.extend.500.bedL2_1
Lambda GC: 1.0895
Mean Chi^2: 1.1093
Intercept: 1.0312 (0.0121)
Ratio: 0.2853 (0.1106)
Reading annot matrix from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.8.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Results printed to AD_celltype_8.results
Analysis finished at Wed Apr  1 01:24:37 2026
Total time elapsed: 1.0m:37.82s

# group 9
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.9.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline. \
--out AD_celltype_9 \
--overlap-annot  \
--pop-prev 0.10 \
--samp-prev 0.344 \
--frqfile-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients  

Beginning analysis at Wed Apr  1 01:24:38 2026
Reading summary statistics from AD.sumstats.gz ...
Read summary statistics for 1207072 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.9.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Read reference panel LD Scores for 1190321 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.[1-22] ...
Read regression weight LD Scores for 1187349 SNPs.
After merging with reference panel LD, 1186816 SNPs remain.
After merging with regression SNP LD, 1183845 SNPs remain.
Removed 35 SNPs with chi^2 > 80 (1183810 SNPs remain)
Total Liability scale h2: 0.0715 (0.0117)
Categories: L2_0 baseL2_1 Coding_UCSC.bedL2_1 Coding_UCSC.extend.500.bedL2_1 Conserved_LindbladToh.bedL2_1 Conserved_LindbladToh.extend.500.bedL2_1 CTCF_Hoffman.bedL2_1 CTCF_Hoffman.extend.500.bedL2_1 DGF_ENCODE.bedL2_1 DGF_ENCODE.extend.500.bedL2_1 DHS_peaks_Trynka.bedL2_1 DHS_Trynka.bedL2_1 DHS_Trynka.extend.500.bedL2_1 Enhancer_Andersson.bedL2_1 Enhancer_Andersson.extend.500.bedL2_1 Enhancer_Hoffman.bedL2_1 Enhancer_Hoffman.extend.500.bedL2_1 FetalDHS_Trynka.bedL2_1 FetalDHS_Trynka.extend.500.bedL2_1 H3K27ac_Hnisz.bedL2_1 H3K27ac_Hnisz.extend.500.bedL2_1 H3K27ac_PGC2.bedL2_1 H3K27ac_PGC2.extend.500.bedL2_1 H3K4me1_peaks_Trynka.bedL2_1 H3K4me1_Trynka.bedL2_1 H3K4me1_Trynka.extend.500.bedL2_1 H3K4me3_peaks_Trynka.bedL2_1 H3K4me3_Trynka.bedL2_1 H3K4me3_Trynka.extend.500.bedL2_1 H3K9ac_peaks_Trynka.bedL2_1 H3K9ac_Trynka.bedL2_1 H3K9ac_Trynka.extend.500.bedL2_1 Intron_UCSC.bedL2_1 Intron_UCSC.extend.500.bedL2_1 PromoterFlanking_Hoffman.bedL2_1 PromoterFlanking_Hoffman.extend.500.bedL2_1 Promoter_UCSC.bedL2_1 Promoter_UCSC.extend.500.bedL2_1 Repressed_Hoffman.bedL2_1 Repressed_Hoffman.extend.500.bedL2_1 SuperEnhancer_Hnisz.bedL2_1 SuperEnhancer_Hnisz.extend.500.bedL2_1 TFBS_ENCODE.bedL2_1 TFBS_ENCODE.extend.500.bedL2_1 Transcribed_Hoffman.bedL2_1 Transcribed_Hoffman.extend.500.bedL2_1 TSS_Hoffman.bedL2_1 TSS_Hoffman.extend.500.bedL2_1 UTR_3_UCSC.bedL2_1 UTR_3_UCSC.extend.500.bedL2_1 UTR_5_UCSC.bedL2_1 UTR_5_UCSC.extend.500.bedL2_1 WeakEnhancer_Hoffman.bedL2_1 WeakEnhancer_Hoffman.extend.500.bedL2_1
Lambda GC: 1.0895
Mean Chi^2: 1.1093
Intercept: 1.0314 (0.0121)
Ratio: 0.287 (0.111)
Reading annot matrix from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.9.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Results printed to AD_celltype_9.results
Analysis finished at Wed Apr  1 01:26:19 2026
Total time elapsed: 1.0m:40.52s

# group 10
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.10.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline. \
--out AD_celltype_10 \
--overlap-annot  \
--pop-prev 0.10 \
--samp-prev 0.344 \
--frqfile-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients  

Beginning analysis at Wed Apr  1 01:26:20 2026
Reading summary statistics from AD.sumstats.gz ...
Read summary statistics for 1207072 SNPs.
Reading reference panel LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.10.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Read reference panel LD Scores for 1190321 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from /projectnb/bs859/data/ldscore_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.[1-22] ...
Read regression weight LD Scores for 1187349 SNPs.
After merging with reference panel LD, 1186816 SNPs remain.
After merging with regression SNP LD, 1183845 SNPs remain.
Removed 35 SNPs with chi^2 > 80 (1183810 SNPs remain)
Total Liability scale h2: 0.0724 (0.0118)
Categories: L2_0 baseL2_1 Coding_UCSC.bedL2_1 Coding_UCSC.extend.500.bedL2_1 Conserved_LindbladToh.bedL2_1 Conserved_LindbladToh.extend.500.bedL2_1 CTCF_Hoffman.bedL2_1 CTCF_Hoffman.extend.500.bedL2_1 DGF_ENCODE.bedL2_1 DGF_ENCODE.extend.500.bedL2_1 DHS_peaks_Trynka.bedL2_1 DHS_Trynka.bedL2_1 DHS_Trynka.extend.500.bedL2_1 Enhancer_Andersson.bedL2_1 Enhancer_Andersson.extend.500.bedL2_1 Enhancer_Hoffman.bedL2_1 Enhancer_Hoffman.extend.500.bedL2_1 FetalDHS_Trynka.bedL2_1 FetalDHS_Trynka.extend.500.bedL2_1 H3K27ac_Hnisz.bedL2_1 H3K27ac_Hnisz.extend.500.bedL2_1 H3K27ac_PGC2.bedL2_1 H3K27ac_PGC2.extend.500.bedL2_1 H3K4me1_peaks_Trynka.bedL2_1 H3K4me1_Trynka.bedL2_1 H3K4me1_Trynka.extend.500.bedL2_1 H3K4me3_peaks_Trynka.bedL2_1 H3K4me3_Trynka.bedL2_1 H3K4me3_Trynka.extend.500.bedL2_1 H3K9ac_peaks_Trynka.bedL2_1 H3K9ac_Trynka.bedL2_1 H3K9ac_Trynka.extend.500.bedL2_1 Intron_UCSC.bedL2_1 Intron_UCSC.extend.500.bedL2_1 PromoterFlanking_Hoffman.bedL2_1 PromoterFlanking_Hoffman.extend.500.bedL2_1 Promoter_UCSC.bedL2_1 Promoter_UCSC.extend.500.bedL2_1 Repressed_Hoffman.bedL2_1 Repressed_Hoffman.extend.500.bedL2_1 SuperEnhancer_Hnisz.bedL2_1 SuperEnhancer_Hnisz.extend.500.bedL2_1 TFBS_ENCODE.bedL2_1 TFBS_ENCODE.extend.500.bedL2_1 Transcribed_Hoffman.bedL2_1 Transcribed_Hoffman.extend.500.bedL2_1 TSS_Hoffman.bedL2_1 TSS_Hoffman.extend.500.bedL2_1 UTR_3_UCSC.bedL2_1 UTR_3_UCSC.extend.500.bedL2_1 UTR_5_UCSC.bedL2_1 UTR_5_UCSC.extend.500.bedL2_1 WeakEnhancer_Hoffman.bedL2_1 WeakEnhancer_Hoffman.extend.500.bedL2_1
Lambda GC: 1.0895
Mean Chi^2: 1.1093
Intercept: 1.0309 (0.0121)
Ratio: 0.2827 (0.111)
Reading annot matrix from /projectnb/bs859/data/ldscore_files/1000G_Phase3_cell_type_groups/cell_type_group.10.,/projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.[1-22] ...
Results printed to AD_celltype_10.results
Analysis finished at Wed Apr  1 01:28:00 2026
Total time elapsed: 1.0m:40.04s

# print out the enrichment scores, std, and p-values
for N in {1..10}; do
  echo -n "Group $N: "
  awk '$1=="L2_0" {print $5, $6, $7}' AD_celltype_${N}.results
done
Group 1: 5.73752919162 1.70698512449 0.0021881685758397842
Group 2: 2.47683904836 1.25954237272 0.26250189310893168
Group 3: 4.45221852312 0.957346487366 0.00020583325421419109
Group 4: 5.72367480954 1.45249272873 4.8374550100772444e-05
Group 5: 4.10062787731 0.845104692174 0.00066623764205058823
Group 6: 4.96955736469 0.853198334309 2.8765959145279837e-07
Group 7: 6.92913011059 2.45754590367 0.0099291283932174731
Group 8: 9.96176003472 2.60554633952 1.7871904506649469e-05
Group 9: 4.85626047693 1.00445465504 2.5687901929833607e-05
Group 10: 5.67939535601 1.35565187875 0.00032979073414727102
```