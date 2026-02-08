## Homework 3: Common Variant GWAS - Association Analysis & Mixed Models
### BS859 Applied Genetic Analysis
### Addison Yam
### February 11, 2026

We will use the TGEN Alzheimer data located in /projectnb/bs859/data/tgen/cleaned/ to do this assignment. For PC covariates, use the file use the file TGEN_pcs.txt in the same directory. In your submitted homework, include the commands you used to get your results, and the parts of the results that support your response. Use tables and/or plots to present results clearly and succinctly.

```bash
# load the necessary modules
module load R 
module load plink/1.90b6.27
#make an alias for the directory with the data
export DATADIR=/projectnb/bs859/data/tgen/cleaned
```

1. Use the file TGEN_pcs.txt to determine which PCs are associated with case status. (TGEN_pcs.txt is the output from smartpca with a column header added). **HINT: “--logistic no-snp” produces results for the specified covariates, with no SNPs in the model.**

a. Show the commands you used for the analysis.
```bash
> head $DATADIR/TGEN_pcs.txt
FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
MAYO_10139 MAYO_10139 -0.0118 -0.0040 0.0100 -0.0080 0.0201 0.0615 -0.0219 -0.0301 0.0293 0.0737
MAYO_10198 MAYO_10198 -0.0083 -0.0074 0.0133 -0.0025 0.0624 -0.0283 -0.0572 -0.0150 0.0341 -0.0103
MAYO_102246 MAYO_102246 -0.0053 0.1093 -0.0200 0.0049 0.0073 0.0272 -0.0010 0.0474 -0.1290 0.0441
MAYO_10249 MAYO_10249 -0.0156 -0.0078 0.0181 0.0116 -0.0007 -0.0200 0.0069 0.0623 -0.0060 0.0203
MAYO_10278 MAYO_10278 -0.0127 0.0043 0.0139 -0.0173 0.0274 0.1298 0.0159 0.0068 0.0140 0.0518
MAYO_10304 MAYO_10304 -0.0049 -0.0004 0.0036 0.0040 0.0150 -0.0051 0.0108 0.0083 -0.0375 0.0445
MAYO_10316 MAYO_10316 -0.0113 -0.0052 0.0020 -0.0089 -0.0102 0.0184 0.0040 -0.0451 0.0274 0.0523
MAYO_10367 MAYO_10367 -0.0128 -0.0125 -0.0011 -0.0045 -0.0041 0.0003 -0.0079 -0.0297 -0.0026 -0.0026
MAYO_10403 MAYO_10403 -0.0118 -0.0051 -0.0075 -0.0123 0.0037 0.0092 -0.0105 0.0013 0.0180 0.0276

# Test to determine which PCs are associated with case status
# --bifle   input binary genotype files
# --covar   covariate file about PCs
# --logistic no snp     Logistic regression testing the covariates and no SNPs
> plink --bfile $DATADIR/TGEN_cleaned --covar $DATADIR/TGEN_pcs.txt --logistic no-snp --allow-no-sex --out PC_case_test
PLINK v1.90b6.27 64-bit (10 Dec 2022)          www.cog-genomics.org/plink/1.9/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to PC_case_test.log.
Options in effect:
  --allow-no-sex
  --bfile /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned
  --covar /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt
  --logistic no-snp
  --out PC_case_test

257325 MB RAM detected; reserving 128662 MB for main workspace.
308462 variants loaded from .bim file.
1237 people (0 males, 0 females, 1237 ambiguous) loaded from .fam.
Ambiguous sex IDs written to PC_case_test.nosex .
1237 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
--covar: 10 covariates loaded.
Before main variant filters, 1237 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.988302.
308462 variants and 1237 people pass filters and QC.
Among remaining phenotypes, 769 are cases and 468 are controls.
Writing logistic model association results to PC_case_test.assoc.logistic ...
done.

# View the contents of the association test
> cat PC_case_test.assoc.logistic
      TEST    NMISS         OR         STAT            P 
       PC1     1237     0.3603      -0.4574       0.6474
       PC2     1237      788.7        2.498      0.01249
       PC3     1237   0.003487       -2.023      0.04311
       PC4     1237      5.832        0.622        0.534
       PC5     1237     0.2724      -0.5759       0.5647
       PC6     1237   0.001533       -2.901     0.003718
       PC7     1237      125.1        2.226      0.02604
       PC8     1237    0.00382       -2.611     0.009026
       PC9     1237     0.2916      -0.5751       0.5652
      PC10     1237     0.1353      -0.9387       0.3479
```

b. Where is plink looking to find case status for your analysis?
- Answer: In the sixth column of the`TGEN_cleaned.fam` file is the Phenotype column where plink is looking for the case status this analysis. Values of 1 represent unaffected and 2 present affected. 

c. How many individuals are included in your analysis?
- Answer: There are 1237 individuals, which I got from the NMISS column in the `PC_case_test.assoc.logistic` file. 

d. Which PCs are associated with case status at p<0.01?
- Answer: PC6 and PC8 are associated with a case status at <0.01

2. Perform a GWAS of Alzheimer Disease case status using logistic regression, adjusting for the PCs associated with case status at p<0.01 that you found in question 1. *Use a plink command that hides the covariates so that your plink output has only the association results for the SNPs, and that gives the effect estimates as regression coefficients rather than odds ratios.*

a. Show your plink command
```bash
# Perform GWAS and adjusting for PCs (PC6, PC8) with a p-value of <0.01
# -- covar-name    Use PC6 and PC8 as covariates
# hide-covar         Only print the results for SNPs, not covariates
# --logistic beta   Gets the regression coefficients
# --ci 0.95           95% confidence internval and standard errors
> plink --bfile $DATADIR/TGEN_cleaned --covar $DATADIR/TGEN_pcs.txt --covar-name PC6, PC8 --logistic beta hide-covar --allow-no-sex --ci 0.95 --out GWAS_PC68
PLINK v1.90b6.27 64-bit (10 Dec 2022)          www.cog-genomics.org/plink/1.9/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to GWAS_PC68.log.
Options in effect:
  --allow-no-sex
  --bfile /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned
  --ci 0.95
  --covar /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt
  --covar-name PC6, PC8
  --logistic beta hide-covar
  --out GWAS_PC68

257325 MB RAM detected; reserving 128662 MB for main workspace.
308462 variants loaded from .bim file.
1237 people (0 males, 0 females, 1237 ambiguous) loaded from .fam.
Ambiguous sex IDs written to GWAS_PC68.nosex .
1237 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
--covar: 2 out of 10 covariates loaded.
Before main variant filters, 1237 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.988302.
308462 variants and 1237 people pass filters and QC.
Among remaining phenotypes, 769 are cases and 468 are controls.
Writing logistic model association results to GWAS_PC68.assoc.logistic ...
done.

> head -6 GWAS_PC68.assoc.logistic
 CHR         SNP         BP   A1       TEST    NMISS       BETA       SE      L95      U95         STAT            P 
   1   rs3094315     742429    C        ADD     1233     0.0525   0.1133  -0.1696   0.2746       0.4633       0.6431
   1   rs4040617     769185    G        ADD     1237    0.07299   0.1211  -0.1644   0.3104       0.6025       0.5468
   1   rs2980300     775852    T        ADD     1182   0.003232   0.1159  -0.2239   0.2304      0.02789       0.9778
   1   rs4075116     993492    G        ADD     1228    -0.0843  0.09289  -0.2664  0.09776      -0.9076       0.3641
   1  rs10907175    1120590    C        ADD     1227    0.01932   0.1511  -0.2769   0.3156       0.1279       0.8983
   1   rs6603781    1148494    A        ADD     1196    0.06635   0.1269  -0.1824   0.3151       0.5229        0.601
   1  rs11260562    1155173    A        ADD     1234    0.05847   0.1935  -0.3207   0.4376       0.3022       0.7625
   1   rs6685064    1201155    T        ADD     1229   -0.06814   0.1747  -0.4105   0.2742      -0.3901       0.6965
   1   rs3766180    1468016    G        ADD     1236    0.01147  0.09626  -0.1772   0.2001       0.1191       0.9052

# Gets a count of the SNPs with a p-value > 0.0001
> awk '$12 < 0.01 { count++ } END { print count+0 }' GWAS_PC68.assoc.logistic
33
```

b. Show the results for the first 5 SNPs in the output file.
```bash
> head -6 GWAS_PC68.assoc.logistic
 CHR         SNP         BP   A1       TEST    NMISS       BETA       SE      L95      U95         STAT            P 
   1   rs3094315     742429    C        ADD     1233     0.0525   0.1133  -0.1696   0.2746       0.4633       0.6431
   1   rs4040617     769185    G        ADD     1237    0.07299   0.1211  -0.1644   0.3104       0.6025       0.5468
   1   rs2980300     775852    T        ADD     1182   0.003232   0.1159  -0.2239   0.2304      0.02789       0.9778
   1   rs4075116     993492    G        ADD     1228    -0.0843  0.09289  -0.2664  0.09776      -0.9076       0.3641
   1  rs10907175    1120590    C        ADD     1227    0.01932   0.1511  -0.2769   0.3156       0.1279       0.8983
```
c. How many SNPs in this GWAS have p-value < 0.0001?
- Answer: There are 33 SNPs with a p-value less than 0.0001

d. Show the plink results for the most significant SNP. Show the 3 genotypes for the SNP, and how they are coded in the plink logistic regression. Compute the odds ratio. Which allele increases the risk of Alzheimer Disease?
- Answer: The most significant SNP is rs41377151. The genotypes for the rs41377151 SNP are AA, GA, GG and they're coded as 0, 1, 2. The β for this SNP is 1.174, using this value I can calculate the odds ratio. OR = e<sup>1.174</sup> = 3.2349. A positive β and OR > 1 tells us that that allele that increases the risk of AD is the G allele.
```bash
# Skip the header, sort the p-values in the 12th column by least to greatest, show the SNP with the most significant p-value
> awk 'NR>1' GWAS_PC68.assoc.logistic | sort -k12,12g | head -1
  19  rs41377151   50114786    G        ADD     1237      1.174   0.1069   0.9643    1.384        10.98    4.869e-28

# Shows the PLINK results for the most significant SNP
> awk 'NR==1 || $2=="rs41377151" {print $0}' GWAS_PC68.assoc.logistic
 CHR         SNP         BP   A1       TEST    NMISS       BETA       SE      L95      U95         STAT            P 
  19  rs41377151   50114786    G        ADD     1237      1.174   0.1069   0.9643    1.384        10.98    4.869e-28

# Find the SNP in the .bim 
grep "rs41377151" $DATADIR/TGEN_cleaned.bim
19	rs41377151	0	50114786	G	A
```

3. Create a GRM from the TGEN data so that you can perform a GWAS using a logistic mixed model. Use only chromosomes 1-22 (exclude X). Don’t forget to prune before you make the GRM! Use pruning parameters: --indep-pairwise 10000kb 1 0.15.

a. Show your code.

b. How many SNPs were left after pruning and removing the chromosome X variants?

c. Provide the first 10 lines of the id grm id file (the file with the ids for the grm)

d. What is the maximum relationship coefficient for MAYO_10139, other than with itself?

4. Perform a GWAS using a logistic mixed model score test, with no covariates. Do a second GWAS using this model and the same PCs you included in the Q2 analysis. You can edit the GMMAT.R program we used in class to do this.

a. How many SNPs have p-value < 0.0001 in each of these two GWAS?

b. Compare the results of the two GLMMs and the logistic model from question 2 for the SNP that is most significantly associated with AD case status in the question 2d analysis. Which analysis yields the most significant association? Can you compare the effect estimates? Which allele increases risk in these models?

5. Produce QQ plots for the three GWAS analyses, and show them here. Looking at the Genomic
Control Lambda and the overall qq plot, which analysis model do you think is best? Explain your
choice.

6. Produce a Manhattan plot for the GWAS that you chose in 5. Based on this plot, are there any regions of the genome OTHER than chromosome 19, that may be of interest? Why or why not?