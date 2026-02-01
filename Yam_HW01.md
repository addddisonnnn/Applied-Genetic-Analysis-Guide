## Homework 1: Genotype Data QC & PLINK
### BS859 Applied Genetic Analysis
### January 21, 2026

```bash
# Set the working directory and load PLINK
export DATADIR=/projectnb/bs859/data/tgen
module load plink/1.90b6.27 
```

1. **How many individuals are in the TGEN dataset? How many variants?** 
- Answer: There are 1411 individuals and 312,315 variants in the TGEN dataset.
```bash
# Check the number of lines in the TGEN.fam file, the first number is the number of individuals
wc $DATADIR/TGEN.fam
1411  8466 40704 /projectnb/bs859/data/tgen/TGEN.fam
# Count the number of lines in the .bim file
wc $DATADIR/TGEN.bim
 312316 1873896 8692055 /projectnb/bs859/data/tgen/TGEN.bim
```

2. **Create a summary file of allele frequencies for all variants. How many variants have no minor (A1) allele homozygotes? There are 9809 variants.**
- Answer: There are 9809 variants.

```bash
# Generate allele frequency file
plink --bfile $DATADIR/TGEN --freq --out freq1   
# Get the genotype  counts
plink --bfile $DATADIR/TGEN --freqx --out freqx

# Check the header
head TGEN_freqx.frqx
CHR	SNP	A1	A2	C(HOM A1)	C(HET)	C(HOM A2)	C(HAP A1)	C(HAP A2)	C(MISSING)
0	rs41408748	A	C	256	709	403	0	0	43
0	rs9678914	A	C	276	673	459	0	0	3
0	rs1246399	A	C	1	83	1314	0	0	13
0	rs28364983	C	A	45	382	956	0	0	28
0	rs4603449	C	A	52	442	903	0	0	14
0	rs2255803	C	A	147	609	608	0	0	47
0	rs2668622	A	C	51	464	887	0	0	9
0	rs10900309	G	A	0	459	936	0	0	16
0	rs12402966	A	G	67	461	883	0	0	0

# Count how many ‘o’s are in the 5th column (C(HOM A1)) as it pertains to minor allele homozyogte count
awk 'NR>1 && $5==0 {count++} END {print count}' freqx.frqx
# NR>1		Skips the header row
# $5==0 	Values in 5th column equal to 0
9809
```

3. **Creare a filtered data set where all variants have mino allele frequencies >0.01, fewer than 5% missing genotype, and the HWE p-valur in controls is greater than 0.000001 and where all individuals have fewer than 5% missing genotypes using a single PLINK command.**

    a. **How many individuals were removed from the dataset, and why?** 
    - Answer: 47 individuals were removed because they were missing >5% genotypes, this was done through the --mind 0.05 parameter
    
    b. **How many variants are removed by the HWE filter?**
    - Answer: 884
    
    c. **How many variants are removed by the missing genotype filter?**
    - Answer: 1364
    
    d. **How many variants are removed by the MAF filter?**
    - Answer: 1353
    
    e. **What was the order in which PLINK applied the filters?**
    - Answer: Individual missing rate (--mind), variant missing rate (--geno), HWE filter (--hwe), and MAF filter (--maf)

```bash
# Single PLINK command creates a filtered data set
# --maf 0.01 		minor allele freq > 0.01
# --geno 0.05 		variant missing rate < 0.05
# --mind 0.05		individual missing rate < 0.05
# --hwe 0.000001	HWE p-value > 0.000001 in controls
# --make-bed		create new binary files
# --out TGEN_filtered	output file prefix
plink --bfile $DATADIR/TGEN --maf 0.01 --geno 0.05 --mind 0.05 --hwe 0.000001 --make-bed --out TGEN_filtered

# look at the contents of the log file
cat TGEN_filtered.log
PLINK v1.90b6.27 64-bit (10 Dec 2022)
Options in effect:
  --bfile /projectnb/bs859/data/tgen/TGEN
  --geno 0.05
  --hwe 0.000001
  --maf 0.01
  --make-bed
  --mind 0.05
  --out TGEN_filtered

Hostname: scc-wl3
Working directory: /projectnb/bs859/students/addisony
Start time: Tue Jan 27 19:08:28 2026

Random number seed: 1769558908
257325 MB RAM detected; reserving 128662 MB for main workspace.
312316 variants loaded from .bim file.
1411 people (0 males, 0 females, 1411 ambiguous) loaded from .fam.
Ambiguous sex IDs written to TGEN_filtered.nosex .
1411 phenotype values loaded from .fam.
Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
phenotypes to be ignored, use the --allow-no-sex flag.
47 people removed due to missing genotype data (--mind).
IDs written to TGEN_filtered.irem .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 1364 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate in remaining samples is 0.987367.
1384 variants removed due to missing genotype data (--geno).
--hwe: 884 variants removed due to Hardy-Weinberg exact test.
1353 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
308695 variants and 1364 people pass filters and QC.
Among remaining phenotypes, 829 are cases and 535 are controls.
--make-bed to TGEN_filtered.bed + TGEN_filtered.bim + TGEN_filtered.fam ...
done.

End time: Tue Jan 27 19:08:28 2026
```

4. **Create a second filtered data set where you first remove SNPs with MAF<0.01 and missing rate>0.05 and HWE p<0.000001 all in one step, and THEN you remove individuals with >5% missing genotypes.**

    a. **How many SNPs are removed by the MAF filter?**
    - Answer: 1215

    b. **How many SNPs are removed by the missing genotype filter?**
    - Answer: 3175
    
    c. **How many SNPs are removed by the HWE filter?**
    - Answer: 895
    
    d. **How many individuals were removed from analysis?** Answer: 42
    
    e. **Which individuals (provide the IDs) were excluded in 1) but included in 2)? (hint: you can compare the individuals in the plink “irem” files - you can use the linux command “diff”, or SAS or R or another method of your choice to find the difference)** 
    - Answer: SHRI_BBBB, WGAAD_134, WGAAD_164, WGACON_144, WGACON_224

5. **How does the data set in part 3 differ from the dataset created in part 4? Are there any advantages or disadvantages to using the default order of operations in PLINK when creating the cleaned dataset?**
    - Answer: The data set in question 3 used the default order of operations in PLINK where it removed individuals first. This led to a final dataset of 308,695 variants and 1,364 individuals by removing 47 individuals and 3,621 variants. For question 4, it performed SNP filtering first, this led to a final dataset of 307,032 variants and 1,369 individuals by removing 42 individuals and 5,285 variants. The main difference seen is, question 4 had kept 5 mote individuals and question 3 kept 1664 more variants. The advantages of using the default PLINK order of operations are it keeps more variants and removes individuals with high missing rates across the genome. The disadvantages of using the PLINK default order of operations are decreased sample size and it removes appropriate individuals that would’ve been kept if we first remove SNPs. 
6. **Perform pairwise LD pruning (using the parameters 10000 kb window size, 1 variant shift, maximum r2=0.2), and remove the pruned variants. How many variants remain after pruning?** 
    - Answer: 79143 variants remain after pruning.
```bash
# First, find SNPs in LD through the prune.in and prune.out files
# 10000 kb (10000 kb window size) 1 (shift by 1 variant) 0.2 (remove SNP from pairs where r^2 > 0.2)
plink --bfile TGEN_filtered --indep-pairwise 10000 kb 1 0.2 --out TGEN_prune

# Second, make dataset of just the prune-in SNPs
plink --bfile TGEN_filtered --extract TGEN_prune.prune.in --make-bed --out TGEN_pruned

# Count the number of lines in .bim file
wc -l TGEN_pruned.bim
79143 TGEN_pruned.bim
```

7. **Compute pairwise IBD estimates for the filtered data set. When you do this, use only chromosomes 1-22 (there are some SNPs with chromosome=0, which are not included in the pruning and should not be trusted)**
```bash
# First, get the pruned SNPs of chromosomes 1-22
plink --bfile TGEN_filtered --chr 1-22 --extract TGEN_prune.prune.in --make-bed --out TGEN_chr1-22_pruned
# Compute pairwise IBD estimates
plink --bfile TGEN_chr1-22_pruned --genome --out TGEN_ibd
```
a. **How many variants remained after pruning?**
- Answer: 79010 variants remain after pruning
```bash
# Count the amount of variants by number of lines
wc -l TGEN_chr1-22_pruned.bim
79010 TGEN_chr1-22_pruned.bim
```

b. **How many pairs share more than 25%? (PI_HAT)**
- Answer: 107 pairs share more than 25%.
```bash
# Skip header and count the lines where the value in the 10th column called PI_HAT are greater than 0.25
awk 'NR>1 && $10>0.25 {count++} END {print count}' TGEN_ibd.genome
107
```
c. **How many unique individuals share more than 25% of alleles IBD with at least one other person? (some individuals are part of multiple pairs)**
- Answer: 127 individuals.
```bash
# Get unique individuals
awk 'NR>1 && $10>0.25 {print $1"\t"$2; print $3"\t"$4}' TGEN_ibd.genome | sort -u | wc -l
127
```
b. **How many of the pairs in b) appear to be genetic duplicates? (Use PI_HAT>0.95 to identify pairs with identical genotypes).**
- Answer: There are 37 pairs. 
```bash
# Get pairs that are genetic duplicates
awk 'NR>1 && $10>0.95 {count++} END {print count}' TGEN_ibd.genome
37
```
e. **Summarize the non-genetic-duplicate relationships you see in the file that are 2nd degree relative or closer (pairs with PI_HAT<0.95 and PI_HAT>0.25). What type of relationship do most of these pairs appear to share? (Compare the estimated proportion of 0, 1, 2 IBD with expected for various relationships).**
- Answer: Summarize the non-genetic-duplicate relationships you see in the file that are 2nd degree relative or closer (pairs with PI_HAT<0.95 and PI_HAT>0.25). What type of relationship do most of these pairs appear to share? (Compare the estimated proportion of 0, 1, 2 IBD with expected for various relationships).
```bash
# Make new file of relatives with 2nd degree or closer
awk 'NR>1 && $10>0.25 && $10<0.95 {print $0}' TGEN_ibd.genome > TGEN_relatives.txt
# check the Z0, Z1, and Z2 in columns 7-9
cat TGEN_relatives.txt
MAYO_10759  MAYO_10759  MAYO_16111  MAYO_16111 UN    NA  0.4903  0.5097  0.0000  0.2549  -1  0.780996  1.0000  4.7374
  MAYO_16111  MAYO_16111  MAYO_16221  MAYO_16221 UN    NA  0.3903  0.6097  0.0000  0.3049  -1  0.790634  1.0000  7.9783
  MAYO_16111  MAYO_16111  MAYO_17623  MAYO_17623 UN    NA  0.4779  0.5221  0.0000  0.2611  -1  0.771620  1.0000  5.1262
  MAYO_16111  MAYO_16111  MAYO_22544  MAYO_22544 UN    NA  0.3656  0.6344  0.0000  0.3172  -1  0.788831  1.0000  6.9454
  MAYO_16111  MAYO_16111   MAYO_2857   MAYO_2857 UN    NA  0.4021  0.5979  0.0000  0.2989  -1  0.784992  1.0000  5.6007
  MAYO_16111  MAYO_16111  MAYO_32401  MAYO_32401 UN    NA  0.4736  0.5264  0.0000  0.2632  -1  0.761921  1.0000  4.7039
  MAYO_16111  MAYO_16111  MAYO_36966  MAYO_36966 UN    NA  0.4534  0.5466  0.0000  0.2733  -1  0.783190  1.0000  5.4306
  MAYO_16111  MAYO_16111  MAYO_41741  MAYO_41741 UN    NA  0.2861  0.7139  0.0000  0.3569  -1  0.800751  1.0000 10.3969
  MAYO_16111  MAYO_16111   MAYO_8397   MAYO_8397 UN    NA  0.4173  0.5827  0.0000  0.2914  -1  0.784585  1.0000  7.6972
  MAYO_16111  MAYO_16111   MAYO_8727   MAYO_8727 UN    NA  0.4866  0.5134  0.0000  0.2567  -1  0.782473  1.0000  5.1180
  MAYO_16111  MAYO_16111   MAYO_9009   MAYO_9009 UN    NA  0.4845  0.5155  0.0000  0.2577  -1  0.769690  1.0000  4.9498
  MAYO_16111  MAYO_16111  NBB_S92304  NBB_S92304 UN    NA  0.4857  0.5143  0.0000  0.2571  -1  0.766779  1.0000  4.5508
  MAYO_16111  MAYO_16111  NBB_S94340  NBB_S94340 UN    NA  0.4248  0.5752  0.0000  0.2876  -1  0.785711  1.0000  6.1508
  MAYO_16111  MAYO_16111  SHRI_99_02  SHRI_99_02 UN    NA  0.4834  0.5166  0.0000  0.2583  -1  0.765768  1.0000  5.0214
  MAYO_16111  MAYO_16111  SHRI_99_06  SHRI_99_06 UN    NA  0.4858  0.5142  0.0000  0.2571  -1  0.767753  1.0000  5.2508
  MAYO_16111  MAYO_16111   SHRI_DDDD   SHRI_DDDD UN    NA  0.4093  0.5907  0.0000  0.2953  -1  0.756712  1.0000  6.3745
  MAYO_16111  MAYO_16111      SHRI_R      SHRI_R UN    NA  0.4690  0.5310  0.0000  0.2655  -1  0.761752  1.0000  5.1175
  MAYO_16111  MAYO_16111  WGACON_108  WGACON_108 UN    NA  0.3887  0.6113  0.0000  0.3057  -1  0.791057  1.0000  7.8225
  MAYO_16111  MAYO_16111  WGACON_114  WGACON_114 UN    NA  0.4939  0.5061  0.0000  0.2530  -1  0.768310  1.0000  5.1050
  MAYO_16111  MAYO_16111   WGACON_38   WGACON_38 UN    NA  0.4758  0.5242  0.0000  0.2621  -1  0.770775  1.0000  5.0774
```
g. **Are there any individuals that appear to be related to a lot of people? If yes, give a few examples.**
 - Answer: Yes, there are 20 individuals who are related to a lot of people. For example, MAYO_16111 is related to 21 people while MAYO_41741 is related to 15 people.

```bash
# Count number of relations where PI_HAT > 0.25
awk 'NR>1 && $10>0.25 {print $2; print $4}' TGEN_ibd.genome | sort | uniq -c | sort -rn | head -20
     21 MAYO_16111
     15 MAYO_41741
      8 MAYO_8397
      8 MAYO_16221
      6 WGACON_108
      5 SHRI_DDDD
      5 MAYO_36966
      5 MAYO_22544
      4 NBB_S94340
      3 WGACON_9
      3 WGAAD_8
      3 SHRI_R
      3 SHRI_QQ
      2 WGACON_39
      2 WGACON_38
      2 WGAAD_157
      2 SHRI_99_06
      2 SHRI_99_02
      2 NBB_S96057
      2 NBB_S91225
```
8. **Remove the individuals you identified as sharing more than 25% of alleles IBD with at least one other person in question 7.b using the --remove function in plink, and then calculate the heterozygote deficit statistic F for all individuals using the pruned genotypes. (Note, in general when cleaning data, one would not remove all such individuals, but only one from each pair if we wanted a data set with no related individuals. However in this case, when there are many pairs and a lot of poor data, we are removing all individuals in the pairs.) Use SAS, R, or any other method that ensures reproducibility to answer the following questions:**

    a. **What is the median, mean, standard deviation, minimum, and maximum F over the individuals?**
    - Answer: 
        Median: 0.002305
    
        Mean: -0.010285
        
        SD: 0.04261
        
        Miniumum: -0.2121

        Maxiumum: 0.1875
    
    b. **Which individuals have the lowest and highest F? Use F< -0.20 and F>0.1 to identify these individuals. How many are there?**
    - Answer: The three individuals with the lowest F is MAYO_8170 (-0.2124), MAYO_92854 (-0.2053), and WGAAD_381 (-0.2005).  The two individual with the highest F is NBB_S98014 (0.1875) and NBB_S90008 (0.1008). 
    c. **How does the distribution of F compare to the in-class example in the class notes? (Compare the values computed in a with the in class example values)**
    - Answer: This distribution of F has a larger standard deviation and a larger range and spread of data. For TGEN, 52% of  the variance is explained while for class, 7% of the variance is explained.

```bash
# Create file of individuals to be removed
awk 'NR>1 && $10>0.25 {print $1"\t"$2; print $3"\t"$4}' TGEN_ibd.genome | sort -u > individuals_to_remove.txt
# Remove related individuals, calculate F stats
plink --bfile TGEN_chr1-22_pruned --remove individuals_to_remove.txt --make-bed --out TGEN_unrelated
plink --bfile TGEN_unrelated --het --out TGEN_Fstat
# Create R script to calculate for F analysis
vi Fstat_analysis.R
```

Contents of Fstat_analysis.R script
```bash
het <- read.table("TGEN_Fstat.het", header=T, as.is=T)

summary(het$F)
mean(het$F)
sd(het$F)

# Count extreme F values
sum(het$F < -0.20)
sum(het$F > 0.1)

# Show individuals with extreme F
het[het$F < -0.20, c("FID", "IID", "F")]
het[het$F > 0.1, c("FID", "IID", "F")]

# Plot N non-missing genotypes vs F statistic
jpeg("hetplot.jpeg")
plot(het$N.NM., het$F, xlab="N Non-Missing", ylab="F")
dev.off()

# Linear regression
summary(lm(het$N.NM. ~ het$F))
```
```bash
# Load R and run script
module load R
Rscript Fstat_analysis.R > Fstat_analysis.log
# Examine contents of log
cat Fstat_analysis.log
    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
-0.212400 -0.022670  0.002305 -0.010285  0.015900  0.187500 
[1] -0.01028515
[1] 0.04260987
[1] 3
[1] 2
           FID        IID       F
284  MAYO_8170  MAYO_8170 -0.2124
310 MAYO_92854 MAYO_92854 -0.2053
909  WGAAD_381  WGAAD_381 -0.2005
           FID        IID      F
368 NBB_S90008 NBB_S90008 0.1008
573 NBB_S98014 NBB_S98014 0.1875
null device 
          1 
Call:
lm(formula = het$N.NM. ~ het$F)
Residuals:
    Min      1Q  Median      3Q     Max 
-4557.9  -209.6   164.2   402.7  1388.2 
Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) 78147.34      17.81  4387.5   <2e-16 ***
het$F       14754.77     406.49    36.3   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Residual standard error: 608.9 on 1235 degrees of freedom
Multiple R-squared:  0.5162,    Adjusted R-squared:  0.5158 
F-statistic:  1318 on 1 and 1235 DF,  p-value: < 2.2e-16
```