## Homework 2: Population Structure Detection
### BS859 Applied Genetic Analysis
### January 28, 2026

```bash
# Set the working directory and load PLINK
export DATADIR=/projectnb/bs859/data/tgen
module load plink/1.90b6.27 
```

1. **In our class example PCA using TGEN merged with hapmap samples, what are the IDs of the two TGEN samples that cluster with the YRI and JPT+CHB on the first two PCs? Which individual appears to be clustering with the CHB+JPT samples, and which close to the Yoruban?** 
    - Answer: 
```bash
```
2. **What are the IDs and case status of the other individuals who did not cluster closely with the CEU on PCs 1 and 2? Explain what filters you use to identify them.**
    - Answer: 
```bash
```

3. **Use plink to remove variants with minor allele frequency <0.02, and with genotype missingness >0.02 and LD prune the variants using the parameters: --indep-pairwise 10000kb 1 0.2 Create a new data set that keeps only the pruned-in variants and has only variants from chromosomes 1-22**

    a. **Fully explain what the pruning parameters we are using mean**

    b. **Print the first 10 SNPs in the pruned data set**
    
    c. **Report the number of SNPs that remain after the filtering on minor allele frequency and pruning**
    
    d. **Run smartpca with this pruned data set. Show your smartpca parameter file and plot PC1 vs PC2.**

4. **Use plink to remove variants with minor allele frequency <0.02, and with genotype missingness >0.01 and LD prune the variants using the parameters: --indep-pairwise 10000kb 1 0.15 Again, create a new data set that keeps only the pruned-in variants and has only variants from chromosomes 1-22. Make sure you give this data set a different name than the pruned data set in 3) so that you do not overwrite it.**

    a. **Fully explain what the pruning parameters we are using mean**

    b. **Print the first 10 SNPs in the pruned data set**
    
    c. **Report the number of SNPs that remain after the filtering on minor allele frequency and pruning**
    
    d. **Run smartpca with this 2nd pruned data set. Show your smartpca parameter file and Plot PC1 vs PC2.**

5. Using the output from the 2 PCAs performed in questions 3 and 4,
    
    a. **Compare: which PCs are significantly different between cases and controls for the two different PCAs in questions 3 and 4? Use p<0.005 to define significance (0.05/10, since we are testing 10 PCs for association with case status), and use the ANOVA output from smartpca.**

    b. **Visually compare the PC1 vs PC2 plots for the two PCAs. Do you see any major differences in the clusters of points? Ignore differences in sign (positive vs negative), and focus more on the general picture.**

6. **Find the two individuals who clustered with the YRI and CHB+JPT populations (see Question 1) in your PC1 vs PC2 plots. Do they appear to be far from most of the other individuals in the plot?**

7. **Individual MAYO_8170 had the lowest F statistic (-0.2124) in Homework 1. Where is this person on the PC1 vs PC2 plot? Does the PCA suggest that this person’s genotype data is extreme compared to the other individuals in the study?**

8. **We will use the PCs from question 3 as covariates for future association analyses. We will need to re-format the “.evec” file to be in the format plink requires, where the first two columns are Family id and individual ID, and the rest of the columns are covariates, and all columns are whitespace delimited.**

    a. **What changes to the evec file are needed so that plink can read it?**
    
    b. **Write code (linux shell, R, SAS, python, your choice!) to convert the file to the needed format. Share your code, and show its results (e.g.,the first 10 lines of the final covariate file).**
    - Answer: 
```bash
```