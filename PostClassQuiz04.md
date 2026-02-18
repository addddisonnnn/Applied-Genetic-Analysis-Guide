## Post-Class Quiz 3: 
### BS859 Applied Genetic Analysis
### February 4 2026

1. Which reference panel attribute is generally most beneficial for imputation accuracy, especially for low-frequency variants?
- Smaller size but single-population ancestry match
- Using only HapMap sample
- Only including unrelated individuals with missing phenotype data
- **Larger size and diverse ancestries (including the study’s ancestry where possible)**

2. For very rare variants (very low MAF), the imputation R2 is a highly reliable indicator of the true imputation accuracy.
- True
- **False**

3. Which of these statements about the imputation R2 (or INFO) metric is correct? (choose all that are correct)
- R2 is a per-individual measure of imputation certainty.
- **R2 estimates the squared correlation between imputed dosages and true genotypes and reflects the loss of information due to imputation; it becomes less reliable at very low MAF.**
- R2 is the exact squared correlation between imputed and true genotypes and is unbiased at all allele frequencies.
- R2 is only reported for genotyped (typed) variants, not imputed ones.

4. The imputed value commonly used as the predictor in regression (e.g., GWAS) is:
- The hard-called (best guess) genotype
- The imputation R2 metric
- The haplotype switch error rate
- **The allele dosage (expected allele count) derived from posterior genotype probabilities**

5. Suppose you imputed a cohort of N = 5,000 samples and observe an imputed SNP with R2 = 0.2 and MAF = 0.02. Approximately what is the effective sample size for association at this SNP, and what are the practical consequences?
- Effective N ≈ 5,000; include the SNP without concern. 
- Effective N ≈ 100; variant is essentially uninformative.
- **Effective N ≈ 1,000; power is reduced and the SNP may be excluded depending on study aims.**
- Effective N cannot be estimated from R2.

6. Which of the following most directly improves the ability to meta-analyze studies that used different genotyping chips?
- Increasing per-sample sequencing depth
- Using principal components as covariates
- **Genotype imputation to a common dense reference panel**
- Applying stricter Hardy–Weinberg filters