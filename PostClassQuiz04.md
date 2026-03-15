## Post-Class Quiz 3: 
### BS859 Applied Genetic Analysis
### February 11, 2026

1. Which reference panel attribute is generally most beneficial for imputation accuracy, especially for low-frequency variants?
- Smaller size but single-population ancestry match
- Using only HapMap sample
- Only including unrelated individuals with missing phenotype data
- **Larger size and diverse ancestries (including the study’s ancestry where possible)**

Feedback: Larger panels increase the chance of observing haplotypes that match study samples; including diverse ancestries helps, particularly for admixed samples and low-frequency variants. Small or older panels (HapMap) are less informative.

2. For very rare variants (very low MAF), the imputation R2 is a highly reliable indicator of the true imputation accuracy.
- True
- **False**

Feedback: At very low MAF, R2 is less reliable and can become noisy or biased; imputation accuracy is harder to estimate and often lower for rare variants.

3. Which of these statements about the imputation R2 (or INFO) metric is correct? (choose all that are correct)
- R2 is a per-individual measure of imputation certainty.
- **R2 estimates the squared correlation between imputed dosages and true genotypes and reflects the loss of information due to imputation; it becomes less reliable at very low MAF.**
- R2 is the exact squared correlation between imputed and true genotypes and is unbiased at all allele frequencies.
- R2 is only reported for genotyped (typed) variants, not imputed ones.

Feedback: R2 is an estimate of the squared correlation (imputed vs true) and indicates effective sample size (N * R2). It is less reliable at very low MAF. It is a per-variant metric, not per-individual, and is reported for imputed variants.

4. The imputed value commonly used as the predictor in regression (e.g., GWAS) is:
- The hard-called (best guess) genotype
- The imputation R2 metric
- The haplotype switch error rate
- **The allele dosage (expected allele count) derived from posterior genotype probabilities**

Feedback: Dosage (expected number of alternate alleles, e.g., 0–2) preserves imputation uncertainty and is typically used in regression. Best-guess GT loses uncertainty; R2 is a quality metric, not a predictor.

5. Suppose you imputed a cohort of N = 5,000 samples and observe an imputed SNP with R2 = 0.2 and MAF = 0.02. Approximately what is the effective sample size for association at this SNP, and what are the practical consequences?
- Effective N ≈ 5,000; include the SNP without concern. 
- Effective N ≈ 100; variant is essentially uninformative.
- **Effective N ≈ 1,000; power is reduced and the SNP may be excluded depending on study aims.**
- Effective N cannot be estimated from R2.

Feedback: Effective N ≈ N * R2 = 5,000 * 0.2 = 1,000. This reduced effective sample size means lower power to detect association; many analyses use R2 thresholds (e.g., 0.3) to exclude such variants, though for specific hypotheses one might keep them with caution.

6. Which of the following most directly improves the ability to meta-analyze studies that used different genotyping chips?
- Increasing per-sample sequencing depth
- Using principal components as covariates
- **Genotype imputation to a common dense reference panel**
- Applying stricter Hardy–Weinberg filters

Feedback: Different chips assay different variant sets. Imputing all studies to the same reference panel fills in variants and creates a common variant set for meta-analysis. PCs and QC are important but do not bridge differing chip content.
