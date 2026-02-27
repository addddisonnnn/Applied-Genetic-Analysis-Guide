## Post-Class Quiz 6: Meta-Analysis for genetic studies
### BS859 Applied Genetic Analysis
### February 25, 2026

Question 1
Which statement about the weighted (signed) Z-score meta-analysis is TRUE?

It always gives the same result as inverse-variance pooling, regardless of study design

It produces a pooled effect size estimate with standard error

It ignores the direction of effect

*It does not require effect size estimates on the same scale across studies and combines signed Z-values weighted by sample size*

2. For meta-analysis of gene-based rare-variant SKAT tests, each study must provide:

Gene-level p-values and sample sizes for each gene

*Allele frequencies, single-variant score statistics, and the genotype covariance (Φ) matrix capturing LD structure
*
Single-variant effect estimates (β) and standard errors for each variant in the gene

Only burden test statistics for each gene from each study

3. Which of the following actions is required before combining effect estimates across studies?

Ensure that all studies use identical genotyping platforms

Standardize all effect sizes to have mean zero and variance one

*Re-align effect estimates so they correspond to the same effect (counted) allele and correctly handle strand orientation (especially for A/T and C/G SNPs)*

Exclude all SNPs showing heterogeneity before meta-analysis

4. By default, which meta-analysis method does METAL use when no scheme is specified in the control file?

Fixed-effect inverse-variance weighting

*Sample-size weighted signed Z-score meta-analysis
*
Random-effects REML meta-analysis

Gene-based SKAT meta-analysis

5. Under a fixed-effects meta-analysis model, as study sample sizes increase, the study-specific effect estimates would converge to the same underlying true effect.

*True*

False

6. Which of the following statements about genomic-control (GC) correction in METAL is TRUE?

*GC correction rescales test statistics using a genomic inflation factor (λ) and generally reduces inflated statistics, making p-values less significant*

GC correction always increases the significance of strongly associated SNPs

GC correction is unnecessary in meta-analysis because combining studies removes population stratification

GC correction alters allele frequencies and effect estimates in the output