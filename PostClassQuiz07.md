## Post-Class Quiz 7: Post-GWAS Analyses - Conditional Analysis and Fine Mapping
### BS859 Applied Genetic Analysis
### March 4, 2026

Question 1 - To approximate joint SNP effects from single-SNP marginal effects and a reference LD panel, which matrix is required or estimated?
- The environmental covariance matrix
- **X'X (the genotype variance-covariance matrix) estimated using allele frequencies and LD from a reference sample**
- The transcriptome co-expression matrix
- The inverse Fisher information matrix from unrelated traits

Question 2 - Which statement correctly contrasts the heuristic LD-threshold approach with penalized regression (e.g., LASSO) for fine mapping?
- Heuristic approaches estimate joint effects while penalized regression ignores multicollinearity
- Both methods are optimized for finding multiple causal variants in high-LD regions and always return identical sets
- **Heuristic uses pairwise LD thresholds around lead SNP but does not model joint effects; penalized regression models SNPs jointly and performs shrinkage/selection**
- Penalized regression never excludes causal variants

Question 3 - In Bayesian fine-mapping, what does PIP represent?
- **The posterior probability that a SNP is included in a causal model**
- The genome-wide adjusted p-value for a SNP
- The probability that a SNP is in high LD with the lead variant
- The estimated effect size after shrinkage

Question 4 - GCTA COJO conditional analysis using GWAS summary statistics does not require an external LD reference panel.
- True
- **False**

Question 5 - How can functional genomic annotation be incorporated into statistical fine-mapping?
- **By incorporating annotation into the prior model or using it to prioritize variants within a credible set**
- By directly estimating SNP–SNP correlations within the associated region
- By replacing the need for conditional or joint modeling
- By restricting analysis only to protein-coding variants

Question 6 - When annotating a list of SNPs with Ensembl VEP, which of the following is true?
- VEP reports a single summary consequence per variant and ignores transcript-specific differences
- VEP estimates linkage disequilibrium between nearby variants during annotation
- VEP performs Bayesian fine-mapping and calculates posterior inclusion probabilities
- **VEP can return multiple transcript-specific consequences per variant, including predicted functional impacts**