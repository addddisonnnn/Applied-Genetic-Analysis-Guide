## Post-Class Quiz 7: Post-GWAS Analyses - Conditional Analysis and Fine Mapping
### BS859 Applied Genetic Analysis
### March 4, 2026

Question 1 - Which of the following statements correctly describe the relationship among polygenic risk scores (PRS), LD score regression heritability estimates, and linkage disequilibrium (LD)? (Select all that apply.)
- **LD score regression estimates SNP heritability from GWAS summary statistics using LD information; this heritability represents an upper limit on how much variance a PRS can explain in an independent sample**
- **When the discovery GWAS sample size is small, LD score regression will generally overestimate heritability relative to the predictive performance achievable by PRS**
- **Differences in LD patterns and allele frequencies across populations can reduce the performance of PRS when applied across ancestries, even if the same causal variants are involved**
- **PRS methods that explicitly model LD (e.g., using an LD matrix from a well-matched reference) can improve prediction accuracy compared to simple clumping and thresholding**
- LD score regression requires individual-level genotype data from the discovery GWAS and therefore cannot be applied using summary statistics alone

Question 2 - If an independent validation dataset is not available, which approach is most appropriate for estimating the out-of-sample predictive performance of a polygenic risk score (PRS) while minimizing overfitting?
- **Use K-fold cross-validation (e.g., 10-fold) within the target dataset to obtain out-of-sample performance estimates**
- Use the base GWAS dataset as the validation sample to ensure independence from the target dataset
- Report the best-performing model from the training data as an unbiased estimate of predictive performance
- Increase the effective sample size by duplicating observations to stabilize performance estimates

Question 3 - In PRSice, empirical p-values are estimated using permutation testing. Suppose you run PRSice with --perm 1000 and obtain an Empirical-P value of 0.000999001. What does this result imply?
- Exactly 999 out of the 1000 permutations produced a p-value smaller than the observed p-value
- **None of the 1000 permutations produced a p-value smaller than the observed p-value, yielding an empirical p-value of 1/(1000+1)​**
- The observed (unadjusted) p-value is effectively treated as zero due to its magnitude
- The permutation procedure becomes unreliable when the observed p-value is smaller than 10^(-10)

Question 4 - Which of the following is NOT a standard step performed by PRSice when constructing a polygenic risk score using the clumping and thresholding (C&T) approach?
- Removing SNPs that are present in only one of the datasets (base or target)
- Excluding strand-ambiguous SNPs (e.g., A/T or C/G) to ensure allele alignment
- **Imputing missing genotype or summary statistic values in the base GWAS dataset**
- Performing LD clumping by retaining the SNP with the smallest p-value within each LD block

Question 5 - Which of the following statements are true regarding the impact of population structure on polygenic risk score (PRS) analyses? (Select all that apply.)
- **Population structure–related bias can be amplified in PRS because a large number of null variants may be included in the score**
- **When base and target samples are drawn from genetically similar populations, adjusting for principal components (PCs) can often adequately control for population stratification**
- PRSs developed in one ancestry generally transfer well to other ancestries as long as principal components are included in the model
- Differences in linkage disequilibrium (LD) patterns and allele frequencies across populations have minimal impact on PRS transferability

Question 6 - Which of the following statements correctly describe Bayesian polygenic risk score (PRS) methods (e.g., LDpred, sBayesR, PRS-CS)? (Select all that apply.)
- **They apply prior distributions to SNP effect sizes to induce shrinkage and may model mixtures of causal and non-causal variants**
- They require individual-level genotype data from the base GWAS to estimate SNP effects
- **They explicitly incorporate linkage disequilibrium (LD) structure (e.g., via an LD matrix or reference panel) when estimating effect sizes**
- They typically yield PRS estimates that are identical to those obtained from the best-performing clumping and thresholding (C&T) model