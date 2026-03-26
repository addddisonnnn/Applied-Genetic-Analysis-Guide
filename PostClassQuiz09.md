## Post-Class Quiz 7: Post-GWAS Analyses - Conditional Analysis and Fine Mapping
### BS859 Applied Genetic Analysis
### March 4, 2026

Question 1 - Which of these is a key limitation of applying LD Score regression to small GWAS (e.g., N < 5,000)?
- It cannot use summary statistics
- It always overestimates liability-scale heritability
- **It has low power and imprecise heritability estimates**
- It requires whole-genome sequencing

Question 2 - What is the primary purpose of using a Genetic Relationship Matrix (GRM) in GRM-based heritability estimation (e.g., in GCTA)?
- To identify causal variants for fine-mapping
- **To measure empirical genetic similarity among individuals for variance-component modeling**
- To replace phenotype measurements with genetic PCs
- To compute LD scores for summary statistics

Question 3 - Under the threshold model for a binary trait, an individual is affected when:
- The genotype contains at least one risk allele
- The observed phenotype exceeds the population mean
- **The underlying liability exceeds a threshold value**
- The environmental exposure exceeds a fixed level


Question 4 - A SNP annotation group or category enrichment score greater than 1 (e.g., in LD score regression) indicates that:
- The total heritability is zero
- The annotation category has no SNPs
- SNPs in the category contribute less heritability than expected
- **SNPs in the category contribute more heritability than expected**

Question 5 - Select all true statements about using cross-trait LD score regression for estimating genetic correlation:
- Sample overlap between two GWAS affects the intercept but not the slope in cross-trait LD score regression.
- **Sample overlap for the GWAS of the two traits biases the genetic correlation estimates in cross-trait LD score regression.**
- **The intercept is a measure of inflation due to cryptic relatedness/population structure for both traits combined.**
- If the two studies and traits are identical, the genetic correlation and heritability are identical.

Question 6 - Which of the following are reasons to remove strand-ambiguous SNPs during summary-statistics preprocessing for LD score regression? (choose all that are true)
- **Their alleles (A/T or C/G) can look the same after strand flipping, making allele alignment across datasets uncertain**
- They often have higher LD scores than other SNPs and may bias the regression if kept in the analysis
- They usually fail quality control because their imputation INFO scores are lower than average
- They are more likely to be rare variants that are not included in the LD reference panel