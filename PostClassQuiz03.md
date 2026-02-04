## Post-Class Quiz 3: 
### BS859 Applied Genetic Analysis
### February 4 2026
1. In a standard GWAS using an additive genotype coding for a biallelic SNP (alleles A and C), with coded allele A and reference allele C, what does the regression coefficient β for the coded genotype represent in a logistic regression model?
    - The difference in mean phenotype between AA and CC
    - **The log of the odds ratio for disease for a genotype with one additional copy of the A allele compared to the reference genotype**
    - The odds ratio comparing heterozygotes to homozygous reference only
    - The variance in phenotype explained by the SNP

In logistic regression with additive coding (0,1,2), the genotype coefficient β is the change in log odds per additional copy of the coded (effect) allele. exp(β) is the odds ratio for a one-allele increase.

2. When deciding whether to include principal components (PCs) of genetic data as covariates in a GWAS, which statement is correct?
    - You must always include the top 10 PCs to avoid any confounding
    - PCs are irrelevant if you use PLINK for association testing 
    - Including more PCs always increases power in logistic regression
    - **PCs need to be included if they are associated with the phenotype of interest.**

If PCs are unrelated to the phenotype, they are unnecessary for avoiding confounding. One common practice is to adjust for PCs that are associated with the trait (or to include a set number of top PCs for conservative control). Including many PCs unnecessarily can reduce power in logistic regression.

3. The genomic control inflation factor λ<sub>GC</sub> is computed from the ratio of:
    - **Observed median test statistic to expected median test statistic under the null**
    - Observed mean p-value to 0.05
    - Number of significant SNPs to total SNPs tested
    - Observed genomic heritability to total phenotypic variance

λ<sub>GC</sub>= (median observed χ2) / (median expected χ2 under null). It quantifies overall inflation of test statistics relative to null expectation.

4. In GWAS with millions of correlated tests, the conventional genome-wide significance threshold ~5×10^-8 is based on a false discovery rate of 0.05
    - True
    - **False**

The conventional genome-wide threshold (~5×10^-8) comes from Bonferroni correction for roughly 1e6 independent tests (0.05/1e6). This is a pragmatic approximation given linkage disequilibrium between SNPs.

5. Which of the following best describes the purpose of QQ plots in GWAS?
    - To display genomic positions vs -log10(p) across the genome
    - To estimate per-SNP effect sizes
    - **To compare observed p-value quantiles to expected under the null to detect inflation and/or enrichment of associations**
    - To compare observed p-value quantiles to expected under the null to detect inflation or enrichment of associations

QQ plots show observed vs expected p-value quantiles; deviations from the diagonal in the middle of the distribution can indicate inflation (population structure, bias) or true signal enrichment at the tail.

6. Which advantage of generalized linear mixed models (GLMMs) like those implemented in GMMAT or GENESIS is most relevant when analyzing binary traits in GWAS with population structure and distantly related individuals?
    - GLMMs avoid the need to compute principal components entirely.
    - **GLMMs model a random effect based on the genetic relationship matrix (GRM) to account for both cryptic relatedness and population structure while allowing a logistic link for binary outcomes.**
    - GLMMs always produce smaller p-values than standard logistic regression, increasing power in all settings.
    - GLMMs require no assumptions about the distribution of the outcome.

GLMMs incorporate a random effect whose covariance is proportional to the GRM, which captures genetic similarity between individuals. This allows proper control for both relatedness and some forms of population structure in analyses of binary traits by combining the logistic (or other) link with random effects. They do not make PCs unnecessary, do not universally reduce p-values, and still rely on appropriate outcome distribution assumptions (e.g., Bernoulli for binary traits with a logistic link).
