## Post-Class Quiz 3: 
### BS859 Applied Genetic Analysis
### February 4 2026
1. In a standard GWAS using an additive genotype coding for a biallelic SNP (alleles A and C), with coded allele A and reference allele C, what does the regression coefficient β for the coded genotype represent in a logistic regression model?
    - The difference in mean phenotype between AA and CC
    - The log of the odds ratio for disease for a genotype with one additional copy of the A allele compared to the reference genotype
    - The odds ratio comparing heterozygotes to homozygous reference only
    - The variance in phenotype explained by the SNP

2. When deciding whether to include principal components (PCs) of genetic data as covariates in a GWAS, which statement is correct?
    - You must always include the top 10 PCs to avoid any confounding
    - PCs are irrelevant if you use PLINK for association testing 
    - Including more PCs always increases power in logistic regression
    - PCs need to  be included if they are associated with the phenotype of interest.

3. The genomic control inflation factor λ<sub>GC</sub> is computed from the ratio of:
    - Observed median test statistic to expected median test statistic under the null
    - Observed mean p-value to 0.05
    - Number of significant SNPs to total SNPs tested
    - Observed genomic heritability to total phenotypic variance

4. In GWAS with millions of correlated tests, the conventional genome-wide significance threshold ~5×10^-8 is based on a false discovery rate of 0.05
    - True
    - False

5. Which of the following best describes the purpose of QQ plots in GWAS?
- To display genomic positions vs -log10(p) across the genome
- To estimate per-SNP effect sizes
- To compare observed p-value quantiles to expected under the null to detect inflation and/or enrichment of associations
- To compare observed p-value quantiles to expected under the null to detect inflation or enrichment of associations

6. Which advantage of generalized linear mixed models (GLMMs) like those implemented in GMMAT or GENESIS is most relevant when analyzing binary traits in GWAS with population structure and distantly related individuals?
- GLMMs avoid the need to compute principal components entirely.
- GLMMs model a random effect based on the genetic relationship matrix (GRM) to account for both cryptic relatedness and population structure while allowing a logistic link for binary outcomes.
- GLMMs always produce smaller p-values than standard logistic regression, increasing power in all settings.
- GLMMs require no assumptions about the distribution of the outcome.