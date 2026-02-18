## Session 2: Genotype imputation and analysis of imputed genotypes
### BS859 Applied Genetic Analysis
### Feburary 11, 2026

### Section 1 - What is Genotype Imputation?

### Section 2 - Why do we Impute? (The "So What? Question)

### Section 3 - How Imputation Works (Conceptual, not Mathematical)

### Section 4 - Hidden Markov Model (HMM) - What you actually need to know

### Section 5 - What Imputed Data Looks Like
### Section 6 - Imputation Quality - R² Metric

### Section 7 - Factors Affacting Imputation Accuracy

### Section 8 - Reference Panels - A Historical Perspective

### Section 9 - The Imputation Workflow
#### Step 1: Pre-imputation QC
More stringent than regular GWAS QC:

- Remove SNPs with low call rate (<95-98%)
- Remove SNPs with Hardy-Weinberg p < 1e-6
- Remove SNPs with MAF < 0.01 (for chip data)
- Check strand alignment with reference
- Update genome build if needed (to build 38)

#### Step 2: Phasing
Why phase? Imputation algorithms work on haplotypes, not genotypes

Software:
- EAGLE2: Default for imputation servers, very fast
- SHAPEIT2: Traditional choice, slower but accurate
- HAPI-UR: For very large datasets

What phasing does: Takes genotype data (0,1,2) and estimates which alleles are on which chromosome copy

#### Step 3: Imputation
Two main approaches:
1. Local imputation (run on your computer): IMPUTE2, Minimac3
2. Server-based imputation: Upload to Michigan/TOPMed/Sanger server

Why servers are better:
- No access to individual-level WGS reference data (privacy!)
- No need for massive computational resources
- Automatic QC and liftover
- Free!

#### Step 4: Post-imputation QC
Must-do filters:
1. R² < 0.3 → Remove low frequency and poorly imputed variants
2. MAF < 0.005 → Remove (unless very large sample size)
3. Frequency mismatches vs reference → Investigate/remove

Optional but recommended:
1. Check that effect alleles match between datasets
2. Compare allele frequencies to external reference
3. Visual inspection of top association signals

#### Section 10 - Association Analysis with Imputed Data

#### Section 11 - Case Study - TGEN Data on Chromosome 19

For homework, chromosome 2 was imputed and analysis on BIN1 variant
chromosomes are ordered from largest to smallest, chromosome 2 is the next largest. Make sure to delete files. Things will take long about four times as long as chromosome 19. Want to identify if there's evidence of replication. Repliation in genetics terms - association with that vairanty in our data set in the same direction, i.e. T allele, we want to see if the T allele has an OR>1. doesn't have to have a genome-wide signficant. significance level of 0.05 since this is a single association (chromosome 2). 