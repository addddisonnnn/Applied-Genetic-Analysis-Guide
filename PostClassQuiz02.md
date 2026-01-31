## Post-Class Quiz 2: 
### BS859 Applied Genetic Analysis
### January 28, 2026
1. Which definition best describes population structure in genetic studies?
    - A population where every individual has identical genotypes
    - **Deviation from random mating, producing subpopulations with different allele frequencies and correlations among variants**
    - Complete panmaxia across all geographic regions
    - A dataset with no missing genotypes

Feedback: Population structure occurs when mating is non-random (geographic separation, selection, inbreeding, assortative mating), producing subgroups with different genotype frequencies.

2. Population structure can cause spurious associations in GWAS when which of the following is true? Select up to 2 options
    - The study uses only unrelated individuals from a single homogenous population
    - Subpopulations differ in trait distributions but have identical genotype frequencies
    - **Subpopulations differ in both allele frequencies and trait distributions**
    - Subpopulations differ in allele frequencies but have identical trait distributions

Spurious associations arise when ancestry (population structure) is associated with both genotype and phenotype. If both vary across subgroups, genotype-phenotype associations can be confounded.

3. Principal Components Analsysis (PCA) on the full set of genome-wide markers is a common and effective method to detect population structure in GWAS data. 
    - True
    - **False**

PCA summarizes genetic variation into orthogonal axes (PCs) that often reflect ancestry; However, we do not use the full set of markers.  LD pruning is used to remove highly correlated markers and reduce correlation among markers before PCA.

4. When using PCA to infer broad ancestry (e.g., European vs. East Asian vs. African), which approach is recommended?
    - Run PCA on only you study damples without external references
    - Use only A/T and G/C SNPs for the PCA
    - Remove all SNPs with MAF > 0.05
    - **Include reference samples with known population labels (e.g., HapMap/1000G) in the PCA so study samples can be compared to known populations**

Including reference samples lets you place your samples relative to known population clusters and identify ancestry or outliers

5. Why are A/T and G/C SNPs often removed before merging datasets from different sources for PCA? 
    - They have the highest allele frequencies
    - They are always multiallelic
    - **Their strand orientation is ambiguous (complements are the same), risking incorrect matching across datasets**
    - They are always on the sex chromosomes

A/T and G/C SNPs look identical after strand flip (A<->T, G<->C), making it hard to determine correct allele orientation when merging datasets; excluding them avoids strand-mismatch errors.