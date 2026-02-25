## Week 5: Association analysis of Rare Variants and Sequence Data
### BS859 Applied Genetic Analysis
### Feburary 18, 2026

### Part 1: Lecture Content - Conceptual Foundations

### 1. Why Rare Variants Matter

#### 1.1 The Missing Heritability Problem
For decades, family studies have estimated the heritability of traits - the proportion of variance explained by genetic factors. With GWAS, we can estimate SNP-based heritability - variance explained by all SNPs in the study

The puzzle: Family-based heritability is almost always LARGER than SNP-based heritability. This gap is called "missed heritability."

Possible explanations for missed heritability:
- Many common variants with very small effects (too small to detect with current sample sizes)
- Structural variants not captured by genotyping chips
- Epigenetic factors
- Over-estimation of family heritability (shared environment, gene-gene interactions)
- Rare variants with larger effect ← focus of today's class

#### 1.2 Why Rare Variants Might Have Larger Effects
Natural selection: if a variant disrups gene function and is deleterious, natural selection will prevent it from becoming common. Carriers produce less, so the variants stay rare.

Evidence from multiple species:
- Fruit flies: negative correlation between allele frequency and effect size
- Yeast: rare variants have larger effects on traits
- Humans: same patterns observed in UK biobank data - rare variants tend to have larger effect sizes

#### 1.3 Not All Rare Variants Are Bad
Important counter-example: PCSK9 gene mutations
- Some rare variants in PCSK9 cause lower LDL cholestorol
- Individuals with these mutations have protection against heart disease
- This discovery led to new drug development (PCSK9 inhibitors)

Key takeaway: Rare variants can be protective too - direction matters!

### 2. Approaches to Study Rare Variants
#### 2.1 Options for Obtaining Rare Variant Data
| Method                       | What It Does                                   | Advantages                                          | Disadvantages                  |
|------------------------------|------------------------------------------------|-----------------------------------------------------|--------------------------------|
| High-depth WGS               | Sequence entire genome at 30X coverage         | Identifies nearly all variants with high confidence | Very expensive                 |
| Low-depth WGS                | Sequence entire genome at 5X coverage          | More cost-effective                                 | Limited accuracy               |
| Whole Exome Sequencing (WES) | Sequence only the protein-coding ~1% of genome | Identifies coding variation at reasonable cost      | Misses non-coding variants     |
| GWAS chip + imputation       | Genotype common variants, impute using TOPMed  | Inexpensive                                         | Low accuracy for rare variants |

#### 2.2 Rare Variation is Abundant!
From TOPMed (160,000 deeply sequenced individuals):

Key Statistic:
- Total SNPs: 781 million
- Singletons (seen once): 46.4& of all variants
- Doubletons (seen twice): 15.7% of all variants
- Functional variants are even rarer:
    - Missense: 46.4% are singletons
    - Stop gain: 53.3% are singletons
    - Frameshift: 60.0% are singletons

The pattern: The most functionally important variants are the rarest. This makes statistical analysis challenging.

### 3. Why Standard GWAS Fails for Rare Variants
#### 3.1 The Power Problem
Example calculation:
- Variant with MAF = 0.0001 (1 in 10,000)
- Significance level α = 5×10⁻⁶ (correcting for 100,000 tests)
- Disease prevalence ≈ 10%
- To detect GRR = 2 with 80% power:
    - Need 342,000 cases AND 342,000 controls

The problem: Even with biobank-scale samples, most rare variants are too rare to test individually.

#### 3.2 The Solution: Aggregation
Instead of testing each rare variant individually, we group them together based on biological units (genes) or genomic regions.

Core idea: Test whether the set of rare variants in a gene is associated with the phenotype, not each variant individually.

### 4. Types of Rare Variant Tests
#### 4.1 General Statistical Framework
For individual i, phenotype yi, covariates Xi, and genotypes Gij for m variants in a set:

Null model: h(μᵢ) = α₀ + α'Xᵢ (no genetic effects)

Score statistic for variant j: Sⱼ = Σᵢ Gᵢⱼ(yᵢ - μ̂ᵢ)
- Sⱼ > 0: variant associated with increased trait/disease risk
- Sⱼ < 0: variant associated with decreased trait/disease risk

#### 4.2 Burden Tests
What they do: Collapse multiple rare variants into a single genetic score

Test statistic: Q_burden = (Σⱼ wⱼSⱼ)²

Equivalent formulation: Create genetic score Cᵢ = Σⱼ wⱼGᵢⱼ, then test Cᵢ in regression

When they're powerful:
- Most variants in the set are causal
- All effects are in the SAME direction
- Example: all variants disrupt gene function (all harmful)

When they fail:
- Effects in opposite directions (some protective, some harmful)
- Only a small fraction of variants are causal

#### 4.3 Variance Component Tests (SKAT)
What they do: Treat variant effects as random and test whether the variance of effects is >0

Test statistic: Q_SKAT = Σⱼ wⱼ²Sⱼ²

Key property: Uses Sⱼ² instead of Sⱼ, so it's robust to direction of effects

When they're powerful:
- Effects in opposite directions
- Many variants are null (no effect)
- Example: a gene where some variants increase risk, others decrease risk

When they're less powerful:
- Most variants are causal with same direction (burden test would be better)

#### 4.4 Hybrid/Omnibus Tests
Combine the best of both approaches:

SKAT-O: Q_ρ = ρ×Q_burden + (1-ρ)×Q_SKAT
- ρ estimated adaptively from data
- ρ near 1 → use burden-like test
- ρ near 0 → use SKAT-like test

SMMAT-E (Efficient hybrid test in GMMAT package):
- Combines two asymptotically independent tests
- Computationally efficient
- Robust across different scenarios

#### 4.5 Weights in Rare Variant Tests
Common weighting schemes:
| Weight Type     | Formula                     | Effect                                   |
|-----------------|-----------------------------|------------------------------------------|
| CAST            | wⱼ = 1 if rare, 0 otherwise | Simple presence/absence                  |
| Madsen-Browning | wⱼ = 1/√[pⱼ(1-pⱼ)]          | Up-weights rarer variants                |
| SKAT Beta       | Beta(pⱼ, a₁=1, a₂=25)       | Smooth down-weighting of common variants |

Visual interpretation: Beta(1,25) gives highest weight to variants with MAF ~0.01-0.03, very low weight to MAF >0.05

### 5. Power Comparisons: Which Test to Use?
#### 5.1 Simulation Results from Chen et al. 2019
Key findings:
| Scenario                    | Best Test      | Why                            |
|-----------------------------|----------------|--------------------------------|
| 100% causal, same direction | Burden         | All variants contribute signal |
| 80% causal, 20% opposite    | SKAT/SMMAT-E   | Handles mixed directions       |
| 50% causal, 50% opposite    | SKAT/SMMAT-E   | Robust to opposite effects     |
| 10-20% causal, mixed        | SMMAT-E/SKAT-O | Adapts to sparse signals       |

Take-home: SMMAT-E and SKAT-O are generally the safest choices because they adapt to the underlying genetic architecture.

#### 5.2 Type I Error Control
All SMMAT tests have well-controlled Type I error:
- Empirical rates close to nominal α (0.05, 0.0001, 2.5×10⁻⁶)
- Works for both continuous and binary traits
- Handles related individuals via kinship matrix

### 6. Practical Considerations for Rare Variant Analysis
#### 6.1 Which Variants to Include?
Gene-based testing options:
- All variants in gene
- Only rare variants (MAF < 0.05 or 0.01)
- Only functional variants (non-synonymous, stop-gain, frameshift)
- Only variants predicted deleterious by multiple algorithms

Beyond genes:
- Sliding windows across genome
- Regulatory regions
- Evolutionarily conserved regions

Caveat: Annotation tools are imperfect. Functional predictions can be wrong.

#### 6.2 Multiple Testing Correction
Key principle: Number of tests = number of genes/regions tested, NOT number of variants

Bonferroni adjustment: α = 0.05 / (# genes with cMAC ≥ threshold)

But what's cMAC?

#### 6.3 Cumulative Minor Allele Count (cMAC)
Formula: cMAC = n_variants × mean_MAF × 2 × N

Where:
- n_variants = number of variants in the gene/region
- mean_MAF = average MAF of variants in the set
- N = sample size

Why cMAC matters:
- If cMAC is very small (e.g., <10), the test has no power
- Even if all rare alleles were in cases only, you couldn't reach significance
- Best practice: Exclude genes with cMAC < 10 from consideration AND from multiple testing correction

Example from class:
- All exonic variants: 454 genes with cMAC ≥ 10, 43 genes with cMAC < 10
- Non-synonymous only: 377 genes with cMAC ≥ 10, 75 with cMAC < 10

#### 6.4 Significance Threshold Calculation
Step 1: Count genes with cMAC ≥ 10 → let this be M
Step 2: Bonferroni threshold = 0.05 / M
Step 3: Only consider genes with p-value < this threshold as significant

If you run multiple tests per gene (e.g., Burden, SKAT, SKAT-O, SMMAT-E), you need to account for that too. Common approach: use ACAT to combine p-values into one per gene.

#### 6.5 ACAT Test (Aggregated Cauchy Association Test)
Combines multiple p-values from the same gene into a single p-value:

T = Σ wₖ × tan[(0.5 - pₖ)π]

Advantages:
- Robust to correlation between tests
- Fast p-value computation (Cauchy distribution)
- Can combine results from different variant filters

### 7. Population Structure Adjustment
Same principles as common variant GWAS:
- Principal components as fixed effects in regression
- Genetic relationship matrix (GRM) as random effect in mixed models
- SMMAT handles both approaches

### 8. VCF Format - Essential for Sequence Data
#### 8.1 VCF Structure
Two main sections:
1. Header (lines starting with ##)
- File format version
- INFO field definitions
- FORMAT field definitions
- Contig information

2. Data (one line per variant, 8 fixed columns + sample columns)
#### 8.2 Fixed Columns (Required)
| Column  | Description               | Example        |
|---------|---------------------------|----------------|
| #CHROM  | Chromosome                | 20             |
| POS     | Position (1-based)        | 68396          |
| ID      | Variant identifier        | rs123 or .     |
| REF     | Reference allele          | C              |
| ALT     | Alternate allele(s)       | T              |
| QUAL    | Phred-scaled quality      | 100            |
| FILTER  | QC filters passed         | PASS           |
| INFO    | Variant-level info        | AF=0.01;DP=100 |
| FORMAT  | Format of genotype fields | GT:DP:GQ       |
| SAMPLE1 | First sample's data       | 0/1:13:99      |
| SAMPLE2 | Second sample's data      | 0/0:20:99      |

#### 8.3 Genotype Encoding
GT field:
- `0/0` = homozygous reference
- `0/1` = heterozygous (reference/ALT1)
- `1/1` = homozygous ALT1
- `0/2` = heterozygous (reference/ALT2) for multi-allelic
- `/` = unphased, `|` = phased

#### 8.4 INFO Field Annotations (ANNOVAR example)
From the example VCF:
```text
##INFO=<ID=Func.refGene,Number=.,Type=String,Description="Func.refGene annotation">
##INFO=<ID=Gene.refGene,Number=.,Type=String,Description="Gene.refGene annotation">
##INFO=<ID=ExonicFunc.refGene,Number=.,Type=String,Description="ExonicFunc.refGene annotation">
```

Example data:
- Variant 1: `Func.refGene=upstream;Gene.refGene=DEFB125` (upstream of gene)
- Variant 3: `Func.refGene=exonic;Gene.refGene=DEFB125;ExonicFunc.refGene=nonsynonymous_SNV` (in gene, changes protein)

These annotations are used to create group files for SMMAT.

#### 8.5 Working with VCF Files: bcftools
Useful commands:
```bash
# Extract specific fields
bcftools query -f '%CHROM %POS %REF %ALT %INFO/Func.refGene\n' file.vcf.gz

# Filter by annotation
bcftools view -i 'INFO/Func.refGene="exonic"' file.vcf.gz

# Count variants
bcftools index --stats file.vcf.gz
```

### Part 2: Code Walkthrough - Step-by-Step Analysis

### 9. Example 1: 1000 Genomes Exome Data
#### 9.1 Data Files
Located in `/projectnb/bs859/data/rarevariant/`
| File                                                          | Description                                                  |
|---------------------------------------------------------------|--------------------------------------------------------------|
| 1000G_exome_chr20_example_softFiltered.calls.hg19_anno.vcf.gz | Genotypes + annotations for 19,101 variants, 822 individuals |
| example.ped                                                   | Phenotype data (simulated quantitative trait)                |
| example.gds                                                   | GDS format of the VCF (pre-converted)                        |
| example.smmat.exonic.groups                                   | Group file with all exonic SNPs per gene                     |
| example.smmat.nonsyn.groups                                   | Group file with non-synonymous SNPs per gene                 |

#### 9.2 Examine the Phenotype File
```bash
head /projectnb/bs859/data/rarevariant/example.ped
```
Output:

```text
#FID IID PAT MAT SEX QPHENO pop.EU pop.AF pop.HI
HG00096 HG00096 0 0 2 1.75284949112253 1 0 0
HG00100 HG00100 0 0 1 -0.696558932247186 1 0 0
...
```

Understanding the data:
- QPHENO: Quantitative phenotype (simulated)
- pop.EU, pop.AF, pop.HI: Dummy-coded population indicators
    - pop.EU=1, pop.AF=0, pop.HI=0 → European ancestry
    - pop.EU=0, pop.AF=1, pop.HI=0 → African ancestry
    - pop.EU=0, pop.AF=0, pop.HI=1 → Hispanic ancestry
    - All zero → Asian ancestry

#### 9.3 Examine Group Files
```bash
head example.smmat.exonic.groups
```
Output:

```text
DEFB125 20 68396 C T 1
DEFB125 20 76689 T C 1
DEFB125 20 76690 T C 1
...
```
Format: Gene Chromosome Position Ref Alt Weight

Note: Weight = 1 means use default SKAT Beta(1,25) weights, not additional weighting.

#### 9.4 R Script: SMMAT_example.R
Step 1: Load Libraries and Data
```r
library(GMMAT)

# Read phenotype file
pheno <- read.table("/projectnb/bs859/data/rarevariant/example.ped", 
                    header=TRUE, comment.char="")

# Fix column names (the # in header causes issues)
colnames(pheno)[1:5] <- c("FID","IID","fa","mo","sex")
```
Step 2: Fit Null Model
```r
# Null model: phenotype ~ covariates, no SNPs
# Quantitative trait → family=gaussian, kins=NULL (no relatedness)
model.0 <- glmmkin(QPHENO ~ sex + pop.EU + pop.AF + pop.HI, 
                   data = pheno, 
                   id = "IID",
                   kins = NULL, 
                   family = gaussian(link = "identity"))
```
What's happening:
- glmmkin fits generalized linear mixed model
- No random effect (kins=NULL) because samples are unrelated
- Covariates: sex + population indicators
- This model gives us μ̂ᵢ for each person under null hypothesis

Step 3: Run SMMAT on All Exonic Variants
```r
example.exonic.groups <- SMMAT(
  null.obj = model.0,
  geno.file = "/projectnb/bs859/data/rarevariant/example.gds",
  group.file = "example.smmat.exonic.groups",
  group.file.se = " ",           # space-separated group file
  MAF.range = c(1e-7, 0.05),     # only variants with MAF ≤ 5%
  MAF.weights.beta = c(1, 25),   # Beta(1,25) weights for SKAT
  miss.cutoff = 1,                # no missing cutoff
  missing.method = "impute2mean", # impute missing to mean
  method = "davies",              # p-value computation method
  tests = c("O", "E"),            # return SKAT-O and SMMAT-E
  rho = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),  # grid for SKAT-O
  use.minor.allele = TRUE,        # use minor allele, not ALT allele
  auto.flip = FALSE,
  is.dosage = FALSE,              # hard calls, not dosages
  ncores = 1,
  verbose = TRUE
)
```
Key parameters explained:
| Parameter        | Value          | Meaning                                                |
|------------------|----------------|--------------------------------------------------------|
| `MAF.range`        | c(1e-7, 0.05)  | Include variants with MAF between essentially 0 and 5% |
| `MAF.weights.beta` | c(1, 25)       | SKAT weights: up-weight rare variants                  |
| `tests`            | c("O", "E")    | Return SKAT-O (optimal) and SMMAT-E (efficient)        |
| `rho`              | grid of values | Correlation parameters for SKAT-O                      |
| `use.minor.allele` | TRUE           | Count minor allele (not necessarily ALT)               |
| `is.dosage`        | FALSE          | Using hard genotypes (0/1/2), not probabilities        |

Step 4: Process Results
```r
# Remove genes with only one variant (can't do gene-based test)
exonic.nvargt1 <- subset(example.exonic.groups, n.variants > 1)

# Sort by SMMAT-E p-value (smallest first)
exonic.nvargt1 <- exonic.nvargt1[order(exonic.nvargt1$E.pval),]

# View top 5 genes
exonic.nvargt1[1:5,]
```
Output columns:
| Column                 | Description                           |
|------------------------|---------------------------------------|
| group                  | Gene name                             |
| n.variants             | Number of variants in gene            |
| freq.min/mean/max      | MAF statistics                        |
| B.score, B.var, B.pval | Burden test results                   |
| S.pval                 | SKAT p-value                          |
| O.pval                 | SKAT-O p-value                        |
| O.minp.rho             | ρ value that gave minimum p in SKAT-O |
| E.pval                 | SMMAT-E p-value                       |

Interpreting O.minp.rho:
- Near 1 → burden-like test was best
- Near 0 → SKAT-like test was best
- Intermediate → mixture

Step 5: Calculate cMAC
```r
# Sample size from null model
N <- length(model.0$id_include)  # 822

# cMAC = n_variants × mean_MAF × 2N
exonic.nvargt1$cMAC <- exonic.nvargt1$n.variants * 
                       exonic.nvargt1$freq.mean * 
                       N * 2

# How many genes have sufficient power?
table(exonic.nvargt1$cMAC < 10)
# FALSE  TRUE 
#   454    43
```

Step 6: Repeat for Non-synonymous Variants
```r
example.nonsyn.groups <- SMMAT(model.0, 
                               geno.file = "example.gds",
                               group.file = "example.smmat.nonsyn.groups",
                               # same parameters as before
                               )

nonsyn.nvargt1 <- subset(example.nonsyn.groups, n.variants > 1)
nonsyn.nvargt1 <- nonsyn.nvargt1[order(nonsyn.nvargt1$E.pval),]
nonsyn.nvargt1$cMAC <- nonsyn.nvargt1$n.variants * 
                       nonsyn.nvargt1$freq.mean * 
                       N * 2

table(nonsyn.nvargt1$cMAC < 10)
# FALSE TRUE 
#   377   75
Observation: More genes excluded with cMAC<10 when using only non-synonymous variants because they're rarer.
```

Step 7: Save Results
```r
write.csv(exonic.nvargt1, "exonic.nvargt1.csv", row.names=FALSE, quote=FALSE)
write.csv(nonsyn.nvargt1, "nonsyn.nvargt1.csv", row.names=FALSE, quote=FALSE)
```

#### 9.5 QQ Plot Script: qqplot_smmat.R
```r
args <- commandArgs(trailingOnly = TRUE)
myfile <- args[1]

# Function to compute expected -log10(p) from uniform distribution
get_expected <- function(xx) {
  yy <- sort(xx[!is.na(xx)])
  exp <- -log10(ppoints(length(yy)))
  obs <- -log10(yy)
  return(cbind(exp=exp, obs=obs))
}

# Read data
mydat <- read.csv(myfile, as.is=TRUE)

# Create plot
bitmap(paste0("qq.", myfile, ".jpeg"), type="jpeg", res=720)

# Plot all four tests with different colors
plot(get_expected(mydat$B.pval), col="cyan", type="p",
     xlab="-log10(exp(P))", ylab="-log10(obs(P))",
     xlim=c(0, max(-log10(mydat$B.pval), na.rm=TRUE)),
     ylim=c(0, max(-log10(mydat$B.pval), na.rm=TRUE)))
points(get_expected(mydat$S.pval), col="green")
points(get_expected(mydat$O.pval), col="blue")
points(get_expected(mydat$E.pval), col="red", lty=2)
abline(0, 1, col="black")
legend("topleft", 
       legend=c("SMMAT-B","SMMAT-S","SMMAT-O","SMMAT-E"),
       col=c("cyan","green","blue","red"), 
       lty=c(1,1,1,2), bty="n")
title(myfile)
dev.off()
```

Run:
```bash
Rscript qqplot_smmat.R exonic.nvargt1.csv
Rscript qqplot_smmat.R nonsyn.nvargt1.csv
```
Interpretation:
- Points along diagonal = no inflation
- Deviation at top right = possible true associations
- Different colors show which test is most significant

### 10. Example 2: TGEN Imputed Data (Chromosome 19)
#### 10.1 Data Files
Located in `/projectnb/bs859/data/tgen/annotated_imputed_vcfs/`
| File                      | Description                                 |
|---------------------------|---------------------------------------------|
| `anno1.chr19.vcf.gz `       | Imputed genotypes + ANNOVAR annotations     |
| `chr19.gds `                | GDS format (pre-converted)                  |
| `tgen.psam`                 | Sample file with FID, IID, sex, case status |
| `grm.rel`, `grm.rel.id`       | Genetic relationship matrix and IDs         |
| `chr19.exonic.smmat.groups` | All exonic variants per gene                |
| `chr19.nonsyn.smmat.groups` | Non-synonymous variants per gene            |
| `TGEN_pcs.txt`              | Principal components (from earlier)         |

#### 10.2 R Script: SMMAT_tgenChr19Example.R
Step 1: Load Data
```r
library(GMMAT)

# Read phenotype
pheno <- read.table("/projectnb/bs859/data/tgen/annotated_imputed_vcfs/tgen.psam",
                    header=TRUE, comment.char="")
colnames(pheno) <- c("FID","IID","sex","case")

# Read PCs
pcs <- read.table("/projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt",
                  header=TRUE, as.is=TRUE)

# Merge phenotype and PCs
pheno1 <- merge(pheno, pcs, by=c("FID","IID"), all.x=TRUE)
```

Step 2: Load Genetic Relationship Matrix
```r
# Read GRM
grm <- as.matrix(read.table("/projectnb/bs859/data/tgen/annotated_imputed_vcfs/grm.rel", 
                             header=FALSE))

# Read GRM IDs
grm.ids <- read.table("/projectnb/bs859/data/tgen/annotated_imputed_vcfs/grm.rel.id",
                       header=FALSE)

# Assign IDs to row/column names
dimnames(grm)[[1]] <- dimnames(grm)[[2]] <- grm.ids[,2]
Why GRM matters: TGEN samples may have cryptic relatedness. GRM accounts for this.
```
Step 3: Fit Null Model (Binary Trait)
```r
# case ~ covariates, with GRM random effect
# Note: case-1 because SMMAT expects 0/1, but data has 1/2
model.0 <- glmmkin(case-1 ~ PC6 + PC8, 
                   data = pheno1, 
                   id = "IID",
                   kins = grm, 
                   family = binomial("logit"))

# Check covariate significance
coef <- model.0$coef
sd <- sqrt(diag(model.0$cov))
coef.pval <- pchisq((coef/sd)^2, 1, lower.tail=FALSE)
coef.table <- cbind(coef, sd, coef.pval)
print(coef.table)

# Expected output:
#               coef        sd     coef.pval
# (Intercept) 0.5151321 0.05941159 4.299806e-18
# PC6        -6.2375138 2.32367994 7.267626e-03
# PC8        -5.0337197 2.19797293 2.201178e-02
```

Key differences from Example 1:
- `family=binomial("logit")` → logistic regression for binary trait
- `kins=grm` → includes random effect for relatedness
- `case-1` → converts 1/2 to 0/1

Step 4: Run SMMAT on Imputed Data
```r
chr19.exonic05 <- SMMAT(
  model.0,
  geno.file = "/projectnb/bs859/data/tgen/annotated_imputed_vcfs/chr19.gds",
  group.file = "chr19.exonic.smmat.groups",
  group.file.se = " ",
  MAF.range = c(1e-7, 0.05),
  MAF.weights.beta = c(1, 25),
  miss.cutoff = 1,
  missing.method = "impute2mean",
  method = "davies",
  tests = c("O", "E"),
  rho = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),
  use.minor.allele = FALSE,
  auto.flip = FALSE,
  Garbage.Collection = FALSE,
  is.dosage = TRUE,        # CRITICAL: imputed data = dosages
  ncores = 1,
  verbose = TRUE
)
```

What's different from Example 1:
- `is.dosage` = TRUE → uses dosage values (0-2) (continous), not hard calls (discrete)
- `use.minor.allele` = FALSE → using ALT allele as coded (common for imputed data)

Step 5: Process Results
```r
# Remove single-variant genes
chr19.exonic05.nvargt1 <- subset(chr19.exonic05, n.variants > 1)

# Sort by SMMAT-E p-value
chr19.exonic05.nvargt1 <- chr19.exonic05.nvargt1[order(chr19.exonic05.nvargt1$E.pval),]

# Calculate cMAC
N <- length(model.0$id_include)
chr19.exonic05.nvargt1$cMAC <- chr19.exonic05.nvargt1$n.variants * 
                                chr19.exonic05.nvargt1$freq.mean * 
                                N * 2

# View top 5 genes
chr19.exonic05.nvargt1[1:5,]

# Check how many genes have sufficient power
table(chr19.exonic05.nvargt1$cMAC >= 10)

# Calculate significance threshold
sig.level <- 0.05 / table(chr19.exonic05.nvargt1$cMAC >= 10)["TRUE"]

# Find significant genes
subset(chr19.exonic05.nvargt1, E.pval <= sig.level)
```

Expected top genes (from class):
```text
     group n.variants freq.mean    E.pval    cMAC
1270 TOMM40   12       0.00724   2.02e-05   213.2
1363 ZNF154   11       0.00566   2.20e-03   152.7
1016 PSG7     18       0.00184   2.47e-03    81.2
...
```
Note: TOMM40 is near APOE on chromosome 19 - biologically plausible for Alzheimer's.

Step 6: Save Results
```r
write.csv(chr19.exonic05.nvargt1, "chr19.exonic05.csv", 
          row.names=FALSE, quote=FALSE)
```

#### 10.3 QQ Plot for TGEN Results
```bash
Rscript qqplot_smmat.R chr19.exonic05.csv
```
What to look for:

Which test (color) shows smallest p-values?

Is there deviation at the top right (possible true associations)?

Is the overall inflation controlled?

### 11. Creating Group Files with bcftools/awk
The README and class5.sh show how group files were created:

#### 11.1 Extract Annotations from VCF
```bash
# Load required modules
module load htslib/1.18
module load bcftools/1.18

# Extract fields: Gene, CHROM, POS, REF, ALT, Func, ExonicFunc
bcftools query -f '%Gene.refGene %CHROM %POS %REF %ALT %Func.refGene %ExonicFunc.refGene\n' \
  /projectnb/bs859/data/tgen/annotated_imputed_vcfs/anno1.chr19.vcf.gz > temp1
```

#### 11.2 Filter for Exonic Variants
```bash
# Keep only lines where Func contains "exonic"
awk '$6~"exonic"{print $0}' temp1 > temp2
```

#### 11.3 Create Exonic Group File
```bash
# Format: Gene Chr Pos Ref Alt Weight
awk '{print $1,$2,$3,$4,$5,1}' temp2 > chr19.exonic.group
```
#### 11.4 Create Non-synonymous Group File
```bash
# Filter for functional variants
awk '$7~"nonsyn"||$7~"stop"||$7~"start"||$7~"frameshift"{print $1,$2,$3,$4,$5,1}' temp2 > chr19.nonsyn.group
```

What the patterns match:
- `nonsyn` = non-synonymous (amino acid changing)
- `stop` = stop-gain or stop-loss
- `start` = start-loss
- `frameshift` = frameshift indel

### 12. Key Takeaways and Best Practices
#### 12.1 When to Use Gene-Based Tests
- Rare variants (MAF < 1-5%) that can't be tested individually
- Biological units (genes) have prior plausibility
- You have functional annotations to prioritize variants
#### 12.2 Test Selection Guide
| Scenario                                | Recommended Test      |
|-----------------------------------------|-----------------------|
| Expect all variants in same direction   | Burden test (SMMAT-B) |
| Expect mixed directions or many nulls   | SKAT (SMMAT-S)        |
| Unknown/unpredictable                   | SMMAT-E or SKAT-O     |
| Want to combine multiple tests per gene | ACAT                  |

#### 12.3 Essential QC Steps
1. Filter imputation quality (R² > 0.3 for imputed data)
2. Set MAF threshold (usually 0.01-0.05 for rare variant tests)
3. Create biologically-informed groups (exonic, non-synonymous, etc.)
4. Remove genes with single variants (can't do set-based test)
5. Calculate cMAC and exclude underpowered genes
6. Apply multiple testing correction based on number of tests with cMAC ≥ 10
7. Check QQ plots for inflation
#### 12.4 Common Pitfalls
| Pitfall                                         | Consequence                        | Solution                    |
|-------------------------------------------------|------------------------------------|-----------------------------|
| Including common variants in rare variant tests | Dilutes signal                     | Set MAF.range appropriately |
| Not filtering by cMAC                           | Wastes multiple testing correction | Exclude cMAC < 10           |
| Using hard calls for imputed data               | Loses uncertainty information      | Set is.dosage=TRUE          |
| Ignoring relatedness                            | Inflated type I error              | Include GRM                 |
| Testing genes with 1 variant                    | Invalid test                       | Remove n.variants=1         |
| Not accounting for multiple tests               | False positives                    | Bonferroni on # tests w     |

#### 12.5 Multiple Testing Workflow
1. Run SMMAT on all genes
2. Filter to n.variants > 1
3. Calculate cMAC
4. Count M = number of genes with cMAC ≥ 10
5. Bonferroni threshold = 0.05 / M
6. Report genes with p < threshold AND cMAC ≥ 10

### 13. Homework Preview
For chromosome 19 (and chromosome 2 in homework), you'll:
1. Run gene-based tests with two different MAF thresholds
2. Try two different weighting schemes
3. Compare results between exonic and non-synonymous variant sets
4. Identify significant genes near
5. Interpret the ρ parameter (O.minp.rho) to understand genetic architecture

Key genes to watch:
- Chromosome 19: APOE region (TOMM40, APOE, APOC1)
- Chromosome 2: BIN1 (known Alzheimer's gene from IGAP)

### 14. Quick Reference: SMMAT Parameters
| Parameter        | Example        | Purpose                            |
|------------------|----------------|------------------------------------|
| `null.obj `        | model.0        | Null model from glmmkin            |
| `geno.file`        | "chr19.gds"    | Genotype file in GDS format        |
| `group.file`       | "groups.txt"   | File defining variant sets         |
| `MAF.range`        | c(1e-7, 0.05)  | MAF range to include               |
| `MAF.weights.beta` | c(1, 25)       | Beta distribution weights for SKAT |
| `tests`            | c("O", "E")    | Which tests to return              |
| `rho`              | grid of values | Correlation grid for SKAT-O        |
| `use.minor.allele` | TRUE/FALSE     | Whether to use minor allele        |
| `is.dosage`        | TRUE/FALSE     | Whether genotypes are dosages      |

With these notes, you should be able to:
- Understand why rare variant analysis requires different methods
- Explain the difference between burden tests and SKAT
- Run SMMAT on both sequenced and imputed data
- Process results, calculate cMAC, and apply multiple testing correction
- Create QQ plots to visualize results
- Interpret output including O.minp.rho to understand genetic architecture

