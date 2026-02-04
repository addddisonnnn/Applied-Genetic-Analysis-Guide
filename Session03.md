## Session 3: Common Varitant GWAS - Association Analysis & Mixed Models
### BS859 Applied Genetic Analysis
### Feburary 4, 2026
The main topics of this week's session are about (1) conducting genetic association tests using PLINK, (2) comparing different association models (unadjusted, covariate-adjusted, mixed models), (3) interpreting GWAS results including QQ plots and Manhattan plots, and (4) using R packages for mixed models and visualization.

### I. Genetic Assocations Studies & GWAS
**A. GWAS Overview**
- Goal: Determine if a phenotype is associated with genetic variants
- Scale: Tests millions of variants individually (separate models)
- Models: Standard regression/logistic regression
- Software: Specialized tools like PLINK for handling large-scale data
- Challenge: Multiple testing requires careful interpretation

**B. Genotype Coding for Regression**
For a biallelic SNP with alleles 1 and 2:
1. Additive Model (Most Common in GWAS)
- Coding: Count number of "2" (effect) alleles
    - Genotype 11 → 0
    - Genotype 12 → 1
    - Genotype 22 → 2
- Test: 1 degree of freedom (H₀: β₁ = 0)
- Assumption: Linear effect per allele copy
2. Other Coding Schemes
- Dominant (2): 12 or 22 = 1, 11 = 0
- Recessive (2): 22 = 1, 11 or 12 = 0
- 2df Models: Additive + dominance deviation
3. Why Use Additive Model?
- Ubiquitous in GWAS
- Reasonable power regardless of true underlying model
- Simple interpretation: per-allele effect

C. Model Equations
For Quantitative Traits (Linear Regression):
​
  = additive-coded genotype (0,1,2)

For Binary Traits (Logistic Regression):
logit

​
 
Interpretation: 
exp(β1)= odds ratio per additional effect allele

-  Odds ratio used in a logistic model

Assumption: 

  with additive coding

### II. Covariate Adjustment
A. Why Adjust for Covariates?
1. Avoid confounding (especially population structure, by ancestry)
2. Increase precision (for linear models)
3. Required covariates: Any variable associated with BOTH genotype and phenotype
- Potential Confounders: Genetic Ancestry, Sec, and Age

B. Common Covariates in GWAS
- Population structure: Principal Components (PCs) of genotype matrix
- Sex (especially for sex chromosomes)
- Age (can proxy for ancestry/environment)

C. Model with Covariates

Quantitative: 

Yi =β<sub>0</sub> +β<sub>1</sub>X<sub>i</sub>​ + β<sub>2</sub>​cov<sub>i</sub> + e<sub>i</sub>
​

Binary: 
logit(p<sub>i</sub>) = β <sub>0</sub> + β<sub>1</sub>X<sub>i</sub>​ + β<sub>2</sub>​cov<sub>i</sub>

Important: For logistic regression, adding covariates can decrease precision/power. Include only necessary covariates.

D. Choosing PCs for Adjustment
Two approaches:
1. Adjust for top n PCs (n=1 to 10): Accounts for all structure
2. Adjust only for PCs associated with phenotype: Better for binary outcomes (fewer covariates)

**III. Population Structure Adjustment Methods**

A. Fixed Effects (PC Adjustment)
Include PCs as covariates in regression

Model: 

B. Mixed Effects Models
- Use genetic similarity as random effect
- Model
- K = kinship/genetic relationship matrix (GRM)

C. Software for Mixed Models
1. GMMAT R package (handles binary traits)
2. GCTA (quantitative traits)
3. Genesis R package
4. SAIGE/REGENIE (for large biobanks)

Important: Linear mixed models (LME) are for quantitative traits only. For binary traits, need specialized methods (GMMAT, Genesis).

**IV. PLINK Association Analysis**

A. Basic Commands

Logistic Regression (Binary Trait):
```bash
# Unadjusted
plink --bfile $DATADIR/wgas3 --logistic beta --ci 0.95 --out logistnoadj

# Adjusted for covariate
plink --bfile $DATADIR/wgas3 --covar $DATADIR/newcov.txt --covar-name pop --logistic beta --ci 0.95 --out logistadjpop

# Hide covariate output lines
plink --bfile $DATADIR/wgas3 --covar $DATADIR/newcov.txt --covar-name pop --logistic beta hide-covar --ci 0.95 --out logistadjpop-hidecov
```

Linear Regression (Quantitative Trait):
Replace `--logistic` with `--linear`

B. PLINK Output Interpretation
Key columns:
- `CHR`: Chromosome
- `SNP`: SNP identifier
- `BP`: Base pair position
- `A1`: Effect allele (coded allele)
- `A2`: Reference allele
- `BETA`: Regression coefficient (log OR for logistic)
- `SE`: Standard error
- `P`: P-value
- `TEST`: Type of test (ADD = additive genotype effect)

Important: By default, PLINK codes A1 as effect allele, A2 as reference. OR > 1 (β > 0) means A1 increases disease odds.

C. Parsing PLINK Output
With covariates, PLINK produces multiple lines per SNP:

- One line for SNP effect (TEST=ADD)
- One line per covariate effect

Extract only SNP lines:
```bash
awk 'NR==1||$5=="ADD"{print $0}' logistadjpop.assoc.logistic > popadj.add.txt
```

***V. Mixed Models with GMMAT***
A. Preparing Genetic Relationship Matrix (GRM)
1. Prune SNPs for LD
```bash
plink --bfile $DATADIR/wgas3 --indep-pairwise 10000kb 1 0.2
# Creates: plink.prune.in (SNPs to keep), plink.prune.out (SNPs to remove)
```

2. Compute GRM
```bash
plink --bfile $DATADIR/wgas3 --chr 1-22 --extract plink.prune.in --make-rel square --out grm
```

Output files:
- `grm.rel.id`: ID order for GRM matrix
- `grm.rel`: N×N genetic relationship matrix

B. GMMAT Analysis in R
See `GMMAT.R` script for details:

    1. Read and Prepare Data
```r
pheno <- read.table("wgas3.fam", header=F)
colnames(pheno) <- c("FID","IID","fa","mo","sex","case")
pcs <- read.table("newcov.txt", header=T, as.is=T)
pheno1 <- merge(pheno, pcs, by=c("FID","IID"), all.x=TRUE)

# Read GRM
grm <- as.matrix(read.table("grm.rel", header=F))
grm.ids <- read.table("grm.rel.id", header=F)
dimnames(grm)[[1]] <- dimnames(grm)[[2]] <- grm.ids[,2]
```

    2. Fit Null Models (without SNPs)
```r
# Note: recode case-1 because PLINK codes: 2=affected, 1=unaffected
# GMMAT expects: 1=affected, 0=unaffected

# No covariates
model1.0 <- glmmkin(case-1 ~ 1, data=pheno1, id="IID", kins=grm, family=binomial("logit"))

# With PC covariates
model2.0 <- glmmkin(case-1 ~ PC1+PC2, data=pheno1, id="IID", kins=grm, family=binomial("logit"))
```
    3. Run Score Test
```r
geno.file <- "/projectnb/bs859/data/plink_tutorial/wgas3/wgas3"
glmm.score(model1.0, infile=geno.file, outfile="test.glmm.score.nocov")
glmm.score(model2.0, infile=geno.file, outfile="test.glmm.score.PC1PC2cov")
```

C. Important Difference: Allele Coding
- PLINK: Beta/OR for A1 allele (A2 is reference)
- GMMAT: Score/effect for A2 allele (A1 is reference)

To make GMMAT output match PLINK column names:
```bash
awk 'NR==1{$4="BP";$11="P"};{print $0}' test.glmm.score.nocov > glmm.score.nocov.txt
```
**VI. Interpreting GWAS Results**
A. Multiple Testing Correction
- Problem
- FWER
- Bonferroni correction
- GWAS standard
B. Genomic Control (λ)
- Purpose: Adjust for modest population structure
- Assumption: Most SNPs are null, test statistics inflated equally
- Calculation:
- Adjustment

In R:
```r
pvalues <- mydat$P
chi1 <- qchisq(pvalues, 1, lower.tail=FALSE)
lambda <- median(chi1, na.rm=TRUE) / 0.455
gc_pvalues <- pchisq(chi1/lambda, 1, lower.tail=FALSE)
```

VII. Visualizing Results
A. QQ Plots
Purpose: Compare observed vs expected p-value distribution under null

Expected under null:
- Most SNPs are unassociated → p-values follow Uniform(0,1)
- Plot: observed -log₁₀(p) vs expected -log₁₀(p)

Interpretation:
- Points on diagonal: Good fit to null distribution
- Deviation along entire line: Systematic inflation (population structure, relatedness)
- Deviation only at tail: Potential true associations

Creating QQ Plots:
```bash
# Using provided scripts
Rscript --vanilla qqplot.R logistnoadj.assoc.logistic crude ADD
Rscript --vanilla qq_umich_gc.R logistnoadj.assoc.logistic noadj ADD
```

Key Parameters in Output:
- λ (lambda): Genomic inflation factor
    - λ ≈ 1: Good calibration
    - λ > 1: Test statistics inflated (possible population structure)
    - λ < 1: Test statistics deflated

B. Manhattan Plots

Purpose: Visualize p-values across genomic positions

Creating Manhattan Plots:
```bash
# Using provided script
Rscript --vanilla gwaplot.R logistnoadj.assoc.logistic "unadjusted analysis" unadj_manhattan

# Using qqman R package (prettier)
library(qqman)
manhattan(mydat, col=c("blue", "grey50"), highlight=top_snps)
```

Interpretation:
- X-axis: Genomic position (chromosomes)
- Y-axis: -log₁₀(p-value)
- Horizontal lines:
    - Suggestive: -log₁₀(1×10⁻⁵) = 5
    - Genome-wide: -log₁₀(5×10⁻⁸) = 7.3
- Peaks above line: Potentially significant associations

**VIII. Comparing Different Models**
A. Example: SNP rs11204005
Compare results from different adjustment methods:

Method	β/Score	SE	P-value
Logistic, PC1+PC2 adj	-2.648	0.639	3.407×10⁻⁵
Logistic, Population adj	-2.708	0.625	1.489×10⁻⁵
GLMM (no cov)	12.638	-	4.291×10⁻⁵
GLMM (PC1+PC2 cov)	12.500	-	1.339×10⁻⁶

Note: GLMM reports score statistics, not β. Also, allele coding differs (GMMAT gives effect for A2 allele).

B. Choosing the Best Model
Use QQ plots to evaluate:
- Which model produces λ closest to 1?
- Which has expected null distribution for majority of SNPs?
- Which shows deviation only at tail (suggesting true signals)?

For our example data:
- All models show reasonable calibration (λ ≈ 1)
- GLMM with PC covariates gives strongest signal for top SNP

IX. Key Takeaways
A. Analysis Pipeline
1. Quality Control (from previous classes)
2. Choose adjustment method based on population structure
3. Run association tests (PLINK for simple models, GMMAT for mixed models)
4. Correct for multiple testing (use 5×10⁻⁸ threshold)
5. Visualize results (QQ and Manhattan plots)
6. Interpret findings considering effect sizes, biological plausibility

B. Model Selection Considerations
1. For unrelated samples without population structure: Simple models sufficient
2. For population structure: Fixed-effect PC adjustment or mixed models
3. For related samples: Must use mixed models
4. For binary traits: Be cautious adding many covariates (reduces power)

C. Software Summary
- PLINK: Standard for basic association tests
- GMMAT: Mixed models for binary/quantitative traits
- R packages: qqman for visualization, various for specialized analyses
- For large biobanks: SAIGE, REGENIE

D. Interpretation Caveats
- Effect sizes in GWAS are typically small (OR < 1.1)
- Sample size needs to be large (10,000s to 100,000s) for adequate power
- Biological interpretation requires follow-up (replication, functional studies)
- Missing heritability: GWAS variants explain only part of genetic contribution

**X. Class Code Execution Summary**
The `class3.sh` script runs through:
1. Unadjusted logistic regression
2. Population-adjusted logistic regression
3. PC-adjusted logistic regression
4. Mixed model analysis with GMMAT
5. QQ plots for all models
6. Manhattan plots for all models

Key files produced:
- `*.assoc.logistic`: PLINK association results
- `test.glmm.score.*`: GMMAT mixed model results
- `qq.*.jpeg`: QQ plots
- `*_manhattan.*`: Manhattan plots

To run the full analysis:

```bash
module load R
module load plink/1.90b6.27
export DATADIR=/projectnb/bs859/data/plink_tutorial/wgas3
bash class3.sh
```

**XI. Homework Preparation**
For the TGEN data homework, you will:
1. Perform GWAS using different adjustment methods
2. Compare results across models
3. Create QQ and Manhattan plots
4. Interpret findings in context of multiple testing
5. Choose the most appropriate model based on diagnostic plots

Remember: Document all steps, save code and output, and be prepared to justify your model choices based on the data characteristics.


