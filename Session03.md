## Session 3: Common Varitant GWAS - Association Analysis & Mixed Models
### BS859 Applied Genetic Analysis
### Feburary 4, 2026
The main topics of this week's session are about (1) conducting genetic association tests using PLINK, (2) comparing different association models (unadjusted, covariate-adjusted, mixed models), (3) interpreting GWAS results including QQ plots and Manhattan plots, and (4) using R packages for mixed models and visualization.

### I. What We're Trying To Do This Week: The Big Picture

Imagine this: We have genetic data from 89 people. Some have a disease (cases), some don't (controls). We want to find which DNA differences (SNPs) are more common in people with the disease.

The challenge:

1. We're testing millions of positions in DNA
2. People might be related or from different populations
3. We need to separate real signals from random noise

This week's goal: Learn how to run these tests properly and understand the results.

### II. Core Concepts Explained Simply
**A. What's a SNP?**
- SNP = Single Nucleotide Polymorphism
- It's a single letter change in DNA at a specific position
- Example: At position 1,000,000 on chromosome 1:
    - Some people have "A"
    - Others have "G"
- We inherit one letter from each parent, so we have three possibilities:
    - AA (from both parents)
    - AG (A from one, G from other)
    - GG (from both parents)

**B. What's the "Additive Model"? (Super Important!)**
This is how we code for genetic data for analysis:

| Genotype | What it means | Additive coding
| --- | --- | --- | 
| AA | 0 copies of G | 0 | 
| AG | 1 copy of G | 1 |
| GG | 2 copies of G | 2 |

The assumption: Each G allele adds the same effect.
- If having one G increases disease risk by some amount...
- Having two Gs increases it by twice that amount

Why we use it: It's simple, works reasonably well, and everyone uses it so we can compare studies.

**C. What's Logistic Regression?**
When our outcome is yes/no (disease/no disease), we use logistic regression.

Simple version: It calculates how much more likely disease is for each G allele.
- Output: β (beta) = log odds ratio
- Odds Ratio = eᵝ (this is easier to understand)
- Example: β = 0.5 → OR = e⁰·⁵ ≈ 1.65
    - Means: Each G allele increases odds of disease by 65%

### III. The Problem: Population Structure & Confounding
**A. What Goes Wrong Without Adjustment
Real-world example:**
- Imagine we're studying diabetes
- Group A: Japanese ancestry, 30% have diabetes
- Group B: Chinese ancestry, 10% have diabetes
- A certain SNP: Common in Japanese (80%), rare in Chinese (20%)

What happens if we don't adjust:

- SNP looks strongly associated with diabetes!
- But really: It's just more common in Japanese, and diabetes is more common in Japanese
- The SNP doesn't cause diabetes - it's a false association

**B. How We Found This Problem (From Last Week)**
Last class we did PCA (Principal Components Analysis) and found:
- PC1 separates Chinese vs Japanese ancestry
- Cases and controls are NOT evenly distributed across populations
- This is BAD - it means population structure could mess up our results

**C. The Covariate File (`newcov.txt`)**
This file contains information to help fix the problem:

```text
FID     IID     PC1     PC2     ... PC10    pop
CH18526 NA18526 -0.1061 0.0520  ... -0.0396 1
CH18524 NA18524 -0.1113 -0.1509 ... -0.2010 1
```

- `PC1-PC10`: Genetic ancestry scores from PCA
- `pop`: Population group (1=CHB/Chinese, 2=JPT/Japanese)

### IV. Step-by-Step: Running Different Types of Analyses
**Step 1: Set Up Your Environment**
```bash
# Load needed software
module load R
module load plink/1.90b6.27

# Set path to data
export DATADIR=/projectnb/bs859/data/plink_tutorial/wgas3

# Go to your working directory
cd /projectnb/bs859/students/YOURUSERNAME/class03
```

**Step 2: Basic Association (No Adjustment)**
```bash
plink --bfile $DATADIR/wgas3 --logistic beta --ci .95 --out logistnoadj
```

What this command does:
- `--bfile $DATADIR/wgas3`: Use the genetic data files (wgas3.bed, wgas3.bim, wgas3.fam)
- `--logistic`: Run logistic regression (for yes/no disease)
- `beta`: Report β coefficients instead of odds ratios
- `--ci .95`: Calculate 95% confidence intervals
- -`-out logistnoadj`: Save output with this name

Output files created:
- `logistnoadj.assoc.logistic`: Main results
- `logistnoadj.log`: Log file with details

**Step 3: Look at the Output**
```bash
head logistnoadj.assoc.logistic
```
You'll see something like:

```text
CHR SNP         BP      A1 TEST NMISS BETA    SE     L95    U95    STAT   P
1   rs3094315   792429  G   ADD  88    0.723  0.5235 -0.3029 1.749  1.381  0.1672
```

What each column means:
- `CHR`: Chromosome (1-22, X, Y)
- `SNP`: SNP name/ID
- `BP`: Base pair position in DNA
- `A1`: Effect allele (the one we're counting)
- `TEST`: Type of test (ADD = additive model)
- `NMISS`: Number of people with data for this SNP
- `BETA`: β coefficient (log odds ratio per A1 allele)
- `SE`: Standard error (uncertainty in β estimate)
- `L95`, `U95`: Lower/Upper bounds of 95% confidence interval
- `STAT`: Test statistic (β/SE)
- `P`: P-value

**Step 4: Adjust for Population Group**
```bash
plink --bfile $DATADIR/wgas3 --covar $DATADIR/newcov.txt --covar-name pop --logistic beta --ci .95 --out logistadjpop
```
What's different:
- --covar $DATADIR/newcov.txt: Use covariate file
- --covar-name pop: Use only the "pop" column from covariate file

Look at output:

```bash
head logistadjpop.assoc.logistic
```
Important: Now you see TWO lines per SNP:

```text
CHR SNP         BP      A1 TEST   BETA    P
1   rs3094315   792429  G   ADD    0.7042  0.2719
1   rs3094315   792429  G   pop    2.762    3.941e-07
```
- First line: SNP effect (what we care about)
- Second line: Population effect (nuisance parameter)

**Step 5: Clean Up the Output**
We don't want the population lines cluttering our results:

Option A: Use `awk` to extract only SNP lines:

```bash
awk 'NR==1||$5=="ADD"{print $0}' logistadjpop.assoc.logistic > popadj.add.txt
```
Option B: Tell PLINK to hide covariate lines:
```bash
plink --bfile $DATADIR/wgas3 --covar $DATADIR/newcov.txt --covar-name pop --logistic beta hide-covar --ci 0.95 --out logistadjpop-hidecov
```

**Step 6: Adjust for Continuous Ancestry (PCs)**
Maybe better than using population groups:

```bash
plink --bfile $DATADIR/wgas3 --covar $DATADIR/newcov.txt --covar-name PC1,PC2 --logistic beta hide-covar --ci .95 --out logistadjPC12-hidecov
```
Why PCs might be better:
- Population groups are crude (Chinese vs Japanese)
- PCs capture continuous ancestry differences
- Some people might be mixed or in-between

### V. Advanced Method: Mixed Models with GMMAT
**Why We Need Something Even Better**

The problem with simple adjustment:
- People aren't just "in a population group" - they have varying degrees of relatedness
- Even within the same population, some people are more genetically similar
- We need to account for ALL genetic similarity, not just ancestry

**Step 1: Prepare Data for Mixed Models**

**A. Prune SNPs (Remove Highly Correlated Ones)**
```bash
plink --bfile $DATADIR/wgas3 --indep-pairwise 10000kb 1 0.2
```
What this does:
- Looks at SNPs that are close together on chromosomes
- If two SNPs are highly correlated (r² > 0.2), remove one
- Why: We want independent SNPs for estimating genetic similarity

Files created:
- `plink.prune.in`: SNPs to KEEP (not too correlated)
- `plink.prune.out`: SNPs to REMOVE (highly correlated)

**B. Create Genetic Relationship Matrix (GRM)**
```bash
plink --bfile $DATADIR/wgas3 --chr 1-22 --extract plink.prune.in --make-rel square --out grm
```
What this does:
- Creates a matrix showing how similar every pair of people is genetically
- Uses only the pruned SNPs (independent ones)

Files created:
- `grm.rel.id`: List of people in order (2 columns: FID, IID)
- `grm.rel`: The actual matrix (89×89 numbers)

**C. Look at the GRM**
```bash
# See who's in the matrix
head grm.rel.id

# See first 5×5 subset of the matrix
cut -f1-5 grm.rel | head -n 5
```
What the numbers mean:
- Diagonal (person with themselves): ~1.0
- Off-diagonal (two different people): Usually close to 0
- Higher numbers mean more genetically similar

**Step 2: Run Mixed Model Analysis in R**
The R script GMMAT.R does this for us:

```bash
Rscript --vanilla GMMAT.R
```
What happens inside GMMAT.R:

**Part 1: Read and Prepare Data**
```r
# Read phenotype data from fam file
pheno <- read.table("wgas3.fam", header=F)
colnames(pheno) <- c("FID","IID","fa","mo","sex","case")
# case: 1=unaffected, 2=affected (PLINK coding)

# Read covariates (PCs)
pcs <- read.table("newcov.txt", header=T, as.is=T)

# Merge them
pheno1 <- merge(pheno, pcs, by=c("FID","IID"), all.x=TRUE)

# Read GRM matrix
grm <- as.matrix(read.table("grm.rel", header=F))
grm.ids <- read.table("grm.rel.id", header=F)

# Label rows/columns of GRM with person IDs
dimnames(grm)[[1]] <- dimnames(grm)[[2]] <- grm.ids[,2]
```

**Part 2: Fit Models**
```r
# Important: GMMAT expects 1=case, 0=control
# PLINK uses: 2=case, 1=control, 0=missing
# So we do: case-1 to convert

# Model 1: No covariates, just GRM
model1.0 <- glmmkin(case-1 ~ 1, data=pheno1, id="IID", 
                   kins=grm, family=binomial("logit"))

# Model 2: With PC1 and PC2 as covariates
model2.0 <- glmmkin(case-1 ~ PC1+PC2, data=pheno1, id="IID",
                   kins=grm, family=binomial("logit"))
```
**Part 3: Run Tests**
```r
# Test all SNPs with both models
geno.file <- "/projectnb/bs859/data/plink_tutorial/wgas3/wgas3"
glmm.score(model1.0, infile=geno.file, outfile="test.glmm.score.nocov")
glmm.score(model2.0, infile=geno.file, outfile="test.glmm.score.PC1PC2cov")
Output files:

test.glmm.score.nocov: Mixed model without covariates

test.glmm.score.PC1PC2cov: Mixed model with PC1, PC2
```

**Step 3: Fix Column Names for Consistency**
GMMAT uses different column names than PLINK:
```bash
# Change "PVAL" to "P" and "POS" to "BP"
awk 'NR==1{$4="BP";$11="P"};{print $0}' test.glmm.score.nocov > glmm.score.nocov.txt
awk 'NR==1{$4="BP";$11="P"};{print $0}' test.glmm.score.PC1PC2cov > glmm.score.PC1PC2cov.txt
```

### VI. Comparing Results Across Methods
**A. Look at One SNP Across All Methods**
```bash
# Check rs11204005 in different analyses
awk '(NR==1||($2=="rs11204005")){print $0}' logistadjpop-hidecov.assoc.logistic
awk '(NR==1||($2=="rs11204005")){print $0}' logistadjPC12-hidecov.assoc.logistic
awk '(NR==1||($2=="rs11204005")){print $0}' glmm.score.nocov.txt
awk '(NR==1||($2=="rs11204005")){print $0}' glmm.score.PC1PC2cov.txt
```
What you'll see:

| Method | β/Score | P-value| 
| --- | --- | --- | 
| Pop-adjusted logistic	 | -2.708 | 	1.489×10⁻⁵
| PC-adjusted logistic | 	-2.648	 | 3.407×10⁻⁵
| Mixed model (no cov) | 	12.638	 | 4.291×10⁻⁵
| Mixed model (PC cov)	 | 12.500	 | 1.339×10⁻⁶
- Important note: Mixed model reports score statistics, not β values. Higher score = more significant.

**B. Understanding the Difference in Effect Directions**

Why is β negative in PLINK but score positive in GMMAT?

Allele coding difference:
- PLINK: Reports effect for A1 allele
- GMMAT: Reports score for A2 allele

If A1 increases risk, β is positive in PLINK.
If A2 increases risk, score is positive in GMMAT.

### VII. Evaluating Quality: QQ Plots
**A. What's a QQ Plot and Why Do We Need It?**

The problem: With millions of tests, some SNPs will have small p-values by chance. We need to check if our p-values look "right."

QQ Plot = Quantile-Quantile Plot

- Compares observed p-values vs expected p-values under null hypothesis (no true associations)
- X-axis: Expected -log₁₀(p) if no SNPs are associated
- Y-axis: Observed -log₁₀(p) from our analysis

**B. Making QQ Plots**
Using the Simple Script:
```bash
# For each analysis method
Rscript --vanilla qqplot.R logistnoadj.assoc.logistic crude ADD
Rscript --vanilla qqplot.R logistadjpop-hidecov.assoc.logistic popadj ADD
Rscript --vanilla qqplot.R logistadjPC12-hidecov.assoc.logistic PC12adj ADD
Rscript --vanilla qqplot.R glmm.score.nocov.txt "glmm nocovariates"
Rscript --vanilla qqplot.R glmm.score.PC1PC2cov.txt "glmm PC1PC2covariates"
```
Files created: `qq.crude.jpeg`, `qq.popadj.jpeg`, etc.

Using the Fancy Script (with confidence bands):
```bash
Rscript --vanilla qq_umich_gc.R logistnoadj.assoc.logistic noadj ADD
Rscript --vanilla qq_umich_gc.R logistadjpop-hidecov.assoc.logistic popadj ADD
Rscript --vanilla qq_umich_gc.R logistadjPC12-hidecov.assoc.logistic pcadj ADD
Rscript --vanilla qq_umich_gc.R glmm.score.nocov.txt glmmnocov
Rscript --vanilla qq_umich_gc.R glmm.score.PC1PC2cov.txt glmmPC1PC2cov
```

C. Interpreting QQ Plots
Perfect null distribution:

```text
Observed
  |
  |   •   (points follow the diagonal line)
  |  •
  | •
  |•
  |________________
      Expected
```
Inflation (λ > 1) - Too many small p-values:

```text
Observed
  |
  |       •
  |      •
  |     •
  |    •
  |   •
  |  •
  | •
  |•
  |________________
      Expected
```
Cause: Population structure, relatedness, other confounding

Deflation (λ < 1) - Too few small p-values:

```text
Observed
  |
  |•
  | •
  |  •
  |   •
  |    •
  |     •
  |      •
  |       •
  |________________
      Expected
```
Cause: Over-correction, wrong model

Good QQ plot with true associations:

```text
Observed
  |
  |           •
  |          •
  |         •
  |        •
  |       •
  |      •
  |     •
  |    •
  |   •
  |  •
  | •
  |•
  |________________
      Expected
```
- Most points follow diagonal (null SNPs)
- Some points deviate upward at the end (true associations)

**D. The λ (lambda) Value**
Shown on QQ plots:
- λ = 1.0: Perfect calibration
- λ > 1.0: Inflation (need better adjustment)
- λ < 1.0: Deflation (maybe over-adjusted)

Rule of thumb:
- λ < 1.05: Good
- 1.05 < λ < 1.10: Acceptable, but could be better
- λ > 1.10: Problem - need better adjustment

### VIII. Visualizing Results: Manhattan Plots
**A. What's a Manhattan Plot?**
Shows p-values across all chromosomes:

- X-axis: Chromosome and position
- Y-axis: -log₁₀(p-value)
- Horizontal line: Genome-wide significance threshold (usually -log₁₀(5×10⁻⁸) ≈ 7.3)

Why "Manhattan": Looks like a city skyline with buildings (peaks) of different heights.

**B. Making Manhattan Plots**
```bash
Rscript --vanilla gwaplot.R logistnoadj.assoc.logistic "unadjusted analysis" unadj_manhattan
Rscript --vanilla gwaplot.R logistadjpop-hidecov.assoc.logistic "Pop adj analysis" popadj_manhattan
Rscript --vanilla gwaplot.R logistadjPC12-hidecov.assoc.logistic "PC adj analysis" pc1pc2adj_manhattan
Rscript --vanilla gwaplot.R glmm.score.nocov.txt "GLMM no cov analysis" GLMM_nocov_manhattan
Rscript --vanilla gwaplot.R glmm.score.PC1PC2cov.txt "GLMM PC1 and PC2 adjusted analysis" GLMM_PC12_manhattan
Files created: unadj_manhattan.png, popadj_manhattan.png, etc.
```
**C. What to Look For**
1. Genome-wide significant hits: Points above the horizontal line
2. Patterns:
    - Clusters of signals in same region = likely real
    - Single isolated peak = might be false positive
3. Chromosome patterns: Some chromosomes might have more signals (biological or artifact?)

### IX. Multiple Testing: The 5×10⁻⁸ Threshold
**A. Why This Crazy-Small Number?**
Simple math: If we test 1,000,000 SNPs at α=0.05:
- Expected false positives: 1,000,000 × 0.05 = 50,000!
- That's way too many.

Bonferroni correction: Divide α by number of tests:
- α = 0.05 / 1,000,000 = 5×10⁻⁸
- This keeps overall false positive rate at 5%

In GWAS: We use p < 5×10⁻⁸ as "genome-wide significant"

**B. What About 1,000,000 Tests?**
Where does this number come from?

- Human genome has ~3 billion base pairs
- But SNPs are correlated (in linkage disequilibrium)
- Approximately 1,000,000 independent tests in Europeans
- Different populations have different numbers

### X. Putting It All Together: Analysis Pipeline
**Complete Workflow:**
1. Data QC (previous weeks) → Clean data
2. PCA (last week) → Check population structure
3. Association testing (this week):
    - Simple logistic (baseline)
    - Adjusted for covariates (population or PCs)
    - Mixed models (best for population structure)
4. Diagnostic checks:
    - QQ plots: Is λ close to 1?
    - Manhattan plots: Any interesting signals?

5. Interpretation:
    - Which SNPs are significant (p < 5×10⁻⁸)?
    - What are their effect sizes?
    - Do they make biological sense?

**Which Method to Choose?**
Based on our QQ plots:
1. If λ ≈ 1.0: Method is working well
2. If λ > 1.1: Need better adjustment
3. Choose method with λ closest to 1

For our data: Mixed model with PCs gave λ closest to 1.

### XI. Key Files and What They Contain
**Input Files:**
1. `wgas3.bed`/`.bim`/`.fam`: Genetic data in PLINK format
2. `newcov.txt`: Covariates (PCs, population)
3. `GMMAT.R`: R script for mixed models
4. `qqplot.R`, `gwaplot.R`, `qq_umich_gc.R`: Plotting scripts

**Output Files:**
1. `logist*.assoc.logistic`: PLINK association results
2. `test.glmm.score.*`: GMMAT mixed model results
3. `grm.rel`, `grm.rel.id`: Genetic relationship matrix
4. `plink.prune.in/out`: Lists of SNPs after pruning
4. `qq.*.jpeg`: QQ plots
5. `_manhattan.`: Manhattan plots

**Intermediate Files:**
1. `ibd.genome`: Pairwise IBD estimates (not used directly)
2. `popadj.add.txt`: Cleaned PLINK output (only SNP lines)

### XII. Common Issues and Solutions
**Problem 1: QQ Plot Shows Inflation (λ > 1.1)**
Solution:
- Try adjusting for more PCs
- Use mixed models instead of simple adjustment
- Check for relatedness in the sample

**Problem 2: Different Programs Give Different Results**
Check:
- Are alleles coded the same way? (A1 vs A2)
- Are the same people included? (NMISS column)
- Are covariates the same?

**Problem 3: No Genome-Wide Significant Hits**
Possible reasons:
- Sample size too small (very common in class examples)
- Trait not strongly genetic
- Need better phenotyping
- Not necessarily a problem - many real studies find nothing

**Problem 4: Too Many "Significant" Hits**
Check:
- QQ plot: Is λ very high?
- Population structure not properly adjusted
- Relatedness not accounted for
- P-value threshold too lenient

### XIII. For Your Homework
You'll need to:

1. Run similar analyses on different data
2. Create QQ and Manhattan plots
3. Compare different adjustment methods
4. Write a report explaining:
    - Which method worked best (based on λ)
    - Any interesting SNPs found
    - Limitations of the analysis

Tips:
- Save all your commands in a script (like class3.sh)
- Keep notes on what each command does
- Look at the log files for errors/warnings
- Always check QQ plots first!

### XIV. Summary: The Three Key Ideas
**1. Adjustment is Crucial**
Population structure causes false positives. You MUST adjust using:

- Population groups
- Principal Components (PCs)
- Mixed models (best)

**2. Visualization is Your Friend**
- QQ plots: Tells you if your analysis is working
- Manhattan plots: Shows where signals are
- Always look at these before interpreting results!

**3. Multiple Testing is Real**
- Testing millions of SNPs = need very small p-values
- Use p < 5×10⁻⁸ as threshold for "genome-wide significant"
- Everything else is "suggestive" at best

### XV. Quick Reference Commands
```bash
# Basic association
plink --bfile data --logistic beta --out results

# Adjust for population
plink --bfile data --covar covariates.txt --covar-name pop --logistic --out adj_results

# Make GRM for mixed models
plink --bfile data --indep-pairwise 10000kb 1 0.2
plink --bfile data --extract plink.prune.in --make-rel square --out grm

# QQ plot
Rscript qqplot.R results.assoc.logistic "My Analysis" ADD

# Manhattan plot
Rscript gwaplot.R results.assoc.logistic "My Results" myplot
```

Golden Rule: If your QQ plot doesn't look right, your results probably aren't right either. Fix the QQ plot first!

