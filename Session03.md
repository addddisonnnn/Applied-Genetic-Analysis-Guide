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