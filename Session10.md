## Week 10: Polygenic (Risk) Scores
### BS859 Applied Genetic Analysis
### April 1, 2026

### 1. What is a Polygenic (Risk) Score?
A polygenic risk score (PRS) quantifies an individual’s genetic predisposition to a trait by summing the effects of many genetic variants (usually SNPs) across the genome.

Formula
PRSj=∑(above N)(below i=1) βi Gji 
- N = number of variants included in the score
- βi = weight for variant i (often from a discovery GWAS)
- Gji = genotype of individual j at variant i (typically coded as 0,1,2 copies of the effect allele)

**Why are PRS useful?
**- Risk stratification: Identify individuals at high genetic risk for a disease, potentially as high as rare monogenic mutations (Khera et al. 2018).
- Disease subtype differentiation: PRS can distinguish between subtypes (e.g., cardioembolic stroke vs other strokes).
- Causal inference: Use genetic risk as a tool to infer causality (Mendelian randomization).
- Clinical utility: Predict benefit from treatments (e.g., evolocumab in patients with high CAD genetic risk).

**History**
Early PRS work (2009) used schizophrenia GWAS; scores built from up to half of all SNPs were associated with disease in independent samples.

Showed shared polygenic basis between schizophrenia and bipolar disorder.

### 2. The PRS Analysis Process
A PRS is developed using two datasets:
1. Base dataset – a GWAS (summary statistics) that provides SNP effect estimates.
2. Target dataset – individual‑level genotype and phenotype data on which the PRS will be tested.

The typical pipeline:
- QC in both datasets (remove SNPs present only in one, remove ambiguous A/T or C/G SNPs, harmonize alleles).
- PRS calculation using methods that select SNPs and estimate weights (clumping+thresholding, penalized regression, Bayesian).
-Validation in an independent sample (or cross‑validation if no independent sample).
- Interpretation (association testing, variance explained, odds ratios, etc.).

### 3. Methods for Building a PRS
#### 3.1 Clumping and Thresholding (C&T)
This is the simplest and most common method.

- Thresholding: Keep only SNPs with a GWAS 
p-value below a certain threshold T (e.g., T=1×10−6,1×10−4,0.01,…,0.5
- Clumping: Remove SNPs in high LD (e.g.,  r2 >0.1) within a window (e.g., 250 kb) by keeping the SNP with the smallest p-value in each LD block.
- Optimization: Try many T values and optionally different LD clumping parameters. Choose the combination that yields the strongest association (smallest p-value) with the target phenotype.
- LD reference: LD can be estimated from a reference panel (e.g., 1000 Genomes) or from the target data.

How to assess the “best” score?
- Model R2 (variance explained by the PRS)
- Area under the curve (AUC) for binary traits
- p-value of the PRS association

Accounting for multiple testing
- The p-value from the best threshold does not account for the many thresholds tested.
- Bonferroni is too conservative because scores across thresholds are highly correlated.
- Permutation is used to obtain an empirical p-value:
    1. Permute the phenotype in the target sample (breaks any true association).
    2. Re‑run the entire C&T optimization and record the best p-value.
    3. Repeat many times (e.g., 10,000).
    4. Empirical p-value = (# permutations with best p-value < observed p-value + 1) / (total permutations + 1).
#### 3.2 Penalized Regression Methods
- Use Lasso, elastic net, or Bayesian approaches to shrink effect estimates.
- Account for LD by converting marginal effects to approximate joint effects using a reference LD matrix.
- Often improve prediction slightly compared to C&T, but are computationally more intensive.
#### 3.3 Bayesian Methods (e.g., LDpred, PRS‑CS)
These place a prior on the joint effect sizes and incorporate LD structure.
- Example: PRS‑CS uses a global‑local scale mixture of normals prior:
βj​ ∼N(0, σ2/N ϕψ j), ψj ∼g
    - One global shrinkage parameter ϕ applied to all variants.
    - A local, marker‑specific parameter ψj modifies shrinkage per SNP.
    - This shrinks small effects toward zero while preserving larger ones.

#### 3.4 Which Method is Best?
- C&T is simplest and most common.
- Penalized/Bayesian methods can yield small improvements in prediction accuracy.
- The choice often depends on sample size, computational resources, and whether you have a well‑matched LD reference.

### 4. Example: PRS for Alzheimer Disease (AD) with PRSice
#### 4.1 Data
- Base GWAS: IGAP AD meta‑analysis (Kunkle et al. 2019). Summary statistics include SNP, effect allele, beta, SE, p-value.
- Target data: TGEN (individual‑level genotypes and AD status).
- Covariates: first few principal components (PCs) to adjust for population structure.

#### 4.2 PRSice Command (simplified)
```bash
Rscript $SCC_PRSICE_BIN/PRSice.R \
    --base Kunkle_etal_Stage1_results2019.txt.gz \
    --target TGEN_cleaned \
    --stat Beta \
    --snp MarkerName \
    --A1 Effect_allele \
    --A2 Non_Effect_allele \
    --pvalue Pvalue \
    --binary-target T \
    --cov-file TGEN_pcs.txt \
    --cov-col PC6,PC8 \
    --clump-r2 0.10 \
    --out IGAP_base_TGEN_target
```
- `--stat Beta`: the column in base containing effect size.
- `--beta`: tells PRSice that effect sizes are betas (not OR).
- `--clump-r2`: LD threshold for clumping (default 0.1).
- `--clump-kb`: window size (default 250 kb).
- `--binary-target T`: outcome is binary.

#### 4.3 Output Files
- `.log`: screen output, shows commands and default parameters.
- `.mismatch`: SNPs omitted due to allele mismatches.
- `.valid`: SNPs that passed QC (useful to extract for re‑runs).
- `.summary`: information for the best model (across thresholds).
- .`prsice`: results for every threshold tested.
- Bar plot and high‑resolution plot (p‑value vs threshold, R2 vs threshold).
- Quantile plot (optional) – divides target sample into quantiles by PRS and plots odds ratios.

#### 4.4 Interpreting the .summary file

| Column      | Meaning                                                             |
|-------------|---------------------------------------------------------------------|
| Threshold   | Best pp-value threshold (e.g., 0.447)                               |
| PRS.R2      | Variance explained by PRS alone (after covariate adjustment)        |
| Full.R2     | Variance explained by full model (covariates + PRS)                 |
| Null.R2     | Variance explained by covariates alone                              |
| Coefficient | Regression coefficient of PRS (positive = higher PRS → higher risk) |
| P           | pp-value for PRS association                                        |
| Num_SNP     | Number of SNPs in the best PRS                                      |

#### 4.5 The .prsice File
- One row per tested threshold.
- Shows threshold, R2 , coefficient, SE, number of SNPs.
- Example: 5995 thresholds tested in the class example.
- The R2 in this file is the additional variance explained by the PRS (Full.R2 – Null.R2). However, this is an in‑sample estimate and may be overfitted – an independent validation set is needed for unbiased estimate.

4.6 Permutation for Empirical 
p-value `Add --perm 1000 --seed 542386` to the command.
- Empirical p-value = (# permuted best p < observed best p + 1) / (N_perm + 1).
- In class example, with observed 
p=5.2×10^−33, no permutation gave a smaller p=1/1001≈0.001.
- Permutation is most useful when the uncorrected p is near significance.

#### 4.7 Quantile Plots
- Divide target into quantiles based on PRS.
- Regress phenotype on quantile indicator, with a reference quantile (e.g., lowest).
- Plot OR per quantile to check for non‑linearity (e.g., sharp increase in top quantile).

### 5. Example: Using a Different Base Phenotype – Cognitive Function
- Base: GWAS of general cognitive function (Davies et al. 2018).
- Target: still AD in TGEN.
- Goal: test if cognitive function genetic variants also influence AD risk.
*Differences in base file:*
- Provided Z‑score, not beta.
- PRSice can use Z‑score with --stat Zscore --beta (the sign of Z is used as the weight).
- The best PRS explained only 0.08% of variance in AD and was not significant (p=0.394).
- Interpretation: either cognitive function variants do not affect AD, or the sample size/power was insufficient to detect a small effect.

### 6. Limitations of PRS
#### 6.1 Population Specificity & Health Disparities
- Most GWAS are European‑ancestry → PRS built from them perform poorly in non‑European populations (Martin et al. 2019).
- Reasons:
    - Differences in LD structure and allele frequencies across populations.
    - Causal variants may be common in one population but rare in another.
    - Gene‑by‑environment and gene‑by‑gene interactions may differ.
- Consequence: Clinical use of current PRS could worsen health disparities.

#### 6.2 Upper Bound of Prediction
- The maximum variance a PRS can explain is the SNP‑based heritability (h^2 SNP) estimated by LD score regression.
- For traits with low heritability, PRS will have limited predictive power.

#### 6.3 Need for Independent Validation
- The  R2 and p-value from the same data used to choose the optimal threshold are over‑optimistic.
- Gold standard: develop PRS in base + target (or training set), then evaluate in a completely independent validation set.
- If no independent sample, use cross‑validation within the target.

### 7. Homework Preparation
#### 7.1 Tuning Clumping Parameters
You will modify the PRSice command to change --clump-r2 and --clump-kb and observe the effect on:
- Number of SNPs after clumping (before thresholding).
- Number of SNPs in the optimal score.
- Proportion of variance explained.
- Optimal p-value threshold.
- p-value of best PRS.

#### 7.2 Hippocampal Volume PRS
- Base: hippocampal volume GWAS (vanderMeer et al.).
- Target: TGEN AD.
- Important: Use default PRSice settings except `--clump-kb 500`.
- Include -`-perm 1000 --seed 1443`.
- You will answer:
    - Is the best PRS associated with AD?
    - Present plots, optimal threshold, R2 , number of SNPs.
    - Interpret the coefficient (direction of effect).
    - Sample size adequacy and alternative approaches.

### 8. Key Concepts for Post‑Class Quiz
#### 8.1 LD Score Regression and PRS
- LD score regression estimates SNP heritability from GWAS summary statistics using LD information.
- This heritability is an upper bound on the variance that a PRS can explain.
- With small discovery sample size, LD score regression can overestimate heritability relative to achievable PRS performance.
- LD score regression does not require individual‑level genotype data; it uses summary statistics and an LD reference.

#### 8.2 Cross‑Validation vs. Independent Validation
- If no independent dataset is available, K‑fold cross‑validation (e.g., 10‑fold) within the target data gives a less biased estimate of out‑of‑sample performance.
- Reporting the best‑performing model from the training data is not unbiased.

#### 8.3 Permutation in PRSice
- Empirical p-value = (count of permuted best p < observed best p + 1) / (N_perm + 1).
- If observed best p is extremely small, it may be smaller than any permuted best p → empirical p = 1/(N_perm+1).
- The procedure does not become unreliable for very small p; it simply cannot give a p lower than 1/(N_perm+1).

#### 8.4 PRSice Steps (C&T)
- Remove SNPs present in only one dataset.
- Exclude strand‑ambiguous SNPs.
- Perform clumping: within each LD block, keep the SNP with the smallest p-value.
- Impute missing genotype or summary statistics? No, PRSice does not impute missing values in the base GWAS; it simply omits SNPs that are not present or mismatched.

#### 8.5 Population Structure and PRS
- Population structure bias can be amplified in PRS because many null variants may be included.
- If base and target are from similar populations, adjusting for principal components is usually adequate.
- PRS developed in one ancestry do not transfer well to others, even with PC adjustment, because of differences in LD and allele frequencies.
- LD differences do impact transferability; they are a major reason for poor performance across ancestries.

#### 8.6 Bayesian PRS Methods
- They apply prior distributions to SNP effect sizes to induce shrinkage and may model mixtures of causal and non‑causal variants.
- They explicitly incorporate LD (via an LD matrix or reference panel).
- They do not require individual‑level genotype data from the base GWAS; they work with summary statistics and an LD reference.
- They typically do not yield identical results to the best C&T model; they often give slightly better prediction.

### 9. Additional Commands and Notes from Class Materials
#### 9.1 Checking Summary Statistics Headers
```bash
zcat $IGAPDIR/Kunkle_etal_Stage1_results2019.txt.gz | head
```
#### 9.2 Running LD Score Regression for Genetic Correlation
```bash
# munge summary stats for cognitive GWAS
munge_sumstats.py \
    --sumstats Davies2018_OPEN_DATASET_summary_results.txt \
    --snp MarkerName \
    --N 250000 \
    --a1 Effect_allele \
    --a2 Other_allele \
    --signed-sumstats Zscore,0 \
    --merge-alleles w_hm3.snplist \
    --out Cog

# compute heritability of cognitive function
ldsc.py \
    --h2 Cog.sumstats.gz \
    --ref-ld-chr eur_w_ld_chr/ \
    --w-ld-chr eur_w_ld_chr/ \
    --out Cog.h2

# genetic correlation between AD and cognitive function
ldsc.py \
    --rg AD.sumstats.gz,Cog.sumstats.gz \
    --ref-ld-chr eur_w_ld_chr/ \
    --w-ld-chr eur_w_ld_chr/ \
    --out AD_COG_rg
```
#### 9.3 PRSice with Hippocampal Volume GWAS (Homework)
- Base file: `/projectnb/bs859/data/hippo/Meta_Whole_hippocampus.txt.gz`
- Headers: CHR, POS, SNP, A1, A2, N, P, Z
- Use `--stat Zscore --beta` because Z is provided.
- Add `--clump-kb 500` (override default 250).
- Include `--perm 1000 --seed 1443`.
- The `N` column is not sample size – ignore.

#### 10. Summary
- PRS aggregate many small‑effect variants to quantify genetic risk.
- C&T is the standard method; penalized/Bayesian approaches can slightly improve prediction.
- Always validate PRS in independent samples or use cross‑validation.
- Population differences (LD, allele frequencies) severely limit PRS transferability, raising concerns about health equity.
- PRSice is a user‑friendly tool for C&T PRS analysis, providing extensive output and permutation‑based significance testing.

Use these notes to complete the homework and answer the quiz questions. For the homework, pay close attention to how changing clumping parameters affects the number of SNPs retained and the performance of the PRS, and be prepared to interpret results from the hippocampal volume GWAS.