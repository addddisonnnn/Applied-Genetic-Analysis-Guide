## Week 6: Meta-Analysis for Genetic Studies
### BS859 Applied Genetic Analysis
### Feburary 25, 2026

### 1. Why Meta-Analysis in Genetics?

#### 1.1 The Power Problem in GWAS
Genome-wide association studies test millions of SNPs, requiring a very stringent significance threshold to control false positives:
- Genome-wide significance: p < 5 × 10⁻⁸

This stringent threshold means:
- A single study may not have enough power to detect associated SNPs
- True associations with small-to-moderate effects can be missed
- Solution: Combine information across multiple studies to increase sample size and power

1.2 Two Approaches to Combine Data
Option 1: Pool individual-level data
- Combine raw data from all studies, analyze as one large sample
- Advantages: Maximum flexibility, consistent analysis methods
- Disadvantages:
    - Confidentiality restrictions may prevent data sharing
    - May not be desirable due to heterogeneity:
        - Different study designs requiring different statistical methods
        - Different population backgrounds (diet, environment, etc.)
        - Different trait measurements or definitions

Option 2: Meta-analysis ← Focus of this week
- Analyze each study separately, then combine summary statistics
- Advantages:
    - Can combine studies even without access to individual-level data
    - Preserves study-specific analyses that account for design differences
    - Widely used and accepted in genetics

### 2. The Meta-Analysis Framework
#### 2.1 Basic Setup
We have k studies (i = 1 to k), each testing the same null hypothesis:
- H₀: β = 0 (no association between SNP and trait)
Each study i provides:
- Effect estimate β̂ᵢ (e.g., regression coefficient from linear or logistic regression)
- Standard error SE(β̂ᵢ)
- P-value pᵢ
- Direction of effect (+ or -)
Goal: Obtain a single pooled estimate and test of H₀ using all k studies.

Key assumption: Independence of study results (each subject contributes to only one study)

2.2 First Step: Harmonize Effect Measures
All studies must report a consistent effect measure:
- For continuous traits: β coefficient from linear regression
- For binary traits: log odds ratio from logistic regression

These can be combined across studies because they share the same interpretation: change in trait per copy of the effect allele.

### 3. Fixed Effects Meta-Analysis
#### 3.1 The Fixed Effects Assumption
Core assumption: All studies are estimating the same underlying population effect β
- If any study had an infinitely large sample size, its effect estimate would exactly equal β
- The only reason study estimates differ is random sampling error
Visual representation: All study estimates scatter around the same true value β, with spread determined by each study's precision.

#### 3.2 Inverse Variance Weighting
The most common fixed effects approach weights each study by its precision:

Weight for study i: wᵢ ∝ 1 / SE(β̂ᵢ)²

Interpretation:
- Larger studies → smaller SE → larger weight
- Smaller studies → larger SE → smaller weight

Scaled weights (sum to 1):
- wᵢ = [1/SE(β̂ᵢ)²] / [Σⱼ 1/SE(β̂ⱼ)²]

#### 3.3 Pooled Effect Estimate and Standard Error
Combined effect estimate:
- β̂ = Σ wᵢ β̂ᵢ
Standard error of combined estimate:
- SE(β̂) = √[1 / Σ 1/SE(β̂ᵢ)²]

#### 3.4 Test Statistic
- Z_metaβ = β̂ / SE(β̂)
- Under H₀: β = 0, Z_metaβ ~ N(0,1) (or Z_metaβ² ~ χ²₁)

### 4. Weighted Z-Score Approach
#### 4.1 When to Use This Method
Sometimes β estimates are not directly comparable across studies because:
- Traits measured on different scales (e.g., different assays)
- Different transformations applied (log-transformed vs. rank-normalized)
- Different instruments or technologies used

In these cases, we can combine test statistics instead of effect estimates.
#### 4.2 The Method
For each study i:
- Obtain p-value pᵢ from association test
- Convert to Z-score Zᵢ, with sign indicating direction of effect:
    - Positive Z = effect allele increases trait/risk
    - Negative Z = effect allele decreases trait/risk
- Weight by √nᵢ (square root of sample size)

Meta-analysis test statistic:
- Z_meta = (Σ √nᵢ Zᵢ) / √(Σ nᵢ)
- Under H₀: Z_meta ~ N(0,1)
Important: This method does NOT produce a pooled effect size estimate, only a p-value and direction.

#### 4.3 Comparison with Inverse Variance Weighting
From Lee et al. 2016:
| Method                     | Optimizes                | Produces                                   |
|----------------------------|--------------------------|--------------------------------------------|
| Inverse variance (Z_metaβ) | Likelihood function      | Pooled β estimate with minimum variance    |
| Weighted Z-score (Z_meta)  | Non-centrality parameter | Maximum shift of Z-score under alternative |

Equivalence: The two methods are equivalent when SE(β̂ᵢ)⁻¹ ∝ √nᵢ. This holds when allele frequencies are similar across studies.

### 5. Heterogeneity
#### 5.1 What is Heterogeneity?
Heterogeneity occurs when studies are not estimating the same underlying effect:
- True effect sizes differ across populations
- Different environmental backgrounds modify genetic effects
- Different phenotype definitions or measurements

Fixed effects assumption violated → need different approach.

#### 5.2 Testing for Heterogeneity: Cochran's Q Test
Q statistic: 
- Q = Σ [1/SE(β̂ᵢ)²] × (β̂ᵢ - β̂)²
- Where β̂ is the fixed effects combined estimate.

Interpretation:
- Q is a weighted sum of squared deviations from the pooled estimate
- Under H₀ of homogeneity, Q ~ χ² with (k-1) degrees of freedom
- Large Q → evidence of heterogeneity

Caveats:
- Low power with few studies
- Can be overly significant with many studies (100+)
- When significant with few studies, strong evidence of heterogeneity

#### 5.3 Measuring Heterogeneity: I²
I² = [Q - (k-1)] / Q
- Proportion of total variation across studies due to heterogeneity (not chance)
- Ranges from 0 to 1 (or 0% to 100%)
- Negative values set to 0
Guidelines:
- I² = 0%: No heterogeneity
- I² = 25%: Low heterogeneity
- I² = 50%: Moderate heterogeneity
- I² = 75%: High heterogeneity

### 6. Random Effects Meta-Analysis
#### 6.1 The Random Effects Model
When heterogeneity exists, we use a two-level model:

Level 1 (within-study): β̂ᵢ | βᵢ, sᵢ² ~ N(βᵢ, sᵢ²)
- sᵢ² = within-study variance (SE² from regression)

Level 2 (between-study): βᵢ ~ N(β, τ²)
- β = mean effect across studies
- τ² = between-study variance (heterogeneity)

#### 6.2 Random Effects Weights
Weights now incorporate both sources of variation:
- wᵢ* ∝ 1 / [SE(β̂ᵢ)² + τ̂²]
Where τ̂² is an estimate of between-study variance.

Scaled weights: wᵢ* = [1/(SEᵢ² + τ̂²)] / [Σⱼ 1/(SEⱼ² + τ̂²)]

#### 6.3 Pooled Estimate and Standard Error
Random effects combined estimate:
- β̂* = Σ wᵢ* β̂ᵢ

Standard error:
- SE(β̂*) = √[1 / Σ 1/(SEᵢ² + τ̂²)]

#### 6.4 Test Statistic
Z* = β̂* / SE(β̂*)

Key relationship: Since τ² ≥ 0, SE(β̂*) ≥ SE(β̂). Therefore:
- Z* ≤ Z almost always
- The more heterogeneity, the larger the random effects p-value compared to fixed effects

#### 6.5 Estimating τ²: DerSimonian and Laird Method
τ̂² = [Q - (k-1)] / [Σ wᵢ - (Σ wᵢ² / Σ wᵢ)]

Where wᵢ = 1/SE(β̂ᵢ)² (fixed effects weights)
- If Q < k-1, set τ̂² = 0 (use fixed effects)
- Otherwise, τ̂² > 0

### 7. Advanced Meta-Analysis Methods
#### 7.1 METASOFT (Han and Eskin)
Alternative Random Effects Model (RE2):
- Tests combined null: H₀: β = 0 AND τ = 0
- More powerful than standard random effects when heterogeneity exists
- Rejection implies either average effect ≠ 0 OR heterogeneity > 0 OR both

Binary Effects Model:
- Assumes effect is either present or absent in each study
- If present, effect size is similar across studies with effect
- Useful when heterogeneity comes from some studies having null effects

#### 7.2 MR-MEGA (Meta-Regression of Multi-Ethnic Genetic Association)
- Partitions heterogeneity into components due to ancestry vs. residual
- Detects SNP associations while quantifying ancestry-correlated heterogeneity
- Useful for multi-ethnic meta-analyses

### 8. Meta-Analysis of Gene-Based Tests
#### 8.1 Review: Score Statistics for Rare Variants
Recall from Week 5:

Score statistic for variant j: Sⱼ = Σᵢ Gᵢⱼ(yᵢ - μ̂ᵢ)

Between-variant relationship matrix: Φ = G'PG
- Measures LD structure among variants
- P is projection matrix accounting for estimated covariates

#### 8.2 Meta-Analysis Framework
For k studies, each provides:
- Allele frequency pⱼ for each variant
- Score statistic Sₖⱼ
- Covariance matrix Φₖ

Meta-analysis burden test:

Q_meta-burden = (Σⱼ Σₖ wₖⱼ Sₖⱼ)² ~ χ²₁

Meta-analysis SKAT:

Q_meta-SKAT = Σⱼ (Σₖ wₖⱼ Sₖⱼ)² ~ mixture of χ²

#### 8.3 Important Considerations
Weights: Should use combined study allele frequency for consistency across studies

Variant selection: Need pre-specified variant subsets (e.g., genes) with consistent definitions across studies

Software options:
- SMMAT.meta in GMMAT R package (linear/logistic, some approximations)
- Raremetal/Raremetalworker (quantitative traits only)

### 9. Special Issues in Genetic Meta-Analysis
#### 9.1 Choice of Effect Allele
The problem: Different studies may choose different alleles as the "effect" (coded) allele.
- P-value and SE are not affected by allele choice
- But sign of β flips depending on which allele is coded

Example:
- Study 1 codes C allele: β = +0.18 (C increases trait)
- Study 2 codes A allele: β = -0.18 (A decreases trait, equivalent to C increasing)

Solution: Before meta-analysis, align all studies to the same effect allele:
- Choose a reference allele for meta-analysis
- For studies using the opposite allele, flip the sign of β
- Update allele frequency to 1-freq if needed

#### 9.2 DNA Strand Issues
Complementary base pairs:
- A ↔ T
- C ↔ G

Easy cases (A/C, A/G, C/T, G/T): Strand can be inferred and flipped unambiguously
- Substitute T for A, G for C, C for G, A for T

Ambiguous cases (A/T, C/G): Cannot determine strand without additional information
- Need manufacturer's strand information
- Must document naming convention (forward/reverse, top/bottom)

Impact on meta-analysis: If strand is inconsistent across studies, alleles may be flipped, causing:
- Apparent allele frequency differences
- Sign inconsistencies
- Heterogeneity

### 10. METAL Software
#### 10.1 What is METAL?
METAL is a commonly used tool for GWAS meta-analysis:
- Fixed effects models (inverse variance or weighted Z-score)
- Heterogeneity testing (Q and I²)
- Automatic allele alignment across studies
- Strand checking for most SNP types

Documentation: http://genome.sph.umich.edu/wiki/METAL_Documentation

#### 10.2 METAL Control File Structure
METAL uses a command file specifying:
- Input files for each study
- Column names for required fields
- Analysis options

Basic commands:
| Command | Purpose                                      |
|---------|----------------------------------------------|
| MARKER  | Column with SNP identifier                   |
| ALLELE  | Columns with effect and non-effect alleles   |
| FREQ    | Column with effect allele frequency          |
| WEIGHT  | Column with sample size (for Z-score method) |
| EFFECT  | Column with β estimate                       |
| STDERR  | Column with standard error                   |
| PVAL    | Column with p-value                          |
| PROCESS | Input file to analyze                        |
| ANALYZE | Run the meta-analysis                        |
| OUTFILE | Output file prefix                           |

#### 10.3 Example Control File
```text
# First study
MARKER   SNP
WEIGHT   N
ALLELE   EFFECT_ALLELE NON_EFFECT_ALLELE
FREQ     EFFECT_ALLELE_FREQ
EFFECT   BETA
STDERR   SE
PVAL     P_VAL
PROCESS study1_results.txt

# Second study
MARKER   SNP
ALLELE   EFFECT_ALLELE NON_EFFECT_ALLELE
FREQ     FREQ_EFFECT
WEIGHT   N
EFFECT   BETA
STDERR   SE
PVAL     PVALUE
PROCESS study2_results.txt.gz

OUTFILE MetaResults .tbl
ANALYZE
```

#### 10.4 METAL Output Files
Main output file (.tbl):
| Column     | Description                              |
|------------|------------------------------------------|
| MarkerName | SNP identifier                           |
| Allele1    | Effect allele (coded)                    |
| Allele2    | Non-effect allele                        |
| Weight     | Sum of weights (e.g., total N)           |
| Zscore     | Combined Z-statistic                     |
| P-value    | Meta-analysis p-value                    |
| Direction  | Effect direction for each study (+ or -) |

Info file (.tbl.info):
- Describes columns in output
- Lists input files in order (crucial for interpreting Direction field)

#### 10.5 Genomic Control in METAL
What is genomic control?
- Adjustment for inflation in test statistics due to population stratification or other confounding
- Lambda (λ) = median(observed χ²) / 0.455
- λ > 1 indicates inflation

In METAL:
```text
# Option 1: Provide known λ from each study
GENOMICCONTROL 1.044   # for first study
GENOMICCONTROL 1.008   # for second study

# Option 2: Let METAL compute λ from the data
GENOMICCONTROL ON
```

Important: For full GWAS meta-analysis, use GENOMICCONTROL ON. For focused analyses (like the class example with only 3 regions), provide λ from published full GWAS.

#### 10.6 Inverse Variance Approach in METAL
To switch from default Z-score method to inverse variance:

```text
SCHEME STDERR
```
Required fields:
- STDERR column specified for each study
- EFFECT column with β estimates

Additional options for allele frequency tracking:

```text
AVERAGEFREQ ON
MINMAXFREQ ON
```

These provide:
- `Freq1`: Average effect allele frequency across studies
- `FreqSE`: Standard error of frequency (large values may indicate allele coding errors)
- `MinFreq`, `MaxFreq`: Range of frequencies across studies

#### 10.7 Heterogeneity Testing in METAL
Replace ANALYZE with:

```text
ANALYZE HETEROGENEITY
```
Output adds:
- HetISq: I² heterogeneity measure
- HetChiSq: Q statistic
- HetDf: Degrees of freedom
- HetPVal: P-value for heterogeneity test

### 11. Practical Example: Glucose Level Meta-Analysis

#### 11.1 The Three Studies
| Study    | File                        | Key Features                                   |
|----------|-----------------------------|------------------------------------------------|
| DGI      | `DGI_three_regions.txt`       | Alleles coded as 1=A,2=C,3=G,4=T               |
| FUSION   | `MAGIC_FUSION_Results.txt.gz` | Has STRAND column (+), sample size 1233        |
| Sardinia | `magic_SARDINIA.tbl`          | Has Rsq (imputation quality), trait glucose_ND |

#### 11.2 Top Finding: rs560887
Meta-analysis results (Z-score method):
- Zscore = -6.853
- P-value = 7.236 × 10⁻¹²
- Direction = `---` (negative effect in all three studies)

Study-specific results:
| Study    | Effect Allele  | β       | SE      | P-value    |
|----------|----------------|---------|---------|------------|
| FUSION   | C (coded as 4) | -0.139  | 0.044   | 0.00169    |
| DGI      | C (coded as 4) | -0.0457 | 0.03945 | 0.257      |
| Sardinia | C              | +0.180  | 0.028   | 1.36×10⁻¹⁰ |

Wait - Sardinia shows positive β? This appears inconsistent with the `---` direction.

Explanation: Sardinia reports β for the C allele as the effect allele. The positive β means C increases glucose. FUSION and DGI report negative β for the C allele? Actually, check FUSION and DGI carefully:
- FUSION: EFFECT_ALLELE = 4 (C), β = -0.139 → C decreases glucose?
- DGI: EFFECT_ALLELE = 4 (C), β = -0.0457 → C decreases glucose?
- Sardinia: AL1 = C, EFFECT = +0.180 → C increases glucose

This is actually consistent if we look at the allele frequencies:
- FUSION: FREQ_EFFECT = 0.314 (C allele frequency)
- DGI: EFFECT_ALLELE_FREQ = 0.298 (C allele frequency)
- Sardinia: FREQ1 = 0.627 (C allele frequency)

Wait - Sardinia's C allele frequency is much higher (0.627 vs 0.3). Something is off.

The key insight: In Sardinia, AL1 = C, AL2 = T, FREQ1 = 0.627. This means:
- If C is the effect allele, frequency = 0.627
- If T is the effect allele, frequency = 0.373

In FUSION and DGI, the C allele frequency is ~0.3, so T allele frequency is ~0.7.

This suggests: The glucose-raising allele is actually the T allele (frequency ~0.7 in most populations, but C in Sardinia if they're using opposite strand). This is a strand issue!

Resolution from the paper: rs560887 is an A/T SNP (strand ambiguous). Sardinia likely used the opposite strand, so their C is actually T on forward strand. When aligned, all studies show the T allele increases glucose.

#### 11.3 Effect of Genomic Control Correction
- Without GC: p = 7.236 × 10⁻¹²
- With GC: p = 2.219 × 10⁻¹¹

Why less significant? GC correction inflates standard errors based on λ, reducing test statistics. This is conservative but necessary when there's inflation.

#### 11.4 Heterogeneity at rs560887
With inverse variance + heterogeneity testing:
- HetPVal = 0.02536
- I² = 72.8%

Interpretation: Moderate-to-high heterogeneity. Looking at study-specific estimates:
- Sardinia: β = 0.180, SE = 0.028
- FUSION: β = -0.139, SE = 0.044 (after alignment, should be positive)
- DGI: β = -0.0457, SE = 0.039 (after alignment, should be positive)

The heterogeneity likely comes from:
1. Different allele frequencies across populations
2. Possible residual strand issues
3. True effect size differences across populations

### 12. Homework Preview: Uganda BMI Meta-Analysis
#### 12.1 The Four Studies
| Study  | Population             | N     | Description                              |
|--------|------------------------|-------|------------------------------------------|
| Uganda | Ugandan                | 6,197 | Uganda Genome Resource                   |
| DCC    | South African (Durban) | 1,478 | Durban Case Control (diabetes study)     |
| DDS    | South African (Durban) | 1,114 | Durban Diabetes Study (population-based) |
| AADM   | Multi-site African     | 5,187 | Nigeria, Ghana, Kenya                    |

Total N: up to 14,126 individuals

#### 12.2 Data Format
One file with all four studies: `Uganda_BMI_GWAS_MAF01.txt.gz`

Columns for each study:
- beta_uganda, se_uganda, pval_uganda, af_uganda, no_uganda
- beta_DCC, se_DCC, pval_DCC, af_DCC, no_DCC
- beta_DDS, se_DDS, pval_DDS, af_DDS, no_DDS
- beta_AADM, se_AADM, pval_AADM, af_AADM, no_AADM

SNP identifier: snpid in format CHR:POS:REF:ALT
- Example: 10:100000012:G:A

#### 12.3 Key Challenge
All four studies are in the same file. Your METAL control file must:
- Read the same file four times
- Extract different columns for each study
- Use the correct column names for each

This requires careful column specification for each PROCESS statement.

#### 12.4 Practical Considerations
- File size: 1.3 GB, 13.7 million markers
- Storage: Delete results files after homework to save space
- Manhattan plot: Need chromosome and position from snpid column
- Efficiency: Use methods that don't require loading full dataset into memory

### 13. Meta-Analysis Workflow Summary
#### 13.1 Before Meta-Analysis
1. Harmonize effect alleles across studies
2. Check strand for ambiguous SNPs
3. Verify allele frequencies are consistent with expectations
4. Generate QQ plots and Manhattan plots for each study
5. Calculate genomic control λ for each study (if full GWAS)

#### 13.2 During Meta-Analysis
1. Create METAL control file with correct column mappings
2. Run METAL and save log file
3. Check for warnings/errors in log
4. Review output for obvious problems

#### 13.3 After Meta-Analysis
1. Generate QQ plot for meta-analysis results
2. Create Manhattan plot to visualize genome-wide results
3. Calculate meta-analysis λ to check for inflation
4. Examine top signals for consistency across studies
5. Check heterogeneity at significant loci
6. Investigate any anomalies (frequency mismatches, direction conflicts)

### 14. Key Takeaways
#### 14.1 When to Use Each Method
| Scenario                               | Recommended Method               |
|----------------------------------------|----------------------------------|
| Comparable β estimates across studies  | Inverse variance (fixed effects) |
| Different trait scales/transformations | Weighted Z-score                 |
| Evidence of heterogeneity              | Random effects or METASOFT RE2   |
| Some studies may have null effects     | METASOFT binary effects          |
| Multi-ethnic samples                   | MR-MEGA                          |
| Rare variant gene-based tests          | SMMAT.meta, Raremetal            |

#### 14.2 Critical Quality Control Steps
1. Allele alignment: Ensure all studies use same effect allele
2. Strand checking: Resolve ambiguous A/T, C/G SNPs
3. Frequency comparison: Large frequency differences may indicate allele mismatch
4. Direction consistency: Check sign of effects across studies
5. Heterogeneity assessment: Identify loci with significant heterogeneity

#### 14.3 Common Pitfalls
| Pitfall                            | Consequence                    | Prevention                              |
|------------------------------------|--------------------------------|-----------------------------------------|
| Ignoring allele alignment          | Sign errors, null results      | Always align before meta-analysis       |
| Strand ambiguity in A/T, C/G SNPs  | Apparent frequency differences | Use manufacturer strand info            |
| Applying GC correction to non-GWAS | Incorrect λ estimates          | Only use GC ON for full GWAS            |
| Ignoring heterogeneity             | Overly significant results     | Check Q and I², consider random effects |
| Not checking frequencies           | Missed allele coding errors    | Use AVERAGEFREQ and MINMAXFREQ          |

### 15. METAL Command Reference
#### 15.1 Essential Commands
| Command | Purpose                       | Example                                |
|---------|-------------------------------|----------------------------------------|
| MARKER  | SNP identifier column         | MARKER SNP                             |
| ALLELE  | Effect and non-effect alleles | ALLELE EFFECT_ALLELE NON_EFFECT_ALLELE |
| FREQ    | Effect allele frequency       | FREQ AF_UGANDA                         |
| WEIGHT  | Sample size or weight         | WEIGHT N or WEIGHT NO_UGANDA           |
| EFFECT  | Beta coefficient              | EFFECT BETA or EFFECT BETA_UGANDA      |
| STDERR  | Standard error                | STDERR SE or STDERR SE_UGANDA          |
| PVAL    | P-value                       | PVAL PVALUE or PVAL PVAL_UGANDA        |
| PROCESS | Input file                    | PROCESS study_file.txt                 |
| ANALYZE | Run analysis                  | ANALYZE or ANALYZE HETEROGENEITY       |
| OUTFILE | Output prefix                 | OUTFILE MyMeta .tbl                    |

#### 15.2 Optional but Useful
| Command              | Purpose                                     |
|----------------------|---------------------------------------------|
| SCHEME STDERR        | Use inverse variance instead of Z-score     |
| GENOMICCONTROL ON    | Compute and apply GC from data              |
| GENOMICCONTROL 1.044 | Apply known GC λ                            |
| AVERAGEFREQ ON       | Calculate average frequency across studies  |
| MINMAXFREQ ON        | Report min and max frequency                |
| DEFAULTWEIGHT        | Set constant weight for all SNPs in a study |

#### 15.3 Example for Homework (Single File, Multiple Studies)
```text
# Study 1: Uganda
MARKER snpid
ALLELE CODED NONCODED
FREQ af_uganda
WEIGHT no_uganda
EFFECT beta_uganda
STDERR se_uganda
PVAL pval_uganda
PROCESS Uganda_BMI_GWAS_MAF01.txt.gz

# Study 2: DCC
MARKER snpid
ALLELE CODED NONCODED
FREQ af_DCC
WEIGHT no_DCC
EFFECT beta_DCC
STDERR se_DCC
PVAL pval_DCC
PROCESS Uganda_BMI_GWAS_MAF01.txt.gz

# Continue for DDS and AADM...

OUTFILE BMI_Meta .tbl
ANALYZE HETEROGENEITY
```

With these notes, you should be able to:
- Understand the rationale for meta-analysis in genetics
- Explain fixed effects vs. random effects models
- Implement both inverse variance and weighted Z-score approaches
- Recognize and handle allele alignment and strand issues
- Use METAL to perform meta-analysis
- Interpret heterogeneity statistics (Q, I²)
- Identify and troubleshoot common meta-analysis problems