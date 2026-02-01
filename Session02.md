## Session 2: Population Structure Detection
### BS859 Applied Genetic Analysis
### January 28, 2026
The main topics of this week's session are about (1) what population structure studies are, (2) how do we know if we have it, (3) what to do if we detect it, and (4) exploring one method (Eigensoft smartpca) for identifying structure and accounting for structure in association analyses

## Part 1: Lecture content - Theory & Concepts
1. **What is Population Structure?**
- Definition: When a population isn’t randomly mating → there are subgroups with different genetic backgrounds
- Types:
    - Discrete: Clear groups (e.g., Yoruba vs. Han Chinese)
    - Admixture: Mixed ancestry (e.g., African American)
    - Continuous: Gradient across geography (e.g., North to South Europe)
    - Fine-scale: Small differences within a region (e.g., Irish vs. English)

**Why does it happen?** Geography, culture, historical migrations, assortative mating.

2. **Why Do We Care?**

    Problem: SPURIOUS ASSOCIATIONS
If:
    1. Genotype frequencies differ across subpopulations, AND
    2. Phenotype/trait differs across subpopulations
→ You get a false genetic association.

    Real Example from Lecture: LCT Gene & Height
- North Europeans: taller, more lactose tolerant
- South Europeans: shorter, more lactose intolerant
- LCT gene varies across Europe
- Result: False association between LCT SNPs and height
- Solution: Account for population structure → association disappears

**Key Insight:** Structure only matters if BOTH genotype AND phenotype differ across groups.

3. **How to Detect Structure?**
    
    A. Heterozygote Deficit (F-statistic)
        
    - HWE says: In a random-mating population: p² + 2pq + q² = 1
    - With structure: Combined population violates HWE → fewer heterozygotes than expected
    - Formula: F = 1 – (Het_observed / Het_expected)
        - F > 0: Heterozygote deficit → suggests structure
        - F = 0: No deficit
        - F < 0: Excess heterozygotes (could be contamination)
    - Example from slides:
        - Subpop1: p=1.0 (all AA), Subpop2: p=0.0 (all aa)
        - Combined: No heterozygotes! F = 1.0

    B. Principal Components Analysis (PCA) – MAIN METHOD
    - What it does: Finds major axes of genetic variation
    - Input: Genotype matrix (individuals × SNPs)
    - Math: Normalize genotypes → covariance matrix → eigenvectors = PCs
    - PC1 explains most variation, PC2 second most, etc.

    Why PCA for genetics?
    - Large eigenvalues = population structure
    - PCs can be used as covariates to “adjust” for ancestry

4. **Reference Populations: HapMap & 1000 Genomes**
- HapMap: 4 populations (YRI, CEU, CHB, JPT) – genotyped SNPs
- 1000 Genomes: Extended HapMap – whole genome sequences
- Use: Merge your data with reference populations to see where your samples cluster
- Expectation: European-ancestry samples should cluster with CEU, etc.

5. **Accounting for Structure in Analysis**
    1. Include top PCs as covariates:
    - Phenotype ~ Genotype + PC1 + PC2 + ... + PC10
    2. How many PCs? Usually top 10, or those significantly associated with phenotype
    3. Random effects models: Use genetic relationship matrix (GRM) – handles both structure and relatedness

## Part 2: Coding Session
### Module 1: PCA on HapMap Example Data

**Step 1: Set up and Load Data**
```bash
# Load necessary software
module load eigensoft
module load plink
module load R

# Copy class files
cp /projectnb/bs859/materials/class02/* .
```
Files you get:
- `test.par` – parameter file for smartpca
- `class2.sh` – main script
- `plotPCs.R` – R script for plotting
- `hapmap.par` – for later

Step 2: Understand the Parameter File (`test.par`)
```bash
genotypename:    /projectnb/bs859/data/plink_tutorial/wgas2/wgas2_pruned.bed
snpname:         /projectnb/bs859/data/plink_tutorial/wgas2/wgas2_pruned.bim
indivname:       /projectnb/bs859/data/plink_tutorial/wgas2/wgas2_pruned.fam
evecoutname:     test.evec      # Output: eigenvectors (PCs)
evaloutname:     test.eval      # Output: eigenvalues
altnormstyle:    NO             # Normalization method
numoutevec:      10             # Save first 10 PCs
numoutlieriter:  0              # Don't remove outliers automatically
```
What this means: Tells smartpca where your data is and what to output.

**Step 3: Run PCA**
```bash
smartpca -p test.par > test.out
```
Check the output:
- `test.out`– log file (ALWAYS READ THIS!)
- `test.evec` – PCs for each individual
- `test.eval` – eigenvalues

**Step 4: Examine Results**
From `test.out`:
```
eigenvector 1:means
Control -0.071
Case 0.060
p-value eigenvector_1_Control_Case_ 1.43468e-10 +++
```
Interpretation: PC1 is significantly different between cases and controls (p < 0.001). This means population structure is confounded with case status!

Why? The data is from HapMap CHB (Chinese) and JPT (Japanese). Cases were simulated more from JPT, controls from CHB.

**Step 5: Plot the PCs**
```bash
Rscript --vanilla plotPCs.R test.evec 1 2 10
```
What this does: Runs the R script to plot PC1 (x) vs PC2 (y) from 10 PCs.

The R script (plotPCs.R):
- Reads test.evec
- Plots specified PCs
- Colors points by phenotype (case/control)
- Saves as JPEG

Key insight from plot: Cases and controls separate along PC1 → structure matters!

Step 6: Re-plot by Population (not case status)
```bash
# Add population info (CH or JA) based on sample ID
awk '{print $0, substr($1,1,2)}' test.evec > test.evec2

# Re-plot
Rscript --vanilla plotPCs.R test.evec2 1 2 10
```
What this shows: The separation is actually by population (CHB vs JPT), not disease. This is why we need to adjust for PCs!
### Module 2: PCA with TGEN + HapMap (More Complex)

**Step 1: Remove Problematic SNPs**
Problem: Strand issues with A/T and G/C SNPs
- DNA has two strands (complementary: A↔T, C↔G)
- If a SNP is A/T or C/G, you can't tell which strand was genotyped
- Causes mismatches when merging datasets
```bash
# Find all A/T and G/C SNPs
awk '($5=="A"&&$6=="T")||($5=="T"&&$6=="A")||($5=="G"&&$6=="C")||($5=="C"&&$6=="G"){print $2}' TGEN_cleaned.bim > atcgSNPs_omit.txt

# Remove them
plink --bfile TGEN_cleaned --exclude atcgSNPs_omit.txt --make-bed --out TGEN_cleaned1
```
Count: 48,561 SNPs removed → 259,901 remain

**Step 2: Merge with HapMap CEU Data**
```bash
# First attempt (will fail due to strand mismatches)
plink --bfile TGEN_cleaned1 --bmerge ceu_tgen --geno 0.03 --make-bed --out merge1
```
Error: "44,808 variants with 3+ alleles present"

Translation: Strand mismatches for 44,808 SNPs

**Step 3: Fix Strand Issues**
```bash
# Flip strands for mismatched SNPs
plink --bfile TGEN_cleaned1 --flip merge1-merge.missnp --make-bed --out TGEN_cleaned2

# Try merge again
plink --bfile TGEN_cleaned2 --bmerge ceu_tgen --make-bed --out merge2
```
Still errors: 7 SNPs still problematic → just exclude them

```bash
# Exclude remaining problematic SNPs
plink --bfile TGEN_cleaned2 --exclude merge2-merge.missnp --make-bed --out TGEN_cleaned3

# Keep only SNPs present in TGEN data
cut -f2 TGEN_cleaned3.bim > keepsnps.txt
```
**Step 4: Merge All HapMap Populations**
```bash
# Merge TGEN + CEU
plink --bfile TGEN_cleaned3 --bmerge ceu_tgen --extract keepsnps.txt --make-bed --out merge3

# Add CHB+JPT
plink --bfile merge3 --bmerge chbjpt_tgen --extract keepsnps.txt --geno 0.03 --make-bed --out merge4

# Add YRI
plink --bfile merge4 --bmerge yri_tgen --extract keepsnps.txt --geno 0.03 --make-bed --out merge5
```
Final dataset: 1,447 individuals (TGEN + CEU + CHB/JPT + YRI), 189,640 SNPs

**Step 5: LD Prune for PCA**
```bash
# Prune to get independent SNPs (r² < 0.2)
plink --bfile merge5 --indep-pairwise 10000 1 0.2 --out prune0.2

# Keep only pruned SNPs, autosomes only
plink --bfile merge5 --extract prune0.2.prune.in --chr 1-22 --make-bed --out TGEN_hapmap_pruned
```
Result: 58,873 LD-pruned SNPs for PCA

**Step 6: Run PCA on Combined Data** 
```bash
# Use hapmap.par parameter file
smartpca -p hapmap.par > hapmaptgen.out
```
**Step 7: Interpret Results**
From `hapmaptgen.out`:

```text
p-value eigenvector_1_Case_Control_ 0.198411
```
Interpretation: PC1 is NOT significantly different between cases and controls (p=0.198) → Good! No confounding.

But check the plot:
```bash
Rscript --vanilla plotPCs.R TGEN_hapmap_pruned.evec 1 2 10
```
What you see:
- CEU (blue) cluster together
- CHB/JPT (another cluster)
- YRI (separate cluster)
- Most TGEN samples cluster with CEU (expected for European ancestry)
- A few outliers: Some TGEN samples near CHB/JPT or YRI

**Step 8: Identify Outliers**
```bash
# Find TGEN samples with PC1 < -0.01
awk '($2 < -0.01) && ($12=="Case" || $12=="Control") {print $0}' TGEN_hapmap_pruned.evec > outliers.txt
```
Result: 6 TGEN samples are outliers → should consider removing them

## Part 3: What it all means - conntecing lecture and code
Key Connections:

1. Theory → Code
- Lecture: "Structure causes false associations"
- Code: We test this by checking if PCs differ between cases/controls (test.out ANOVA)
- Action: If significant, include PCs as covariates

2. HapMap Reference Populations
Lecture: Use known populations to interpret your samples
- Code: Merge your data with HapMap, run PCA, see where your samples cluster
- Result: TGEN samples should cluster with CEU (they mostly do, except outliers)

3. Strand Issues
- Lecture: A/T and G/C SNPs are problematic for merging
- Code: We identify and remove them before merging
- Why? Prevents allele mismatches that would ruin PCA

4. LD Pruning
- Lecture: Need independent SNPs for PCA
- Code: --indep-pairwise 10000 1 0.2
- Translation: In 10,000 kb windows, shift 1 SNP at a time, remove if r² > 0.2

5. Outlier Detection
- Lecture: Remove samples that don't match expected ancestry
- Code: Use awk to find samples with extreme PC values
- Next step: Remove outliers, re-run PCA on just your cleaned data

Common Pitfalls & Solutions:
|Problem	|Cause|	Solution
| -------- | ------- | -------- | 
|Merge fails with "3+ alleles"	|Strand mismatch	|Use `--flip` or `--exclude`
|PCA shows case-control separation	|Population stratification	|Include PCs as covariates
|Samples cluster wrong population	|Wrong ancestry or outliers	|Remove outliers
|Too few SNPs after pruning	|Too strict threshold	|Use r² < 0.5 or larger window

### Summary
1. Population structure = non-random mating → subgroups
2. Detect with PCA on LD-pruned SNPs
3. Use HapMap to verify ancestry
4. Remove A/T, G/C SNPs before merging datasets
5. Check if PCs differ between cases/controls (ANOVA in .out)
6. If significant: include PCs as covariates
7. Remove outliers that don't cluster with expected population
8. Final analysis: Clean data + PCs as covariates = unbiased associations

### Command Cheat Sheet
```bash
# PCA workflow:
smartpca -p parameter_file.par > output.log

# Plot PCs:
Rscript plotPCs.R file.evec PCx PCy numPCs

# Find outliers:
awk '($col < threshold) && phenotype_condition {print}' file.evec

# LD prune:
plink --bfile data --indep-pairwise 10000 1 0.2 --out pruned

# Remove A/T, G/C SNPs:
awk '($5=="A"&&$6=="T")||($5=="T"&&$6=="A")||($5=="G"&&$6=="C")||($5=="C"&&$6=="G"){print $2}' file.bim > to_remove.txt
```
Final Thought: Population structure is a confounder. PCA helps detect it, and including PCs as covariates adjusts for it. The code implements the theory: prepare data (remove problematic SNPs, LD prune), run PCA, interpret results, take action (remove outliers, adjust models).
