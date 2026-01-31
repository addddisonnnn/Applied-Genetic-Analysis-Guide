## Session 1: Genotype Data QC & PLINK
### BS859 Applied Genetic Analysis
### January 21, 2026
The main topics of this week's session are about (1) genetic association studies and GWAS (Genome Wide Association Studies), (2) data summaries and QC and the first steps prior to any analysis (the type of summaries and tests, their meaning, and what to do with the results), (3) introduction to the PLINK tool for GWAS data, and (4) using PLINK to implement tests and summaries (and using the Linux tool `awk` to parse large files).
### I. Genetic Assocations Studies & GWAS
**A. Core Concepts**
- Goal: Find SNPs associated with complex traits (non-Mendelian)
- Human genome: ~3 billion base pairs, 23 chromosomes
- Complex traits: Height, BMI, Alzheimer's, heart disease (genes + environment)
- GWAS: Hypothesis-generating, scans whole genome (not candidate gene)

**B. Historical Timeline**
- 1978: RFLPs for linkage (Huntington's chr4)
- 1987: First complete human genetic map
- 1989: Microsatellites enable genome-wide linkage
- 1990-2003: Human Genome Project
- 2005: SNP genotyping arrays
- 2006: First GWAS
- 2010+: Next-gen sequencing affordable

### II. GENOTYPE DATA & QC RATIONALE
**A. Data Structure**
- Format: Individuals (rows) × SNPs (columns)
- SNP/SNV: Single nucleotide variant (A/C/G/T substitution)
- Genotype: Two alleles per individual (11, 12, 22)
- Alleles: "1" and "2" (coded arbitrarily)

**B. Why QC is Critical**
- Systematic errors: Batch effects → spurious associations
- Random errors: Missing data → reduced power
- Sample issues: Mix-ups, contamination, relatedness
- Goal: Remove problematic data BEFORE analysis

### III. QC METRICS: SNP-BASED
**A. Call Rate (Missingness)**
- SNP call rate: Proportion of individuals with non-missing genotype
- Threshold: Usually 95-99% (newer chips = stricter)
- Low call rate indicates: Poor genotyping assay
- PLINK command:
```bash
plink --bfile data --missing --out miss1
```
- Filter:
```bash
plink --bfile data --geno 0.05 --make-bed --out filtered  # removes SNPs with >5% missing
```
**B. Hardy-Weinberg Equilibrium (HWE)**
- HWE Law: p² + 2pq + q² = 1 (random mating, large population)
- Test: χ² or exact test (preferred for low MAF)
- HWE violation causes:
    1. Non-random genotyping errors
    2. Population structure
    3. Assay bias (heterozygotes called missing)
- When to test:
    - Case-control: Use CONTROLS only
    - Population sample: Use all individuals
    - Related individuals: Use unrelated subset
- Threshold: p < 1e-4 (in controls)
- PLINK command:
```bash
plink --bfile data --hardy --out hwe1
```
- Filter:
```bash
plink --bfile data --hwe 1e-4 --make-bed --out filtered
```

### IV. QC METRICS: INDIVIDUAL-BASED
**A. Individual Call Rate**
- Individual call rate: Proportion of SNPs with non-missing genotype
- Threshold: Usually 95-99%
- Low call rate indicates: Poor DNA quality
- Strategy: Remove low call rate SNPs FIRST, then individuals
- PLINK filter:
```bash
plink --bfile data --mind 0.05 --make-bed --out filtered  # removes individuals with >5% missing
```
**B. Heterozygosity (F-statistic)** 
- Formula: F = 1 - (Het<sub>obs</sub>/ Het<sub>exp</sub>)
- Interpretation:
    - F ≈ 0: Expected (random mating)
    - F > 0: Heterozygote deficit (inbreeding, population structure)
    - F < 0: Heterozygote excess (sample contamination)
    - F = 1: No heterozygotes observed
- Requires LD-pruned SNPs (independent markers)
- PLINK command:
```bash
plink --bfile data --het --out het_results
```
**C. Sex Check**
- Biology: Males = XY (hemizygous), Females = XX
- Expectation: Males should be homozygous for X chromosome SNPs (except PAR)
- PLINK sex check:
```bash
plink --bfile data --check-sex --out sexcheck
```
- F-statistic thresholds:
    - F > 0.8 → Male
    - F < 0.2 → Female
    0.2 ≤ F ≤ 0.8 → Ambiguous
- Action: Exclude mismatches (reported vs. genotype sex)

**D. Relatedness (IBD/IBS)**
- IBS: Identical by State (alleles look identical)
- IBD: Identical by Descent (from common ancestor)
- π (pi-hat): Proportion of alleles shared IBD = P(IBD=2) + 0.5×P(IBD=1)
- Expected π values:
    - MZ twins: 1.0
    - Parent-child/Siblings: 0.5
    - Grandparent/grandchild: 0.25
    - Cousins: 0.125
    - Unrelated: 0
- PLINK command (requires LD-pruned SNPs):

```bash
plink --bfile data --genome --out ibd_results
```
- Filter for relatives:

```bash
plink --bfile data --genome --min 0.05 --out ibd05  # keeps pairs with π ≥ 0.05
```

### V. PLINK COMMANDS & WORKFLOW
**A. File Formats**
- Binary set: `.bed` (genotypes), `.bim` (SNP info), `.fam `(sample info)
- .fam columns (6): FID, IID, Father ID, Mother ID, Sex (1=M, 2=F, 0=unknown), Phenotype
- .bim columns (6): Chr, SNP ID, Genetic distance, BP position, Allele1, Allele2

**B. Complete QC Workflow**
```bash
# 1. Load data
export DATADIR=/path/to/data
module load plink/1.90b6.27

# 2. Initial summaries
plink --bfile $DATADIR/data --freq --out freq1
plink --bfile $DATADIR/data --missing --out miss1
plink --bfile $DATADIR/data --hardy --out hwe1

# 3. Filter SNPs first (MAF, missingness)
plink --bfile $DATADIR/data --maf 0.01 --geno 0.05 --make-bed --out step1

# 4. Filter individuals (missingness)
plink --bfile step1 --mind 0.05 --make-bed --out step2

# 5. Filter HWE violations (in controls)
plink --bfile step2 --hwe 1e-4 --make-bed --out cleaned

# 6. LD pruning for independent SNPs
plink --bfile cleaned --indep-pairwise 10000 1 0.2 --out pruned
# 10000 kb window, 1 SNP shift, r² threshold 0.2

# 7. Sex check
plink --bfile cleaned --check-sex --out sexcheck

# 8. Heterozygosity (use pruned SNPs)
plink --bfile cleaned --extract pruned.prune.in --het --out het_results

# 9. IBD estimation (use pruned SNPs)
plink --bfile cleaned --extract pruned.prune.in --genome --out ibd_results
```

**C. Order of Operations MATTERS**
1. Correct order:
2. Remove low call rate SNPs (--geno)
3. Remove low MAF SNPs (--maf)
4. Remove low call rate individuals (--mind)
5. Remove HWE violations (--hwe)
6. LD pruning
7. Sex check, heterozygosity, IBD

**D. AWK Commands for Parsing**
```bash
# Find individuals with >3% missing
awk 'NR==1 || $6>0.03 {print $0}' miss1.imiss > high_missing.txt

# Find SNPs with HWE p<0.0001 in controls
awk '($3=="UNAFF" && $9<0.0001) || NR==1 {print $0}' hwe1.hwe > hwe_fail.txt

# Find IBD pairs with π > 0.05
awk '$10>0.05 {print $0}' ibd.genome > relatives.txt

# Count lines in file (number of individuals/SNPs)
wc -l file.fam  # individuals
wc -l file.bim  # SNPs
```
**E. Special Cases & Flags**
```bash
# Samples with unknown sex (coded as 0)
plink --bfile data --allow-no-sex --maf 0.01 ...  # Prevents exclusion

# Remove specific individuals
plink --bfile data --remove bad_samples.txt --make-bed --out filtered

# Remove specific SNPs
plink --bfile data --exclude bad_snps.txt --make-bed --out filtered

# Keep only specific chromosomes
plink --bfile data --chr 1-22 --make-bed --out autosomes  # Remove X,Y,MT
```

### VI. INTERPRETATION & TROUBLESHOOTING
**A. Common Problems & Solutions**
1. All individuals excluded from analysis:
    - Check: Sex column has 0 (unknown) → Use `--allow-no-sex`
    - Check: All phenotypes missing → Use `--allow-no-pheno`
2. High missingness:
    - Older chips: Use 95% threshold
    - Newer chips: Use 98-99% threshold
3. Population structure affects HWE/IBD:
    - Use methods robust to structure (KING, PC-RELATE)
    - Prune SNPs more aggressively (r² < 0.1)
4. Sample contamination:
    - Signs: Negative F, high heterozygosity
    - Action: Exclude contaminated samples

**B. QC Threshold Guidelines**
| Metric | Typical Threshold |Strict Threshold|
| -------- | ------- | -------- | 
| SNP call | rate	95% (--geno 0.05)	| 99% (--geno 0.01)
| Individual call rate	| 95% (--mind 0.05)	| 99% (--mind 0.01)
| MAF | 1% (--maf 0.01)	| 5% (--maf 0.05)
| HWE (controls) | p < 1e-4	| p < 1e-6
| IBD (π)| > 0.05 (2nd degree+)	| > 0.125 (1st degree)

C. Expected Output Files
- `.log`: Command log, check for errors
- `.frq`: Allele frequencies
- `.imiss`: Individual missingness
- `.lmiss`: SNP missingness
- `.hwe`: HWE test results
- `.het`: Heterozygosity statistics
- `.genome`: IBD estimates
- `.sexcheck`: Sex concordance

### VII. KEY POINTS
**A. Must-Know Concepts** 
1. HWE assumptions: Random mating, no selection/mutation/migration, large population
2. IBD vs IBS: IBD always ≤ IBS
3. X chromosome in males: Should be homozygous (except PAR)
4. LD pruning purpose: Get independent SNPs for F-stat and IBD
5. Case-control HWE: Test CONTROLS only (cases may violate HWE)

**B. Command Pattern Recognition**
```bash
# Summary command:
plink --bfile INPUT --METRIC --out PREFIX

# Filter command:
plink --bfile INPUT --FILTER_THRESHOLD --make-bed --out OUTPUT

# Combined filter:
plink --bfile INPUT --maf X --geno Y --mind Z --hwe W --make-bed --out OUTPUT
```
**C. Commands for reporting and filtering data**
 What? | Report/summary File |Filtering|
| -------- | ------- | -------- | 
|Allele Freq |--freq | --maf
|By-variant missing proportion |--missing |--geno 
|By-individual missing variant population |--missing |--mind 
|Hardy-Weinberg exact test |--hardy |--hwe


**C. Data Sizes & Computation**
- n individuals → `.imiss` has n+1 lines
- m SNPs → `.frq `has m+1 lines
- IBD file size: n×(n-1)/2 pairs
- Example: 1000 samples → ~500,000 IBD pairs

**D. Homework Data (TGEN) Specifics**
- Issue: Sex column all 0s
- Solution: Use --allow-no-sex
- Data: 1411 individuals, 312,316 SNPs
- Paper: Alzheimer's GWAS

### STUDY CHECKLIST:
- Can define all QC metrics and their interpretations
- Can write full QC pipeline from memory
- Know thresholds for each filter
- Understand when to use which subset (controls for HWE, pruned for F/IBD)
- Can troubleshoot common PLINK errors
- Know how to parse output files with awk
- Understand biological rationale for each test

Study Tip: Practice writing commands WITHOUT looking at notes. Time yourself on a mock dataset with different parameters (e.g., --maf 0.05 --geno 0.02 --mind 0.03).