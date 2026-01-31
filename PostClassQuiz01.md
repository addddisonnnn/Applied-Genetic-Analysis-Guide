## Post-Class 1: Genotype Data QC & PLINK
### BS859 Applied Genetic Analysis
### January 21, 2026
1. Which statement about call rate filtering is correct? Select up to 4 options
    - You should always remove low call-rate individuals first to maximize SNP retention.
    - **Low SNP call rate means many individuals lack a genotype for that marker.**
    - Call rate refers only to the number of SNPs on the chip, not to samples(individuals).
    - Call rate filtering is irrelevant for modern chips and can be skipped.

SNP call rate = proportion of samples with a non-missing genotype at that marker. Call rate applies to both markers and samples.  While modern chips tend to have higher call rates, filtering out low call rate is still recommended.

2. In a case-control GWAS, which sample should you generally use when testing SNPs for departure from Hardy–Weinberg equilibrium (HWE) as a QC measure?
All individuals (cases + controls) always
    - Only affected individuals (cases)
    - **Only unaffected individuals (controls)**
    - Only related individuals

Controls are used for HWE QC because cases may deviate due to true association; related individuals complicate HWE testing.

3. Because males have one X chromosome, X-chromosome SNPs outside the pseudo-autosomal region are expected to appear homozygous in males; therefore a high proportion of heterozygous calls on X in a male sample suggests a problem (e.g., assay error or sample mix-up).
    - **True**
    - False

Males are hemizygous for non-PAR X SNPs, so heterozygous genotypes are unexpected and indicate possible issues.

4. Identical-by-state (IBS) is always less than identical-by-descent (IBD) for a pair of individuals.
    - True
    - **False**

IBS counts alleles identical in state regardless of ancestry and is always greater than or equal to IBD (IBS ≥ IBD).

5. Why do we perform LD-pruning (creating a set of variants with low pairwise linkage disequilibrium) before some analyses?
    - To increase the number of variants for association testing
    - To remove monomorphic markers only
    - **To obtain an approximately independent set of markers for analyses like heterozygosity estimation, IBD estimation, and PCA**
    - To force all SNPs to have the same minor allele frequency

LD-pruning reduces correlation among markers so methods that assume independence perform correctly.

6. A pair of samples has PI_HAT ≈ 0.5 in a PLINK IBD output. Which is the most likely relationship?
    - Unrelated
    - Monozygotic twins or duplicate sample
    - **Parent–child or full siblings (first-degree)**
    - Third-degree relatives (first cousins)

PI_HAT ≈ 0.5 corresponds to sharing on average half their alleles IBD, consistent with first-degree relationships.

7. Which are indicators that a DNA sample might be contaminated (not a pure sample from one individual)?  Choose all that are correct.
    - **The individual has high IBD sharing with multiple other individuals in the - sample**
    - **The individual has an F statistic that is negative and large in absolute value**
    - The sample has a low call rate (high missing rate)
    - The sample has 0 IBD sharing with all other individuals in the sample
    - The individual has an F statistics that is positive and large in absolute value

Contamination can induce false heterozygote calls (excess heterozygosity, negative F) and artificially inflate IBS/IBD with multiple samples.
