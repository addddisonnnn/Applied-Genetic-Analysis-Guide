## Post-Class Quiz 5: Association analysis of Rare Variants and Sequence Data
### BS859 Applied Genetic Analysis
### February 18, 2026

1. Which of these is a recommended way to reduce multiple-testing burden when performing many variant-set tests per gene (different variant filters/weights)?
- Apply Bonferroni across every individual test without exception
- Only test single variants and ignore gene-based tests
- Combine multiple test p-values per gene using ACAT to produce a single p-value
- Use Ti/Tv ratio as a correction factor

2. Which of the following is the primary rationale for aggregating rare variants in gene- or region-based tests instead of testing each rare variant individually?

- Sequencing errors make single-variant tests impossible
- Single-variant tests have low power for very low MAF variants
- Aggregation removes the need for covariate adjustment
- Aggregation guarantees causal variants will be identified

3. A gene-set whose cumulative minor allele count (cMAC) is extremely low should still always be included in the multiple-testing correction because it might become significant after Bonferroni adjustment.
- True
- False

4. Which test is most appropriate when you expect many variants in a set to have effects in different directions (some increasing, some decreasing the trait)?
- Burden test
- CAST
- SKAT (variance-component test)
- Madsen-Browning burden

5. When designing gene-based rare-variant association tests, which variant inclusion strategy is most likely to increase statistical power by enriching the tested set for variants that are plausibly functional (i.e., more likely to affect protein function)?
- Include all variants within a 1 Mb window around the gene, without additional filtering
- Include only loss-of-function (LoF) variants and nonsynonymous variants predicted to be deleterious by functional annotation
- Include only intergenic variants with no predicted regulatory function
- Randomly select variants located within the gene region regardless of annotation

6. Which of the following is a proper reason to use a mixed-model framework (GLMM) for rare-variant set tests like SMMAT?
- To account for relatedness or cryptic population structure via a genetic relationship matrix (GRM)
- To remove the need to annotate variants by function
- Because mixed models give exact p-values without asymptotics
- To automatically increase the minor allele counts of rare variants