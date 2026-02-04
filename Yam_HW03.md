## Homework 3: Common Variant GWAS - Association Analysis & Mixed Models
### BS859 Applied Genetic Analysis
### Addison Yam
### February 11, 2026

We will use the TGEN Alzheimer data located in /projectnb/bs859/data/tgen/cleaned/ to do this assignment. For PC covariates, use the file use the file TGEN_pcs.txt in the same directory. In your submitted homework, include the commands you used to get your results, and the parts of the results that support your response. Use tables and/or plots to present results clearly and succinctly.

1. Use the file TGEN_pcs.txt to determine which PCs are associated with case status. (TGEN_pcs.txt is the output from smartpca with a column header added). **HINT: “--logistic no-snp” produces results for the specified covariates, with no SNPs in the model.**

a. Show the commands you used for the analysis.

b. Where is plink looking to find case status for your analysis?

c. How many individuals are included in your analysis?

d. Which PCs are associated with case status at p<0.01?

2. Perform a GWAS of Alzheimer Disease case status using logistic regression, adjusting for the PCs associated with case status at p<0.01 that you found in question 1. *Use a plink command that hides the covariates so that your plink output has only the association results for the SNPs, and that gives the effect estimates as regression coefficients rather than odds ratios.*

a. Show your plink command

b. Show the results for the first 5 SNPs in the output file.

c. How many SNPs in this GWAS have p-value < 0.0001?

d. Show the plink results for the most significant SNP. Show the 3 genotypes for the SNP, and how they are coded in the plink logistic regression. Compute the odds ratio. Which allele increases the risk of Alzheimer Disease?

3. Create a GRM from the TGEN data so that you can perform a GWAS using a logistic mixed model. Use only chromosomes 1-22 (exclude X). Don’t forget to prune before you make the GRM! Use pruning parameters: --indep-pairwise 10000kb 1 0.15.

a. Show your code.

b. How many SNPs were left after pruning and removing the chromosome X variants?

c. Provide the first 10 lines of the id grm id file (the file with the ids for the grm)

d. What is the maximum relationship coefficient for MAYO_10139, other than with itself?

4. Perform a GWAS using a logistic mixed model score test, with no covariates. Do a second GWAS using this model and the same PCs you included in the Q2 analysis. You can edit the GMMAT.R program we used in class to do this.

a. How many SNPs have p-value < 0.0001 in each of these two GWAS?

b. Compare the results of the two GLMMs and the logistic model from question 2 for the SNP that is most significantly associated with AD case status in the question 2d analysis. Which analysis yields the most significant association? Can you compare the effect estimates? Which allele increases risk in these models?

5. Produce QQ plots for the three GWAS analyses, and show them here. Looking at the Genomic
Control Lambda and the overall qq plot, which analysis model do you think is best? Explain your
choice.

6. Produce a Manhattan plot for the GWAS that you chose in 5. Based on this plot, are there any regions of the genome OTHER than chromosome 19, that may be of interest? Why or why not?