## Homework 4: Genotype Imputation
### BS859 Applied Genetic Analysis
### Addison Yam
### February 18, 2026

Use the chromosome 2 TOPMed-reference imputed genotype data and genotyped data in /projectnb/bs859/data/tgen/2023/chr2 to complete this assignment.

Imputed genotypes and information files (using the TOPMed reference panel, build 38 locations):
- chr2.dose.vcf.gz
- chr2.info.gz

The Chromosome 2 genotyped data file (using genome build 37 locations) uploaded to the imputation server are in the same directory:
- chr2.vc

The PC covariate file and fam file with case/control status are:
- /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt
- /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.fam

(the same files that we used in the in-class analysis.)

Chromosome 2 has robust, replicated association at the BIN1 gene in published analyses. The end goal today is to 1) assess the overall imputation quality and 2) determine whether there is evidence for association at or near this locus in the TGEN data set. The top SNP in the IGAP GWAS is rs6733839, which has location
2:127135234 (GRCh38)
2:127892810 (GRCh37)
(see class notes for link to paper).
```bash
# load the necessary modules
module load R 
module load plink/1.90b6.27
#make an alias for the directory with the data
export DATADIR=/projectnb/bs859/data/tgen/cleaned
```

1. The imputation report and QC report for chromosome 2 from the imputation server are in the files chr2_imputationreport.pdf and chr2_imputationQCreport.pdf. Looking at the plot in the QC report, explain what is being plotted. Do you have any concerns, or any suggestions that might improve the imputation?

2) How many imputed SNPs are in the chromosome 2 dosage vcf?
3) a) How many imputed SNPs have both R2>=0.3 and MAF>=0.005?

b) How does this number compare to the number of SNPs that were used to impute these variants.

c) How does this number compare to the total number of variants imputed?

4) Provide a crosstabulation of the imputed variants with MAF>0, broken down into the same MAF and R2 bins as we did for chr19 in class. Which bin has the largest number of SNPs?

5) a) Provide Rsq box plots (or violin plots if you are feeling fancy) for all of the imputed variants, broken down into the same MAF bins as we did for chr19 in class.
b) Explain the plot, including why the distribution of R2 is different by MAF bin.

c) How does this plot compare to the same type of plot for the CHR19 imputed data we looked at in class?

6) a) Plot the R2 boxplots vs BP location for chr2 as we did for chr 19 in class.

b) How do these compare to the plots from chr 19?

c) Based on this, do you think that the chromosome 2 imputation is overall better or worse than the chromosome 19 imputation? Explain your answer.

7) Perform an association analysis of the imputed chromosome 2 variants with MAF>=0.005 and imputation R2>=0.3 with Alzheimer disease in the TGEN data. Use the same PCs we used in class for the chromosome 19 analyses as covariates.
a) Make a QQ and (single chromosome) manhattan plot for the tgen imputed analysis of chromosome 2 to display the results. Interpret them in a few sentences.

b) Find the 4 most significant associations (i.e., the 4 variants with the smallest p-values), and show the full results for these 4. Are any of them within 100,000 bp of the top BIN1 gene variant reported in the IGAP GWAS? (That snp is located at 127135234 on genome build 38).

c) Report the association for the BIN1 gene variant reported in the IGAP paper (rs6733839, which is at location 127135234 on genome build 38) in your imputed chromosome 2 association analyses. In the published paper (see paper or the table in the class notes), the
T allele of this SNP increased odds of Alzheimer Disease by 1.22. Are your results for this SNP consistent with the published results? Explain your answer.

d) What is the allele frequency of the T allele, and imputation R2, for this variant? Is this variant reliably imputed? Explain your answer.

e) Was this SNP among the genotyped variants used to do the imputation? Explain your
answer.