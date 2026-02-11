### BS859 Class 4 in-class analsyes

#plink2 has a lot of features for working with different file types such 
# as vcfs that plink 1.9 does not have.  It does not have all for the 
# analysis options that plink 1.9 has, however.  We need plink2 to work
# with the imputed data vcf files

module load plink2/alpha3.7

#location of the imputed genotypes for TGEN chromosome 19:
export CHR19DIR=/projectnb/bs859/data/tgen/2023/chr19/

#How many variants are in the imputed chromosome 19 file? 
#the number of lines -1 in the "info" file is the number of imputed 
#variants.
#(we can also determine this by looking at the log file for the plink2 run
#below that converts the vcf to plink format.)

zcat $CHR19DIR/chr19.info.gz |wc

##The SNP information file has MAF, imputatation quality R2, and
##other metrics:
zcat $CHR19DIR/chr19.info.gz |head
##There are *a lot* of variants with MAF<<0.01 and R2<<0.3!


##Before we do an association anslysis, we will take a look
##at the R2 and MAF distribution of the SNPs and make sure that
##they are close to what we expect  
module load R
Rscript MAFinfo.R > MAFinfo.log

zcat $CHR19DIR/chr19.dose.vcf.gz|head -n 21 |cut -f 1-10

zcat $CHR19DIR/chr19.info.gz|head -n 2


##Take a look at the imputed data vcf (its a gzipped text file):
zcat $CHR19DIR/chr19.dose.vcf.gz|head -n 21|cut -f 1-10

##Convert the vcf that comes from the imputation server to 
##plink2 "pgen" format, which is much faster to read in/out for analyes
##at the same time, filters the imputed variants so that variants with 
##Minor Allele Frequency (MAF)<0.005 are excluded, and variants with 
##imputation R2<0.3 are excluded.  VCF format does not include phenotypes,
##so this command also tells plink2 where to get the phenotype information:


##In plink2, you must specify the actual column number for the phenotype.
##in the fam file, the phenotype is the 6th column (fid iid mid fid sex pheno)

plink2  --double-id\
  --exclude-if-info "MAF<0.005" \
  --extract-if-info "R2>=0.3" \
  --vcf $CHR19DIR/chr19.dose.vcf.gz dosage=DS \
  --pheno /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.fam\
  --pheno-col-nums 6\
  --make-pgen \
  --out tgen_chr19_imputed


head tgen_chr19_imputed.psam
head /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.fam

head tgen_chr19_imputed.pvar -n 20
head /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.bim


cat tgen_chr19_imputed.log


##In the last homework, we found that PCs 6 and 8 were associated with AD
##with p<0.01.  We will use these 2 PCs as covariates in our analyses.  
plink2 --pfile tgen_chr19_imputed \
--covar /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt \
--covar-name PC6,PC8 --logistic hide-covar --out chr19_01pcs 

head chr19_01pcs.PHENO1.glm.logistic.hybrid

module load R
Rscript  qqplot_pgen.R chr19_01pcs.PHENO1.glm.logistic.hybrid chr19.tgen ADD

Rscript  gwaplot_pgen.R chr19_01pcs.PHENO1.glm.logistic.hybrid "chr19.imputed" "chr19.imputed" 

##let's compare with the chr19 data that we used for the imputation:
##this is the same data as we had in the *.bed/*.bim/*.fam files we used last
##week.  I extracted chromosome 19 and did a little additional filtering, 
##and put it in the vcf format needed when uploading to the imputation server.
plink2  --double-id\
   --vcf $CHR19DIR/chr19.vcf \
  --pheno /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.fam\
  --pheno-col-nums 6\
   --covar /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt\
   --covar-name PC6,PC8\
     --logistic hide-covar --out chr19_GENOTYPED 
 

Rscript  qqplot_pgen.R chr19_GENOTYPED.PHENO1.glm.logistic.hybrid chr19.tgen.genotyped ADD

Rscript  gwaplot_pgen.R chr19_GENOTYPED.PHENO1.glm.logistic.hybrid chr19.tgen.genotyped "chr19.genotyped"

##since the chromosome position is different in the genotyped-only and imputed
##genotypes files, we'd need to look up the variants to determine which 
##which positions match up.
##We can start by just looking at the most significant variants (APOE region)

sort -gk13 chr19_GENOTYPED.PHENO1.glm.logistic.hybrid |head -n 5|cut -f 1-6,10-13


sort -gk13 chr19_01pcs.PHENO1.glm.logistic.hybrid |head -n 5|cut -f 1-6,10-13