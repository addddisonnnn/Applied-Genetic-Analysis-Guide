export IGAP='/projectnb/bs859/data/igap'
export LDSCORES_DIR='/projectnb/bs859/data/ldscore_files'
module load R
module load python2
module load ldsc

## The LD scores directory has 1000G EUR using hapmap3 snps in 
## $LDSCORES/eur_w_ld_chr
## and UK Biobank-based LD scores using hapmap3 snps 
## in $LDSCORES/UKBB.ALL.ldscore
ls $LDSCORES_DIR/eur_w_ld_chr
ls $LDSCORES_DIR/UKBB.ALL.ldscore

### IGAP 2019 meta-analysis 
## 21,982 cases, 41,944 cognitively normal controls

cat $IGAP/Kunkle_README

zcat $IGAP/Kunkle_etal_Stage1_results2019.txt.gz|head
zcat $IGAP/Kunkle_etal_Stage1_results2019.txt.gz|wc

##Step 1:  format 2019 AD summary statistics for ldsc
munge_sumstats.py \
--sumstats $IGAP/Kunkle_etal_Stage1_results2019.txt.gz \
--snp MarkerName \
--N-cas 21982 \
--N-con 41944 \
--a1 Effect_allele \
--a2 Non_Effect_allele \
--signed-sumstats Beta,0 \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out AD

##output:  reformated, gzipped summary statistics (in your own directory)
zcat  AD.sumstats.gz|head
zcat  AD.sumstats.gz|wc


## note: could have done much of what the munge tool did with awk just 
## as easily, but this tool provides
## some useful information about what is kept/removed,and matches up alleles
## between the summary statistics and the snps used in the LD scores

##Note that 1 snp had out-of-bounds p-value.  It would be a good idea to
##investigate this and make sure its not a problem

##what do the LD score files look like?
ls $LDSCORES_DIR/eur_w_ld_chr/
cat /projectnb/bs859/data/ldscore_files/eur_w_ld_chr/README

##what is in these files?
head $LDSCORES_DIR/eur_w_ld_chr/1.l2.M_5_50
zcat $LDSCORES_DIR/eur_w_ld_chr/1.l2.ldscore.gz|head

ls /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/
cat $LDSCORES_DIR/UKBB.ALL.ldscore/README 

##what is in these files?
zcat /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.AFR.l2.ldscore.gz|head 
zcat /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.AFR.rsid.l2.ldscore.gz|head 

##note for the UKBB files:
##UKBB files:  1 file set for all chromosomes, 5 different population sets
##Two different SNP name formats: rsid and chr:pos:ref:alt
##


### AD h2 estimated from 2019 data, on original (not liability) scale
### using 1000G EUR ldscores
## the Kunkle et al summary statistics have rsids for the markers, so
##the 1000G ldscore files will work
#OPTIONS USED:
#--h2 tells ldsc to compute heritability
#--ref-ld-chr tells ldsc which LD score files to use as 
#the independent variable in the LD score regression
#--w-ld-chr tells ldsc the location of the LD scores used 
#for regression weights (we use the same for both)
#remove the -chr if all the chromosomes are in one file, 
#as for UKBB LD scores)

ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--out AD_h2_orig_scale


##how sensitive is the h2 estimate to the reference LD scores used?
##### AD h2 estimated from 2019 data, on original (not liability) scale
### using UKBB EUR ldscores
ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out AD_h2_orig_scale_UKBB



##how do the results differ?

#if this trait were continuous, we would be done.  Here, AD is a dichotomous
#trait, and we want to estimate the heritability on the liability scale.
##  This requires a simple transformation.  We need K (population prevalence)
## and  p (sample prevalence) of cases.  Can re-run ldsc with these parameters
##included, or use R to do the transformation

#In R, K=0.1, p=21982/(21982+41944)=0.344:
#K<-0.1
#p<-21982/(21982+41944)
#H2.liab.K0.1<-0.0713*(K^2)*(1-K)^2/(p*(1-p)*dnorm(qnorm(K,lower=F))^2)
#H2.liab.K0.1
#result is 0.083, the same as re-running ldsc with pop-prev and samp-prev
#options

#Assume AD population prevalence is 0.10, specify the sample prevalence 0.314:
ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--pop-prev 0.10 \
--samp-prev 0.344 \
--out AD_h2_prev.1

## assume AD population prevalence is 0.05
ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--pop-prev 0.05 \
--samp-prev 0.344 \
--out AD_h2_prev.05

## assume AD population prevalence is 0.15
ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--pop-prev 0.15 \
--samp-prev 0.344 \
--out AD_h2_prev.15

###Genetic Correlation
###GWAS summary statistics from a study of Brain Volume in 

ls /projectnb/bs859/data/brainvolume

zcat /projectnb/bs859/data/brainvolume/UKB_bv_height_cor.sumstats.txt.gz|head 

munge_sumstats.py \
--sumstats /projectnb/bs859/data/brainvolume/UKB_bv_height_cor.sumstats.txt.gz \
--snp RSID \
--a1 A1 \
--a2 A2 \
--signed-sumstats BETA,0 \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out BV.HGT

zcat BV.HGT.sumstats.gz|head

ldsc.py \
--h2 BV.HGT.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--out BV.HGT.h2

ldsc.py \
--rg AD.sumstats.gz,BV.HGT.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--out AD_BV.HGT_rg



#####partitioned:
###
###baseline ld score and weight files:
ls  $LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline.*
ls  $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.*

cat /projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.1.l2.M
cat /projectnb/bs859/data/ldscore_files/1000G_EUR_Phase3_baseline/baseline.1.l2.M_5_50

zcat  $LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline.1.l2.ldscore.gz|head -n 2
zcat  $LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline.1.annot.gz|head -n 2


ls  $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.*
zcat  $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.1.l2.ldscore.gz|head -n 2



###compute heritability for each of the 53 categories:
ldsc.py \
 --h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR//1000G_EUR_Phase3_baseline/baseline. \
--w-ld-chr  $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
--pop-prev 0.10 \
--samp-prev 0.344 \
--print-coefficients \
--out AD_h2_partitioned




## We can look at h2 for SNPs that influence different types of cells.
##cell type groups are in the directory
## $LDSCORES_DIR/1000G_Phase3_cell_type_groups
ls  $LDSCORES_DIR/1000G_Phase3_cell_type_groups

##file "names" describes the groups:

cat $LDSCORES_DIR/1000G_Phase3_cell_type_groups/names

##cell type groups: CNS


ldsc.py \
 --h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/1000G_Phase3_cell_type_groups/cell_type_group.3.,$LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline.  \
--w-ld-chr  $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
--pop-prev 0.10 \
--samp-prev 0.344 \
--out AD_CNS \
--print-coefficients


