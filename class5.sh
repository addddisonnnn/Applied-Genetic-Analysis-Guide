module load R  

DATADIR="/projectnb/bs859/data/rarevariant"

ls $DATADIR

##take a look at the example data.  The first 33 lines of the 
##VCF file are "information"
##the 34th line gives the IDs, and the 35th is the first variant.

zcat $DATADIR/1000G_exome_chr20_example_softFiltered.calls.hg19_anno.vcf.gz |head -n 37 |cut -f1-10

## the VCF only stores genotype data (and annotations).  A separate file stores 
##phenotype and covariate data:

head $DATADIR/example.ped

head example.smmat.exonic.groups
head example.smmat.nonsyn.groups

## can run the SMMAT analyses using:
Rscript --vanilla SMMAT_example.R

##I will step through this in class

Rscript --vanilla qqplot_smmat.R exonic.nvargt1.csv
Rscript --vanilla qqplot_smmat.R nonsyn.nvargt1.csv


##example code for doing SMMAT gene-based analyses using the imputed genotype
##data we used last session.  Send screen output to a log file:
Rscript --vanilla SMMAT_tgenChr19Example.R > SMMAT_tgenChr19Example.log

Rscript --vanilla qqplot_smmat.R chr19.exonic05.csv




###using bcftools and awk to create the group files:
module load htslib/1.18
module load bcftools/1.18

bcftools query -f '%Gene.refGene %CHROM %POS  %REF %ALT %Func.refGene %ExonicFunc.refGene\n' /projectnb/bs859/data/tgen/annotated_imputed_vcfs/anno1.chr19.vcf.gz >temp1

awk '$6~"exonic"{print $0}' temp1 > temp2

##all exonic:
awk '{print $1,$2,$3,$4,$5,1}' temp2 > chr19.exonic.group
##just nonsynonymous, frameshift, stop/start gain and stop/start loss:
awk '$7~"nonsyn"||$7~"stop"||$7~"start"||$7~"frameshift"{print $1,$2,$3,$4,$5,1}' temp2 > chr19.nonsyn.group