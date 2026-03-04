module load metal
DATADIR=/projectnb/bs859/data/meta/METAL_example
ls $DATADIR

head $DATADIR/DGI_three_regions.txt
head $DATADIR/magic_SARDINIA.tbl
zcat $DATADIR/MAGIC_FUSION_Results.txt.gz |head

more metal.txt

metal metal.txt > metal.log
cat metal.log
cat ZscoreMeta1.tbl.info

head ZscoreMeta1.tbl

##rs560887 has the smallest p-value (see log file).
##What do the meta-analysis and study-specific results look
##like for this variant?
awk 'NR==1||$1=="rs560887"{print $0}' ZscoreMeta1.tbl

##which allele of rs560887 increases glucose levels?


##get the study specific results for rs560887:
#For the FUSION study and DGI study, the SNP name is in column 3:
zcat $DATADIR/MAGIC_FUSION_Results.txt.gz |awk 'NR==1||$3=="rs560887" {print $0}'

awk 'NR==1||$3=="rs560887" {print $0}' $DATADIR/DGI_three_regions.txt

##for the Sardinia study, the SNP name is in column 1:
awk 'NR==1||$1=="rs560887" {print $0}' $DATADIR/magic_SARDINIA.tbl


##
##Alteration 1:  uncomment the genomic control lines to
##rescale the test statistics accounting for lambda_GC>1.
##change the output file name, and save as a new script.

metal metal_GC.txt > metal_GC.log
cat metal_GC.log

##Alteration 2: change to inverse variance weighted meta-analysis and add
##options to compute min, mean, and max allele frquencies for each variant


