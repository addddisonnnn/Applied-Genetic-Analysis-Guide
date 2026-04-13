module load gcta/1.94.1
module load R

##UKBB BMI summary statistics from Locke et al:
head /projectnb/bs859/data/bmi/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt
## format of summary statistics required by GCTA gsmr is the same as for
##the conditional analysis we did a few weeks ago:  
## SNP A1 A2 freq b se p N
## where SNP=snp name
##A1 is the coded allele, A2 is the reference allele
##freq is the frequency of the A1 allele , b,se, and p are the regression
##estimate, SE, and p-value
##
##for this data set, these are columns 3-10

##reformat the BMI summary statistics for GSMR:
awk 'NR==1{print "SNP A1 A2 freq b se p N"};NR>1{print $3,$4,$5,$6,$7,$8,$9,$10}'  /projectnb/bs859/data/bmi/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt > bmi.ss.txt

##there are some duplicate SNPs in the file.  We don't need to worry 
##about deleting something important, as there are so many associated variants
##to use as instruments.
##this will remove lines that have a duplicate first column (snp name)
sort -t ' ' -k 1,1 -u bmi.ss.txt > bmi.ss.dedup.txt
#-t ' ' fields are separated by spaces
#-k 1,1: only look at the first field
#-u: delete duplicates based on that first field

wc bmi.ss.txt
wc bmi.ss.dedup.txt
##keep the de duplicated version, in original file name
mv bmi.ss.dedup.txt bmi.ss.txt



##Type 2 Diabetes (T2D) summary statistics
##1. Xue et al 2018
head /projectnb/bs859/data/T2D/T2D_Xue_et_al_2018.txt
# reformat Xue summary stats for GSMR
awk 'NR==1{print "SNP A1 A2 freq b se p N"};NR>1{print $3,$4,$5,$6,$7,$8,$9,$10}'  /projectnb/bs859/data/T2D/T2D_Xue_et_al_2018.txt > T2D.ss.txt
sort -t ' ' -k 1,1 -u T2D.ss.txt > T2D.ss.dedup.txt
wc T2D.ss.txt T2D.ss.dedup.txt
mv T2D.ss.dedup.txt T2D.ss.txt

##2. Mahajan et al 2018 Exome variant GWAS:
head /projectnb/bs859/data/T2D/T2D_European.BMIunadjusted.txt
##There are a lot of variants with no rsid (SNP name is ".").  
##We want to skip those when
##reformatting the file:
awk 'NR==1{print "SNP A1 A2 freq b se p N"};(NR>1&&$1!="."){print $1,$5,$6,$7,$9,$10,$11,$8}'  /projectnb/bs859/data/T2D/T2D_European.BMIunadjusted.txt > T2De.ss.txt
sort -t ' ' -k 1,1 -u T2De.ss.txt > T2De.ss.dedup.txt
wc T2De.ss.txt T2De.ss.dedup.txt
mv T2De.ss.dedup.txt T2De.ss.txt



##write the exposure and outcome file names to files for gsmr to read:
echo "bmi bmi.ss.txt" > exposure.txt

echo "T2D T2D.ss.txt" >outcome.txt
echo "T2De T2De.ss.txt" >>outcome.txt  
## note that >> writes to a file without overwriting what is already in it

###run gsmr to estimate causal effect of bmi on T2D using two different 
###T2D summary statistics files (first is GWAS, second is exome-only GWAS)
##Use 1000 Genomes Europeans to estimate LD among the SNPs (this is needed
##to eliminate SNPs that are in high LD)
gcta64 --gsmr-file exposure.txt outcome.txt --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR --gsmr-direction 0 --out BMI-T2D-gsmr --effect-plot

##What output files were created?
ls BMI-T2D-gsmr

##check log file for errors first!
more T2D-BMI-gsmr.log
##SNPs removed due to HEIDI outlier procedure (pleiotropic):
cat T2D-BMI-gsmr.pleio_snps
##results:
cat T2D-BMI-gsmr.gsmr



gcta64 --gsmr-file exposure.txt outcome.txt --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR --gsmr-direction 1 --out T2D-BMI-gsmr --effect-plot




##note:  to run both forward and backward in the same run, you can
##use --gsmr-direction 2

##re-run with bi-directional MR and with HEIDI-outlier
gcta64 --gsmr-file exposure.txt outcome.txt --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR --gsmr-direction 2 --out T2D-BMI-gsmr-bidir-02Heidithresh --effect-plot  --heidi-thresh 0.05

