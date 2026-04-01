export IGAPDIR='/projectnb/bs859/data/igap/'
export TGENDIR='/projectnb/bs859/data/tgen/cleaned/'

ls $IGAPDIR
zcat $IGAPDIR/Kunkle_etal_Stage1_results2019.txt.gz | head

ls $TGENDIR

head /projectnb/bs859/data/cognitive/Davies_NC_2018_OPEN_DATASET/Davies2018_OPEN_DATASET_summary_results.txt




## how are the PRS results related to genetic correlation?
##Ldscore regression for genetic correlation of the two traits:
## we have the reformatted AD summary statistics file from last week:

zcat ../class09/AD_sumstats.gz |head

##create sumstats file for Davies et al cognitive gwas:
module load R
module load python2
module load ldsc
export LDSCORES_DIR='/projectnb/bs859/data/ldscore_files'

munge_sumstats.py \
--sumstats /projectnb/bs859/data/cognitive/Davies_NC_2018_OPEN_DATASET/Davies2018_OPEN_DATASET_summary_results.txt  \
--snp MarkerName \
--N 250000 \
--a1 Effect_allele \
--a2 Other_allele \
--signed-sumstats Zscore,0 \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out Cog



#heritability of the cognitive performance:
ldsc.py \
--h2 Cog.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--out Cog.h2

#genetic correlation of AD and cognitive performance:
ldsc.py \
--rg ../class09/AD.sumstats.gz,Cog.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--out AD_COG_rg

