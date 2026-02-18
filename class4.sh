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
6497881 84472466 476068950
# six million variants

##The SNP information file has MAF, imputatation quality R2, and
##other metrics:
zcat $CHR19DIR/chr19.info.gz |head
SNP	REF(0)	ALT(1)	ALT_Frq	MAF	AvgCall	Rsq	Genotyped	LooRsq	EmpR	EmpRsq	Dose0	Dose1
chr19:60547:T:C	T	C	0.00000	0.00000	1.00000	0.00027	Imputed	-	--	-	-
chr19:60562:C:G	C	G	0.00000	0.00000	1.00000	0.00018	Imputed	-	--	-	-
chr19:60585:G:A	G	A	0.00001	0.00001	0.99999	0.00021	Imputed	-	--	-	-
chr19:60636:GAC:G	GAC	G	0.00003	0.00003	0.99997	0.00054	Imputed--	-	-	-
chr19:60649:C:CT	C	CT	0.00018	0.00018	0.99982	0.00209	Imputed--	-	-	-
chr19:60658:G:A	G	A	0.00001	0.00001	0.99999	0.00056	Imputed	-	--	-	-
chr19:60668:G:A	G	A	0.00001	0.00001	0.99999	0.00236	Imputed	-	--	-	-
chr19:60759:ACT:A	ACT	A	0.00000	0.00000	1.00000	0.00006	Imputed--	-	-	-
chr19:60770:TG:T	TG	T	0.00002	0.00002	0.99998	0.00077	Imputed--	-	-	-
##There are *a lot* of variants with MAF<<0.01 and R2<<0.3!


##Before we do an association anslysis, we will take a look
##at the R2 and MAF distribution of the SNPs and make sure that
##they are close to what we expect  
module load R
Rscript MAFinfo.R > MAFinfo.log

zcat $CHR19DIR/chr19.dose.vcf.gz|head -n 21 |cut -f 1-10
#ouput from running this zcat command
##fileformat=VCFv4.1
##filedate=2022.2.23
##contig=<ID=chr19>
##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated Alternate Allele Frequency">
##INFO=<ID=MAF,Number=1,Type=Float,Description="Estimated Minor Allele Frequency">
##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square)">
##INFO=<ID=ER2,Number=1,Type=Float,Description="Empirical (Leave-One-Out) R-square (available only for genotyped variants)">
##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="Marker was imputed but NOT genotyped">
##INFO=<ID=TYPED,Number=0,Type=Flag,Description="Marker was genotyped AND imputed">
##INFO=<ID=TYPED_ONLY,Number=0,Type=Flag,Description="Marker was genotyped but NOT imputed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]">
##FORMAT=<ID=HDS,Number=2,Type=Float,Description="Estimated Haploid Alternate Allele Dosage ">
##FORMAT=<ID=GP,Number=3,Type=Float,Description="Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 ">
##pipeline=michigan-imputationserver-1.6.0
##imputation=minimac4-1.0.2
##phasing=eagle-2.4
##panel=apps@topmed-r2@1.0.0
##r2Filter=0.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	MAYO_10139
chr19	60547	chr19:60547:T:C	T	C	.	PASS	AF=0.00000;MAF=0.00000;R2=0.00027;IMPUTED	GT:DS:HDS:GP	0|0:0:0,0:1,0,0

zcat $CHR19DIR/chr19.info.gz|head -n 2
SNP	REF(0)	ALT(1)	ALT_Frq	MAF	AvgCall	Rsq	Genotyped	LooRsq	EmpR	EmpRsq	Dose0	Dose1
chr19:60547:T:C	T	C	0.00000	0.00000	1.00000	0.00027	Imputed	-	--	-	-


##Take a look at the imputed data vcf (its a gzipped text file):
zcat $CHR19DIR/chr19.dose.vcf.gz|head -n 21|cut -f 1-10
##fileformat=VCFv4.1
##filedate=2022.2.23
##contig=<ID=chr19>
##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated Alternate Allele Frequency">
##INFO=<ID=MAF,Number=1,Type=Float,Description="Estimated Minor Allele Frequency">
##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square)">
##INFO=<ID=ER2,Number=1,Type=Float,Description="Empirical (Leave-One-Out) R-square (available only for genotyped variants)">
##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="Marker was imputed but NOT genotyped">
##INFO=<ID=TYPED,Number=0,Type=Flag,Description="Marker was genotyped AND imputed">
##INFO=<ID=TYPED_ONLY,Number=0,Type=Flag,Description="Marker was genotyped but NOT imputed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]">
##FORMAT=<ID=HDS,Number=2,Type=Float,Description="Estimated Haploid Alternate Allele Dosage ">
##FORMAT=<ID=GP,Number=3,Type=Float,Description="Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 ">
##pipeline=michigan-imputationserver-1.6.0
##imputation=minimac4-1.0.2
##phasing=eagle-2.4
##panel=apps@topmed-r2@1.0.0
##r2Filter=0.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	MAYO_10139
chr19	60547	chr19:60547:T:C	T	C	.	PASS	AF=0.00000;MAF=0.00000;R2=0.00027;IMPUTED	GT:DS:HDS:GP	0|0:0:0,0:1,0,0

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
PLINK v2.00a3.7LM 64-bit Intel (24 Oct 2022)   www.cog-genomics.org/plink/2.0/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to tgen_chr19_imputed.log.
Options in effect:
  --double-id
  --exclude-if-info MAF<0.005
  --extract-if-info R2>=0.3
  --make-pgen
  --out tgen_chr19_imputed
  --pheno /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.fam
  --pheno-col-nums 6
  --vcf /projectnb/bs859/data/tgen/2023/chr19//chr19.dose.vcf.gz dosage=DS

Start time: Wed Feb 11 15:36:22 2026
257325 MiB RAM detected; reserving 128662 MiB for main workspace.
Using up to 28 threads (change this with --threads).
--vcf: 6497881 variants scanned.
--vcf: tgen_chr19_imputed-temporary.pgen +
tgen_chr19_imputed-temporary.pvar.zst + tgen_chr19_imputed-temporary.psam
written.
1228 samples (0 females, 0 males, 1228 ambiguous; 1228 founders) loaded from
tgen_chr19_imputed-temporary.psam.
218897 out of 6497881 variants loaded from
tgen_chr19_imputed-temporary.pvar.zst.
1 binary phenotype loaded (767 cases, 460 controls).
218897 variants remaining after main filters.
Writing tgen_chr19_imputed.psam ... done.
Writing tgen_chr19_imputed.pvar ... done.
Writing tgen_chr19_imputed.pgen ... done.
End time: Wed Feb 11 15:38:17 2026

head tgen_chr19_imputed.psam
#FID	IID	SEX	PHENO1
MAYO_10139	MAYO_10139	NA	1
MAYO_10198	MAYO_10198	NA	2
MAYO_102246	MAYO_102246	NA	2
MAYO_10249	MAYO_10249	NA	1
MAYO_10278	MAYO_10278	NA	2
MAYO_10304	MAYO_10304	NA	1
MAYO_10316	MAYO_10316	NA	1
MAYO_10367	MAYO_10367	NA	2
MAYO_10403	MAYO_10403	NA	2

head /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.fam
MAYO_10139 MAYO_10139 0 0 0 1
MAYO_10198 MAYO_10198 0 0 0 2
MAYO_102246 MAYO_102246 0 0 0 2
MAYO_10249 MAYO_10249 0 0 0 1
MAYO_10278 MAYO_10278 0 0 0 2
MAYO_10304 MAYO_10304 0 0 0 1
MAYO_10316 MAYO_10316 0 0 0 1
MAYO_10367 MAYO_10367 0 0 0 2
MAYO_10403 MAYO_10403 0 0 0 2
MAYO_10529 MAYO_10529 0 0 0 2

head tgen_chr19_imputed.pvar -n 20
##filedate=2022.2.23
##contig=<ID=chr19>
##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated Alternate Allele Frequency">
##INFO=<ID=MAF,Number=1,Type=Float,Description="Estimated Minor Allele Frequency">
##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square)">
##INFO=<ID=ER2,Number=1,Type=Float,Description="Empirical (Leave-One-Out) R-square (available only for genotyped variants)">
##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="Marker was imputed but NOT genotyped">
##INFO=<ID=TYPED,Number=0,Type=Flag,Description="Marker was genotyped AND imputed">
##INFO=<ID=TYPED_ONLY,Number=0,Type=Flag,Description="Marker was genotyped but NOT imputed">
##pipeline=michigan-imputationserver-1.6.0
##imputation=minimac4-1.0.2
##phasing=eagle-2.4
##panel=apps@topmed-r2@1.0.0
##r2Filter=0.0
#CHROM	POS	ID	REF	ALT	FILTER	INFO
19	225487	chr19:225487:C:T	C	T	PASS	AF=0.00590;MAF=0.00590;R2=0.37044;IMPUTED
19	229783	chr19:229783:G:A	G	A	PASS	AF=0.23020;MAF=0.23020;R2=0.52562;IMPUTED
19	230529	chr19:230529:G:A	G	A	PASS	AF=0.00908;MAF=0.00908;R2=0.37094;IMPUTED
19	231076	chr19:231076:TAAAAAG:T	TAAAAAG	T	PASS	AF=0.03152;MAF=0.03152;R2=0.45045;IMPUTED
19	238434	chr19:238434:C:G	C	G	PASS	AF=0.00804;MAF=0.00804;R2=0.42738;IMPUTED

head /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.bim
1	rs3094315	0	742429	C	T
1	rs4040617	0	769185	G	A
1	rs2980300	0	775852	T	C
1	rs4075116	0	993492	G	A
1	rs10907175	0	1120590	C	A
1	rs6603781	0	1148494	A	G
1	rs11260562	0	1155173	A	G
1	rs6685064	0	1201155	T	C
1	rs3766180	0	1468016	G	A
1	rs6603791	0	1490804	G	A

cat tgen_chr19_imputed.log
PLINK v2.00a3.7LM 64-bit Intel (24 Oct 2022)
Options in effect:
  --double-id
  --exclude-if-info MAF<0.005
  --extract-if-info R2>=0.3
  --make-pgen
  --out tgen_chr19_imputed
  --pheno /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.fam
  --pheno-col-nums 6
  --vcf /projectnb/bs859/data/tgen/2023/chr19//chr19.dose.vcf.gz dosage=DS

Hostname: scc-wl4
Working directory: /projectnb/bs859/students/addisony/Class04
Start time: Wed Feb 11 15:36:22 2026

Random number seed: 1770842182
257325 MiB RAM detected; reserving 128662 MiB for main workspace.
Using up to 28 threads (change this with --threads).
--vcf: 6497881 variants scanned.
--vcf: tgen_chr19_imputed-temporary.pgen +
tgen_chr19_imputed-temporary.pvar.zst + tgen_chr19_imputed-temporary.psam
written.
1228 samples (0 females, 0 males, 1228 ambiguous; 1228 founders) loaded from
tgen_chr19_imputed-temporary.psam.
218897 out of 6497881 variants loaded from
tgen_chr19_imputed-temporary.pvar.zst.
1 binary phenotype loaded (767 cases, 460 controls).
218897 variants remaining after main filters.
Writing tgen_chr19_imputed.psam ... done.
Writing tgen_chr19_imputed.pvar ... done.
Writing tgen_chr19_imputed.pgen ... done.

End time: Wed Feb 11 15:38:17 2026


##In the last homework, we found that PCs 6 and 8 were associated with AD
##with p<0.01.  We will use these 2 PCs as covariates in our analyses.  
plink2 --pfile tgen_chr19_imputed \
--covar /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt \
--covar-name PC6,PC8 --logistic hide-covar --out chr19_01pcs 
PLINK v2.00a3.7LM 64-bit Intel (24 Oct 2022)   www.cog-genomics.org/plink/2.0/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to chr19_01pcs.log.
Options in effect:
  --covar /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt
  --covar-name PC6,PC8
  --glm hide-covar
  --out chr19_01pcs
  --pfile tgen_chr19_imputed

Start time: Wed Feb 11 16:02:23 2026
257325 MiB RAM detected; reserving 128662 MiB for main workspace.
Using up to 28 threads (change this with --threads).
1228 samples (0 females, 0 males, 1228 ambiguous; 1228 founders) loaded from
tgen_chr19_imputed.psam.
218897 variants loaded from tgen_chr19_imputed.pvar.
1 binary phenotype loaded (767 cases, 460 controls).
2 covariates loaded from /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt.
Calculating allele frequencies... done.
--glm logistic-Firth hybrid regression on phenotype 'PHENO1': done.
Results written to chr19_01pcs.PHENO1.glm.logistic.hybrid .
End time: Wed Feb 11 16:02:24 2026

head chr19_01pcs.PHENO1.glm.logistic.hybrid
#CHROM	POS	ID	REF	ALT	A1	FIRTH?	TEST	OBS_CT	OR	LOG(OR)_SE	Z_STAT	P	ERRCODE
19	225487	chr19:225487:C:T	C	T	T	N	ADD	1227	0.846543	0.898871	-0.185337	0.852964	.
19	229783	chr19:229783:G:A	G	A	A	N	ADD	1227	0.948903	0.138485	-0.378729	0.704889	.
19	230529	chr19:230529:G:A	G	A	A	N	ADD	1227	0.401842	0.704174	-1.2947	0.195423	.
19	231076	chr19:231076:TAAAAAG:T	TAAAAAG	T	T	N	ADD	1227	0.678	0.349437	-1.1121	0.266096	.
19	238434	chr19:238434:C:G	C	G	G	N	ADD	1227	0.395376	0.696712	-1.33185	0.182909	.
19	240963	chr19:240963:G:A	G	A	A	N	ADD	1227	1.72325	0.568867	0.956657	0.33874	.
19	241380	chr19:241380:C:T	C	T	T	N	ADD	1227	1.03808	0.139354	0.268153	0.788582	.
19	243192	chr19:243192:A:G	A	G	G	N	ADD	1227	0.836628	0.808702	-0.220571	0.825427	.
19	244421	chr19:244421:A:G	A	G	G	N	ADD	1227	0.99056	0.148902	-0.0636999	0.949209	.

module load R
Rscript  qqplot_pgen.R chr19_01pcs.PHENO1.glm.logistic.hybrid chr19.tgen ADD
[1] "chr19_01pcs.PHENO1.glm.logistic.hybrid"
[2] "chr19.tgen"                            
[3] "ADD"                                   
null device 
          1 

Rscript  gwaplot_pgen.R chr19_01pcs.PHENO1.glm.logistic.hybrid "chr19.imputed" "chr19.imputed" 
[1] "chr19_01pcs.PHENO1.glm.logistic.hybrid"
[2] "chr19.imputed"                         
[3] "chr19.imputed"                         
null device 
          1 

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
PLINK v2.00a3.7LM 64-bit Intel (24 Oct 2022)   www.cog-genomics.org/plink/2.0/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to chr19_GENOTYPED.log.
Options in effect:
  --covar /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt
  --covar-name PC6,PC8
  --double-id
  --glm hide-covar
  --out chr19_GENOTYPED
  --pheno /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.fam
  --pheno-col-nums 6
  --vcf /projectnb/bs859/data/tgen/2023/chr19//chr19.vcf

Start time: Wed Feb 11 16:12:34 2026
257325 MiB RAM detected; reserving 128662 MiB for main workspace.
Using up to 28 threads (change this with --threads).
--vcf: 4069 variants scanned.
--vcf: chr19_GENOTYPED-temporary.pgen + chr19_GENOTYPED-temporary.pvar.zst +
chr19_GENOTYPED-temporary.psam written.
1228 samples (0 females, 0 males, 1228 ambiguous; 1228 founders) loaded from
chr19_GENOTYPED-temporary.psam.
4069 variants loaded from chr19_GENOTYPED-temporary.pvar.zst.
1 binary phenotype loaded (767 cases, 460 controls).
2 covariates loaded from /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt.
Calculating allele frequencies... done.
--glm logistic-Firth hybrid regression on phenotype 'PHENO1': done.
Results written to chr19_GENOTYPED.PHENO1.glm.logistic.hybrid .
End time: Wed Feb 11 16:12:34 2026

Rscript  qqplot_pgen.R chr19_GENOTYPED.PHENO1.glm.logistic.hybrid chr19.tgen.genotyped ADD
[1] "chr19_GENOTYPED.PHENO1.glm.logistic.hybrid"
[2] "chr19.tgen.genotyped"                      
[3] "ADD"                                       
null device 
          1 

Rscript  gwaplot_pgen.R chr19_GENOTYPED.PHENO1.glm.logistic.hybrid chr19.tgen.genotyped "chr19.genotyped"
[1] "chr19_GENOTYPED.PHENO1.glm.logistic.hybrid"
[2] "chr19.tgen.genotyped"                      
[3] "chr19.genotyped"                           
null device 
          1 
#the genotyped is not as crowded as the imputed plots

##since the chromosome position is different in the genotyped-only and imputed
##genotypes files, we'd need to look up the variants to determine which 
##which positions match up.
##We can start by just looking at the most significant variants (APOE region)

sort -gk13 chr19_GENOTYPED.PHENO1.glm.logistic.hybrid |head -n 5|cut -f 1-6,10-13
#CHROM	POS	ID	REF	ALT	A1	OR	LOG(OR)_SE	Z_STAT	P
19	45422946	rs4420638	A	G	G	3.22347	0.107503	10.8877	1.31997e-27
19	38313102	rs8110293	A	G	G	0.711935	0.0919062	-3.6969	0.000218247
19	21703545	rs2650757	A	G	G	0.667809	0.109466	-3.68838	0.000225686
19	54815093	rs741584	G	C	C	1.3924	0.0952339	3.47599	0.000508967
# get the top five variances and columns of interest, strongest variant has a beta of 3.22347

sort -gk13 chr19_01pcs.PHENO1.glm.logistic.hybrid |head -n 5|cut -f 1-6,10-13
#CHROM	POS	ID	REF	ALT	A1	OR	LOG(OR)_SE	Z_STAT	P
19	44908684	chr19:44908684:T:C	T	C	C	4.78809	0.141816	11.0434	2.36026e-28
19	44919689	chr19:44919689:A:G	A	G	G	3.23839	0.10802	10.8783	1.46303e-27
19	44888997	chr19:44888997:C:T	C	T	T	4.66369	0.142942	10.7723	4.65426e-27
19	44919589	chr19:44919589:G:A	G	A	A	3.63223	0.122104	10.5635	4.40049e-26
# we have four variants small p-values (which are in the same region)
# can look up variant rs4420638 by going to dbsnp database website, we can find the location 44919689
# But why do they have different ORs? theyre similar but the imputed is filled in with the dosages, so imputed reflects more uncertainty, bigger standard error

# if you look up 19:44908684, showing it's in the APOE, go to Clinical significance, the c coded allele has an OR>1 --> increased odds of AD