###Class 3 in-class analyses 
### February 2026

module load R 
module load plink/1.90b6.27
##make an alias for the directory with the data we are using today:
export DATADIR=/projectnb/bs859/data/plink_tutorial/wgas3

# Simple regression with no covariate
> plink --bfile $DATADIR/wgas3  --logistic  beta --ci .95 --out logistnoadj
PLINK v1.90b6.27 64-bit (10 Dec 2022)          www.cog-genomics.org/plink/1.9/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to logistnoadj.log.
Options in effect:
  --bfile /projectnb/bs859/data/plink_tutorial/wgas3/wgas3
  --ci .95
  --logistic beta
  --out logistnoadj

257325 MB RAM detected; reserving 128662 MB for main workspace.
179493 variants loaded from .bim file.
89 people (44 males, 45 females) loaded from .fam.
89 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 89 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.996308.
179493 variants and 89 people pass filters and QC.
Among remaining phenotypes, 48 are cases and 41 are controls.
Writing logistic model association results to logistnoadj.assoc.logistic ...
done.

> head  logistnoadj.assoc.logistic
CHR         SNP         BP   A1       TEST    NMISS       BETA       SE      L95      U95         STAT            P 
   1   rs3094315     792429    G        ADD       88      0.723   0.5235  -0.3029    1.749        1.381       0.1672
   1   rs4040617     819185    G        ADD       89     0.5901    0.527  -0.4429    1.623         1.12       0.2629
   1   rs4075116    1043552    C        ADD       89    -0.6343   0.6841   -1.975   0.7065      -0.9272       0.3538
   1   rs9442385    1137258    T        ADD       88    -0.1988   0.2894   -0.766   0.3683       -0.687       0.4921
   1  rs11260562    1205233    A        ADD       87    -0.5521   0.9393   -2.393    1.289      -0.5877       0.5567
   1   rs6685064    1251215    C        ADD       89    -0.1733   0.2708   -0.704   0.3575      -0.6398       0.5223
   1   rs3766180    1563420    T        ADD       89     0.6822   0.4615  -0.2224    1.587        1.478       0.1394
   1   rs6603791    1586208    A        ADD       89     0.8213   0.4806  -0.1206    1.763        1.709      0.08745
   1   rs7519837    1596068    C        ADD       88     0.7637    0.481   -0.179    1.706        1.588       0.1123
# in this we're getting the phenotype coded allele
# what are we missing? This won't be able to tell us the other allele, only one seen in A1. This would be stored in the bim file

# see what is in the covar txt
> head /projectnb/bs859/data/plink_tutorial/wgas3/newcov.txt
FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 pop
CH18526	NA18526 -0.1061 0.0520 0.0853 0.0530 0.1469 -0.2413 0.0869 -0.0844 -0.0277 -0.0396  1
CH18524	NA18524 -0.1113 -0.1509 0.1828 0.0854 -0.0094 0.0079 0.0486 -0.2673 0.2114 -0.2010  1
CH18529	NA18529 -0.0970 -0.0232 -0.0743 0.1261 -0.0483 -0.0803 0.1161 -0.1977 -0.1706 0.0351  1
CH18558	NA18558 -0.0967 -0.1336 -0.0604 0.0008 0.2042 -0.0620 -0.1064 -0.1895 -0.0023 0.1040  1
CH18532	NA18532 -0.0889 0.1528 -0.0739 0.1461 0.0167 -0.0366 -0.1674 -0.1467 0.2217 -0.0785  1
CH18561	NA18561 -0.1054 -0.2970 -0.0717 -0.0170 -0.2613 -0.2161 0.1400 -0.0881 -0.1362 0.0828  1
CH18562	NA18562 -0.1000 0.0151 -0.1185 -0.0186 0.0808 0.1350 -0.1409 -0.2131 -0.1136 0.1655  1
CH18537	NA18537 -0.1039 -0.1011 0.0710 0.2528 -0.0688 0.1604 0.0862 0.0676 -0.0570 0.0326  1
CH18603	NA18603 -0.1091 0.0504 -0.0796 -0.0238 0.0039 -0.0666 -0.0160 0.1325 -0.0268 -0.0577  1

# adjusting for the covariate name to pop (chinese versus japanese)
> plink   --bfile $DATADIR/wgas3 --covar $DATADIR/newcov.txt --covar-name pop --logistic beta  --ci .95 --out logistadjpop
PLINK v1.90b6.27 64-bit (10 Dec 2022)          www.cog-genomics.org/plink/1.9/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to logistadjpop.log.
Options in effect:
  --bfile /projectnb/bs859/data/plink_tutorial/wgas3/wgas3
  --ci .95
  --covar /projectnb/bs859/data/plink_tutorial/wgas3/newcov.txt
  --covar-name pop
  --logistic beta
  --out logistadjpop

257325 MB RAM detected; reserving 128662 MB for main workspace.
179493 variants loaded from .bim file.
89 people (44 males, 45 females) loaded from .fam.
89 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
--covar: 1 out of 11 covariates loaded.
Before main variant filters, 89 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.996308.
179493 variants and 89 people pass filters and QC.
Among remaining phenotypes, 48 are cases and 41 are controls.
Writing logistic model association results to logistadjpop.assoc.logistic ...
done.

> head logistadjpop.assoc.logistic
 CHR         SNP         BP   A1       TEST    NMISS       BETA       SE      L95      U95         STAT            P 
   1   rs3094315     792429    G        ADD       88     0.7042    0.641  -0.5521    1.961        1.099       0.2719
   1   rs3094315     792429    G        pop       88      2.762   0.5446    1.695     3.83        5.072    3.941e-07
   1   rs4040617     819185    G        ADD       89      0.642   0.6499  -0.6318    1.916       0.9878       0.3233
   1   rs4040617     819185    G        pop       89      2.803   0.5438    1.737    3.869        5.154    2.553e-07
   1   rs4075116    1043552    C        ADD       89     -1.029   0.8637   -2.722   0.6637       -1.191       0.2335
   1   rs4075116    1043552    C        pop       89       2.86   0.5544    1.774    3.947        5.159    2.478e-07
   1   rs9442385    1137258    T        ADD       88    -0.2138   0.3593   -0.918   0.4903      -0.5952       0.5517
   1   rs9442385    1137258    T        pop       88      2.769   0.5413    1.708     3.83        5.116    3.124e-07
   1  rs11260562    1205233    A        ADD       87     0.2275    1.084   -1.897    2.352       0.2099       0.8338
# the difference we see pop in the test column, this gives the two variables but not the intercept (logistic regression)
# we don't need the variable info, how do we get the results we want

#two methods to do so, professor didn't run it but also an option
awk 'NR==1||$5=="ADD"{print $0}' logistadjpop.assoc.logistic>popadj.add.txt
head popadj.add.txt

#did this one in class, telling no covariates
> plink   --bfile $DATADIR/wgas3 --covar $DATADIR/newcov.txt --covar-name pop --logistic beta hide-covar  --ci 0.95 --out logistadjpop-hidecov
PLINK v1.90b6.27 64-bit (10 Dec 2022)          www.cog-genomics.org/plink/1.9/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to logistadjpop-hidecov.log.
Options in effect:
  --bfile /projectnb/bs859/data/plink_tutorial/wgas3/wgas3
  --ci 0.95
  --covar /projectnb/bs859/data/plink_tutorial/wgas3/newcov.txt
  --covar-name pop
  --logistic beta hide-covar
  --out logistadjpop-hidecov

257325 MB RAM detected; reserving 128662 MB for main workspace.
179493 variants loaded from .bim file.
89 people (44 males, 45 females) loaded from .fam.
89 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
--covar: 1 out of 11 covariates loaded.
Before main variant filters, 89 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.996308.
179493 variants and 89 people pass filters and QC.
Among remaining phenotypes, 48 are cases and 41 are controls.
Writing logistic model association results to
logistadjpop-hidecov.assoc.logistic ... done.

> head logistadjpop-hidecov.assoc.logistic
CHR         SNP         BP   A1       TEST    NMISS       BETA       SE      L95      U95         STAT            P 
   1   rs3094315     792429    G        ADD       88     0.7042    0.641  -0.5521    1.961        1.099       0.2719
   1   rs4040617     819185    G        ADD       89      0.642   0.6499  -0.6318    1.916       0.9878       0.3233
   1   rs4075116    1043552    C        ADD       89     -1.029   0.8637   -2.722   0.6637       -1.191       0.2335
   1   rs9442385    1137258    T        ADD       88    -0.2138   0.3593   -0.918   0.4903      -0.5952       0.5517
   1  rs11260562    1205233    A        ADD       87     0.2275    1.084   -1.897    2.352       0.2099       0.8338
   1   rs6685064    1251215    C        ADD       89     0.4822   0.3801  -0.2627    1.227        1.269       0.2045
   1   rs3766180    1563420    T        ADD       89      1.121   0.5813 -0.01834     2.26        1.928       0.0538
   1   rs6603791    1586208    A        ADD       89      1.233   0.5978  0.06124    2.405        2.062      0.03917
   1   rs7519837    1596068    C        ADD       88      1.212   0.6003  0.03534    2.388        2.019       0.0435


> plink   --bfile $DATADIR/wgas3 --covar $DATADIR/newcov.txt --covar-name PC1,PC2 --logistic beta  --ci .95 --out logistadjPC12
PLINK v1.90b6.27 64-bit (10 Dec 2022)          www.cog-genomics.org/plink/1.9/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to logistadjPC12.log.
Options in effect:
  --bfile /projectnb/bs859/data/plink_tutorial/wgas3/wgas3
  --ci .95
  --covar /projectnb/bs859/data/plink_tutorial/wgas3/newcov.txt
  --covar-name PC1,PC2
  --logistic beta
  --out logistadjPC12

257325 MB RAM detected; reserving 128662 MB for main workspace.
179493 variants loaded from .bim file.
89 people (44 males, 45 females) loaded from .fam.
89 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
--covar: 2 out of 11 covariates loaded.
Before main variant filters, 89 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.996308.
179493 variants and 89 people pass filters and QC.
Among remaining phenotypes, 48 are cases and 41 are controls.
Writing logistic model association results to logistadjPC12.assoc.logistic ...
done.

> head logistadjPC12.assoc.logistic
CHR         SNP         BP   A1       TEST    NMISS       BETA       SE      L95      U95         STAT            P 
   1   rs3094315     792429    G        ADD       88     0.6249   0.6496  -0.6483    1.898        0.962       0.3361
   1   rs3094315     792429    G        PC1       88      13.89    2.713    8.575    19.21         5.12     3.05e-07
   1   rs3094315     792429    G        PC2       88     -1.776    2.589    -6.85    3.298       -0.686       0.4927
   1   rs4040617     819185    G        ADD       89     0.5399   0.6612  -0.7559    1.836       0.8166       0.4141
   1   rs4040617     819185    G        PC1       89      14.09    2.715    8.772    19.42        5.191    2.094e-07
   1   rs4040617     819185    G        PC2       89     -1.704     2.58   -6.762    3.353      -0.6605       0.5089
   1   rs4075116    1043552    C        ADD       89     -1.279   0.9115   -3.065   0.5079       -1.403       0.1607
   1   rs4075116    1043552    C        PC1       89      14.59      2.8    9.103    20.08        5.211    1.877e-07
   1   rs4075116    1043552    C        PC2       89     -2.042    2.575    -7.09    3.006       -0.793       0.4278

###when you run a logistic or linear regression in SAS or R, you get much more information than what the PLINK output provides
###What information is not included in the PLINK output that is provided in a full regression output?

###What information is missing from the PLINK output that you might need in order to fully interpret your results?  Where can you look
###for that information?


##normally, we'd only run as below, hiding the covariate lines. 
##the previous run was just to show you what we are supressing when we use
## hide-covar below:
plink   --bfile $DATADIR/wgas3 --covar $DATADIR/newcov.txt --covar-name PC1,PC2 --logistic beta hide-covar --ci .95 --out logistadjPC12-hidecov

head logistadjPC12-hidecov.assoc.logistic


## now we will try a mixed model using GMMAT with the plink tutorial data
##Get set up to do glmm model for binary trait (mixed model with genetic relationship matrix)
## Recall that 2 weeks ago, we computed pairwise IBD sharing using the --genome command, 
## after pruning the SNPs to get a set that was not highly correlated with each other.
## --genome created a file with each pair of individuals as a line, and the IBD 
## sharing estimates for each pair.

##prune to remove variants in high LD:
##note that in a real-world situation, we'd probably want to prune
##with a more stringent r2 (0.2 means that we are allowing correlation of
##as high as sqrt(0.2)=0.447, which is pretty high 
> plink --bfile $DATADIR/wgas3  --indep-pairwise 10000kb 1  0.2 
PLINK v1.90b6.27 64-bit (10 Dec 2022)          www.cog-genomics.org/plink/1.9/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to plink.log.
Options in effect:
  --bfile /projectnb/bs859/data/plink_tutorial/wgas3/wgas3
  --indep-pairwise 10000kb 1 0.2

257325 MB RAM detected; reserving 128662 MB for main workspace.
179493 variants loaded from .bim file.
89 people (44 males, 45 females) loaded from .fam.
89 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 89 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.996308.
179493 variants and 89 people pass filters and QC.
Among remaining phenotypes, 48 are cases and 41 are controls.
Pruned 10692 variants from chromosome 1, leaving 4558.
Pruned 10244 variants from chromosome 2, leaving 4272.
Pruned 8247 variants from chromosome 3, leaving 3647.
Pruned 6601 variants from chromosome 4, leaving 3204.
Pruned 8158 variants from chromosome 5, leaving 3450.
Pruned 8256 variants from chromosome 6, leaving 3227.
Pruned 6320 variants from chromosome 7, leaving 2910.
Pruned 6767 variants from chromosome 8, leaving 2866.
Pruned 6045 variants from chromosome 9, leaving 2674.
Pruned 7976 variants from chromosome 10, leaving 3076.
Pruned 7248 variants from chromosome 11, leaving 2767.
Pruned 6367 variants from chromosome 12, leaving 2893.
Pruned 4091 variants from chromosome 13, leaving 2063.
Pruned 3971 variants from chromosome 14, leaving 1864.
Pruned 3843 variants from chromosome 15, leaving 1856.
Pruned 4227 variants from chromosome 16, leaving 1973.
Pruned 3043 variants from chromosome 17, leaving 1761.
Pruned 3347 variants from chromosome 18, leaving 1757.
Pruned 1751 variants from chromosome 19, leaving 1165.
Pruned 3479 variants from chromosome 20, leaving 1562.
Pruned 1609 variants from chromosome 21, leaving 865.
Pruned 1834 variants from chromosome 22, leaving 967.
Pruning complete.  124116 of 179493 variants removed.
Marker lists written to plink.prune.in and plink.prune.out .

## We used --genome to get the pairwise IBD last week:
> plink --bfile $DATADIR/wgas3 --chr 1-22 --extract plink.prune.in --genome --out ibd 
PLINK v1.90b6.27 64-bit (10 Dec 2022)          www.cog-genomics.org/plink/1.9/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to ibd.log.
Options in effect:
  --bfile /projectnb/bs859/data/plink_tutorial/wgas3/wgas3
  --chr 1-22
  --extract plink.prune.in
  --genome
  --out ibd

257325 MB RAM detected; reserving 128662 MB for main workspace.
179493 variants loaded from .bim file.
89 people (44 males, 45 females) loaded from .fam.
89 phenotype values loaded from .fam.
--extract: 55377 variants remaining.
Using up to 27 threads (change this with --threads).
Before main variant filters, 89 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.996045.
55377 variants and 89 people pass filters and QC.
Among remaining phenotypes, 48 are cases and 41 are controls.
IBD calculations complete.  
Finished writing ibd.genome .

> wc ibd.genome
  3917  54838 423036 ibd.genome

##first 10 lines:
head ibd.genome
     FID1     IID1     FID2     IID2 RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO
  CH18526  NA18526  CH18524  NA18524 UN    NA  1.0000  0.0000  0.0000  0.0000  -1  0.745392  0.1208  1.9158
  CH18526  NA18526  CH18529  NA18529 UN    NA  0.9804  0.0000  0.0196  0.0196  -1  0.747905  0.3380  1.9695
  CH18526  NA18526  CH18558  NA18558 UN    NA  0.9875  0.0125  0.0000  0.0063  -1  0.745061  0.4218  1.9855
  CH18526  NA18526  CH18532  NA18532 UN    NA  0.9581  0.0419  0.0000  0.0210  -1  0.744858  0.9476  2.1255
  CH18526  NA18526  CH18561  NA18561 UN    NA  0.9917  0.0075  0.0008  0.0045  -1  0.745123  0.6848  2.0360
  CH18526  NA18526  CH18562  NA18562 UN    NA  1.0000  0.0000  0.0000  0.0000  -1  0.744175  0.0368  1.8731
  CH18526  NA18526  CH18537  NA18537 UN    NA  1.0000  0.0000  0.0000  0.0000   0  0.745130  0.1025  1.9093
  CH18526  NA18526  CH18603  NA18603 UN    NA  1.0000  0.0000  0.0000  0.0000   0  0.744447  0.4412  1.9891
  CH18526  NA18526  CH18540  NA18540 UN    NA  1.0000  0.0000  0.0000  0.0000  -1  0.743700  0.0206  1.8556


##We need the "pihat" information, but in a different format -- as a matrix where the 
##diagonal is the variance of each person, and the offdiagonals are the pairwise 
##covariance of genotypes for each pair of individuals.  
## --make-rel square is one way to achieve this:
> plink --bfile $DATADIR/wgas3 --chr 1-22 --extract plink.prune.in --make-rel square --out grm
PLINK v1.90b6.27 64-bit (10 Dec 2022)          www.cog-genomics.org/plink/1.9/
(C) 2005-2022 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to grm.log.
Options in effect:
  --bfile /projectnb/bs859/data/plink_tutorial/wgas3/wgas3
  --chr 1-22
  --extract plink.prune.in
  --make-rel square
  --out grm

257325 MB RAM detected; reserving 128662 MB for main workspace.
179493 variants loaded from .bim file.
89 people (44 males, 45 females) loaded from .fam.
89 phenotype values loaded from .fam.
--extract: 55377 variants remaining.
Using up to 27 threads (change this with --threads).
Before main variant filters, 89 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.996045.
55377 variants and 89 people pass filters and QC.
Among remaining phenotypes, 48 are cases and 41 are controls.
Relationship matrix calculation complete.
Relationship matrix written to grm.rel , and IDs written to grm.rel.id .

##--make-grm-gz square makes the same file, but gzipped (compressed).  This is 
##not necessary for this small tutorial data set, but is very important for larger datasets!
##plink --bfile $DATADIR/wgas3 --exclude plink.prune.out --make-grm-gz square --out grm1
#
#creates 2 files: 
# grm.rel.id -- the order of the ids for the GRM file

> head grm.rel.id
CH18526	NA18526
CH18524	NA18524
CH18529	NA18529
CH18558	NA18558
CH18532	NA18532
CH18561	NA18561
CH18562	NA18562
CH18537	NA18537
CH18603	NA18603
CH18540	NA18540

# grm.rel -- the actual genetic relationship matrix.  For our data, with 89 individuals, 
# it is an 89x89  matrix.  Let's look at the first 5 individuals submatrix:

> cut -f1-5 grm.rel |head -n 5
0.999166	-0.000686896	-1.64415e-05	-0.00294458	-0.00432704
-0.000686896	1.00915	-0.00825227	-0.000538692	-0.00110805
-1.64415e-05	-0.00825227	0.993131	-0.00148868	-0.00436289
-0.00294458	-0.000538692	-0.00148868	1.01764	-0.00396274
-0.00432704	-0.00110805	-0.00436289	-0.00396274	0.993579

##the numbers don't match ibd.genome, because the two methods normalize 
##the estimates differently.

> head -n 1 grm.rel > grm.CH18526
> awk '$1=="CH18526"||$3=="CH18526"{print $0}' ibd.genome> ibd.CH18526
> ibd.CH18526
  CH18526  NA18526  CH18524  NA18524 UN    NA  1.0000  0.0000  0.0000  0.0000  -1  0.745392  0.1208  1.9158
  CH18526  NA18526  CH18529  NA18529 UN    NA  0.9804  0.0000  0.0196  0.0196  -1  0.747905  0.3380  1.9695
  CH18526  NA18526  CH18558  NA18558 UN    NA  0.9875  0.0125  0.0000  0.0063  -1  0.745061  0.4218  1.9855
  CH18526  NA18526  CH18532  NA18532 UN    NA  0.9581  0.0419  0.0000  0.0210  -1  0.744858  0.9476  2.1255
  CH18526  NA18526  CH18561  NA18561 UN    NA  0.9917  0.0075  0.0008  0.0045  -1  0.745123  0.6848  2.0360
  CH18526  NA18526  CH18562  NA18562 UN    NA  1.0000  0.0000  0.0000  0.0000  -1  0.744175  0.0368  1.8731
  CH18526  NA18526  CH18537  NA18537 UN    NA  1.0000  0.0000  0.0000  0.0000   0  0.745130  0.1025  1.9093
  CH18526  NA18526  CH18603  NA18603 UN    NA  1.0000  0.0000  0.0000  0.0000   0  0.744447  0.4412  1.9891
  CH18526  NA18526  CH18540  NA18540 UN    NA  1.0000  0.0000  0.0000  0.0000  -1  0.743700  0.0206  1.8556
  CH18526  NA18526  CH18605  NA18605 UN    NA  1.0000  0.0000  0.0000  0.0000  -1  0.742606  0.2942  1.9606


> Rscript --vanilla GMMAT.R

##creates 2 output files --  test.glmm.score.nocov, test.glmm.score.PC1PC2cov
## the first is the mixed model with no covariates, the second had PC1 and PC2 covariates
## let's look at the first few lines of each.
> head test.glmm.score.nocov
CHR	SNP	cM	POS	A1	A2	N	AF	SCORE	VAR	PVAL
1	rs3094315	0.792429	792429	G	A	88	0.880682-2.32331	3.31703	0.202078
1	rs4040617	0.819185	819185	G	A	89	0.88764-1.85531	3.21946	0.301131
1	rs4075116	1.04355	1043552	C	T	89	0.94382	1.220131.85665	0.370546
1	rs9442385	1.13726	1137258	T	G	88	0.602273	1.99713	10.0874	0.529477
1	rs11260562	1.20523	1205233	A	G	87	0.971264	0.510051	0.98758	0.607778
1	rs6685064	1.25121	1251215	C	T	89	0.589888	1.70415	11.4393	0.614361
1	rs3766180	1.56342	1563420	T	C	89	0.859551	-3.04227	4.56626	0.154534
1	rs6603791	1.58621	1586208	A	G	89	0.865169	-3.48533	4.47123	0.0992957
1	rs7519837	1.59607	1596068	C	T	88	0.869318	-3.2068	4.36158	0.124661

> head test.glmm.score.PC1PC2cov
CHR	SNP	cM	POS	A1	A2	N	AF	SCORE	VAR	PVAL
1	rs3094315	0.792429	792429	G	A	88	0.880682-1.46658	2.33439	0.337115
1	rs4040617	0.819185	819185	G	A	89	0.88764-1.22619	2.26829	0.415555
1	rs4075116	1.04355	1043552	C	T	89	0.94382	1.589931.28739	0.161133
1	rs9442385	1.13726	1137258	T	G	88	0.602273	1.2661	6.99343	0.632106
1	rs11260562	1.20523	1205233	A	G	87	0.971264	-0.0905551	0.792047	0.918955
1	rs6685064	1.25121	1251215	C	T	89	0.589888	-2.60625	6.86213	0.319777
1	rs3766180	1.56342	1563420	T	C	89	0.859551	-3.20116	3.29468	0.0777989
1	rs6603791	1.58621	1586208	A	G	89	0.865169	-3.48847	3.26659	0.0535899
1	rs7519837	1.59607	1596068	C	T	88	0.869318	-3.38653	3.19449	0.0581243

#these awk commands replace the "PVAL" column header with "P", and the "POS" with "BP" . 
#This will make it match the P-value and base pair column headers from PLINK output.
##this will be useful when we use the R program qqplot.R to do a qqplot and manhattan plot 
##of the p-values
awk 'NR==1{$4="BP";$11="P"};{print $0}' test.glmm.score.nocov>glmm.score.nocov.txt
awk 'NR==1{$4="BP";$11="P"};{print $0}' test.glmm.score.PC1PC2cov>glmm.score.PC1PC2cov.txt

> head glmm.score.nocov.txt
CHR SNP cM BP A1 A2 N AF SCORE VAR P
1	rs3094315	0.792429	792429	G	A	88	0.880682-2.32331	3.31703	0.202078
1	rs4040617	0.819185	819185	G	A	89	0.88764-1.85531	3.21946	0.301131
1	rs4075116	1.04355	1043552	C	T	89	0.94382	1.220131.85665	0.370546
1	rs9442385	1.13726	1137258	T	G	88	0.602273	1.99713	10.0874	0.529477
1	rs11260562	1.20523	1205233	A	G	87	0.971264	0.510051	0.98758	0.607778
1	rs6685064	1.25121	1251215	C	T	89	0.589888	1.70415	11.4393	0.614361
1	rs3766180	1.56342	1563420	T	C	89	0.859551	-3.04227	4.56626	0.154534
1	rs6603791	1.58621	1586208	A	G	89	0.865169	-3.48533	4.47123	0.0992957
1	rs7519837	1.59607	1596068	C	T	88	0.869318	-3.2068	4.36158	0.124661

> head glmm.score.PC1PC2cov.txt
CHR SNP cM BP A1 A2 N AF SCORE VAR P
1	rs3094315	0.792429	792429	G	A	88	0.880682-1.46658	2.33439	0.337115
1	rs4040617	0.819185	819185	G	A	89	0.88764-1.22619	2.26829	0.415555
1	rs4075116	1.04355	1043552	C	T	89	0.94382	1.589931.28739	0.161133
1	rs9442385	1.13726	1137258	T	G	88	0.602273	1.2661	6.99343	0.632106
1	rs11260562	1.20523	1205233	A	G	87	0.971264	-0.0905551	0.792047	0.918955
1	rs6685064	1.25121	1251215	C	T	89	0.589888	-2.60625	6.86213	0.319777
1	rs3766180	1.56342	1563420	T	C	89	0.859551	-3.20116	3.29468	0.0777989
1	rs6603791	1.58621	1586208	A	G	89	0.865169	-3.48847	3.26659	0.0535899
1	rs7519837	1.59607	1596068	C	T	88	0.869318	-3.38653	3.19449	0.0581243

##Let's compare the effect estimates for one of the most associated variants:  rs11204005 
> awk '(NR==1||($2=="rs11204005")){print $0}' logistadjpop-hidecov.assoc.logistic
 CHR         SNP         BP   A1       TEST    NMISS       BETA       SE      L95      U95         STAT            P 
   8  rs11204005   12895576    A        ADD       89     -2.708   0.6254   -3.934   -1.482        -4.33    1.489e-05
>  awk '(NR==1||($2=="rs11204005")){print $0}' logistadjPC12.assoc.logistic 
 CHR         SNP         BP   A1       TEST    NMISS       BETA       SE      L95      U95         STAT            P 
   8  rs11204005   12895576    A        ADD       89     -2.648   0.6389     -3.9   -1.396       -4.144    3.407e-05
   8  rs11204005   12895576    A        PC1       89      21.32    4.597    12.31    30.33        4.638    3.523e-06
   8  rs11204005   12895576    A        PC2       89     -2.886    3.471   -9.688    3.917      -0.8314       0.4058
> awk '(NR==1||($2=="rs11204005")){print $0}' glmm.score.nocov.txt
CHR SNP cM BP A1 A2 N AF SCORE VAR P
8	rs11204005	12.8956	12895576	A	G	89	0.52247212.6384	9.54287	4.29135e-05
> awk '(NR==1||($2=="rs11204005")){print $0}' glmm.score.PC1PC2cov.txt
CHR SNP cM BP A1 A2 N AF SCORE VAR P
8	rs11204005	12.8956	12895576	A	G	89	0.52247212.4996	6.68642	1.33882e-06


#QQplots of the GWAS for the various models we've tried

Rscript --vanilla qqplot.R logistnoadj.assoc.logistic crude ADD 

Rscript --vanilla qqplot.R   logistadjpop-hidecov.assoc.logistic popadj ADD 

Rscript --vanilla qqplot.R  logistadjPC12-hidecov.assoc.logistic PC12adj ADD 

Rscript --vanilla qqplot.R  glmm.score.nocov.txt "glmm nocovariates" 

Rscript --vanilla qqplot.R glmm.score.PC1PC2cov.txt "glmm PC1PC2covariates" 

# creates the following files: 	qq.crude.jpeg, 	qq.glmm nocovariates.jpeg, 	qq.glmm PC1PC2covariates.jpeg, 	qq.popadj.jpeg

### alternative tool for making prettier qq plots (using Umich qqplot code):

Rscript  --vanilla qq_umich_gc.R  logistnoadj.assoc.logistic noadj ADD 
Rscript --vanilla qq_umich_gc.R  logistadjpop-hidecov.assoc.logistic popadj ADD 
Rscript --vanilla qq_umich_gc.R  logistadjPC12-hidecov.assoc.logistic pcadj ADD 
Rscript --vanilla qq_umich_gc.R  glmm.score.nocov.txt glmmnocov  
Rscript --vanilla qq_umich_gc.R  glmm.score.PC1PC2cov.txt glmmPC1PC2cov  


##Manhattan plots  for all of the GWAS models we tried

Rscript --vanilla gwaplot.R  logistnoadj.assoc.logistic "unadjusted analysis" unadj_manhattan
Rscript --vanilla gwaplot.R  logistadjpop-hidecov.assoc.logistic "Pop adj analysis" popadj_manhattan
Rscript --vanilla gwaplot.R logistadjPC12-hidecov.assoc.logistic "PC adj analysis" pc1pc2adj_manhattan
Rscript --vanilla gwaplot.R  glmm.score.nocov.txt "GLMM no cov analysis" GLMM_nocov_manhattan 
Rscript --vanilla gwaplot.R  glmm.score.PC1PC2cov.txt "GLMM PC1 and PC2 adjusted analysis" GLMM_PC12_manhattan 

