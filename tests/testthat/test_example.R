#PC002566:1KGPLINK vincentcarey$ plink2 --bfile /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//Rtmpy6XqlY/geuv445 --snp chr18_10644_C_G_b38 --geno-counts
#PLINK v2.0.0-a.6.29 M1 (28 Nov 2025)                cog-genomics.org/plink/2.0/
#(C) 2005-2025 Shaun Purcell, Christopher Chang    GNU General Public License v3
#Logging to plink2.log.
#Options in effect:
#  --bfile /var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//Rtmpy6XqlY/geuv445
#  --geno-counts
#  --snp chr18_10644_C_G_b38
#
#Start time: Fri Feb 13 16:20:16 2026
#32768 MiB RAM detected; reserving 16384 MiB for main workspace.
#Using up to 10 threads (change this with --threads).
#445 samples (236 females, 209 males; 445 founders) loaded from
#/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//Rtmpy6XqlY/geuv445.fam.
#367759 variants loaded from
#/var/folders/yw/gfhgh7k565v9w83x_k764wbc0000gp/T//Rtmpy6XqlY/geuv445.bim.
#Note: No phenotype data present.
#--snp: 1 variant remaining.
#Calculating allele frequencies... done.
#--geno-counts: Genotype counts written to plink2.gcount .
#End time: Fri Feb 13 16:20:16 2026
#PC002566:1KGPLINK vincentcarey$ cat plink2.gcount
##CHROM	ID	REF	ALT	PROVISIONAL_REF?	HOM_REF_CT	HET_REF_ALT_CTS	TWO_ALT_GENO_CTS	HAP_REF_CT	HAP_ALT_CTS	MISSING_CT
#18	chr18_10644_C_G_b38	C	G	Y	430	15	0	0	0	0
#PC002566:1KGPLINK vincentcarey$ 


test_that("example is as expected", {
 ex = example_PlinkMatrix()
 expect_true(all(dim(ex) == c(367759L, 445L)))
 tt = table(as.numeric(ex["chr18_10644_C_G_b38",]))
 expect_true(all(as.integer(tt)==c(430L, 15L)))
})

#test_that("subset RangedPlink functions", {
# example(RangedPlink) # creates lit and myGR
# expect_true(max(IRanges::start((slot(IRanges::subsetByOverlaps(lit, myGR), "ranges")))) < 500000)
#})

test_that("subsetting of example with as_RSE is successful", {
 ex = example_PlinkMatrix(as_RSE=TRUE)
 GenomeInfoDb::seqlevelsStyle(rowRanges(ex)) = "Ensembl"
 simp = GenomicRanges::GRanges("18:1-200000")
 rs = subsetByOverlaps(ex, simp)
 expect_true(nrow(rs) == 835L)
})
