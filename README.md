# PlinkMatrix: DelayedArray interface for plink genotypes

Given a collection of .bed, .bim and .fam files produced
for use with Plink, this package facilitates simple
array-like operations for retrieving genotype data.

As an example, the GEUVADIS study genotypes distributed
with tensorQTL were transformed to `geuv445.{bed,bim,fam}`.
In a folder containing these files, we would proceed
as follows to interrogate the genotypes.

```
> library(PlinkMatrix)
> gdemo = example_PlinkMatrix()
> colnames(gdemo) = gsub("0_", "", colnames(gdemo))
> gdemo
<367759 x 445> DelayedMatrix object of type "double":
                        HG00096 HG00097 HG00099 ... NA20826 NA20828
    chr18_10644_C_G_b38       0       0       0   .       0       0
    chr18_10847_C_A_b38       0       0       0   .       0       0
    chr18_11275_G_A_b38       0       0       0   .       0       0
    chr18_11358_G_A_b38       0       0       0   .       0       0
    chr18_11445_G_A_b38       0       0       0   .       0       0
                    ...       .       .       .   .       .       .
chr18_80259028_AG_A_b38       1       2       2   .       0       1
 chr18_80259147_G_C_b38       0       0       0   .       0       0
 chr18_80259181_A_G_b38       0       0       0   .       0       0
   chr18_8025919C_G_b38       1       0       0   .       1       0
 chr18_80259245_C_A_b38       0       0       0   .       0       0
```

The associated data can be obtained in this package with `example_PlinkArray()`.

See "Get started" tab above for more information.



