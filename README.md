# PlinkArray: DelayedArray interface for plink genotypes

Given a collection of .bed, .bim and .fam files produced
for use with Plink, this package facilitates simple
array-like operations for retrieving genotype data.

As an example, the GEUVADIS study genotypes distributed
with tensorQTL were transformed to `geuv445.{bed,bim,fam}`.
In a folder containing these files, we would proceed
as follows to interrogate the genotypes.

```
> library(PlinkArray)
> gdemo = PlinkArray("geuv445")
> gdemo
<445 x 367759> DelayedMatrix object of type "double":
             chr18_10644_C_G_b38 ... chr18_80259245_C_A_b38
0_HG00096                      0   .                      0
0_HG00097                      0   .                      0
0_HG00099                      0   .                      0
0_HG00100                      0   .                      0
0_HG00101                      0   .                      0
      ...                      .   .                      .
0_NA20814                      0   .                      0
0_NA20815                      0   .                      0
0_NA20819                      0   .                      0
0_NA20826                      0   .                      0
0_NA20828                      0   .                      0
> sum(gdemo[4:10, 50:100])
[1] 114
> summary(colSums(gdemo))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    9.0    20.0    77.0   190.1   290.0   881.0 
```

The associated data can be obtained in this package with `example_PlinkArray()`.

See "Get started" tab above for more information.



