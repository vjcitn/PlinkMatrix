# Constructor function for seed for plink bed format

Constructor function for seed for plink bed format

## Usage

``` r
PlinkSeed(filepath)
```

## Arguments

- filepath:

  character string without suffixes

## Value

PlinkSeed instance

## Examples

``` r
PlinkSeed
#> function (filepath) 
#> {
#>     filepath <- sub("\\.bed$", "", filepath)
#>     if (!file.exists(paste0(filepath, ".bed"))) {
#>         stop("BED file not found: ", filepath, ".bed")
#>     }
#>     if (!file.exists(paste0(filepath, ".fam"))) {
#>         stop("FAM file not found: ", filepath, ".fam")
#>     }
#>     if (!file.exists(paste0(filepath, ".bim"))) {
#>         stop("BIM file not found: ", filepath, ".bim")
#>     }
#>     fam <- read.table(paste0(filepath, ".fam"), header = FALSE, 
#>         col.names = c("FID", "IID", "PID", "MID", "Sex", "Phenotype"), 
#>         colClasses = "character")
#>     bim <- read.table(paste0(filepath, ".bim"), header = FALSE, 
#>         col.names = c("CHR", "ID", "cM", "POS", "A1", "A2"), 
#>         colClasses = "character")
#>     n_samples <- nrow(fam)
#>     n_variants <- nrow(bim)
#>     sample_ids <- paste(fam$FID, fam$IID, sep = "_")
#>     variant_ids <- bim$ID
#>     new("PlinkSeed", filepath = filepath, dim = c(n_variants, 
#>         n_samples), dimnames = list(variant_ids, sample_ids), 
#>         fam = fam, bim = bim)
#> }
#> <bytecode: 0x7a574b6238>
#> <environment: namespace:PlinkMatrix>
```
