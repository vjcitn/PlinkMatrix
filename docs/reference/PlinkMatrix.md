# Constructor for DelayedArray

Constructor for DelayedArray

## Usage

``` r
PlinkMatrix(filepath)
```

## Arguments

- filepath:

  path to plink bed, bim, fam resources without suffixes

## Value

An instance of PlinkMatrix

## Examples

``` r
PlinkMatrix
#> function (filepath) 
#> {
#>     filepath <- normalizePath(filepath, mustWork = FALSE)
#>     seed <- PlinkSeed(filepath)
#>     DelayedArray(seed)
#> }
#> <bytecode: 0x7a4cac9818>
#> <environment: namespace:PlinkMatrix>
```
