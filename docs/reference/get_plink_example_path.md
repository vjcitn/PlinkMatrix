# operate with BiocFileCache to retrieve zip file of plink example data

operate with BiocFileCache to retrieve zip file of plink example data

## Usage

``` r
get_plink_example_path(ca = BiocFileCache::BiocFileCache())
```

## Arguments

- ca:

  BiocFileCache instance

## Value

character string with path to example in cache

## Examples

``` r
get_plink_example_path()
#> [1] "/Users/vincentcarey/Library/Caches/org.R-project.R/R/BiocFileCache/e56c6de7e795_geuv445_plink.zip"
```
