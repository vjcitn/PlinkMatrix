# sample characteristics of 445 GEUVADIS samples

sample characteristics of 445 GEUVADIS samples

## Usage

``` r
data("g445samples", package="PlinkMatrix")
```

## Format

data.frame

## Note

Example data are those provided with tensorqtl, see
<https://github.com/broadinstitute/tensorqtl/tree/0c4db65a0cdc47f3b824ae530b89d270ef5e0096/example/data>.

## Examples

``` r
data("g445samples", package="PlinkMatrix")
g445samples[seq_len(4), seq_len(4)]
#>            Sex Population.code Population.name Superpopulation.code
#> HG00096   male             GBR         British                  EUR
#> HG00097 female             GBR         British                  EUR
#> HG00099 female             GBR         British                  EUR
#> HG00100 female             GBR         British                  EUR
```
