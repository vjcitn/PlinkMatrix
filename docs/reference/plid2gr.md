# produce GRanges from variant notation for plink example from geuvadis

produce GRanges from variant notation for plink example from geuvadis

## Usage

``` r
plid2gr(x, sepused = "_")
```

## Arguments

- x:

  character vector of variant names

- sepused:

  single character, defaults to "\_"

## Value

GRanges instance

## Examples

``` r
plid2gr("chr18_80259028_AG_A_b38")
#> GRanges object with 1 range and 0 metadata columns:
#>                           seqnames    ranges strand
#>                              <Rle> <IRanges>  <Rle>
#>   chr18_80259028_AG_A_b38    chr18  80259028      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
