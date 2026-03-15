<div id="main" class="col-md-9" role="main">

# produce GRanges from variant notation for plink example from geuvadis

<div class="ref-description section level2">

produce GRanges from variant notation for plink example from geuvadis

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
plid2gr(x, sepused = "_")
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    character vector of variant names

-   sepused:

    single character, defaults to "\_"

</div>

<div class="section level2">

## Value

GRanges instance

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
plid2gr("chr18_80259028_AG_A_b38")
#> GRanges object with 1 range and 0 metadata columns:
#>                           seqnames    ranges strand
#>                              <Rle> <IRanges>  <Rle>
#>   chr18_80259028_AG_A_b38    chr18  80259028      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

</div>

</div>

</div>
