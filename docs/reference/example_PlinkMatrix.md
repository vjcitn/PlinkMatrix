<div id="main" class="col-md-9" role="main">

# produce PlinkMatrix from example data

<div class="ref-description section level2">

produce PlinkMatrix from example data

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
example_PlinkMatrix(folder = tempdir(), as_RSE = FALSE)
```

</div>

</div>

<div class="section level2">

## Arguments

-   folder:

    a path where unzipped example data will be managed

-   as\_RSE:

    logical(1) if TRUE (default is FALSE) a RangedSummarizedExperiment
    is returned with rowRanges calculated from SNP ids

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
example_PlinkMatrix()
#> <367759 x 445> PlinkMatrix object of type "double":
#>                         0_HG00096 0_HG00097 0_HG00099 ... 0_NA20826 0_NA20828
#>     chr18_10644_C_G_b38         0         0         0   .         0         0
#>     chr18_10847_C_A_b38         0         0         0   .         0         0
#>     chr18_11275_G_A_b38         0         0         0   .         0         0
#>     chr18_11358_G_A_b38         0         0         0   .         0         0
#>     chr18_11445_G_A_b38         0         0         0   .         0         0
#>                     ...         .         .         .   .         .         .
#> chr18_80259028_AG_A_b38         1         2         2   .         0         1
#>  chr18_80259147_G_C_b38         0         0         0   .         0         0
#>  chr18_80259181_A_G_b38         0         0         0   .         0         0
#>  chr18_80259190_C_G_b38         1         0         0   .         1         0
#>  chr18_80259245_C_A_b38         0         0         0   .         0         0
```

</div>

</div>

</div>
