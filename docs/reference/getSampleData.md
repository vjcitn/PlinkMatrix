<div id="main" class="col-md-9" role="main">

# Get sample metadata

<div class="ref-description section level2">

Get sample metadata

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
getSampleData(x)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    DelayedArray instance

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
tst = example_PlinkMatrix()
head(getSampleData(tst))
#>   FID     IID PID MID Sex Phenotype
#> 1   0 HG00096   0   0   1        -9
#> 2   0 HG00097   0   0   2        -9
#> 3   0 HG00099   0   0   2        -9
#> 4   0 HG00100   0   0   2        -9
#> 5   0 HG00101   0   0   1        -9
#> 6   0 HG00102   0   0   2        -9
```

</div>

</div>

</div>
