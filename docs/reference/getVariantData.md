<div id="main" class="col-md-9" role="main">

# Get variant metadata

<div class="ref-description section level2">

Get variant metadata

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
getVariantData(x)
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
dim(getVariantData(tst))
#> [1] 367759      6
```

</div>

</div>

</div>
