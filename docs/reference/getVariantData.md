# Get variant metadata

Get variant metadata

## Usage

``` r
getVariantData(x)
```

## Arguments

- x:

  DelayedArray instance

## Value

data.frame

## Examples

``` r
tst <- example_PlinkMatrix()
dim(getVariantData(tst))
#> [1] 367759      6
```
