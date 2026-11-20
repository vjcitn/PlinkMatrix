# Get sample metadata

Get sample metadata

## Usage

``` r
getSampleData(x)
```

## Arguments

- x:

  DelayedArray instance

## Value

data.frame

## Examples

``` r
tst <- example_PlinkMatrix()
head(getSampleData(tst))
#>   FID     IID PID MID Sex Phenotype
#> 1   0 HG00096   0   0   1        -9
#> 2   0 HG00097   0   0   2        -9
#> 3   0 HG00099   0   0   2        -9
#> 4   0 HG00100   0   0   2        -9
#> 5   0 HG00101   0   0   1        -9
#> 6   0 HG00102   0   0   2        -9
```
