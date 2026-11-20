# Helper function to decode BED format

Helper function to decode BED format

## Usage

``` r
decode_bed_genotypes(raw_bytes, n_samples)
```

## Arguments

- raw_bytes:

  byte stream from file

- n_samples:

  number of samples

## Value

numeric vector

## Examples

``` r
PlinkMatrix:::decode_bed_genotypes(raw(4), 4L)
#> [1] 2 2 2 2
```
