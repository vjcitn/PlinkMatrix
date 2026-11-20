# read subset for use in delayed matrix approach

read subset for use in delayed matrix approach

## Usage

``` r
read_bed_subset(prefix, snp_indices, sample_indices, n_total_samples = NULL)
```

## Arguments

- prefix:

  character path to bed resources and file prefix

- snp_indices:

  integer vector of snp indices

- sample_indices:

  integer vector of sample indices

- n_total_samples:

  optional

## Value

matrix

## Examples

``` r
read_bed_subset
#> function (prefix, snp_indices, sample_indices, n_total_samples = NULL) 
#> {
#>     .Call(`_PlinkMatrix_read_bed_subset`, prefix, snp_indices, 
#>         sample_indices, n_total_samples)
#> }
#> <bytecode: 0x7a5b604cb8>
#> <environment: namespace:PlinkMatrix>
```
