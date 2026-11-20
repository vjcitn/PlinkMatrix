# Method: dim for delayed plink

Method: dim for delayed plink

## Usage

``` r
# S4 method for class 'PlinkSeed'
dim(x)
```

## Arguments

- x:

  PlinkSeed instance

## Value

2-vector

## Examples

``` r
methods(class="PlinkSeed")
#>  [1] DelayedArray                     [                               
#>  [3] [<-                              [[                              
#>  [5] as.array                         as.character                    
#>  [7] as.complex                       as.data.frame                   
#>  [9] as.integer                       as.logical                      
#> [11] as.matrix                        as.numeric                      
#> [13] as.raw                           as.vector                       
#> [15] bindROWS                         cbind                           
#> [17] chunkdim                         coerce                          
#> [19] dim                              dim<-                           
#> [21] dimnames                         drop                            
#> [23] extract_array                    isEmpty                         
#> [25] kronecker                        length                          
#> [27] rbind                            show                            
#> [29] showAsCell                       subassign_Array_by_Lindex       
#> [31] subassign_Array_by_Mindex        subassign_Array_by_Nindex       
#> [33] subassign_Array_by_logical_array subset_Array_by_Lindex          
#> [35] subset_Array_by_Mindex           subset_Array_by_Nindex          
#> [37] subset_Array_by_logical_array    t                               
#> see '?methods' for accessing help and source code
```
