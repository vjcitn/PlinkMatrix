#' Method: extract_array, internal 
#' @param x seed instance
#' @param index list of suitable values for extracting elements
setMethod("extract_array", "PlinkSeed",
  function(x, index) {
    # index is a list of length 2: list(row_indices, col_indices)
    # If NULL, means select all
    
    n_samples <- x@dim[1]
    n_variants <- x@dim[2]
    
    # Handle NULL indices (select all)
    if (is.null(index[[1]])) {
      row_idx <- seq_len(n_samples)
    } else {
      row_idx <- index[[1]]
    }
    
    if (is.null(index[[2]])) {
      col_idx <- seq_len(n_variants)
    } else {
      col_idx <- index[[2]]
    }
    
    # Allocate result matrix
    result <- matrix(NA_real_, 
                     nrow = length(row_idx), 
                     ncol = length(col_idx))
    
    # Open BED file
    bed_file <- paste0(x@filepath, ".bed")
    con <- file(bed_file, "rb")
    on.exit(close(con))
    
    # Check magic number (first 3 bytes should be 0x6C, 0x1B, 0x01)
    magic <- readBin(con, "raw", n = 3)
    if (!identical(magic, as.raw(c(0x6C, 0x1B, 0x01)))) {
      stop("Invalid BED file format")
    }
    
    # Calculate bytes per variant (SNP-major mode)
    bytes_per_variant <- ceiling(n_samples / 4)
    
    # Read each requested variant
    for (j in seq_along(col_idx)) {
      variant_idx <- col_idx[j]
      
      # Seek to the position of this variant
      # 3 bytes header + (variant_idx - 1) * bytes_per_variant
      seek(con, where = 3 + (variant_idx - 1) * bytes_per_variant, 
           origin = "start")
      
      # Read the bytes for this variant
      raw_data <- readBin(con, "raw", n = bytes_per_variant)
      
      # Decode genotypes
      genotypes <- decode_bed_fast(raw_data, n_samples)
      
      # Extract requested samples
      result[, j] <- genotypes[row_idx]
    }
    
    result
  }
)
