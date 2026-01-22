#' Helper function to decode BED format
#' @param raw_bytes byte stream from file
#' @param n_samples number of samples
decode_bed_genotypes <- function(raw_bytes, n_samples) {
  # Lookup table for decoding (NOW COUNTING A1 ALLELE)
  # 00 = A1A1 -> 2 copies of A1
  # 01 = missing -> NA
  # 10 = A1A2 -> 1 copy of A1
  # 11 = A2A2 -> 0 copies of A1
  lookup <- matrix(c(
    2, NA, 1, 0,  # First genotype in byte
    2, NA, 1, 0,  # Second genotype
    2, NA, 1, 0,  # Third genotype
    2, NA, 1, 0   # Fourth genotype
  ), nrow = 4, byrow = TRUE)
  
  n_bytes <- length(raw_bytes)
  n_genotypes_total <- n_bytes * 4
  
  result <- numeric(n_genotypes_total)
  
  for (i in seq_along(raw_bytes)) {
    byte_val <- as.integer(raw_bytes[i])
    
    idx_base <- (i - 1) * 4
    result[idx_base + 1] <- lookup[1, (byte_val %% 4) + 1]
    result[idx_base + 2] <- lookup[2, ((byte_val %/% 4) %% 4) + 1]
    result[idx_base + 3] <- lookup[3, ((byte_val %/% 16) %% 4) + 1]
    result[idx_base + 4] <- lookup[4, (byte_val %/% 64) + 1]
  }
  
  result[1:n_samples]
}
