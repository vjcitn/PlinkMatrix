
#' Define the PlinkSeed class
#' @importFrom utils read.table unzip
#' @import methods
#' @import DelayedArray
#' @export
setClass("PlinkSeed",
  slots = c(
    filepath = "character",      # Path to .bed file (without extension)
    dim = "integer",            # Dimensions: c(n_samples, n_variants)
    dimnames = "list",          # List of sample IDs and variant IDs
    fam = "data.frame",         # .fam file content
    bim = "data.frame"          # .bim file content
  ),
  contains = "Array"
)

#' Constructor function for seed for plink bed format
#' @param filepath character string without suffixes
#' @export
PlinkSeed <- function(filepath) {
  # Remove .bed extension if present
  filepath <- sub("\\.bed$", "", filepath)  # consider filepath_as_absolute from tools
  
  # Check files exist
  if (!file.exists(paste0(filepath, ".bed"))) {
    stop("BED file not found: ", filepath, ".bed")
  }
  if (!file.exists(paste0(filepath, ".fam"))) {
    stop("FAM file not found: ", filepath, ".fam")
  }
  if (!file.exists(paste0(filepath, ".bim"))) {
    stop("BIM file not found: ", filepath, ".bim")
  }
  
  # Read .fam file (sample information)
  fam <- read.table(paste0(filepath, ".fam"), 
                    header = FALSE,
                    col.names = c("FID", "IID", "PID", "MID", "Sex", "Phenotype"),
                    colClasses = "character")
  
  # Read .bim file (variant information)
  bim <- read.table(paste0(filepath, ".bim"),
                    header = FALSE,
                    col.names = c("CHR", "ID", "cM", "POS", "A1", "A2"),
                    colClasses = "character")
  
  # Dimensions: samples x variants
  n_samples <- nrow(fam)
  n_variants <- nrow(bim)
  
  # Create dimnames
  sample_ids <- paste(fam$FID, fam$IID, sep = "_")
  variant_ids <- bim$ID
  
  new("PlinkSeed",
      filepath = filepath,
      dim = c(n_samples, n_variants),
      dimnames = list(sample_ids, variant_ids),
      fam = fam,
      bim = bim)
}

#' Method: dim for delayed plink
#' @param x PlinkSeed instance
#' @export
setMethod("dim", "PlinkSeed", function(x) x@dim)

#' Method: dimnames for delayed plink
#' @param x PlinkSeed instance
#' @export
setMethod("dimnames", "PlinkSeed", function(x) x@dimnames)


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

# Optional: Define chunkdim for better performance
setMethod("chunkdim", "PlinkSeed",
  function(x) {
    as.integer(c(x@dim[1], min(50000, x@dim[2])))
  }
)

#' Constructor for DelayedArray
#' @param filepath path to plink bed, bim, fam resources without suffixes
#' @export
PlinkArray <- function(filepath) {
  seed <- PlinkSeed(filepath)
  DelayedArray(seed)
}

#' present seed concisely
#' @param object instance of PlinkSeed
#' @export
setMethod("show", "PlinkSeed",
  function(object) {
    cat("PlinkSeed object\n")
    cat("Filepath:", object@filepath, "\n")
    cat("Dimensions:", object@dim[1], "samples x", 
        object@dim[2], "variants\n")
  }
)

