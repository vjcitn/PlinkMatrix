library(PlinkMatrix)
zz = example_PlinkMatrix()
zzs = slot(zz, "seed")
usepath = slot(zzs, "filepath")

test_that("read_bed_subset returns correct dimensions", {
  # Test various dimension combinations
  result <- read_bed_subset(usepath, 
                            snp_indices = 1:10, 
                            sample_indices = 1:20)
  expect_equal(dim(result), c(10, 20))
  
  # Single SNP, multiple samples
  result <- read_bed_subset(usepath, 
                            snp_indices = 1, 
                            sample_indices = 1:445)
  expect_equal(dim(result), c(1, 445))
  
  # Multiple SNPs, single sample
  result <- read_bed_subset(usepath, 
                            snp_indices = 1:100, 
                            sample_indices = 1)
  expect_equal(dim(result), c(100, 1))
})

test_that("read_bed_subset handles edge cases", {
  # First SNP, first sample
  result <- read_bed_subset(usepath, 
                            snp_indices = 1, 
                            sample_indices = 1)
  expect_equal(dim(result), c(1, 1))
  expect_true(is.numeric(result[1, 1]))
  
  # Last SNP, last sample
  result <- read_bed_subset(usepath, 
                            snp_indices = 367759, 
                            sample_indices = 445)
  expect_equal(dim(result), c(1, 1))
  expect_true(is.numeric(result[1, 1]))
  
  # First and last SNPs
  result <- read_bed_subset(usepath, 
                            snp_indices = c(1, 367759), 
                            sample_indices = 1:10)
  expect_equal(dim(result), c(2, 10))
})

test_that("read_bed_subset produces valid genotype values", {
  result <- read_bed_subset(usepath, 
                            snp_indices = 1:1000, 
                            sample_indices = 1:100)
  
  # All values should be 0, 1, 2, or NA
  valid_values <- c(0, 1, 2, NA)
  expect_true(all(result %in% valid_values | is.na(result)))
  
  # Check that we have at least some variation
  unique_vals <- unique(as.vector(result))
  expect_true(length(unique_vals) > 1)
})

test_that("read_bed_subset is consistent across different access patterns", {
  # Read a block of data
  full_block <- read_bed_subset(usepath, 
                                snp_indices = 100:110, 
                                sample_indices = 50:60)
  
  # Read same data SNP by SNP
  snp_by_snp <- lapply(100:110, function(s) {
    read_bed_subset(usepath, 
                    snp_indices = s, 
                    sample_indices = 50:60)
  })
  snp_by_snp <- do.call(rbind, snp_by_snp)
  
  expect_equal(full_block, snp_by_snp)
  
  # Read same data sample by sample
  sample_by_sample <- lapply(50:60, function(samp) {
    read_bed_subset(usepath, 
                    snp_indices = 100:110, 
                    sample_indices = samp)
  })
  sample_by_sample <- do.call(cbind, sample_by_sample)
  
  expect_equal(full_block, sample_by_sample)
})

test_that("read_bed_subset handles non-contiguous indices", {
  # Random SNPs and samples
  set.seed(123)
  random_snps <- sort(sample(1:367759, 50))
  random_samples <- sort(sample(1:445, 30))
  
  result <- read_bed_subset(usepath, 
                            snp_indices = random_snps, 
                            sample_indices = random_samples)
  
  expect_equal(dim(result), c(50, 30))
  expect_true(all(result %in% c(0, 1, 2, NA) | is.na(result)))
})

test_that("read_bed_subset respects index order", {
  # Forward order
  forward <- read_bed_subset(usepath, 
                             snp_indices = c(10, 20, 30), 
                             sample_indices = c(5, 10, 15))
  
  # Reverse order
  reverse <- read_bed_subset(usepath, 
                             snp_indices = c(30, 20, 10), 
                             sample_indices = c(15, 10, 5))
  
  # Should NOT be equal (order matters)
  expect_false(identical(forward, reverse))
  
  # But reverse of reverse should match forward
  expect_equal(forward, reverse[3:1, 3:1])
})

test_that("read_bed_subset auto-detects sample count correctly", {
  # With auto-detection
  result1 <- read_bed_subset(usepath, 
                             snp_indices = 1:10, 
                             sample_indices = 1:10)
  
  # With explicit sample count
  result2 <- read_bed_subset(usepath, 
                             snp_indices = 1:10, 
                             sample_indices = 1:10,
                             n_total_samples = 445)
  
  expect_equal(result1, result2)
})

test_that("read_bed_subset handles boundary samples correctly", {
  # Samples at byte boundaries (every 4th sample)
  # This tests proper bit extraction
  byte_boundary_samples <- seq(1, 445, by = 4)
  
  result <- read_bed_subset(usepath, 
                            snp_indices = 1:100, 
                            sample_indices = byte_boundary_samples)
  
  expect_equal(ncol(result), length(byte_boundary_samples))
  expect_true(all(result %in% c(0, 1, 2, NA) | is.na(result)))
})

test_that("read_bed_subset gives informative errors", {
  # Invalid SNP index
  expect_error(
    read_bed_subset(usepath, 
                    snp_indices = 367760,  # Beyond max
                    sample_indices = 1:10),
    "Error reading SNP"
  )
  
  # Invalid sample index (beyond n_total_samples)
  expect_error(
    read_bed_subset(usepath, 
                    snp_indices = 1:10, 
                    sample_indices = 446,
                    n_total_samples = 445),
    NA  # May not error, but shouldn't crash
  )
  
  # Nonexistent file
  expect_error(
    read_bed_subset("nonexistent_file", 
                    snp_indices = 1, 
                    sample_indices = 1),
    "Cannot open"
  )
})

test_that("read_bed_subset handles duplicate indices", {
  # Duplicate SNPs
  result <- read_bed_subset(usepath, 
                            snp_indices = c(10, 10, 20), 
                            sample_indices = 1:5)
  
  expect_equal(nrow(result), 3)
  expect_equal(result[1, ], result[2, ])  # First two rows should be identical
  
  # Duplicate samples
  result <- read_bed_subset(usepath, 
                            snp_indices = 1:5, 
                            sample_indices = c(10, 10, 20))
  
  expect_equal(ncol(result), 3)
  expect_equal(result[, 1], result[, 2])  # First two columns should be identical
})

