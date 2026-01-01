#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector decode_bed_fast(RawVector raw_bytes, int n_samples) {
  int n_bytes = raw_bytes.size();
  NumericVector result(n_samples);

  // Lookup table - COUNTS A1 ALLELE (matches bed_reader)
  // 00->2 (A1A1), 01->NA, 10->1 (A1A2), 11->0 (A2A2)
  double lookup[4] = {2.0, NA_REAL, 1.0, 0.0};

  int geno_idx = 0;
  for (int i = 0; i < n_bytes && geno_idx < n_samples; i++) {
    unsigned char byte = raw_bytes[i];

    for (int j = 0; j < 4 && geno_idx < n_samples; j++) {
      result[geno_idx++] = lookup[byte & 0x03];
      byte >>= 2;
    }
  }

  return result;
}

