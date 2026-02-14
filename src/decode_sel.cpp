#include <Rcpp.h>
#include <fstream>
#include <string>

using namespace Rcpp;

// Helper function to count samples from .fam file
int count_samples_from_fam(std::string fam_file) {
  std::ifstream file(fam_file.c_str());
  if (!file.is_open()) {
    stop("Cannot open FAM file: " + fam_file);
  }
  
  int count = 0;
  std::string line;
  while (std::getline(file, line)) {
    if (!line.empty()) count++;
  }
  
  file.close();
  return count;
}

//' read subset for use in delayed matrix approach
//' @param prefix character path to bed resources and file prefix 
//' @param sample_indices integer vector of snp indices
//' @param sample_indices integer vector of sample indices
//' @param n_total_samples optional
//' @export
// [[Rcpp::export]]
NumericMatrix read_bed_subset(std::string prefix,
                              IntegerVector snp_indices,    // 1-based
                              IntegerVector sample_indices, // 1-based
                              Nullable<int> n_total_samples = R_NilValue) {
  
  // Construct file paths
  std::string bed_file = prefix + ".bed";
  std::string fam_file = prefix + ".fam";
  std::string bim_file = prefix + ".bim";  // For future use
  
  // Determine total number of samples
  int n_samples_total;
  if (n_total_samples.isNotNull()) {
    n_samples_total = as<int>(n_total_samples);
  } else {
    n_samples_total = count_samples_from_fam(fam_file);
  }
  
  int n_snps = snp_indices.size();
  int n_samples = sample_indices.size();
  NumericMatrix result(n_snps, n_samples);
  
  // Lookup table for BED encoding
  double lookup[4] = {2.0, NA_REAL, 1.0, 0.0};
  
  // Open BED file
  std::ifstream file(bed_file.c_str(), std::ios::binary);
  if (!file.is_open()) {
    stop("Cannot open BED file: " + bed_file);
  }
  
  // Verify BED file header (0x6c 0x1b 0x01 for SNP-major mode)
  unsigned char header[3];
  file.read(reinterpret_cast<char*>(header), 3);
  if (header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
    stop("Invalid BED file header - expected SNP-major mode");
  }
  
  // Bytes per SNP (each SNP stored in ceil(n_total_samples/4) bytes)
  int bytes_per_snp = (n_samples_total + 3) / 4;
  
  // Buffer to hold one SNP's worth of data
  std::vector<unsigned char> buffer(bytes_per_snp);
  
  for (int i = 0; i < n_snps; i++) {
    // Calculate file position for this SNP (convert from 1-based to 0-based)
    int snp_idx = snp_indices[i] - 1;
    std::streampos pos = 3 + static_cast<std::streampos>(snp_idx) * bytes_per_snp;
    
    // Seek and read
    file.seekg(pos);
    file.read(reinterpret_cast<char*>(buffer.data()), bytes_per_snp);
    
    if (!file) {
      stop("Error reading SNP at index " + std::to_string(snp_indices[i]));
    }
    
    // Extract genotypes for requested samples
    for (int j = 0; j < n_samples; j++) {
      int sample_idx = sample_indices[j] - 1; // Convert to 0-based
      
      // Which byte and which 2-bit position within that byte?
      int byte_pos = sample_idx / 4;
      int bit_pos = (sample_idx % 4) * 2;
      
      // Extract 2-bit genotype
      unsigned char genotype = (buffer[byte_pos] >> bit_pos) & 0x03;
      result(i, j) = lookup[genotype];
    }
  }
  
  file.close();
  return result;
}

