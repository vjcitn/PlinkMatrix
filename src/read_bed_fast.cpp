#include <Rcpp.h>

#include <fstream>

using namespace Rcpp;

static void decode_variant(const char *buf, IntegerVector row_idx, int *out)
{
        const int lookup[4] = {2, NA_INTEGER, 1, 0};
        for (int i = 0; i < row_idx.size(); i++) {
                int sample_idx = row_idx[i] - 1;
                unsigned char byte = buf[sample_idx / 4];
                byte >>= 2 * (sample_idx % 4);
                *(out++) = lookup[byte & 0x03];
        }
        return;
}

// [[Rcpp::export]]
IntegerVector read_bed_fast(CharacterVector file_name,
                IntegerVector row_idx, IntegerVector col_idx,
                int n_samples)
{
        std::ifstream file(file_name[0]);
        int bytes_per_variant = (n_samples + 3) / 4;
        char buf[bytes_per_variant];
        file.read(buf, 3);
        if (buf[0] != 0x6C || buf[1] != 0x1B || buf[2] != 0x01)
                throw std::runtime_error("Invalid BED file format");

        int ans_nrow = row_idx.size(), ans_ncol = col_idx.size();
        IntegerMatrix ans(ans_nrow, ans_ncol);
        int *out = reinterpret_cast<int*>(dataptr(ans));

        for (int j = 0; j < ans_ncol; j++) {
                int variant_idx = col_idx[j] - 1;
                std::streampos pos = 3 + variant_idx * bytes_per_variant;
                file.seekg(pos);
                file.read(buf, bytes_per_variant);
                decode_variant(buf, row_idx, out + j * ans_nrow);
        }
        return ans;
}
