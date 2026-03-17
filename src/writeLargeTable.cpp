#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include <vector>

// Cross-platform 64-bit file seek.
// fseek() takes a `long` offset, which is 32-bit on Windows — arithmetic byte
// offsets into large files would overflow. Use _fseeki64 on Windows and
// fseeko (with 64-bit off_t enabled via _FILE_OFFSET_BITS=64) elsewhere.
#ifdef _WIN32
#include <io.h>
static inline int pkg_fseek(FILE* fp, long long offset, int origin) {
  return _fseeki64(fp, (__int64)offset, origin);
}
#else
#include <sys/types.h>
static inline int pkg_fseek(FILE* fp, long long offset, int origin) {
  return fseeko(fp, (off_t)offset, origin);
}
#endif

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Wrapper to suppress -Wunused-result for intentionally discarded fread return values.
// The file format is fixed-width and written by this package; errors are non-recoverable.
static void fread_discard(void* buf, size_t size, FILE* fp) {
  size_t n = fread(buf, 1, size, fp);
  (void)n;
}

// Parse a fixed-width character field as a double. Warns (once per call site)
// if the string cannot be parsed, returning 0.0 in that case.
static double safe_atof(const char* s, int var, int sample) {
  char* endptr;
  double v = strtod(s, &endptr);
  if (endptr == s)
    Rf_warning("calcCorr: could not parse value for variable %d, sample %d; treating as 0",
               var, sample);
  return v;
}

// Replace values with their ranks (average rank for ties). O(n log n).
void rankVec(NumericVector& x, int n) {
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int a, int b) { return x[a] < x[b]; });
  std::vector<double> ranks(n);
  int i = 0;
  while (i < n) {
    int j = i;
    while (j < n - 1 && x[idx[j + 1]] == x[idx[j]]) j++;
    double avg = (i + j) / 2.0 + 1.0;
    for (int k = i; k <= j; k++) ranks[idx[k]] = avg;
    i = j + 1;
  }
  for (int k = 0; k < n; k++) x[k] = ranks[k];
}

void stdvec(NumericVector& x, int n) {
  double mnx = mean(x);
  double sdx = sd(x);
  if (sdx == 0.0) {
    Rf_warning("stdvec: zero standard deviation — vector has constant values and cannot be standardized");
    return;
  }
  for (int i = 0; i < n; i++) {
    x[i] = (x[i] - mnx) / sdx;
  }
}

double corr(const NumericVector& x, const NumericVector& y, int n) {
  if (n <= 1) {
    Rf_warning("corr: sample size n=%d is too small to compute correlation (need n >= 2)", n);
    return 0.0;
  }
  double s = 0.0;
  for (int i = 0; i < n; i++) s += x[i] * y[i];
  return s / (n - 1); // to match R
}

//' Calculate pairwise correlations and print to a file.
//'
//' @param fls A character vector which includes the names of the input and output files.
//' @param MaxNameLen The length of the first column in the input data file (predictor names).
//' @param recSize The total length of input rows (fixed-width format).
//' @param P The number of predictors.
//' @param n The sample size.
//' @param zeroSD If the standard deviation of a column is below this threshold, this column will be considered as uncorrelated with all other columns (Default=1e-3).
//' @param spearman If TRUE, compute Spearman rank correlations instead of Pearson (Default=FALSE).
// [[Rcpp::export]]
void calcCorr(CharacterVector fls, int MaxNameLen, int recSize, int P, int n, double zeroSD = 1e-3, bool spearman = false) {
  if (n <= 0)
    Rf_error("calcCorr: n must be positive; got n = %d", n);
  if (P <= 1)
    Rf_error("calcCorr: P must be >= 2; got P = %d", P);
  if (recSize <= MaxNameLen + 1)
    Rf_error("calcCorr: recSize (%d) must be > MaxNameLen + 1 (%d); check input file parameters",
             recSize, MaxNameLen + 1);
  FILE *finptr = NULL;
  FILE *foutptr = NULL;
  finptr = fopen(fls[0], "rb");
  if (finptr == NULL) {
    Rf_error("Cannot open input file: %s", (const char *)fls[0]);
  }
  foutptr = fopen(fls[1], "wb");
  if (foutptr == NULL) {
    fclose(finptr);
    Rf_error("Cannot open output file: %s", (const char *)fls[1]);
  }
  int numLen = (recSize - 1 - MaxNameLen) / n;
  std::vector<char> pname(MaxNameLen);
  std::vector<char> val(numLen);
  NumericVector x1(n);
  NumericVector x2(n);
  for (int i = 0; i < P - 1; i++) {
    pkg_fseek(finptr, (long long)i * recSize, SEEK_SET);
    fread_discard(pname.data(), MaxNameLen, finptr);
    for (int j = 0; j < n; j++) {
      fread_discard(val.data(), numLen, finptr);
      x1[j] = safe_atof(val.data(), i + 1, j + 1);
    }
    if (spearman) rankVec(x1, n);
    if (sd(x1) < zeroSD) {
      for (int j = i + 1; j < P; j++) {
        fprintf(foutptr, "%d|%d|0|0\n", i + 1, j + 1);
      }
    } else {
      stdvec(x1, n);
      for (int k = i + 1; k < P; k++) {
        pkg_fseek(finptr, (long long)k * recSize, SEEK_SET);
        fread_discard(pname.data(), MaxNameLen, finptr);
        for (int j = 0; j < n; j++) {
          fread_discard(val.data(), numLen, finptr);
          x2[j] = safe_atof(val.data(), k + 1, j + 1);
        }
        if (spearman) rankVec(x2, n);
        if (sd(x2) < zeroSD) {
          fprintf(foutptr, "%d|%d|0|0\n", i + 1, k + 1);
        } else {
          stdvec(x2, n);
          double tmpcor = corr(x1, x2, n);
          fprintf(foutptr, "%d|%d|%lf|%lf\n", i + 1, k + 1, tmpcor, 1.0 - tmpcor * tmpcor);
        }
      }
    }
  }
  fclose(finptr);
  fclose(foutptr);
}
