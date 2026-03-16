#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Wrapper to suppress -Wunused-result for intentionally discarded fread return values.
// The file format is fixed-width and written by this package; errors are non-recoverable.
static void fread_discard(void* buf, size_t size, FILE* fp) {
  size_t n = fread(buf, 1, size, fp);
  (void)n;
}

void stdvec(NumericVector& x, int n) {
  double mnx = mean(x);
  double sdx = sd(x);
  for (int i = 0; i < n; i++) {
    x[i] = (x[i] - mnx) / sdx;
  }
}

double corr(const NumericVector& x, const NumericVector& y, int n) {
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
// [[Rcpp::export]]
void calcCorr(CharacterVector fls, int MaxNameLen, int recSize, int P, int n, double zeroSD = 1e-3) {
  FILE *finptr = NULL;
  FILE *foutptr = NULL;
  finptr = fopen(fls[0], "r");
  if (finptr == NULL) {
    Rf_error("Cannot open input file: %s", (const char *)fls[0]);
  }
  foutptr = fopen(fls[1], "w");
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
    fseek(finptr, i * recSize, SEEK_SET);
    fread_discard(pname.data(), MaxNameLen, finptr);
    for (int j = 0; j < n; j++) {
      fread_discard(val.data(), numLen, finptr);
      x1[j] = atof(val.data());
    }
    if (sd(x1) < zeroSD) {
      for (int j = i + 1; j < P; j++) {
        fprintf(foutptr, "%d|%d|0|0\n", i + 1, j + 1);
      }
    } else {
      stdvec(x1, n);
      for (int k = i + 1; k < P; k++) {
        fseek(finptr, k * recSize, SEEK_SET);
        fread_discard(pname.data(), MaxNameLen, finptr);
        for (int j = 0; j < n; j++) {
          fread_discard(val.data(), numLen, finptr);
          x2[j] = atof(val.data());
        }
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
