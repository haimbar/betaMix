#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

//[[Rcpp::depends(RcppArmadillo)]]                                                             
using namespace Rcpp;


void stdvec(NumericVector& x, int n) {
  double mnx = mean(x);
  double sdx = sd(x);
  for (int i=0; i<n; i++) {
    x[i] = (x[i] - mnx)/sdx;
  }
}

double corr(const NumericVector& x, const NumericVector& y, int n) {
  double s = 0.0;
  for (int i = 0; i < n; i++) s += x[i] * y[i];
  return s / (n-1); // to match R
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
void calcCorr(CharacterVector fls, int MaxNameLen,  int recSize, int P, int n, double zeroSD=1e-3) {
  FILE * finptr = NULL;
  FILE * foutptr = NULL;
  finptr = fopen(fls[0], "r");
  if (finptr == NULL) {
    printf("Cannot open input file\n");
    return;
  }
  foutptr = fopen(fls[1], "w");
  if (foutptr == NULL) {
    printf("Cannot open output file\n");
    return;
  }
  int numLen = (recSize-1-MaxNameLen)/n;
  //printf("%d\n",numLen);
  std::vector<char> pname(MaxNameLen);
  std::vector<char> val(numLen);
  NumericVector x1(n);
  NumericVector x2(n);
  for (int i=0; i < P-1; i++) {
    fseek(finptr, i*recSize, SEEK_SET);
    int acread = fread(pname.data(), MaxNameLen, 1, finptr);
    //printf("%s,%d\n",pname.data(),i*recSize);
    for(int j=0; j<n; j++){
      acread = fread(val.data(), numLen, 1, finptr);
      x1[j] = atof(val.data());
      //printf("%d,%lf\t",j,atof(val.data()));
    }
    if (sd(x1) < zeroSD) {
      for (int j=i+1; j<P; j++) {
        fprintf(foutptr, "%d|%d|0|0\n",i+1,j+1);
      }
    } else {
      stdvec(x1, n);
      for (int k=i+1; k<P; k++) {
        fseek(finptr, k*recSize, SEEK_SET);
        int acread = fread(pname.data(), MaxNameLen, 1, finptr);
        for(int j=0; j<n; j++){
          acread = fread(val.data(), numLen, 1, finptr);
          x2[j] = atof(val.data());
          //printf("%d,%lf\t",j,atof(val.data()));
        }
        if (sd(x2) < zeroSD) {
          fprintf(foutptr,"%d|%d|0|0\n",i+1,k+1);
        } else {
          stdvec(x2, n);
          double tmpcor = corr(x1,x2,n);
          //printf("%f\n",tmpcor);
          fprintf(foutptr,"%d|%d|%lf|%lf\n",i+1,k+1, tmpcor, 1.0 - tmpcor*tmpcor);
        }
      }
    }
  }
  fclose(finptr);
  fclose(foutptr);
}
