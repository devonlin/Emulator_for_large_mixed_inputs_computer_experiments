#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate the sum
double psum(NumericVector x1, NumericVector x2, NumericVector par2) {
  int n = x1.size();
  double result = 0.0;
  for (int i = 0; i < n; ++i) {
    result += -par2[i] * pow((x1[i] - x2[i]), 2);
  }
  return result;
}

// Inner function
double covx(NumericVector w1, NumericVector w2, NumericVector parv, int q, int p, IntegerVector m) {
  int npar = 1 + q + p + p * sum(m);
  NumericVector par1 = parv[Range(0, q)];
  NumericVector par2 = parv[Range(q + 1, q + p)];
  NumericVector par3 = parv[Range(q + 1 + p, npar - 1)];
  NumericVector x1 = w1[Range(0, p - 1)];
  NumericVector z1 = w1[Range(p, p + q - 1)];
  NumericVector x2 = w2[Range(0, p - 1)];
  NumericVector z2 = w2[Range(p, p + q - 1)];
  double res1 = par1[0] * exp(psum(x1, x2, par2));
  for (int i = 0; i < q; ++i) {
    if (z1[i] != z2[i]) {
      res1 += 0;
    } else {
      int l = z1[i];
      int start_index = sum(m[Range(0, i)]) * p - m[i] * p + (l - 1) * p;
      int end_index = sum(m[Range(0, i)]) * p - m[i] * p + (l - 1) * p + p - 1;
      double sum_term = psum(x1, x2, par3[Range(start_index, end_index)]);
      res1 += par1[i + 1] * exp(sum_term);
    }
  }
  return res1;
}

// Rcpp Exporting function
// [[Rcpp::export]]
NumericVector covx_wrapper(NumericVector w1, NumericVector w2, NumericVector parv, int q, int p, IntegerVector m) {
  double result = covx(w1, w2, parv, q, p, m);
  return NumericVector::create(result);
}
