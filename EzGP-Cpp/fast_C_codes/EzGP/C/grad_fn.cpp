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

// 1. Derivative function of sigma_0 for wi, wj
// [[Rcpp::export]]
NumericVector gradf_var_cpp(NumericVector w1, NumericVector w2, NumericVector parv, int q, int p) {
  NumericVector x1 = w1[Range(0, p - 1)];
  NumericVector x2 = w2[Range(0, p - 1)];
  NumericVector par2 = parv[Range(q + 1, q + p)];
  double gx = exp(psum(x1, x2, par2));
  return NumericVector::create(gx);
}

// 2.Derivative function of sigma_h for wi, wj
// [[Rcpp::export]]
NumericVector gradf_var1_cpp(NumericVector w1, NumericVector w2, NumericVector parv, int q, int p, int h, IntegerVector m) {
  NumericVector x1 = w1[Range(0, p - 1)];
  NumericVector x2 = w2[Range(0, p - 1)];
  NumericVector z1 = w1[Range(p, p + q - 1)];
  NumericVector z2 = w2[Range(p, p + q - 1)];
  NumericVector par3 = parv[Range(q + 1 + p, parv.size() - 1)];
  if (z1[h - 1] != z2[h - 1]) {
    return NumericVector::create(0);
  } else {
    int l = z1[h - 1];
    int start_index = sum(m[Range(0, h - 1)]) * p - m[h - 1] * p + (l - 1) * p;
    int end_index = sum(m[Range(0, h - 1)]) * p - m[h - 1] * p + (l - 1) * p + p - 1;
    double gx = exp(psum(x1, x2, par3[Range(start_index, end_index)]));
    return NumericVector::create(gx);
  }
}



// 3.Derivative function of theta_0_s for wi, wj
// [[Rcpp::export]]
NumericVector gradf_cor0_cpp(NumericVector w1, NumericVector w2, int s, NumericVector parv, int q, int p, IntegerVector m) {
  NumericVector x1 = w1[Range(0, p - 1)];
  NumericVector x2 = w2[Range(0, p - 1)];
  NumericVector par2 = parv[Range(q + 1, q + p)];
  double gx = -parv[0] * pow((x1[s - 1] - x2[s - 1]), 2) * exp(psum(x1, x2, par2));
  return NumericVector::create(gx);
}

// 4. Derivative function of theta_0_s for wi, wj
// [[Rcpp::export]]
NumericVector gradf_corhs_cpp(NumericVector w1, NumericVector w2, int h, int l, int s, NumericVector parv, int q, int p, IntegerVector m) {
  NumericVector x1 = w1[Range(0, p - 1)];
  NumericVector z1 = w1[Range(p, p + q - 1)];
  NumericVector x2 = w2[Range(0, p - 1)];
  NumericVector z2 = w2[Range(p, p + q - 1)];
  NumericVector par1 = parv[Range(0, q)];
  NumericVector par3 = parv[Range(q + 1 + p, parv.size() - 1)];
  if (z1[h - 1] != l || z2[h - 1] != l) {
    return NumericVector::create(0);
  } else {
    int start_index = sum(m[Range(0, h - 1)]) * p - m[h - 1] * p + (l - 1) * p;
    int end_index = sum(m[Range(0, h - 1)]) * p - m[h - 1] * p + (l - 1) * p + p - 1;
    double gx = -par1[h] * pow((x1[s - 1] - x2[s - 1]), 2) * exp(psum(x1, x2, par3[Range(start_index, end_index)]));
    return NumericVector::create(gx);
  }
}