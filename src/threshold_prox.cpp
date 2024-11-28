#include <RcppArmadillo.h>
#include <algorithm>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Helper function for soft thresholding
inline double soft_threshold(double x, double lambda) {
  double threshold = std::abs(x) - lambda;
  return (threshold < 0) ? 0 : threshold * ((x > 0) ? 1 : -1);
}

// [[Rcpp::export]]
arma::vec threshold_prox(const arma::vec& y, double vartheta, std::string penalty, double lambda, double gamma) {
  int n = y.n_elem;
  arma::vec x(n);

  // Convert penalty to lowercase
  std::transform(penalty.begin(), penalty.end(), penalty.begin(), ::tolower);

  if (penalty == "lasso") {
    for (int i = 0; i < n; i++) {
      x[i] = soft_threshold(y[i], lambda / vartheta);
    }
  } else if (penalty == "mcp") {
    double const_val = 1 - 1 / (gamma * vartheta - 1);
    for (int i = 0; i < n; i++) {
      if (std::abs(y[i]) <= lambda * gamma) {
        x[i] = soft_threshold(y[i], lambda / vartheta) * const_val;
      } else {
        x[i] = y[i];
      }
    }
  } else if (penalty == "scad") {
    double const_val = vartheta * (gamma - 1);
    for (int i = 0; i < n; i++) {
      if (std::abs(y[i]) <= lambda + lambda / vartheta) {
        x[i] = soft_threshold(y[i], lambda / vartheta);
      } else if (std::abs(y[i]) <= lambda * gamma) {
        x[i] = soft_threshold(y[i], lambda * gamma / const_val) * const_val / (const_val - 1);
      } else {
        x[i] = y[i];
      }
    }
  }

  return x;
}