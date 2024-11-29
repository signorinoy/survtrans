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
arma::vec close_update(const arma::vec& z, const arma::vec& v, std::string penalty, double lambda, double gamma) {
  int n = z.n_elem;
  arma::vec eta(n);

  // Convert penalty to lowercase
  std::transform(penalty.begin(), penalty.end(), penalty.begin(), ::tolower);

  if (penalty == "lasso") {
    for (int i = 0; i < n; i++) {
      eta[i] = soft_threshold(z[i], lambda) / v[i];
    }
  } else if (penalty == "mcp") {
    for (int i = 0; i < n; i++) {
      if (std::abs(z[i]) <= gamma * lambda * v[i]) {
        eta[i] = soft_threshold(z[i], lambda) / (v[i] - 1.0 / gamma);
      } else {
        eta[i] = z[i] / v[i];
      }
    }
  } else if (penalty == "scad") {
    for (int i = 0; i < n; i++) {
      if (std::abs(z[i]) <= (v[i] + 1.0) * lambda) {
        eta[i] = soft_threshold(z[i], lambda) / v[i];
      } else if (std::abs(z[i]) <= gamma * lambda * v[i]) {
        eta[i] = soft_threshold(z[i], gamma * lambda / (gamma - 1.0));
        eta[i] = eta[i] / (v[i] - 1.0 / (gamma - 1.0));
      } else {
        eta[i] = z[i] / v[i];
      }
    }
  }

  return eta;
}