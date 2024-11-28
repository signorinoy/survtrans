#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List approx_likelihood(const arma::vec& offset, const arma::vec& time, const arma::vec& status) {
  int n_samples = time.n_elem;
  if (n_samples == 1) {
    return List::create(Named("weights") = 0, Named("residuals") = 0);
  }

  // Update hazard
  arma::vec haz = arma::exp(offset);

  // Update risk set
  arma::vec risk_set = arma::cumsum(haz);
  for (int i = n_samples - 1; i > 0; --i) {
    if (time[i] == time[i - 1]) {
      risk_set[i - 1] = risk_set[i];
    }
  }

  // Update weights and residuals
  arma::vec weights(n_samples, arma::fill::zeros);
  arma::vec residuals(n_samples, arma::fill::zeros);
  double denominator_risk_set = 0.0;
  double denominator_risk_set_sq = 0.0;
  int i = n_samples - 1;
  while (i >= 0) {
    double ti = time[i];
    int c = 0;
    int k = i;
    while (k >= 0 && ti == time[k]) {
      if (status[k] == 1) {
        c += 1;
      }
      k -= 1;
    }

    double v = risk_set[i];
    denominator_risk_set += c / v;
    denominator_risk_set_sq += c / (v * v);

    while (i > k) {
      v = haz[i];
      weights[i] = v * (denominator_risk_set - v * denominator_risk_set_sq);
      residuals[i] = (status[i] - v * denominator_risk_set) / weights[i];
      i -= 1;
    }
  }

  // Replace NaN residuals with 0
  residuals.elem(arma::find_nonfinite(residuals)).zeros();

  return List::create(Named("weights") = weights, Named("residuals") = residuals);
}