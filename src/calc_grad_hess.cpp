#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Optimized helper function to calculate the maximum value within each group
NumericVector ave_max(NumericVector x, NumericVector group) {
  int n = x.size();
  // Pre-compute maximum values for each group
  for (int i = n - 2; i >= 0; --i) {
    if (group[i] == group[i + 1]) {
      x[i] = std::max(x[i], x[i + 1]);
    }
  }
  return x;
}

// [[Rcpp::export]]
List calc_grad_hess(NumericVector lp, NumericMatrix x, 
                    NumericVector time, NumericVector status) {
  int n_samples = time.size();
  int n_features = x.ncol();

  if (n_samples == 1) {
    return List::create(
      Named("grad") = NumericVector(n_features, 0.0),
      Named("hess") = NumericVector(n_features * n_features, 0.0)
    );
  }

  NumericVector hazard = exp(lp);
  NumericVector cumsum_hazard = cumsum(hazard);
  NumericVector risk_set = ave_max(cumsum_hazard, time);

  // Precompute hazard_x for gradient calculations
  mat arma_x(x.begin(), n_samples, n_features, false);
  mat hazard_x = arma_x.each_col() % vec(hazard.begin(), hazard.size(), false);

  // Cumulative sum for each feature
  mat risk_set_x = hazard_x;
  for (int j = 0; j < n_features; j++) {
    vec col = hazard_x.col(j);
    vec cumsum_col = cumsum(col);
    risk_set_x.col(j) = vec(ave_max(wrap(cumsum_col), time).begin(), n_samples, false);
  }

  // Compute risk_set_x_ratio
  mat risk_set_x_ratio = risk_set_x.each_col() / vec(risk_set.begin(), risk_set.size(), false);

  // Compute outer product for Hessian calculation
  cube hazard_xx(n_features, n_features, n_samples, fill::zeros);
  #pragma omp parallel for
  for (int i = 0; i < n_samples; i++) {
    hazard_xx.slice(i) = hazard[i] * (arma_x.row(i).t() * arma_x.row(i));
  }

  // Compute risk_set_xx_ratio
  cube risk_set_xx_ratio(n_features, n_features, n_samples);
  for (int j = 0; j < n_features; j++) {
    for (int k = 0; k < n_features; k++) {
      vec col = hazard_xx.tube(j, k);
      vec cumsum_col = cumsum(col);
      vec max_col = vec(ave_max(wrap(cumsum_col), time).begin(), n_samples, false);
      risk_set_xx_ratio.tube(j, k) = max_col / vec(risk_set.begin(), risk_set.size(), false);
    }
  }

  // Compute risk_set_x2_ratio
  cube risk_set_x2_ratio(n_features, n_features, n_samples);
  #pragma omp parallel for
  for (int i = 0; i < n_samples; i++) {
    mat outer_product = risk_set_x.row(i).t() * risk_set_x.row(i);
    risk_set_x2_ratio.slice(i) = outer_product / (risk_set[i] * risk_set[i]);
  }

  // Gradient and Hessian initialization
  mat grad(n_samples, n_features);
  mat hess(n_samples, n_features * n_features);

  // Parallel loop for gradient and Hessian computation
  #pragma omp parallel for
  for (int i = 0; i < n_samples; i++) {
    grad.row(i) = (arma_x.row(i) - risk_set_x_ratio.row(i)) * status[i];
    hess.row(i) = vectorise(risk_set_x2_ratio.slice(i) - risk_set_xx_ratio.slice(i)).t() * status[i];
  }

  return List::create(
    Named("grad") = wrap(grad),
    Named("hess") = wrap(hess)
  );
}