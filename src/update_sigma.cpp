#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double cauchy_jacobian(double tau, double sigma_scale) {
  double sigma = pow(tau, -0.5);
  int give_log = 1;
  
  double out = Rf_dcauchy(sigma, 0.0, sigma_scale, give_log);
  out = out - M_LN2 - 3.0 / 2.0 * log(tau);
  
  return out;
  
}

// [[Rcpp::export]]
double update_sigma(const arma::vec& r, double sigma_scale, double sigma_old) {
  
  double SSE = dot(r,r);
  double n = r.size();
  
  double shape = 0.5 * n + 1.0;
  double scale = 2.0 / SSE;
  double sigma_prop = pow(Rf_rgamma(shape, scale), -0.5);
  
  double tau_prop = pow(sigma_prop, -2.0);
  double tau_old = pow(sigma_old, -2.0);

  // Compute prior ratio
  double loglik_rat = cauchy_jacobian(tau_prop, sigma_scale) -
    cauchy_jacobian(tau_old, sigma_scale);
  
  // // Compute Gaussian likelihood ratio
  // double loglik_data_prop = -n/2.0 * log(2 * M_PI * pow(sigma_prop, 2)) - (SSE / (2 * pow(sigma_prop, 2)));
  // double loglik_data_old = -n/2.0 * log(2 * M_PI * pow(sigma_old, 2)) - (SSE / (2 * pow(sigma_old, 2)));
  // 
  // // Add the data log-likelihoods to the loglik_rat
  // loglik_rat += (loglik_data_prop - loglik_data_old);
  // 
  // // Compute the proposal ratio
  // double prop_prop = Rf_dgamma(tau_prop, shape, scale, 1);
  // double prop_old = Rf_dgamma(tau_old, shape, scale, 1);
  // 
  // // Add the proposal ratio to the loglik
  // loglik_rat += (prop_old - prop_prop);
  
  return log(unif_rand()) < loglik_rat ? sigma_prop : sigma_old;
  
}
