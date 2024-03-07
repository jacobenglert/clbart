#include <RcppArmadillo.h>
// #include <unordered_map>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec clr_loglik_cpp(const arma::vec& beta,
                         const arma::mat& x,
                         const arma::vec& y,
                         const arma::ivec& stratum,
                         const arma::vec& offset,
                         bool total = true) {
  
  int n = x.n_rows;
  int max_stratum = arma::max(stratum);
  arma::vec a = arma::zeros<arma::vec>(max_stratum);
  arma::vec b = arma::zeros<arma::vec>(max_stratum);
  
  // Accumulate within-strata sums needed for log-likelihood computation
  for (int i = 0; i < n; ++i) {
    int stratum_idx = stratum[i] - 1;
    double linear_predictor = arma::as_scalar(x.row(i) * beta) + offset[i];
    a[stratum_idx] += y[i] * linear_predictor;
    b[stratum_idx] += exp(linear_predictor);
  }
  
  // Compute log-likelihood for each stratum
  arma::vec loglik_per_stratum = a - arma::log(b);
  
  // If total requested (default) return total log-likelihood
  if (total) {
    return arma::vec(1, arma::fill::value(arma::sum(loglik_per_stratum)));
  } else {
    return loglik_per_stratum; // Otherwise return vector of log-likelihoods
  }
  
}

// [[Rcpp::export]]
arma::vec clr_grd_cpp(const arma::colvec& beta,
                      const arma::mat& x,
                      const arma::vec& y,
                      const arma::ivec& stratum,
                      const arma::colvec& offset) {
  int n = x.n_rows;
  int p = x.n_cols;
  int max_stratum = arma::max(stratum);
  arma::mat a = arma::zeros<arma::mat>(max_stratum, p);
  arma::mat b = arma::zeros<arma::mat>(max_stratum, p);
  arma::vec c = arma::zeros<arma::vec>(max_stratum);
  arma::vec exp_sc = exp(x * beta + offset);
  arma::rowvec grd = arma::zeros<arma::rowvec>(p);
         
  // // Accumulate within-strata sums needed for gradient computation
  // for (int i = 0; i < n; ++i) {
  //   int idx = stratum[i] - 1;
  //   c[idx] += exp_sc[i];
  //   a.row(idx) += y[i] * x.row(i);
  //   b.row(idx) += exp_sc[i] * x.row(i);
  // }
  //        
  // // Compute the gradient
  // for (int k = 0; k < max_stratum; ++k) {
  //   grd += a.row(k) - b.row(k) / c[k];
  // }
  
  // For some reason, accumulating the gradient as you go is faster for now
  for (int i = 0; i < n; ++i) {
    int idx = stratum[i] - 1;
    int next_idx = (i + 1 < n) ? stratum[i + 1] - 1 : n - 1;
    
    c[idx] += exp_sc[i];
    
    for (int j = 0; j < p; j++) {
      a(idx, j) += y[i] * x(i,j);
      b(idx, j) += x(i,j) * exp_sc[i];
      
      if(idx != next_idx){
        grd[j] += a(idx, j) - (b(idx, j) / c[idx]);
      }
    }
    
  }
  
  // Transpose gradient to a column vector before returning   
  return grd.t();
}


// [[Rcpp::export]]
arma::mat clr_fish_cpp(const arma::colvec& beta,
                       const arma::mat& x,
                       const arma::vec& y,
                       const arma::ivec& stratum,
                       const arma::colvec& offset) {
  
  int n = x.n_rows;
  int p = x.n_cols;
  int max_stratum = arma::max(stratum);
  arma::vec   a = arma::zeros<arma::vec>(max_stratum);
  arma::cube  b = arma::zeros<arma::cube>(max_stratum, p, p);
  arma::mat   c = arma::zeros<arma::mat>(max_stratum, p);
  arma::vec   exp_sc = exp(x * beta + offset);
  arma::mat   I = arma::zeros<arma::mat>(p, p);
  
  // Accumulate within-strata sums needed for Fisher information computation
  for (int i = 0; i < n; ++i) {
    int idx = stratum[i] - 1;
    a[idx] += exp_sc[i];
    for (int j = 0; j < p; ++j) {
      c(idx, j) += x(i, j) * exp_sc[i];
      for (int k = 0; k <= j; ++k) { // Only compute lower triangle of b
        b(idx, j, k) += x(i, j) * x(i, k) * exp_sc[i];
      }
    }
  }
  
  // Calculate the Fisher Information (only lower triangle)
  for (int i = 0; i < max_stratum; ++i) {
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k <= j; ++k) {
        double term = ((a[i] * b(i, j, k)) - (c(i, j) * c(i, k))) / pow(a[i], 2);
        I(j, k) += term;
        if (j != k) {
          I(k, j) = I(j, k); // Symmetrically fill the upper triangle
        }
      }
    }
  }
  
  return I;
}

// [[Rcpp::export]]
Rcpp::List clr_cpp(const arma::mat& x,
                   const arma::vec& y,
                   const arma::ivec& stratum,
                   const Rcpp::Nullable<arma::vec>& offset = R_NilValue,
                   const Rcpp::Nullable<arma::vec>& start_values = R_NilValue,
                   const int max_iter = 100,
                   const double tol = 1e-7) {
  
  // Initialize beta vector, log-likelihood, gradient, and Fisher information
  int n = x.n_rows;
  int p = x.n_cols;
  arma::vec beta;
  
  // Use starting values if provided, or initialize at the zero vector
  if (start_values.isNotNull()) {
    beta = Rcpp::as<arma::vec>(start_values);
  } else {
    beta = arma::zeros<arma::vec>(p);
  }
  
  // Set offset if provided
  arma::vec offset_vec;
  if (offset.isNotNull()) {
    offset_vec = Rcpp::as<arma::vec>(offset);
  } else {
    offset_vec = arma::zeros<arma::vec>(n);
  }
  
  arma::vec beta_old;
  arma::mat I_inv;
  
  // Use Fisher scoring to obtain estimates of beta
  for (int i = 0; i < max_iter; ++i) {
    beta_old = beta;
    arma::vec grd = clr_grd_cpp(beta, x, y, stratum, offset_vec);
    arma::mat I = clr_fish_cpp(beta, x, y, stratum, offset_vec);
    
    // Invert the fisher information matrix
    I_inv = arma::inv(I);
    
    // Update beta
    beta += I_inv * grd;
    
    // Check for convergence
    double max_change = arma::max(arma::abs(beta - beta_old));
    if (max_change < tol) {
      break;
    }
  }
  
  // Return final estimate of beta
  return Rcpp::List::create(Rcpp::Named("beta") = beta,
                            Rcpp::Named("vcov") = I_inv);
}


// Versions which do not required stratum to be pre-sorted
// // [[Rcpp::export]]
// arma::vec clr_loglik_cpp(const arma::vec& beta,
//                          const arma::mat& x,
//                          const arma::vec& y,
//                          const arma::ivec& stratum,
//                          const arma::vec& offset,
//                          bool return_sum = true) {
//   
//   int n = x.n_rows;
//   std::unordered_map<int, int> stratum_map;
//   int stratum_count = 0;
//   
//   // Map stratum to contiguous indices
//   for (int i = 0; i < n; ++i) {
//     if (stratum_map.find(stratum[i]) == stratum_map.end()) {
//       stratum_map[stratum[i]] = stratum_count++;
//     }
//   }
//   
//   arma::vec sum_linear_predictor = arma::zeros<arma::vec>(stratum_count);
//   arma::vec sum_exp_linear_predictor = arma::zeros<arma::vec>(stratum_count);
//   //arma::vec linear_predictor = x * beta + offset;
//   for (int i = 0; i < n; ++i) {
//     int mapped_idx = stratum_map[stratum[i]];
//     double linear_predictor = arma::as_scalar(x.row(i) * beta) + offset[i];
//     sum_linear_predictor[mapped_idx] += y[i] * linear_predictor;
//     sum_exp_linear_predictor[mapped_idx] += exp(linear_predictor);
//   }
//   
//   arma::vec loglik_per_stratum = sum_linear_predictor - arma::log(sum_exp_linear_predictor);
//   
//   if (return_sum) {
//     return arma::vec(1, arma::fill::value(arma::sum(loglik_per_stratum)));
//   } else {
//     // Create a vector with the log-likelihoods in the order of the original stratum
//     arma::vec loglik_original_order(stratum_count);
//     for (int i = 0; i < stratum_count; ++i) {
//       loglik_original_order[i] = loglik_per_stratum[stratum_map[stratum[i]]];
//     }
//     return loglik_original_order;
//   }
// }
// // [[Rcpp::export]]
// arma::vec clr_grd_cpp2(const arma::colvec& beta,
//                        const arma::mat& x,
//                        const arma::vec& y,
//                        const arma::ivec& stratum,
//                        const arma::colvec& offset) {
//   
//   int n = x.n_rows;
//   int p = x.n_cols;
//   std::unordered_map<int, int> stratum_map;
//   int stratum_count = 0;
//   
//   // Map stratum to contiguous indices
//   for (int i = 0; i < n; ++i) {
//     if (stratum_map.find(stratum[i]) == stratum_map.end()) {
//       stratum_map[stratum[i]] = stratum_count++;
//     }
//   }
//   
//   arma::mat a = arma::zeros<arma::mat>(stratum_count, p);
//   arma::mat b = arma::zeros<arma::mat>(stratum_count, p);
//   arma::vec c = arma::zeros<arma::vec>(stratum_count);
//   arma::vec exp_sc = exp(x * beta + offset);
//   arma::rowvec grd = arma::zeros<arma::rowvec>(p);
//   
//   // Accumulate sums in a, b, and c
//   for (int i = 0; i < n; ++i) {
//     int mapped_idx = stratum_map[stratum[i]];
//     c[mapped_idx] += exp_sc[i];
//     a.row(mapped_idx) += y[i] * x.row(i);
//     b.row(mapped_idx) += exp_sc[i] * x.row(i);
//   }
//   
//   // Compute the gradient
//   for (int k = 0; k < stratum_count; ++k) {
//     grd += a.row(k) - b.row(k) / c[k];
//   }
//   
//   return grd.t();
// }
                                                            