mc_clbart <- function(w, x, y, z, strata, hypers, opts) {

  # Set up cluster
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # Load required packages
  parallel::clusterEvalQ(cl, {
    library(progress)
    library(mvtnorm)
    library(armspp)
    Rcpp::sourceCpp("/Users/jacob/Desktop/clbart3/src/update_sigma.cpp")
    Rcpp::sourceCpp("/Users/jacob/Desktop/clbart3/src/rlgam.cpp")
    Rcpp::sourceCpp("/Users/jacob/Desktop/clbart3/src/clr.cpp")
  })
  
  # Export function input and clbart functions
  parallel::clusterExport(cl, varlist = c("w","x","y","z","strata","hypers","opts",
                                          "clbart",
                                          #"clr_cpp","clr_loglik_cpp","clr_grd_cpp","clr_fish_cpp",
                                          "make_forest","get_node_map",
                                          "get_node_depth","node_depth_prior","p_moves","get_goodvars",
                                          "get_sibling","get_children","get_parent","get_ancestors",
                                          "get_mapped_observations",
                                          "update_beta","update_sigma_beta",
                                          "update_sigma_mu","update_mu_mu",
                                          "update_forest","grow","prune","change",
                                          "compute_grow_ratio","compute_prune_ratio","compute_change_ratio",
                                          "get_tree_predictions","get_forest_predictions","get_forest_posterior_predictions",
                                          "update_mu","slice_sampler","log_sum_exp",
                                          "update_alpha","update_s"))
  
  # Run MCMC in parallel
  results <- parallel::parLapply(cl, 1:opts$num_chains, function (chain_id) {
    opts_chain <- opts
    opts_chain$seed <- opts$seed + chain_id
    clbart(w, x, y, z, strata, hypers, opts_chain)
  })
  


  # if (combine_results) {
  #   
  # }
  # else return(results)
  
  return (results)
}
