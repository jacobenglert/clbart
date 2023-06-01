
# Simulated dataset based on a single tree --------------------------------


# Set Simulation Parameters -----------------------------------------------
set.seed(1234)        # Random Seed
n     <- 10000       # Initial Population Size
t     <- 500          # Initial Observation Period Length
p_c   <- 5          # Number of confounders
p_m   <- 10           # Number of moderators
rho_c <- .5         # Correlation between confounders
rho_m <- 0.6          # Correlation between moderators
corr_mat_c <- rho_c * matrix(1, p_c, p_c) + (1 - rho_c) * diag(p_c)
corr_mat_m <- sapply(1:p_m, \(p) rho_m^(abs(p - (1:p_m))))
beta <- seq(-0.1, 0.1, length.out = p_c) # Confounder Effect Sizes

compute_logOR <- function(x){

  # Method 1: Single Tree
  t1 <- log(0.8) * I(x[,1] <= 0 & x[,2] <= 0) +
    log(1.1) * I(x[,1] <= 0 & x[,2] > 0) +
    log(1.3) * I(x[,1] > 0 & x[,3] <= 0) +
    log(1.5) * I(x[,1] > 0 & x[,3] > 0)

  return(t1)
}

# Simulate Data -----------------------------------------------------------

# Generate (un)correlated confounders from multivariate normal
amps <- sample(seq(0.1, 1, by = 0.1), p_c, replace = TRUE)
freqs <- sample(c(1, 2, 3, 4, 6), p_c, replace = TRUE)

x_c_trend <- mapply(\(A, f) A * sin(seq(0, 2 * pi * f, length.out = t)), A = amps, f = freqs)

x_c <- t(apply(x_c_trend, 1, \(x) round(mvtnorm::rmvnorm(1, mean = x, sigma = corr_mat_c), 1)))

colnames(x_c) <- paste0('XC', 1:p_c)

# Generate exposure
z <- rnorm(t, sin(seq(0, 2 * pi * 2, length.out = t)))

# Generate (un)correlated moderators from multivariate normal
x_m <- (round(mvtnorm::rmvnorm(n, sigma = corr_mat_m), 1) > 0) |>
  as.data.frame() |>
  lapply(as.numeric) |>
  as.data.frame()
colnames(x_m) <- paste0('XM', 1:p_m)

# Randomly select confounders and exposure for each individual in the initial population
z_idx <- sample(4:(t-3), n, replace = TRUE)
z_ord <- replace(z_idx %% 4, which(z_idx %% 4 == 0), 4) # sample(1:4, n, replace = TRUE)
x_c_obs <- x_c[z_idx,]
z_obs <- z[z_idx]

# Compute the heterogeneous exposure effects on the log-odds scale
x_m_used <- sample(1:ncol(x_m), 5)
tau <- compute_logOR(x_m[x_m_used])

# Generate random intercepts for each individual in the initial population
alpha <- rnorm(n, -3, 1)

# Simulate the binary outcome
expit <- function(x) exp(x) / (1 + exp(x))
p <- expit(alpha + x_c_obs %*% beta + tau * z_obs)
y <- rbinom(n, 1, p)

# Identify which observations experienced the event
event_idx <- which(y == 1)

# Data for initial population
data_obs <- data.frame(y = y[event_idx],
                       x_c_obs[event_idx,],
                       x_m[event_idx,],
                       z = z_obs[event_idx],
                       tau = tau[event_idx])


# Create the Case-Crossover Dataset ---------------------------------------

# Duplicate the moderators and exposure effects from the observations that
# experienced the event according to the length of the window
tau_cc <- tau[event_idx][rep(seq_len(length(event_idx)), each = 4)]
x_m_cc <- x_m[event_idx,][rep(seq_len(length(event_idx)), each = 4),] |>
  as.data.frame()

# Identify the windows indexes surrounding the observed exposure and confounders
z_window <- mapply(\(idx, ord){
  if(ord == 1) return(idx:(idx + 3))
  else if(ord == 2) return((idx-1):(idx+2))
  else if(ord == 3) return((idx-2):(idx+1))
  else if(ord == 4) return((idx-3):idx)},
  idx = z_idx, ord = z_ord, SIMPLIFY = FALSE)

# Create "control" exposures and confounders for each event
x_c_cc <- lapply(z_window[event_idx], \(w) x_c[w,]) |> do.call(what = rbind)
z_cc <- lapply(z_window[event_idx], \(w) z[w]) |> unlist()
y_cc <- mapply(\(w, idx) w == idx,
               w = z_window[event_idx], idx = z_idx[event_idx]) |>
  as.numeric()

# Create a strata variable indexing which observations correspond to each original event
strata_cc <- rep(1:length(event_idx), each = 4)

# Combine these to create the case-crossover dataset
data_cc <- data.frame(y = y_cc,
                      x_c_cc,
                      x_m_cc,
                      z = z_cc,
                      strata = strata_cc,
                      tau = tau_cc)


# Export ------------------------------------------------------------------
readr::write_csv(data_cc, 'data-raw/sim_data.csv')

usethis::use_data(sim_data, overwrite = TRUE)
