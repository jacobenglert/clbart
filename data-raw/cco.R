

# Set Simulation Parameters -----------------------------------------------
set.seed(2187)
t <- 5
p_x <- 9
p_w <- 9
n <- 100
s <- rep(1:n, each = t)

# Simulate Data -----------------------------------------------------------
x <- round(matrix(stats::rnorm(n * t * p_x), nrow = n * t, ncol = p_x), 2)
colnames(x) <- paste0('X', 1:p_x)
w <- round(matrix(stats::rnorm(n * p_w), nrow = n, ncol = p_w), 2) |>
  kronecker(Y = rep(1, t)) |>
  data.frame()
colnames(w) <- paste0('W', 1:p_w)
z <- round(rnorm(n * t), 2)

beta <- sample(-3:3, size = p_x, replace = TRUE)
x_eff <- x %*% beta
x_eff <- (x_eff - mean(x_eff)) / stats::sd(x_eff)
w_eff <- (10 * sin(pi * w[,1] * w[,2]) + 20 * (w[,3] - .5)^2 + 10 * w[,4] + 5 * w[,5]) * z
w_eff <- (w_eff - mean(w_eff)) / stats::sd(w_eff)

log_odds <- x_eff + w_eff
p <- (exp(log_odds) / (1 + exp(log_odds)))
y <- tapply(p, s, \(p) stats::rmultinom(1, size = 1, prob = p)) |> unlist(use.names = FALSE)


# Export Data -------------------------------------------------------------
cco <- data.frame(x, w, Z = z, Strata = s, Y = y)
write.csv(cco, "data-raw/cco.csv")
usethis::use_data(cco, overwrite = TRUE)
