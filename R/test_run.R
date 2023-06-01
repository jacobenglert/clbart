
sim_data <- readr::read_csv('data-raw/sim_data.csv')

w <- dplyr::select(sim_data, XM1:XM10)
x <- dplyr::select(sim_data, XC1:XC5)
y <- sim_data$y
z <- sim_data$z
strata <- sim_data$strata

fit <- clbart(w = dplyr::select(sim_data, XM1:XM10),
              x = dplyr::select(sim_data, XC1:XC5),
              y = sim_data$y, z = sim_data$z, strata = sim_data$strata,
              num_trees = 1, iter = 5000)

survival::clogit(y ~ as.matrix(x) + z + strata(strata))
