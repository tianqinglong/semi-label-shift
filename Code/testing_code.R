source("g_c_1.R")
source("rej-sample.R")
source("fit_model.R")
source("evaluate.R")

n <- 1000
m <- 500
mu <- c(0, 0, 0)
sigma <- diag(3)*0.8
beta <- c(1, 0.4, 0.4, 0.4, 0.5, 0.5)
sd <- 0.75

beta_rho <- c(-0.1, 0.05)

xn <- gen_x_marginal(1000, mu, sigma)
yn <- gen_y(xn, beta, sd)
max_out <- max_rho(yn, beta_rho)
ac_list <- rej_sample(yn, max_out)

t_samp <- make_target_sample(m, mu, sigma, beta, sd, beta_rho)
s_samp <- make_source_sample(n, mu, sigma, beta, sd, beta_rho)

dat <- s_samp
xn_t <- t_samp[,-1]
yx_all <- Generate_sampler(dat, xn_t, B = 1000)

foptim <- optim(c(0, 0.05), Compute_S_eff, dat = dat, xn_t = xn_s, yx_all = yx_all)
ESS <- Compute_beta_var(foptim$par, dat, xn_t, yx_all)
solve(ESS)
