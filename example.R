source("lmFBST.R")

# Water droplet dataset
elapsed_time <- seq(from = 0, to = 7, by = .5)
s <- c(4.7, 8.1, 11.5, 14.5, 17.4, 19.9, 22.4, 24.6, 26.6, 28.5, 30.1, 31.5,
       32.8, 33.8, 34.6)
delta_s <- c(NA, diff(s))
Ks <- 8.446
drop_rad <- sqrt(2*delta_s*Ks)

x <- elapsed_time
y <- drop_rad
Y <- y[-1]
X <- X_new <- t(t(x))
X <- X[-1,,drop = F]

# Thresholds
eps2 <- 0.621835/sqrt(length(Y)+1)
eps1 <- eps2/3 # more restrictive

# Prior specification
sigma2 <- 0.01
mu.fun <- function(X) rep(6, dim(X)[1])
my_ker <- rbfdot(sigma = .5)
ker.fun <- function(X1, X2)  kernelMatrix(my_ker, X1, X2)

# Linear model test on finite covariate space
lmFBST(Y, X, mu.fun, ker.fun, sigma2, alpha = .05, eps = 0, X_new=X_new)
lmFBST(Y, X, mu.fun, ker.fun, sigma2, alpha = .05, eps = eps1, X_new=X_new)
lmFBST(Y, X, mu.fun, ker.fun, sigma2, alpha = .05, eps = eps2, X_new=X_new)

# Variable selection on finite covariate space
lmFBST(Y, X, mu.fun, ker.fun, sigma2, alpha = .05, eps = 0,
       add_intercept = T, vars = 1, X_new=X_new)
lmFBST(Y, X, mu.fun, ker.fun, sigma2, alpha = .05, eps = eps1,
       add_intercept = T, vars = 1, X_new=X_new)
lmFBST(Y, X, mu.fun, ker.fun, sigma2, alpha = .05, eps = eps2,
       add_intercept = T, vars = 1, X_new=X_new)

# Linear model test on infinite covariate space
lmFBST(Y, X, mu.fun, ker.fun, sigma2, alpha = 0.05, eps = 0,
       finite_X = F, B = 10000, seed = 42)
lmFBST(Y, X, mu.fun, ker.fun, sigma2, alpha = 0.05, eps = eps1,
       finite_X = F, B = 10000, B.X = 1000, seed = 42)
lmFBST(Y, X, mu.fun, ker.fun, sigma2, alpha = 0.05, eps = eps2,
       finite_X = F, B = 10000, B.X = 1000, seed = 42)

# Variable selection on infinite covariate space
lmFBST(Y, X, mu.fun, ker.fun, sigma2, alpha = 0.05, eps = 0,
       add_intercept = T, vars = 1,
       finite_X = F, B = 10000, seed = 42)
lmFBST(Y, X, mu.fun, ker.fun, sigma2, alpha = 0.05, eps = eps1,
       add_intercept = T, vars = 1,
       finite_X = F, B = 10000, B.X = 1000, seed = 42)
lmFBST(Y, X, mu.fun, ker.fun, sigma2, alpha = 0.05, eps = eps2,
       add_intercept = T, vars = 1,
       finite_X = F, B = 10000, B.X = 1000, seed = 42)
