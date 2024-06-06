library(kernlab)
library(MASS)

my_inv <- function(mat){
  mat_inv <- try(chol2inv(chol(mat)), silent = T)
  if(any(class(mat_inv) == 'try-error')) mat_inv <- ginv(mat)
  mat_inv
}

Ks.fun <- function(x, a, b, A_inv, B_inv){
  part.A <- A_inv/(1-x)
  part.B <- B_inv/x
  (1 - t(a - b)%*%my_inv(part.A + part.B)%*%(a - b))[1]
}

inter_ellipsoid <- function(a, b, A_inv, B_inv){
  optim(.5, Ks.fun, a = a, b = b, A_inv = A_inv, B_inv = B_inv,
        method = "Brent", lower = 0, upper = 1)$value > 0
}

lm_ellipsoid <- function(X, mu, Sig_inv){
  beta.hat <- my_inv(t(X)%*%Sig_inv%*%X)%*%t(X)%*%Sig_inv%*%mu
  fitted.values <- X%*%beta.hat
  t(fitted.values - mu)%*%Sig_inv%*%(fitted.values - mu)
}

GP_pred <- function(Y, X, X_new, sigma2, mu.fun, ker.fun){
  mu.X <- mu.fun(X)
  mu.Xn <- mu.fun(X_new)
  ker.XX <- ker.fun(X,X)
  mat_inv <- my_inv(ker.XX + diag(sigma2, length(Y)))
  ker.XnX <- ker.fun(X_new,X)
  ker.XnXn <- ker.fun(X_new,X_new)
  f.mean <- mu.Xn + ker.XnX%*%mat_inv%*%(Y - mu.X)
  f.covmat <- ker.XnXn - ker.XnX%*%mat_inv%*%t(ker.XnX)
  list(f.mean = f.mean, f.covmat = f.covmat)
}

gp_hpd <- function(Y, X, sigma2, mu.fun, ker.fun, alpha = 0.05, B = 2000, seed = NULL){
  post.parm <- GP_pred(Y, X, X, sigma2, mu.fun, ker.fun)
  dim.X <- dim(X)[1]
  f.X <- numeric()
  set.seed(seed)
  f.X <- mvrnorm(n = B, mu = post.parm$f.mean, Sigma = post.parm$f.covmat)
  res_mat <- matrix(Y - f.X, ncol = B)
  hpd_cond <- apply(res_mat, 2, function(x) sum(x^2))
  quantile(hpd_cond, prob = 1-alpha)
}

lmFBST <- function(Y, X, mu.fun, ker.fun, sigma2 = NULL, alpha = 0.05, eps = 0,
                   add_intercept = T, vars = "all",
                   finite_X = T, X_new = X, probs.X = NULL,
                   B = 2000, X_sampler = NULL, B.X = B, Ab = NULL, seed = NULL){
  
  # Stopping if the HPD contains less than 100 functions
  if(B*(1-alpha) < 100 & !finite_X){
    stop("Choose a higher B so the HPD contains at least 100 functions.")
  }
  
  # If sigma2 was not informed, we assume that sigma2 ~ IG(alpha, beta)
  ## and plug the MAP for the GP, assuming (alpha, beta) converging to 0
  if(is.null(sigma2)) sigma2 <- sum((Y - mu.fun(X))^2)/(dim(X)[1] + 2)
  
  # Is the covariate space finite?
  if(finite_X){
    N <- dim(X_new)[1]
    chi <- qchisq(p = 1 - alpha, df = N)
    
    # Posterior parameters
    post.parm <- GP_pred(Y, X, X_new, sigma2, mu.fun, ker.fun)
    mu <- post.parm$f.mean
    Sig <- post.parm$f.covmat
    Sig_inv <- my_inv(Sig)
    
    # Setting X_new
    if(add_intercept) X_new <- cbind(1, X_new)
    if (vars != "all") {
      X_new <- X_new[, vars, drop = F]
    }
    
    # Decision for eps = 0
    min.val <- lm_ellipsoid(X_new, mu, Sig_inv)
    rej <- min.val[1] > chi
    dec <- ifelse(rej, "Reject H0", "Do not reject H0")
    
    # Decision for eps != 0
    if(rej & eps != 0){
      # If prob vector for X is not informed, assume discrete uniform
      if(is.null(probs.X)){
        probs.X <- rep(1/N, N)
        M <- diag(1, N) - X_new%*%my_inv(t(X_new)%*%X_new)%*%t(X_new)
      } else{
        Dprob <- diag(probs.X, N)
        M <- diag(1, N) - X_new%*%my_inv(t(X_new)%*%Dprob%*%X_new)%*%t(X_new)%*%Dprob
      }
      
      A_inv <- my_inv(M)%*%diag(1/probs.X, N)*eps^2
      dec <- ifelse(inter_ellipsoid(a = rep(0, N), A_inv = A_inv,
                                    b = mu, B_inv = Sig*chi),
                    "Do not reject Pg(H0)",
                    "Reject Pg(H0)")
    }
    
  } else{
    if(eps == 0){
      # Posterior parameters
      post.parm <- GP_pred(Y, X, X, sigma2, mu.fun, ker.fun)
      
      # Finding the HPD
      set.seed(seed)
      f.X <- mvrnorm(n = B, mu = post.parm$f.mean, Sigma = post.parm$f.covmat)
      hpd_cond <- apply(Y - t(f.X), 2, function(x) sum(x^2))
      hpd_check <- quantile(hpd_cond, prob = 1-alpha)
      
      if(add_intercept) X_new <- cbind(1, X)
      if (vars != "all") {
        X_new <- X_new[, vars, drop = F]
      }
      
      # Reaching a decision
      res <-lm(Y~X_new)$residuals
      dec <- ifelse(sum(res^2) >= hpd_check, "Reject H0", "Do not reject H0")
      
    } else{
      N <- length(Y)
      p <- dim(X)[2]
      
      # Was the sampler of X defined?
      if(is.null(X_sampler)){
        X.max <- X[max.col(t(X)), 1:p]
        X.min <- X[max.col(-t(X)), 1:p]
        X_sampler <- function(n){
          matrix(runif(n*p)*(X.max-X.min) + X.min, ncol = p, byrow = T)
          }
      }
      set.seed(seed)
      X_new <- rbind(X, X_sampler(n = B.X))
      
      # Finding the HPD
      post.parm <- GP_pred(Y, X, X_new, sigma2, mu.fun, ker.fun)
      f.X <- mvrnorm(n = B, mu = post.parm$f.mean,
                     Sigma = post.parm$f.covmat)
      hpd_cond <- apply(Y - t(f.X[,1:N]), 2, function(x) sum(x^2))
      hpd_check <- quantile(hpd_cond, prob = 1-alpha)
      
      # Functions on the HPD
      f.X <- f.X[hpd_cond <= hpd_check,-(1:N)]
      
      # Setting Ab
      if(add_intercept) X_new <- cbind(1, X_new[-(1:N),])
      if (vars != "all") {
        X_new <- X_new[, vars, drop = F]
      }
      p <- dim(X_new)[2] # Updating p
      if(is.null(Ab)){
        Ab <- matrix(NA, nrow = p, p)
        for(i in 1:p) Ab[i,i:p] <- colMeans(
          apply(X_new[,i:p, drop = F], 2, function(x) X_new[,i]*x)
        )
        Ab[lower.tri(Ab)] <- t(Ab[upper.tri(Ab)])
      }
      
      # Obtaining beta and checking the dissimilarity
      b <- 1
      cond <- F
      while(!cond & b <= dim(f.X)[1]){
        cond <- sqrt(mean((X_new%*%solve(Ab, colMeans(f.X[b,]*X_new)) - f.X[b,])^2)) <= eps
        b <- b + 1
      }
      dec <- ifelse(cond, "Do not reject Pg(H0)", "Reject Pg(H0)")
    }
    
  }
  
  dec
}
