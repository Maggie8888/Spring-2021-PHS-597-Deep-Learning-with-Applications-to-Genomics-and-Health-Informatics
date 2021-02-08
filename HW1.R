rm(list=ls())
dir<-"/Users/mengluliang/Dropbox/biostatisitcs_PSU/PHS597-2/Ex"
setwd(dir)
require(Rcpp)
require(RcppArmadillo)
sourceCpp("HW1.cpp")
fk_sum <- function(x, omega, h, x_eval = NULL, beta = c(.25, .25), nbin = NULL, type = "ksum"){
  n <- length(x)
  
  # center x
  mn <- mean(x)
  x <- x - mn
  
  if(is.null(nbin)){
    
    o <- order(x)
    xo <- x[o]
    omo <- omega[o]
    if(is.null(x_eval)){
     
      if(type == "ksum") ksum(xo, omo, x, h, beta, match(x, xo))
      else if(type=="dksum") dksum(xo, omo, x, h, beta, match(x, xo))
      else if(type=="both") kndksum(xo, omo, x, h, beta, match(x, xo))
    }
    else{
      x_eval <- x_eval - mn
      o_eval <- order(x_eval)
      x_eo <- x_eval[o_eval]
      if(type == "ksum") ksum(xo, omo, x_eo, h, beta)[match(x_eval, x_eo)]
      else if(type == "dksum") dksum(xo, omo, x_eo, h, beta)[match(x_eval, x_eo)]
      else if(type == "both") kndksum(xo, omo,x_eo, h, beta)[match(x_eval, x_eo),]
    }
  }
  else{ 
    xo <- seq(min(x) - 1e-5, max(x) + 1e-5, length = nbin)
    omo <- sm_bin_wts(x, omega, nbin, xo[1], xo[nbin])
    if(is.null(x_eval)){
      if(type == "ksum") ksum(xo, omo, x, h, beta, cbin_alloc(x, nbin, xo[1], xo[nbin]))
      else if(type == "dksum") dksum(xo, omo, x, h, beta, cbin_alloc(x, nbin, xo[1], xo[nbin]))
      else if(type == "both") kndksum(xo, omo, x, h, beta, cbin_alloc(x, nbin, xo[1], xo[nbin]))
    }
    else{
     
      x_eval <- x_eval - mn
      if(type == "ksum") ksum(xo, omo, x_eval, h, beta, cbin_alloc(x_eval, nbin, xo[1], xo[nbin]))
      else if(type == "dksum") dksum(xo, omo, x_eval, h, beta, cbin_alloc(x_eval, nbin, xo[1], xo[nbin]))
      else if(type == "both") kndksum(xo, omo, x_eval, h, beta, cbin_alloc(x_eval, nbin, xo[1], xo[nbin]))
    }
  }
}

phi_ppr <- function(w, X, r, h, beta){
  n <- nrow(X)
  p <- X %*% w / sqrt(sum(w^2))
  Sr <- fk_sum(p, r, h, beta = beta) - beta[1] * r
  S1 <- fk_sum(p, rep(1, n), h, beta = beta) - beta[1]
  S1[S1 < 1e-20] <- 1e-20
  r_hat <- Sr / S1
  sum((r - r_hat)^2)
}


dphi_ppr <- function(w, X, r, h, beta){
  n <- nrow(X)
  nw <- sqrt(sum(w^2))
  p <- X %*% w / nw
  S1 <- fk_sum(p, rep(1, n), h, beta = beta, type = "both")
  S1[, 1] <- S1[, 1] - beta[1]
  S1[S1[, 1] < 1e-20, 1] <- 1e-20
  Sr <- fk_sum(p, r, h, beta = beta, type = "both")
  Sr[, 1] <- Sr[, 1] - beta[1] * r
  r_hat <- Sr[, 1] / S1[, 1]
  T1 <- fk_sum(p, r_hat * (r_hat - r) / S1[, 1], h,
               beta = beta, type = "dksum")
  T2 <- r * fk_sum(p, (r_hat - r) / S1[, 1], h,
                   beta = beta, type = "dksum")
  T3 <- (r_hat - r) / S1[, 1] * (r_hat * S1[, 2] - Sr[, 2])
  dphi_dp <- (T1 - T2 + T3) * 2 / h
  dp_dw <- (X / nw - p %*% t(w) / nw^2)
  c(t(dphi_dp) %*% dp_dw)
}

set.seed(1234567)
n_dat <- 1000
n_dim <- 10
X <- matrix(rnorm(n_dat * n_dim), n_dat, n_dim) %*%
  matrix(2 * runif(n_dim^2) - 1, n_dim, n_dim)
wtrue1 <- rnorm(n_dim)
wtrue2 <- rnorm(n_dim)
y <- (X %*% wtrue1 > 1) * (X %*% wtrue1 - 1) + tanh(X %*% wtrue2 / 2) *
  (X %*% wtrue1) + (X %*% (wtrue1 - wtrue2) / 5)^2 + rnorm(n_dat)
w <- rnorm(n_dim)
h <- runif(1)
beta <- c(0.25, 0.25)
Eh <- diag(n_dim) * 1e-5
dphi_approx <- apply(Eh, 1, function(eh) (phi_ppr(w + eh, X, y, h, beta)
- phi_ppr(w - eh, X, y, h, beta)) / 2e-5)
dphi <- dphi_ppr(w, X, y, h, beta)
max(abs(dphi / dphi_approx - 1))
max(abs(dphi - dphi_approx))


ppr_nw <- function(X, y, w = NULL){
  n <- nrow(X)
  d <- ncol(X)
  if(is.null(w)) w <- solve(t(X) %*% X + .01 * diag(d)) %*% t(X) %*% y
  h <- sqrt(eigen(cov(X))$values[1]) / n^.2
  w <- optim(w, phi_ppr, dphi_ppr, X, y, h, c(.25, .25),
             method = "L-BFGS-B")$par
  w <- w / sqrt(sum(w^2))
  loo_sse <- function(h) phi_ppr(w, X, y, h, c(.25, .25))
  h <- optimise(loo_sse, c(h/50, h))$minimum
  list(w = w, h = h)
}


w_opt <- ppr_nw(X, y, w = w)
p <- X %*% w_opt$w
S1 <- fk_sum(p, rep(1, n_dat), w_opt$h)
Sy <- fk_sum(p, y, w_opt$h)
fitted <- Sy / S1
par(mfrow = c(1, 3))
plot(X %*% w, y, main = "Response against initial projection",
         xlab = "optimal projection", ylab = "y")
plot(X %*% w_opt$w, y, main = "Response against optimal projection",
       xlab = "optimal projection", ylab = "y")
points(X %*% w_opt$w, fitted, col = 2)
plot(fitted, y, main = "Response against fitted values",xlab = "fitted", ylab = "y")
abline(0, 1, lwd = 2, col = 2)

n_rep <- 100
t_stats <- numeric(n_rep)
t_nw <- numeric(n_rep)
R2_stats <- numeric(n_rep)
R2_nw <- numeric(n_rep)

for(rep in 1:n_rep){
  set.seed(rep)
  X <- matrix(rnorm(n_dat * n_dim), n_dat, n_dim) %*%
    matrix(2 * runif(n_dim * n_dim) - 1, n_dim, n_dim)
  wtrue1 <- rnorm(n_dim)
  wtrue2 <- rnorm(n_dim)
  y <- (X %*% wtrue1 > 1) * (X %*% wtrue1 - 1) + tanh(X %*% wtrue2 / 2) *
    (X %*% wtrue1) + (X %*% (wtrue1 - wtrue2) / 5)^2 + rnorm(n_dat)
  t_stats[rep] <- system.time(model <- ppr(X[1:(n_dat / 2),],
                                           y[1:(n_dat / 2)], nterms = 1))[1]
  yhat <- predict(model, X[(n_dat / 2 + 1):n_dat,])
  R2_stats[rep] <- 1 - mean((yhat - y[(n_dat / 2 + 1):n_dat])^2) /
    var(y[(n_dat / 2 + 1):n_dat])
  t_nw[rep] <- system.time(model <- ppr_nw(X[1:(n_dat / 2),],
                                           y[1:(n_dat / 2)]))[1]
  p <- X[1:(n_dat / 2),] %*% model$w
  ptest <- X[(n_dat / 2 + 1):n_dat,] %*% model$w
  S1 <- fk_sum(p, rep(1, n_dat / 2), model$h, x_eval = ptest)
  Sy <- fk_sum(p, y[1:(n_dat / 2)], model$h, x_eval = ptest)
  yhat <- Sy / S1
  R2_nw[rep] <- 1 - mean((yhat - y[(n_dat / 2 + 1):n_dat])^2) /
    var(y[(n_dat / 2 + 1):n_dat])
}

#### comparison with PPR 
colMeans(cbind(t_stats, t_nw))
#t_stats    t_nw 
#0.00354 0.03574
colMeans(cbind(R2_stats, R2_nw))
#R2_stats     R2_nw 
#0.7120016 0.6697174 

