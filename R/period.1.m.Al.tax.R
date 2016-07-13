#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

period.1.m.Al.tax <- function(i, r, pc = FALSE){
  library(alabama)
  
  # Distributing parameters
  par <- param()
  N <- par$N
  reps1 <- par$reps1
  reps2 <- par$reps2
  seed <- par$seed
  
  theta <- par$theta
  mu <- par$mu
  beta <- par$beta
  ybar <- par$ybar
  delta <- par$delta
  
  wt.h <- par$wt.h
  wt.w <- par$wt.w
  u <- par$u
  U <- par$U
  u.grad <- par$u.grad
  
  tau <- 0.05
  
  # Cartesian product
  index <- expand.grid(1:N, 1:reps1)
  index2 <- expand.grid(1:N, 1:reps1, 1:reps2)
  
  fn <- function(c, i, r) - U(c[1], c[2])- beta * E.1.u.m(c[3], type = "U")[i,r]
  gr <- function(c, i, r) c(-wt.h * u.grad(c[1]), -wt.w * u.grad(c[2]), (-beta) * E.1.u.m.der(c[3], type = "U")[i,r] )
  
  heq <- function(c, i, r) c(c[1] + c[2] + c[3] - y1[i,1,r] - y1[i,2,r])
  heq.jac <- function(c, i, r) c(1, 1, 1)
  
  # without participation constraints
  hin <- function(c, i, r) c( c[1], c[2])
  hin.jac <- function(c, i, r) diag(2)
  
  if (pc){
    hin <- function(c, i, r) c(E.1.u.m(c[3], type = "u", spouse = "h")[i,r] - (1 - tau) * E.1.u.d(c[3], type = "u", spouse = "h")[i,r],
                               E.1.u.m(c[3], type = "u", spouse = "w")[i,r] - (1 - tau) * E.1.u.d(c[3], type = "u", spouse = "w")[i,r], 
                               c[1],
                               c[2])
    hin.jac <- function(c, i, r) matrix(c(0, 0, 1, 0,
                                          0, 0, 0, 1, 
                                          E.1.u.m.der(c[3], type = "u", spouse = "h")[i,r] - (1 - tau) * E.1.u.d.der(c[3], type = "u", spouse = "h")[i,r],
                                          E.1.u.m.der(c[3], type = "u", spouse = "w")[i,r] - (1 - tau) * E.1.u.d.der(c[3], type = "u", spouse = "w")[i,r],
                                          0,
                                          0), ncol = 3)
  }
  
  
  equal <- (y1[i,1,r] + y1[i,2,r])/2
  if (equal < 0) {
    y1.expand <- array(rep(y1,reps2), dim = c(N,2,reps1, reps2))
    equal.neg.expand <- (y1.expand[i,1,r,] + y1.expand[i,2,r,] + y2[i,1,r,] + y2[i,2,r,])/4
    equal.neg <- mean(equal.neg.expand)
    c0 <- c(equal.neg, equal.neg, 0)
  } else{
    c0 <- c(equal, equal, 0)
  }
  
  res <- alabama::auglag(par = c0, fn = fn, gr = gr, heq = heq, heq.jac = heq.jac, hin = hin, hin.jac = hin.jac, i = i, r = r)
  
  return(list("sol" = res$par, "Lag" = res$lambda))
}