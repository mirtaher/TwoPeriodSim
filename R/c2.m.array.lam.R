#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export


c2.m.array.lam <- function(lambda, S){

  # just for testing the git

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
  chi <- par$chi

  c2 <- mapply(function(yh, yw) period.2.lam(lambda, S, yh, yw, Analytical = T), y2[,1,,], y2[,2,,])
  c2.array <- array(c2, dim = c(2,N,reps1,reps2) )
  res <- aperm(c2.array, c(2,1,3,4))
  return(res)
}
