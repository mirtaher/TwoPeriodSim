#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export


c2.m.array.lam <- function(lambda, S, i, r1, r2){

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

  res <- period.2.lam(lambda, S, y2[i, 1, r1, r2], y2[i, 2, r1, r2], Analytical = T)
  return(res)
}
