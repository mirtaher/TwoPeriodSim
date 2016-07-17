#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export


period.2.d <- function(S, i, r1, r2){

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


  res.h <- 1/beta * delta * S + y2[i, 1, r1, r2]
  res.w <- 1/beta * (1 - delta) * S + y2[i, 2, r1, r2]
  return(c(res.h, res.w))
}
