#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

global.initial <- function(i, r){
  library(nloptr)

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

  # Cartesian product
  index <- expand.grid(1:N, 1:reps1)
  index2 <- expand.grid(1:N, 1:reps1, 1:reps2)

  start.points <- feasible.region(i, r)
  values <- mapply(function(j) obj.lam( c.h = start.points[j, 1],
                              c.w <- start.points[j, 2],
                              S <- start.points[j, 3],
                              lam <- start.points[j, 4],  i, r ), 1: nrow(start.points))

  sol <- start.points[which.min(values), ]
  res <- list("S0" = sol[3], "lam0" = sol[4], "c0" = sol[1], "status" = "Stay Married with New Terms")
  return(res)

}
