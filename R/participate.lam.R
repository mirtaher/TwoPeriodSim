#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

participate.lam <- function(lam, S, i, r){
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
  tau.d <- par$tau.d

  wt.h <- par$wt.h
  wt.w <- par$wt.w
  u <- par$u
  U <- par$U
  u.grad <- par$u.grad
  chi <- par$chi

  con.h <- ifelse(E.1.u.d(S, type = "u", spouse = "h", i, r) <= E.1.u.m.lam(lam, S, type = "u", spouse = "h", i, r), 1, 0)
  con.w <- ifelse(E.1.u.d(S, type = "u", spouse = "w", i, r) <= E.1.u.m.lam(lam, S, type = "u", spouse = "w", i, r), 1, 0)
  con <- con.h & con.w
  return(con)

  }
