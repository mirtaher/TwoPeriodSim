#' Simulating the income process
#'
#' This function simulate the income process for two spaouses
#' @param Phi Rho
#' @keywords cats
#' @export


E.1.m.ana <- function(S, i, r){

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
  b <- par$b

  wt.h <- par$wt.h
  wt.w <- par$wt.w
  u <- par$u
  U <- par$U
  u.grad <- par$u.grad


  B <- 1/2 * S + (1-theta) * ybar
  term.1 <- B^2/beta^2
  term.2 <- 1/4 * (mu^2 * (epsilon1[i,1,r] + epsilon1[i,2,r])^2 + sigma_eta_h + sigma_eta_w + 2 * cov_eta)
  term.3 <- B/beta * mu * (epsilon1[i,1,r] + epsilon1[i,2,r])
  res2 <- (term.1 + term.2 + term.3) * b/2
  res1 <- B/beta + mu/2 * (epsilon1[i,1,r] + epsilon1[i,2,r])
  res <- res1 - res2
  return(res)
}
