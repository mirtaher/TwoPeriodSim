#' Joint Utility in the Collective Framework
#'
#' This function returns the joint utility function for a given set of sharing rule, husband, and wife consumptions.
#'
#' @param lambda Sharing rule
#' @param ch Husband consumption
#' @param cw Wife consumption
#' @export

U2.lam <- function(lambda, ch, cw){

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

  res <- lambda * u(ch) + (1 - lambda) * u(cw)
  return(res)


}
