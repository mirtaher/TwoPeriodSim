#' Joint Objective Function Given Period One Info
#'
#' This function returns the (negative) value of joint objective function for a set of given first period varaibles as of period one information set.
#' @param c.h First period husband consumption
#' @param c.w First period wife consumption
#' @param S Saving
#' @param lam Sharing rule for the second period
#' @param i The marriage index
#' @param r First period repetition. It is not necessary to be greater than one for the first period. It is needed for taking expectaions, which is required in the second period

#' @export

obj.lam <- function(c.h, c.w, S, lam, i, r){
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

  res <- - U(c.h, c.w) - beta * E.1.u.m.lam(lam, S, type = "U", i = i, r = r)
  return(res)

}
