#' Participation Constraints
#'
#' This function checks whether the participation constraints of a specific family satisfy for a given value of saving and sharing rule
#' @param lam Sharing rule
#' @param S Saving
#' @param i The marriage index
#' @param r First period repetition. It is not necessary to be greater than one for the first period. It is needed for taking expectaions, which is required in the second period
#' @export


participate.lam <- function(lam, S, i, r){

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

  sol.d <- period.1.d(i, r)
  s.d <- sol.d[3]

  con.h <- ifelse(E.1.u.d(s.d, type = "u", spouse = "h", i, r) <= E.1.u.m.lam(lam, S, type = "u", spouse = "h", i, r), 1, 0)
  con.w <- ifelse(E.1.u.d(s.d, type = "u", spouse = "w", i, r) <= E.1.u.m.lam(lam, S, type = "u", spouse = "w", i, r), 1, 0)
  con <- con.h & con.w
  return(con)

  }
