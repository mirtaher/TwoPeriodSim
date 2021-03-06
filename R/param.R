#' Parameters Distribution
#'
#' This function distributes the parameters values
#' @export

param <- function(){

  # size pars
  N <- 5 #number of couples
  reps1 <- 2
  reps2 <- 5 # simulation repetitions for
  seed <- 123

  # income process parameters
  theta <-  0.3
  mu <-  0.9
  beta <-  1/(1 + 0.04)
  ybar <-  1 # the same for husband and wife

  # income shocks parameter baseline
  sigma_eta_h <- 0.01
  rho <- 0.35
  phi <- 1.92

  # Divorce parameters
  delta <- 0.5 # divorce dividing rule upon divorce
  tau.d <- 0.02 # the foregone utility of singlehood by divorce
  # note that when using quadratic utility or any utility function with positive
  # values use negative values for tau.d

  # Quadratic Utility
    b <- 0.5
    #u <- function(c) c - 1/2 * b * c^2
    #u.grad <-  function(c) 1 - b * c

  # CRRA Utility
    chi <- 1.5
    u <- function(c) 1/(1 - chi) * (c^(1 - chi))
    u.grad <- function(c) c^(- chi)

  # Common Utility parameters
  wt.h <- 0.5
  wt.w <- 0.5
  U <- function(ch, cw) wt.h * u(ch) + wt.w * u(cw)

  res <- list("N" = N , "reps1" = reps1 , "reps2" = reps2, "seed" = seed, "theta" = theta,
          "mu" = mu, "beta" = beta, "ybar" = ybar , "delta" = delta, "tau.d" = tau.d,
          "wt.h" = wt.h, "wt.w" = wt.w, "u" = u, "U" = U, "u.grad" = u.grad, "chi" = chi,
          "sigma_eta_h" = sigma_eta_h, "phi" = phi, "rho" = rho)


}
