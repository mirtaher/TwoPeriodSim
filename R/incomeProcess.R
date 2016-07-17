#' Simulating the income process
#'
#' This function simulate the income process for two spaouses
#' @param Phi Rho
#' @keywords cats
#' @export


incomeProcess <- function(sigma_eta_h, Rho, Phi){
  suppressWarnings(library(MASS))

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


  # Transitory shocks parameters
  sigma_eta_w = sigma_eta_h * Phi^2
  cov_eta = sqrt(sigma_eta_h * sigma_eta_w)  * Rho
  Omega <- (1 + 2 * Rho * Phi + Phi^2)/4
  sigma_eta = matrix(c(sigma_eta_h, cov_eta, cov_eta, sigma_eta_w),2,2)
  if (det(sigma_eta) <= 0){
    warning("covariance matrix is not S.P.D")
    res <- NA
  }
  else {
    # shocks
    set.seed(seed); eta1 <- replicate(n = reps1, mvrnorm(n = N, mu = rep(0,2), Sigma = sigma_eta))
    eta1.expand <- array(rep(eta1, reps2), dim = c(N,2,reps1, reps2))
    set.seed(seed *100); eta2 <- replicate(n = reps2, replicate(n = reps1, mvrnorm(n = N, mu = rep(0,2), Sigma = sigma_eta)))
    epsilon1 <- eta1
    epsilon2 <- mu * eta1.expand + eta2
    ybar1 <- array(rep(theta* ybar, 2 * N * reps1), dim = c(N,2,reps1))
    y1 <- ybar1 + epsilon1
    ybar2 <- array(rep(1/beta * (1-theta) * ybar, 2 * N * reps1 *reps2), dim = c(N,2,reps1,reps2))
    y2 <- ybar2 + epsilon2
    y1.expand <- array(rep(y1,reps2), dim = c(N,2,reps1, reps2))
    y.life <- y1.expand + beta * y2
    if (any(y.life <0)){
      res <- NA
      warning ("life time earnings are negative")
    }
    else{
      res <- list("y1" = y1, "y2" = y2, "Omega" = Omega, "epsilon1" = epsilon1, "epsilon2" = epsilon2)
    }
  }

  return(res)

}




























