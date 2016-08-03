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
    y2.min <- apply(y2, MARGIN = c(1,2,3), FUN = min)
    y.life.min <- y1 + beta * y2.min
    ind <- which(y.life.min < 0, arr.ind = TRUE)
    y1[ind] <- NA
    cases <- length(ind[,1])
    for (ii in 1:cases){
      y2[ind[ii,1], ind[ii,2], ind[ii,3], ] <- NA
    }
    y1.tot <- y1[, 1,] + y1[, 2, ]
    ind.1.neg <- which(y1.tot < 0, arr.ind = TRUE)

    #y1.expand <- array(rep(y1,reps2), dim = c(N,2,reps1, reps2))
    #y.life <- y1.expand + beta * y2

    res <- list("y1" = y1, "y2" = y2, "Omega" = Omega, "epsilon1" = epsilon1, "epsilon2" = epsilon2, "ind.1.neg" = ind.1.neg)
  }

  return(res)

}




























