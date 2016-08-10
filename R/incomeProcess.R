#' Simulate the Income Process
#'
#' This function simulates the income process for two spaouses for a given set of parameters of income shocks process. First, it checks the
#' covariance matrix determinant for the provided shock parameters. If the cov.det is negative, it exits and assigns FALSE to cov.det outcome.
#' Otherwise, it simulates the shocks. It is possible that the lifetime earnings of some individuals turn out to be negative, which is problematic computationally
#' and theoretecally. I have dropped those observations (in fact dropping the families for whom at least one spouse has negative lifetime earnings)
#' . However, this in fact messes the randomization. I assigend NA for both couples in both periods if the lifetime earninga of one of the
#' spouses is negative. In addition, we have identified those fammilies for who the total first period income is negative despite the fact that their lifetime earning is not.
#'
#' @param sigma_eta_h The husband's variance of transitory shock. If not specified the default is the baseline value specified in the param()
#' @param Rho The contemporaneous correlation coefficient of the husband and wife income shocks. If not specified the default is the baseline value specified in the param()
#' @param Phi The ratio of the wife's standard deviation of the transitory shock to that of the husband. If not specified the default is the baseline value specified in the param()
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
    cov.det <- FALSE
  } else {
    cov.det <- TRUE
  }
  if (cov.det){
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
    if (cases > 0){
      for (ii in 1:cases){
        y2[ind[ii,1], ind[ii,2], ind[ii,3], ] <- NA
      }
    }

    y1.tot <- y1[, 1,] + y1[, 2, ]
    ind.1.neg <- which(y1.tot < 0, arr.ind = TRUE)

    res <- list("y1" = y1, "y2" = y2, "Omega" = Omega, "epsilon1" = epsilon1,
                "epsilon2" = epsilon2, "ind" = ind, "ind.1.neg" = ind.1.neg, "cov.det" = cov.det)
  }

  return(res)

}




























