#' The Second Period Problem Given Decision to Divorce
#'
#' This function find the optimal spouses' consumption levels in the second period given their decision to get divorced.
#' There is no optimization involved here; they consume all their income in the second period.
#'
#' @param i The marriage index
#' @param r1 First period repetition. It is not necessary to be greater than one for the first period. It is needed for taking expectaions, which is required in the second period
#' @param r2 Second period repetition. It is needed for taking expectaions, which is required in the second period
#' @param sigma_eta_h The husband's variance of transitory shock. If not specified the default is the baseline value specified in the param()
#' @param rho The contemporaneous correlation coefficient of the husband and wife income shocks. If not specified the default is the baseline value specified in the param()
#' @param phi The ratio of the wife's standard deviation of the transitory shock to that of the husband. If not specified the default is the baseline value specified in the param()
#' @export

period.2.d <- function(S, i, r1, r2, sigma_eta_h = param()$sigma_eta_h, rho = param()$rho, phi = param()$phi){

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

  income <- TwoPeriodSim::incomeProcess(sigma_eta_h = sigma_eta_h, Rho = rho, Phi = phi)
  if (!income$cov.det) {
    res <- rep(NA, 2)
    return(res)
  } else {
    y1 <- income$y1
    y2 <- income$y2

    miss <- is.na(y2[i, 1, r1, r2]) | is.na(y2[i, 2, r1, r2])
    if (miss){
      res.h <- NA
      res.w <- NA
    } else{
      res.h <- 1/beta * delta * S + y2[i, 1, r1, r2]
      res.w <- 1/beta * (1 - delta) * S + y2[i, 2, r1, r2]
    }
    return(c(res.h, res.w))
  }
}
