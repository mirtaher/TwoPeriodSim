#' The Second Period Problem Given Decision to Stay Married
#'
#' This function finds the optimal spouses' consumption levels in the second period given their decision to stay married.
#' We can use the optimization procedure or just use the analytical solution. The default is to use the analytical solution
#'  and it is necessary to be worked in some other parts of the code.
#'
#' @param S Saving
#' @param i The marriage index
#' @param r1 First period repetition. It is not necessary to be greater than one for the first period. It is needed for taking expectaions, which is required in the second period
#' @param r2 Second period repetition. It is needed for taking expectaions, which is required in the second period
#' @param Analytical If it is on it uses the analytical solution instead of using the optimization procedure.
#' @param sigma_eta_h The husband's variance of transitory shock. If not specified the default is the baseline value specified in the param()
#' @param rho The contemporaneous correlation coefficient of the husband and wife income shocks. If not specified the default is the baseline value specified in the param()
#' @param phi The ratio of the wife's standard deviation of the transitory shock to that of the husband. If not specified the default is the baseline value specified in the param()
#' @export

period.2.m <- function(S, i, r1, r2,  Analytical = TRUE, sigma_eta_h = param()$sigma_eta_h, rho = param()$rho, phi = param()$phi){

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

  income <- TwoPeriodSim::incomeProcess(sigma_eta_h = sigma_eta_h, Rho = rho, Phi = phi)
  if (!income$cov.det){
    res <- rep(NA, 2)
    return(res)
  } else {
    y1 <- income$y1
    y2 <- income$y2

    max.c <- 1/beta * S + y2[i, 1, r1, r2] + y2[i, 2, r1, r2]

    if (!Analytical) {
      f <- function(x) U(x, (max.c -x))
      res <- optimize(f, interval = c(0, max.c), maximum = T)
      return(c(res$maximum, (max.c -res$maximum)))
    } else {
      res.h <- 1/2 * max.c
      res.w <- max.c - res.h
      return(c(res.h, res.w))
    }
  }
}
