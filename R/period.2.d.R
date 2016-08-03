#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export


period.2.d <- function(S, i, r1, r2){

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

  sigma_eta_h <- par$sigma_eta_h
  rho <- par$rho
  phi <- par$phi

  income <- TwoPeriodSim::incomeProcess(sigma_eta_h = sigma_eta_h, Rho = rho, Phi = phi)
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
