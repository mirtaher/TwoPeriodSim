#' Borrowing Constraints
#'
#' This function finds the max and min amount of savings given a specific income shocks
#' parameters and individual family (i) and specific repetition (r). If income sgocks pars are not specified
#' by defualt the values in the param function would be adopted.
#'  First, it checks whether the incomes are NAs, whuch in this case just returns two NAs.
#'  Otherwise, the minimum savings (max borrowing) in the first period is the minimum realized
#'  second period income in simulation in the state of marriage and divorce for each of spouses.
#'  The minimum of these three anounts, i.e. min second period income in the state of divorce for husband, wife, and
#'  family if tehy decide to stay married, would constitute as the min savings. In addition, an starvation consumption in the second
#'  period (which specifies by tol variable) sets aside.
#'  The maximum saving is the whole total first period income after setting aside the starvation consumption in the first period.
#' @export


borrow.const <- function(i, r, sigma_eta_h = param()$sigma_eta_h, rho = param()$rho, phi = param()$phi){


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

  income <- TwoPeriodSim::incomeProcess(sigma_eta_h = sigma_eta_h, Rho = rho, Phi = phi)
  y1 <- income$y1
  y2 <- income$y2

  miss <- is.na(y1[i,1,r]) | is.na(y1[i,2,r])
  if (miss) {
    s.min <- NA
    s.max <- NA
  } else {
    equal <- (y1[i,1,r] + y1[i,2,r])/2
    y2.h.min <- min(y2[i,1,r,])
    y2.w.min <- min(y2[i,2,r,])
    y2.min <- min(y2[i,1,r,] + y2[i,2,r,])
    tol.cons <- 0.01

    b.mar <- (-beta) * y2.min  * (1 - tol.cons)
    b.h <- (-beta)/delta * y2.h.min * (1 - tol.cons)
    b.w <- (-beta)/(1 - delta) * y2.w.min * (1 - tol.cons)

    s.min <- max(b.mar, b.h, b.w)

    if (equal > 0){
      s.max <- (y1[i,1,r] + y1[i,2,r]) * (1 - tol.cons)
    } else {
      s.max <- (y1[i,1,r] + y1[i,2,r]) * (1 + tol.cons)
      if (s.max < s.min){
        s.min <- NA
        s.max <- NA
      }
    }
  }

  return(list("lower" = s.min, "upper" = s.max))

}





