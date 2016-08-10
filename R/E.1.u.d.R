#' Expected Second Period Divorce Utility Given The First Period Info
#'
#' This function calculates the expected second period divorce utility given the information revealed at the first period for a given level of savings for an individual family.
#' Essentially, it calculates the second period optimal consumption for the given level of saving in the state of divorce for different
#' relizations of the shocks in the second period and then substiute them in the individual utility functions (u) or the joint utility function(U).
#' Finally, it takes the average to return the expectations.
#' @param S Saving
#' @param type Could be individual utility (u) or the joint utility of the couple (U) in a unitary framework
#' @param spouse If the type is individual utility (u), we should specify we mean husband or wife utility
#' @param i The marriage index
#' @param r First period repetition. It is not necessary to be greater than one for the first period. It is needed for taking expectaions, which is required in the second period
#' @param sigma_eta_h The husband's variance of transitory shock. If not specified the default is the baseline value specified in the param()
#' @param rho The contemporaneous correlation coefficient of the husband and wife income shocks. If not specified the default is the baseline value specified in the param()
#' @param phi The ratio of the wife's standard deviation of the transitory shock to that of the husband. If not specified the default is the baseline value specified in the param()
#' @export

E.1.u.d <-function(S, type = c("u", "U"), spouse = c("h", "w"), i, r, sigma_eta_h = param()$sigma_eta_h, rho = param()$rho, phi = param()$phi){

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
  y1 <- income$y1
  y2 <- income$y2

  # Cartesian product
  index <- expand.grid(1:N, 1:reps1)
  index2 <- expand.grid(1:N, 1:reps1, 1:reps2)

  if (!missing(spouse)){
    if (type == "u"){
      if (spouse == "h"){
        U.2.m.vec <- mapply(function(h) u(h),
                            matrix(period.2.d(S, i, r, 1:reps2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), ncol = 2)[,1])
      }
      if (spouse == "w"){
        U.2.m.vec <-  mapply(function(w) u(w),
                             matrix(period.2.d(S, i, r, 1:reps2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), ncol = 2)[,2])
      }

    }

    if (type == "U"){
      U.2.m.vec <- mapply(function(h,w) U(h, w),
                          matrix(period.2.d(S, i, r, 1:reps2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), ncol = 2)[,1],
                          matrix(period.2.d(S, i, r, 1:reps2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), ncol = 2)[,2])

    }
  } else {
    U.2.m.vec <- mapply(function(h,w) U(h, w),
                        matrix(period.2.d(S, i, r, 1:reps2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), ncol = 2)[,1],
                        matrix(period.2.d(S, i, r, 1:reps2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), ncol = 2)[,2])
  }


  res <- mean(U.2.m.vec, na.rm = TRUE)
  return(res)
}
