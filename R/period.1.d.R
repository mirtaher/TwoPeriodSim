#' The First Period Problem Given Decision to Divorce in the Second Period
#'
#' This function find the optimal spouses' consumption levels and saving given the fact that they know that they are going to be divorced in the following period.
#' The decision to divorce shows itself in the expected next period utility in the objective function; we use E.1.u.d as opposed to E.1.u.m
#' @param i The marriage index
#' @param r First period repetition. It is not necessary to be greater than one for the first period. It is needed for taking expectaions, which is required in the second period
#' @param sigma_eta_h The husband's variance of transitory shock. If not specified the default is the baseline value specified in the param()
#' @param rho The contemporaneous correlation coefficient of the husband and wife income shocks. If not specified the default is the baseline value specified in the param()
#' @param phi The ratio of the wife's standard deviation of the transitory shock to that of the husband. If not specified the default is the baseline value specified in the param()
#' @export


period.1.d <- function(i, r, sigma_eta_h = param()$sigma_eta_h, rho = param()$rho, phi = param()$phi){
  library(nloptr)

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

  # Cartesian product
  index <- expand.grid(1:N, 1:reps1)
  index2 <- expand.grid(1:N, 1:reps1, 1:reps2)

  income <- TwoPeriodSim::incomeProcess(sigma_eta_h = sigma_eta_h, Rho = rho, Phi = phi)
  if (!income$cov.det) {
    res <- rep(NA, 3)
    return(res)
  }else {
    y1 <- income$y1
    y2 <- income$y2

    miss <- is.na(y1[i,1,r]) | is.na(y1[i,2,r])
    if (miss){
      res <- rep(NA, 3)
      return(res)
    } else {
      eval_f_1_d <- function(c, i, r){
        return(list("objective" = - U(c[1], c[2])- beta * E.1.u.d(c[3], type = "U", i = i, r = r, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi),
                    "gradient" = c(-wt.h * u.grad(c[1]),
                                   -wt.w * u.grad(c[2]),
                                   (-beta) * E.1.u.d.der(c[3], type = "U", i = i, r = r, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi) )))
      }

      eval_g_eq_1 <- function(c, i, r){
        constr <- c( c[1] + c[2] + c[3] - y1[i,1,r] - y1[i,2,r])
        grad <-  c(1, 1, 1)
        return( list( "constraints"=constr, "jacobian"=grad ) )
      }

      local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                          "xtol_rel"  = 1.0e-7)
      opts <- list( "algorithm" = "NLOPT_LD_SLSQP",
                    "xtol_rel"  = 1.0e-7,
                    "maxeval"   = 10000,
                    "local_opts" = local_opts )

      ini <- borrow.const(i, r, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi)
      S.lower <- ini$lower
      S.upper <- ini$upper

      c.upper <- y1[i,1,r] + y1[i,2,r] - S.lower
      c.lower <- 0

      ub <- c(c.upper, c.upper, S.upper)
      lb <- c(c.lower, c.lower, S.lower)

      S0 <- (S.lower + S.upper)/2
      c0 <- (y1[i,1,r] + y1[i,2,r] - S0)/2
      x0 <- c(c0, c0, S0)

      res <- nloptr(x0 = x0, eval_f = eval_f_1_d, lb = lb, ub = ub,
                    eval_g_eq = eval_g_eq_1, opts = opts, i = i, r = r)
      return(res$solution)
    }
  }
}
