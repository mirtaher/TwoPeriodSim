#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

period.1.m <- function(i, r){
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

  sigma_eta_h <- par$sigma_eta_h
  rho <- par$rho
  phi <- par$phi

  income <- TwoPeriodSim::incomeProcess(sigma_eta_h = sigma_eta_h, Rho = rho, Phi = phi)
  y1 <- income$y1
  y2 <- income$y2

  # Cartesian product
  index <- expand.grid(1:N, 1:reps1)
  index2 <- expand.grid(1:N, 1:reps1, 1:reps2)

  eval_f_1_m <- function(c, i, r){
    return(list("objective" = - U(c[1], c[2])- beta * E.1.u.m(c[3], type = "U", i = i, r = r),
                "gradient" = c(-wt.h * u.grad(c[1]), -wt.w * u.grad(c[2]), (-beta) * E.1.u.m.der(c[3], type = "U", i = i, r = r) )))
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

  equal <- (y1[i,1,r] + y1[i,2,r])/2
  if (equal < 0) {
    y1.expand <- array(rep(y1,reps2), dim = c(N,2,reps1, reps2))
    equal.neg.expand <- (y1.expand[i,1,r,] + y1.expand[i,2,r,] + y2[i,1,r,] + y2[i,2,r,])/4
    equal.neg <- mean(equal.neg.expand)
    c0 <- c(equal.neg, equal.neg, 0)
  }
  else{
    c0 <- c(equal, equal, 0)
  }

  ini <- borrow.const(i, r)
  S.lower <- ini$lower
  S.upper <- ini$upper

  ub <- c(equal*2, equal*2, S.upper)
  lb <- c(0, 0, S.lower)

  res <- nloptr(x0 = c0, eval_f = eval_f_1_m, lb = lb, ub = ub,
                eval_g_eq = eval_g_eq_1, opts = opts, i = i, r = r)
  return(res$solution)
}
