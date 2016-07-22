#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

period.1.m.lam <- function(i, r){
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
  chi <- par$chi

  # Cartesian product
  index <- expand.grid(1:N, 1:reps1)
  index2 <- expand.grid(1:N, 1:reps1, 1:reps2)

  ini <- feasible.region(i, r)

  if (ini$status == "Stay Married with New Terms") {
    eval_f_1_m <- function(c, i, r){
      return(list("objective" = - U(c[1], c[2])- beta * E.1.u.m.lam(c[4], c[3], type = "U", i = i, r = r),
                  "gradient" = c(-wt.h * u.grad(c[1]),
                                 -wt.w * u.grad(c[2]),
                                 (-beta) * E.1.u.m.der.lam(c[4], c[3], var =  "S", type = "U", i = i, r = r),
                                 (-beta) * E.1.u.m.der.lam(c[4], c[3], var = "lambda", type = "U", i = i, r = r))))
    }

    eval_g_eq_1 <- function(c, i, r){
      constr <- c( c[1] + c[2] + c[3] - y1[i,1,r] - y1[i,2,r])
      grad <-  c(1, 1, 1, 0)
      return( list( "constraints"=constr, "jacobian"=grad ) )
    }


    eval_g_ineq <- function(c, i, r){
      c(E.1.u.d(c[3], type = "u", spouse = "h", i, r) - E.1.u.m.lam(c[4], c[3], type = "u", spouse = "h", i, r),
        E.1.u.d(c[3], type = "u", spouse = "w", i, r) - E.1.u.m.lam(c[4], c[3], type = "u", spouse = "w", i, r))
    }

    eval_jac_g_ineq <- function(c, i, r){
      matrix(c(0,
               0,
               E.1.u.d.der(c[3], type = "u", spouse = "h", i, r) - E.1.u.m.der.lam(c[4], c[3], var = "S", type = "u", spouse = "h", i, r),
               - E.1.u.m.der.lam(c[4], c[3], var = "lambda", type = "u", spouse = "h", i, r),
               0,
               0,
               E.1.u.d.der(c[3], type = "u", spouse = "w", i, r) - E.1.u.m.der.lam(c[4], c[3], var = "S", type = "u", spouse = "w", i, r),
               - E.1.u.m.der.lam(c[4], c[3], var = "lambda", type = "u", spouse = "w", i, r)) ,
             ncol = 2)
    }

    local_opts <- list( "algorithm" = "NLOPT_LN_COBYLA",
                        "xtol_rel"  = 1.0e-3)
    opts <- list( "algorithm" = "NLOPT_GN_ISRES",
                  "xtol_rel"  = 1.0e-3,
                  "maxeval"   = 100000,
                  "local_opts" = local_opts )

    # initial values
    c0 <- ini$c0
    S0 <- ini$S0
    lam0 <- ini$lam0
    x0 <- c(c0, c0, S0, lam0)

    # upper and lower bounds
    BC <- borrow.const(i, r)
    S.upper <- BC$upper
    S.lower <- BC$lower

    lam.upper <- 1
    lam.lower <- 0

    c.upper <- y1[i,1,r] + y1[i,2,r] - S.lower
    c.lower <- 0

    ub <- c(c.upper, c.upper, S.upper, lam.upper)
    lb <- c(c.lower, c.lower, S.lower, lam.lower)

    res <- nloptr(x0 = x0, eval_f = eval_f_1_m, lb = lb, ub = ub,
                  eval_g_eq = eval_g_eq_1, eval_g_ineq = eval_g_ineq,
                  eval_jac_g_ineq = eval_jac_g_ineq, opts = opts, i = i, r = r)
    #return(res)
    return(list("sol" = res$solution, "status" = "Stay Married with New Terms"))
  }

  if (ini$status == "Stay Married with Old Terms"){
    sol <- c(ini$c.h.uncon, ini$c.w.uncon, ini$s.uncon, 0.5)
    return(list("sol" = sol, "status" = "Stay Married with Old Terms"))
  }

  if (ini$status == "Divorce"){
    sol <- period.1.d(i, r)
    sol <- c(sol, NA)
    return("sol" = sol, "status" = "Divorce")

  }


}
