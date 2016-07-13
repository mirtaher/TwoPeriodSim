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
  
  eval_f_1_m <- function(c, i, r){
    return(list("objective" = - U(c[1], c[2])- beta * E.1.u.m.lam(c[4], c[3], type = "U")[i,r], 
                "gradient" = c(-wt.h * u.grad(c[1]), 
                               -wt.w * u.grad(c[2]),
                               (-beta) * E.1.u.m.der.lam(c[4], c[3], var =  "S", type = "U")[i,r], 
                               (-beta) * E.1.u.m.der.lam(c[4], c[3], var = "lambda", type = "U")[i,r])))
  }
  
  eval_g_eq_1 <- function(c, i, r){
    constr <- c( c[1] + c[2] + c[3] - y1[i,1,r] - y1[i,2,r]) 
    grad <-  c(1, 1, 1, 0)
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
    c0 <- c(equal.neg, equal.neg, 0, 0.7)
  }
  else{
    c0 <- c(equal, equal, 0, 0.7)
  }
  
  # borrowing constraint
  y2.h.min <- apply(y2[,1,,], MARGIN = c(1,2), FUN = min)
  y2.w.min <- apply(y2[,2,,], MARGIN = c(1,2), FUN = min)
  min.consumption <- 0.01
  S.bor.const <- (-beta) * (y2.h.min[i, r]  + y2.w.min[i, r]) * (1 - min.consumption) 
  
  
  ub <- c(equal*2, equal*2, equal*2, 1)
  lb <- c(0, 0, S.bor.const, 0)
  res <- nloptr(x0 = c0, eval_f = eval_f_1_m, lb = lb, ub = ub, 
                eval_g_eq = eval_g_eq_1, opts = opts, i = i, r = r)
  return(res$solution)
}