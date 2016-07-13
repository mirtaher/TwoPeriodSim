#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

period.1.m.KT.01 <- function(i, r){
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
  tau.d <- par$tau.d
  
  wt.h <- par$wt.h
  wt.w <- par$wt.w
  u <- par$u
  U <- par$U
  u.grad <- par$u.grad
  
  # Cartesian product
  index <- expand.grid(1:N, 1:reps1)
  index2 <- expand.grid(1:N, 1:reps1, 1:reps2)
  
  eval_f_1_m <- function(c, i, r){
    return(list("objective" = - U(c[1], c[2])- beta * E.1.u.m(c[3], type = "U")[i,r] -
                  c[4] * log((-1) * ((1 + tau.d) * E.1.u.d(c[3], type = "u", spouse = "h")[i,r] - E.1.u.m(c[3], type = "u", spouse = "h")[i,r]))  -
                  c[5] * ((1 + tau.d) * E.1.u.d(c[3], type = "u", spouse = "w")[i,r] - E.1.u.m(c[3], type = "u", spouse = "w")[i,r]), 
                "gradient" = c(-wt.h * u.grad(c[1]), 
                               -wt.w * u.grad(c[2]), 
                               (-beta) * E.1.u.m.der(c[3], type = "U")[i,r] - c[4] * 1/((1 + tau.d) * E.1.u.d(c[3], type = "u", spouse = "h")[i,r] - E.1.u.m(c[3], type = "u", spouse = "h")[i,r]) * ((1 + tau.d) * E.1.u.d.der(c[3], type = "u", spouse = "h")[i,r] - E.1.u.m.der(c[3], type = "u", spouse = "h")[i,r]) - c[5] * ((1 + tau.d) * E.1.u.d.der(c[3], type = "u", spouse = "w")[i,r] - E.1.u.m.der(c[3], type = "u", spouse = "w")[i,r]),
                               (-1) * log((-1) * ((1 + tau.d) * E.1.u.d(c[3], type = "u", spouse = "h")[i,r] - E.1.u.m(c[3], type = "u", spouse = "h")[i,r])),
                               (-1) * ((1 + tau.d) * E.1.u.d(c[3], type = "u", spouse = "w")[i,r] - E.1.u.m(c[3], type = "u", spouse = "w")[i,r]) )
    )
    )
    
  }
  
  eval_g_eq_1 <- function(c, i, r){
    constr <- c( c[1] + c[2] + c[3] - y1[i,1,r] - y1[i,2,r]) 
    grad <-  c(1, 1, 1, 0, 0)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel"  = 1.0e-7)
  opts <- list( "algorithm" = "NLOPT_LD_SLSQP",
                "xtol_rel"  = 1.0e-7,
                "maxeval"   = 10000,
                "local_opts" = local_opts )
  
  # finding a feasible saving level 
  s.length <- 100
  b.const <- -ybar/2
  s.const <-  ybar/2
  s.seq <-seq(from = b.const, to = s.const, length.out = s.length)
  si <- 0
  cond <- TRUE
  while(cond) {
    si <- si + 1
    if (si > s.length) stop("No feasible saving found")
    else{
      S0 <- s.seq[si]
      cond.h <- ifelse((1 + tau.d) * E.1.u.d(S0, type = "u", spouse = "h")[i,r] > E.1.u.m(S0, type = "u", spouse = "h")[i,r], 1, 0)
      cond <- cond.h
    }
  }
  
  equal <- (y1[i,1,r] + y1[i,2,r])/2
  if (equal < 0) {
    y1.expand <- array(rep(y1,reps2), dim = c(N,2,reps1, reps2))
    equal.neg.expand <- (y1.expand[i,1,r,] + y1.expand[i,2,r,] + y2[i,1,r,] + y2[i,2,r,])/4
    equal.neg <- mean(equal.neg.expand)
    c0 <- c(equal.neg, equal.neg, S0, 0, 0)
  }
  else{
    c0 <- c(equal, equal, S0, 0, 0)
  }
  
  # borrowing constraint
  y2.h.min <- apply(y2[,1,,], MARGIN = c(1,2), FUN = min)
  y2.w.min <- apply(y2[,2,,], MARGIN = c(1,2), FUN = min)
  S.bor.const <- (-beta) * (y2.h.min[i, r]  + y2.w.min[i, r])
  
  ub <- c(equal*2, equal*2, equal*2, Inf, Inf)
  lb <- c(0, 0, S.bor.const, 0, 0)
  res <- nloptr(x0 = c0, eval_f = eval_f_1_m, lb = lb, ub = ub, 
                eval_g_eq = eval_g_eq_1, opts = opts, i = i, r = r)
  return(res$solution)
}