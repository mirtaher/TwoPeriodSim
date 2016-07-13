#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

feasible.lam <- function(i, r){
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
  chi <- par$chi
  
  S0 <- Borrow.const.lam(i, r)
  S0.lower <- S0$lower
  S0.upper <- S0$upper
  uncon <- period.1.m(i, r)
  S.uncon <- uncon[3]
  S.length <- 20
  lam.length <- 20

  
  cond.h <- ifelse(E.1.u.d(S.uncon, type = "u", spouse = "h")[i,r] > E.1.u.m(S.uncon, type = "u", spouse = "h")[i,r], 1, 0)
  cond.w <- ifelse(E.1.u.d(S.uncon, type = "u", spouse = "w")[i,r] > E.1.u.m(S.uncon, type = "u", spouse = "w")[i,r], 1, 0)
  cond <- (cond.h | cond.w)
  if (!cond) {
    res <- list("status" = "Stay Married with Old Terms")
  } 
  if (cond) {
    if (cond.h){
      lam.seq <- seq(from = 0.51, to = 0.99, length.out = lam.length) 
    }
    if(cond.w){
      lam.seq <- seq(from = 0.49, to = 0.01, length.out = lam.length)
    }
    s.seq <- seq(from = S0.lower, to = S0.upper, length.out = S.length )
    ind <- 0
    con <- TRUE
    while(con){
      ind <- ind + 1
      if (ind > S.length) break
      else{
        S0.lam <- s.seq[ind]
        for (l in lam.seq){
          con.h <- ifelse(E.1.u.d(S0.lam, type = "u", spouse = "h")[i,r] > E.1.u.m.lam(l, S0.lam, type = "u", spouse = "h")[i,r], 1, 0)
          con.w <- ifelse(E.1.u.d(S0.lam, type = "u", spouse = "w")[i,r] > E.1.u.m.lam(l, S0.lam, type = "u", spouse = "w")[i,r], 1, 0)
          con <- con.h | con.w
          if (!con) break
        }
      }
    }
    if (con) res <- list("status" = "Divorce") 
    if (!con) res <- list("S0" = S0.lam, "lam0" = l, "status" = "Stay Married with New Terms")
  }
  
  return(res)
}
  
  
  
  
  
  
  
  