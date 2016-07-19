#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

feasible.region <- function(i, r){
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

  S0 <- borrow.const(i, r)
  S0.lower <- S0$lower
  S0.upper <- S0$upper
  uncon <- period.1.m(i, r)
  S.uncon <- uncon[3]
  S.length <- 20
  lam.length <- 20

  cond.h <- ifelse(E.1.u.d(S.uncon, type = "u", spouse = "h", i, r) > E.1.u.m(S.uncon, type = "u", spouse = "h", i, r), 1, 0)
  cond.w <- ifelse(E.1.u.d(S.uncon, type = "u", spouse = "w", i, r) > E.1.u.m(S.uncon, type = "u", spouse = "w", i, r), 1, 0)
  cond <- (cond.h | cond.w)

  if (!cond) {
    res <- list("status" = "Stay Married with Old Terms", "c.h.uncon" = uncon[1], "c.w.uncon" = uncon[2], "s.uncon" = uncon[3])

  }
  if (cond) {
    if (cond.h){
      lam.seq <- seq(from = 0.99, to = 0.51, length.out = lam.length)
    }
    if(cond.w){
      lam.seq <- seq(from = 0.01, to = 0.49, length.out = lam.length)
    }
    s.seq <- seq(from = S0.lower, to = S0.upper, length.out = S.length )
    reg <- vector()
    temp <- c()
    for (is in 1:S.length){
      for(l in lam.seq){
        con.h <- ifelse(E.1.u.d(s.seq[is], type = "u", spouse = "h", i, r) <= E.1.u.m.lam(l, s.seq[is], type = "u", spouse = "h", i, r), 1, 0)
        con.w <- ifelse(E.1.u.d(s.seq[is], type = "u", spouse = "w", i, r) <= E.1.u.m.lam(l, s.seq[is], type = "u", spouse = "w", i, r), 1, 0)
        con <- con.h & con.w
        if (con) temp <- c( temp, c(s.seq[is], l) )
      }
    }
    sl <-matrix(temp, nrow = 2)
    c <- (y1[i,1,r] + y1[i,2,r] - sl[1, ])/2
    res <- cbind(c, c, sl[1, ], sl[2, ])

  }

  return(res)
}







