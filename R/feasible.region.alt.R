#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

feasible.region.alt <- function(i, r, Extensive = FALSE, Optimize = FALSE){
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

  sigma_eta_h <- par$sigma_eta_h
  rho <- par$rho
  phi <- par$phi

  income <- TwoPeriodSim::incomeProcess(sigma_eta_h = sigma_eta_h, Rho = rho, Phi = phi)
  y1 <- income$y1
  y2 <- income$y2

  S0 <- borrow.const(i, r)
  S0.lower <- S0$lower
  S0.upper <- S0$upper
  uncon <- period.1.m(i, r)
  S.uncon <- uncon[3]
  S.length <- 200
  lam.length <- 50

  cond.h <- ifelse(E.1.u.d(S.uncon, type = "u", spouse = "h", i, r) > E.1.u.m(S.uncon, type = "u", spouse = "h", i, r), 1, 0)
  cond.w <- ifelse(E.1.u.d(S.uncon, type = "u", spouse = "w", i, r) > E.1.u.m(S.uncon, type = "u", spouse = "w", i, r), 1, 0)
  cond <- (cond.h | cond.w)

  if (!cond) {
    res <- list("status" = "Stay Married with Old Terms", "c.h.uncon" = uncon[1], "c.w.uncon" = uncon[2], "s.uncon" = uncon[3])
    if (Optimize){
      res <- c(uncon[1], uncon[2], uncon[3], 0.5)
    }

  }
  if (cond) {
    if (cond.h){
      lam.seq <- seq(from = 0.99, to = 0.51, length.out = lam.length)
    }
    if(cond.w){
      lam.seq <- seq(from = 0.01, to = 0.49, length.out = lam.length)
    }
    s.seq <- seq(from = S0.lower, to = S0.upper, length.out = S.length )
    ind <- expand.grid(lam.seq, s.seq)
    consumption <- mapply(function(S) (y1[i,1,r] + y1[i,2,r] - S)/2 , ind[, 2])
    participate <- mapply(function(i1, i2) participate.lam(i1, i2, i, r), ind[, 1], ind[, 2] )
    df <- data.frame( consumption, consumption, ind[,2], ind[,1], participate)
    names(df) <- c( "consumption.h", "consumption.w", "saving", "lambda", "participate")
    start.points <- df[df$participate == TRUE,]

    if (  is.null(start.points) ){
      res <- list("status" = "Divorce")
      sol <- period.1.d(i, r)
      res <- list("sol" = sol, "status" = "Divorce")

      if (Optimize){
        res <- c(sol, NA)
      }

    } else {
      start.points$values <- mapply(function(c.h, c.w, S, lam) obj.lam(c.h, c.w, S, lam, i, r),
                       start.points$consumption.h,
                       start.points$consumption.w,
                       start.points$saving,
                       start.points$lambda)

      sol <- start.points[which.min(start.points$values), ]
      sol <- as.vector(sol)
      res <- res <- list("S0" = sol[3], "lam0" = sol[4], "c0" = sol[1], "status" = "Stay Married with New Terms")

      if (Optimize){
        res <- sol
      }

      if (Extensive){
        res <- list("S0" = sol[3], "lam0" = sol[4], "c0" = sol[1], "status" = "Stay Married with New Terms",
                    "region" = start.points, "c.h.uncon" = uncon[1], "c.w.uncon" = uncon[2], "s.uncon" = uncon[3])
      }
    }

  }

  return(res)
}







