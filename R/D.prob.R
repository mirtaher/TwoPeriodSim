#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

D.prob <- function(sigma_eta_h, Rho, Phi){

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

  cons <- optimal.path(sigma_eta_h, Rho, Phi)
  options(warn = -1)
  if (is.na(cons)){
    res <- NA
  }
  else{
    c1.d <- cons$c1.d
    c1.m <- cons$c1.m
    c2.d <- cons$c2.d
    c2.m <- cons$c2.m

    # non-cooperative decision to divorce
    D.h <- ifelse(E.1.u.d.opt(c2.d, type = "u", spouse = "h") > E.1.u.m.opt(c2.m, type = "u", spouse = "h"), 1, 0)
    D.w <- ifelse(E.1.u.d.opt(c2.d, type = "u", spouse = "w") > E.1.u.m.opt(c2.m, type = "u", spouse = "w"), 1, 0)
    D.n <- ifelse(D.h | D.w, 1 , 0)

    # Cooperative decision to divorce
    U.opt.m.per1 <- matrix(mapply(function(i, r) U(c1.m[1,i,r], c1.m[2,i,r]), index[,1], index[,2]), nrow = N, ncol = reps1)
    U.opt.m <- U.opt.m.per1  + beta * E.1.u.m.opt(c2.m, type = "U")
    U.opt.d <- U.opt.m.per1  + beta * E.1.u.d.opt(c2.d, type = "U")
    D.c <- ifelse(U.opt.d > U.opt.m, 1, 0)

    # Taking expectation
    res.h <- apply(D.h, MARGIN = 2, FUN = mean, na.rm = T)
    res.w <- apply(D.w, MARGIN = 2, FUN = mean, na.rm = T)
    res.n <- apply(D.n, MARGIN = 2, FUN = mean, na.rm = T)
    res.c <- apply(D.c, MARGIN = 2, FUN = mean, na.rm = T)
    res <- list("h" = res.h, "w" = res.w, "n" = res.n, "c" = res.c)
  }

  return(res)
}
















