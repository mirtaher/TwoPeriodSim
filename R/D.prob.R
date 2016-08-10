#' Probability of Divorce and Final Consumption Path for a Given Income Shocks Pars
#'
#' This function calculates the probability of divorce for a given income shock process parameters. In addition, it
#' returns the optimal consumption paths after taking into account the decision to divorce. The decision to divorec could be cooperatively ("c")
#' Non-cooperatively ("n"). In non-cooperative way, each couple compares the expected second period utility between divorec and marriage.
#' The marriage dessolves once one of the couples betetr off to get divorced. In the cooperative way, in a unitary framework, the couple substitute
#' the consumption levels in the state of divore and marraige in the joint utility function with the same sharing rule (i.e. 0.5). They
#' decide to diviorce t=if the collective benefit of divorce outweight the collective utility of staying married. The function also returns
#' the probability that the participations constraints of husband and wife is binding.
#' @param Consumption Returns the optimal consumption path and saving after taking into account the optimal decision divorce. By default it is TRUE
#' @export

# It returns Regime
D.prob <- function(sigma_eta_h = param()$sigma_eta_h, rho = param()$rho, phi = param()$phi, Consumption = TRUE){

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

  cons <- optimal.path(sigma_eta_h = sigma_eta_h, rho = rho, phi = phi)

  c1.d <- cons$c1.d
  c1.m <- cons$c1.m
  c2.d <- cons$c2.d
  c2.m <- cons$c2.m

  # non-cooperative decision to divorce
  D.h <- ifelse(E.1.u.opt(c2.d, type = "u", spouse = "h") > E.1.u.opt(c2.m, type = "u", spouse = "h"), 1, 0)
  D.w <- ifelse(E.1.u.opt(c2.d, type = "u", spouse = "w") > E.1.u.opt(c2.m, type = "u", spouse = "w"), 1, 0)
  D.n <- ifelse(D.h | D.w, 1 , 0)

  # Cooperative decision to divorce
  U.opt.m.per1 <- matrix(mapply(function(i, r) U(c1.m[1,i,r], c1.m[2,i,r]), index[,1], index[,2]), nrow = N, ncol = reps1)
  U.opt.m <- U.opt.m.per1  + beta * E.1.u.opt(c2.m, type = "U")
  U.opt.d <- U.opt.m.per1  + beta * E.1.u.opt(c2.d, type = "U")
  D.c <- ifelse(U.opt.d > U.opt.m, 1, 0)

  # Final decision about consumption and savings after deciding about divorce
  D.c.1 <- array(rep(D.c, each = 3), dim = c(3, N, reps1))
  D.n.1 <- array(rep(D.n, each = 3), dim = c(3, N, reps1))
  D.c.2 <- array(rep(D.c, each = 2), dim = c(2, N, reps1, reps2))
  D.n.2 <- array(rep(D.n, each = 2), dim = c(2, N, reps1, reps2))

  c1.n <- c1.d * D.n.1 + c1.m * (1 - D.n.1)
  c2.n <- c2.d * D.n.2 + c2.m * (1 - D.n.2)
  c1.c <- c1.d * D.c.1 + c1.m * (1 - D.c.1)
  c2.c <- c2.d * D.c.2 + c2.m * (1 - D.c.2)


  # Taking expectation
  res.h <- apply(D.h, MARGIN = 2, FUN = mean, na.rm = T)
  res.w <- apply(D.w, MARGIN = 2, FUN = mean, na.rm = T)
  res.n <- apply(D.n, MARGIN = 2, FUN = mean, na.rm = T)
  res.c <- apply(D.c, MARGIN = 2, FUN = mean, na.rm = T)

  if (Consumption){
    res <- list("h" = res.h, "w" = res.w, "n" = res.n, "c" = res.c,
                "c1.n" = c1.n, "c2.n" = c2.n, "c1.c" = c1.c, "c2.c" = c2.c)
  } else{
    res <- list("h" = res.h, "w" = res.w, "n" = res.n, "c" = res.c)
  }

  return(res)
  }


















