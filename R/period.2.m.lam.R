#' The Second Period Problem in the Collective Framework
#'
#' This function finds the optimal spouses' consumption levels in the second period for agiven level of saving and sharing rule.
#' We can use the optimization procedure or just use the analytical solution. The default is to use the analytical solution
#'  and it is necessary to be worked in some other parts of the code.
#'
#' @param lambda Sharing Rule
#' @param S Saving
#' @param i The marriage index
#' @param r1 First period repetition. It is not necessary to be greater than one for the first period. It is needed for taking expectaions, which is required in the second period
#' @param r2 Second period repetition. It is needed for taking expectaions, which is required in the second period
#' @param Analytical If it is on it uses the analytical solution instead of using the optimization procedure.
#' @export

period.2.m.lam <- function(lambda, S, i, r1, r2, Analytical = TRUE){

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

  sigma_eta_h <- par$sigma_eta_h
  rho <- par$rho
  phi <- par$phi

  income <- TwoPeriodSim::incomeProcess(sigma_eta_h = sigma_eta_h, Rho = rho, Phi = phi)
  y1 <- income$y1
  y2 <- income$y2

  # Cartesian product
  index <- expand.grid(1:N, 1:reps1)
  index2 <- expand.grid(1:N, 1:reps1, 1:reps2)

  max.c <- 1/beta * S + y2[i, 1, r1, r2] + y2[i, 2, r1, r2]

  if (!Analytical) {
    f <- function(x) U2.lam(lambda, x, (max.c -x))
    res <- optimize(f, interval = c(1e-10, max.c), maximum = T)
    return(c(res$maximum, (max.c -res$maximum)))
  } else {
    res.h <- 1/(((1 - lambda)/lambda)^(1/chi) + 1) * max.c
    res.w <- max.c - res.h
    return(c(res.h, res.w))
  }

}
