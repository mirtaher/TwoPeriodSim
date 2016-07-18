#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

optimal.path <- function(sigma_eta_h, Rho, Phi){

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

  y <- incomeProcess(sigma_eta_h, Rho, Phi)
  options(warn = -1)
  if (is.na(y)){
    res <- NA
  }
  else{
    y1 <<- y$y1
    y2 <<- y$y2

    # obtaining optimal consumption paths
    c1.d <- array(mapply(function(i, r) period.1.d(i,r), index[,1], index[,2]), dim = c(3, N, reps1))
    c1.m <- array(mapply(function(i, r) period.1.m(i,r), index[,1], index[,2]), dim = c(3, N, reps1))
    c2.m <- aperm(array(mapply(function(i, r1, r2) period.2.m(c1.m[3,i,r1], i, r1, r2), index2[,1], index2[,2], index2[,3]), dim = c(2, N, reps1, reps2)), c(2,1,3,4))
    c2.d <- aperm(array(mapply(function(i, r1, r2) period.2.d(c1.d[3,i,r1], i, r1, r2), index2[,1], index2[,2], index2[,3]), dim = c(2, N, reps1, reps2)), c(2,1,3,4))

    res <- list("c1.d" = c1.d, "c1.m" = c1.m, "c2.m" = c2.m, "c2.d" = c2.d)
  }

  return(res)
}
















