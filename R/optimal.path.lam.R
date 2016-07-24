#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

optimal.path.lam <- function(Precise = FALSE, local = TRUE){

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

  # Set up the parallel
  library(parallel)

  if (local){
    nc <- detectCores() - 1
    cl <- makeForkCluster(nc)
  } else {
    nc <- 12
    cl <- makeCluster(nc)
  }


  # obtaining optimal consumption paths
  if (!Precise){
    c1.m <- clusterMap(cl, function(i, r) feasible.region(i,r, Optimize = TRUE), index[,1], index[,2])
    c1.m <- array(unlist(c1.m), dim = c(4, N, reps1))

  } else {
    c1.m <- clusterMap(cl, function(i, r) period.1.m.lam(i,r)$sol, index[,1], index[,2])
    c1.m <- array(unlist(c1.m), dim = c(4, N, reps1))
  }

  c1.d <- clusterMap(cl, function(i, r) period.1.d(i,r), index[,1], index[,2])
  c1.d <- array(unlist(c1.d), dim = c(3, N, reps1))

  c2.m <- clusterMap(cl, function(i, r1, r2) period.2.m.lam(c1.m[4,i,r1], c1.m[3,i,r1], i, r1, r2), index2[,1], index2[,2], index2[,3])
  c2.m <- array(unlist(c2.m), dim = c(2, N, reps1, reps2))

  c2.d <- clusterMap(cl, function(i, r1, r2) period.2.d(c1.d[3,i,r1], i, r1, r2), index2[,1], index2[,2], index2[,3])
  c2.d <- array(unlist(c2.d), dim = c(2, N, reps1, reps2))


  res <- list("c1.d" = c1.d, "c1.m" = c1.m, "c2.m" = c2.m, "c2.d" = c2.d)
  return(res)
  }



















