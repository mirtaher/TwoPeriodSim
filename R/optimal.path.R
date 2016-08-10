#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

optimal.path <- function(sigma_eta_h = param()$sigma_eta_h, rho = param()$rho, phi = param()$phi, Parallel = FALSE, local = TRUE){

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

  if (Parallel){
    if (local){
      library(parallel)
      nc <- detectCores() - 1
      cl <- makeForkCluster(nc)
    } else {
      library(snow)
      nc <- 150
      cl <- makeCluster(nc, type = "MPI")
    }
    c1.d.arr <- clusterMap(cl, function(i, r) period.1.d(i,r, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), index[,1], index[,2])
    c1.d <- array(unlist(c1.d.arr), dim = c(3, N, reps1))
    c1.m.arr <- clusterMap(cl, function(i, r) period.1.m(i,r, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), index[,1], index[,2])
    c1.m <- array(unlist(c1.m.arr), dim = c(3, N, reps1))
    c2.m.arr <- clusterMap(cl, function(i, r1, r2) period.2.m(c1.m[3,i,r1], i, r1, r2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), index2[,1], index2[,2], index2[,3])
    c2.m <- array(unlist(c2.m.arr), dim = c(2, N, reps1, reps2))
    c2.d.arr <- clusterMap(cl, function(i, r1, r2) period.2.d(c1.d[3,i,r1], i, r1, r2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), index2[,1], index2[,2], index2[,3])
    c2.d <- array(unlist(c2.d.arr), dim = c(2, N, reps1, reps2))


  } else {
    c1.d <- array(mapply(function(i, r) period.1.d(i,r, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), index[,1], index[,2]), dim = c(3, N, reps1))
    c1.m <- array(mapply(function(i, r) period.1.m(i,r, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), index[,1], index[,2]), dim = c(3, N, reps1))
    c2.m <- array(mapply(function(i, r1, r2) period.2.m(c1.m[3,i,r1], i, r1, r2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), index2[,1], index2[,2], index2[,3]), dim = c(2, N, reps1, reps2))
    c2.d <- array(mapply(function(i, r1, r2) period.2.d(c1.d[3,i,r1], i, r1, r2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), index2[,1], index2[,2], index2[,3]), dim = c(2, N, reps1, reps2))

  }



  res <- list("c1.d" = c1.d, "c1.m" = c1.m, "c2.m" = c2.m, "c2.d" = c2.d)
  return(res)
  }



















