#' Comparative Statics WRT Rho and Phi
#'
#' This function computes comparative statics wrt to rho (correlation coefficient of spouses' shocks)
#' and phi(the ratio of wife's standard deviation shock to her husband's).
#' @param parameter Pick either rho or phi to do comparative statics
#' @param  num the Length of grid for rho or phi. For cluster use recommended to be a multiplier of 150
#' @param  par.max The maximum of parameter value in the grid. The minimums are calculated automatically; for rho, it is negative of the rho max, for phi, it is the inverse of the max value. For sigma, the length of grid for values greater than one is equal to that of the ones less than one.
#' @param local Specifies the parallelization should happen on local machine or Acropolis cluster. Default is the local machine
#' @export

cs.shocks <- function(parameter = c("rho", "phi"), num, par.max, local = TRUE){

  if (local){
    library(parallel)
    nc <- detectCores() - 1
    cl <- makeForkCluster(nc)
  } else {
    library(snow)
    nc <- 150
    cl <- makeCluster(nc, type = "MPI")
  }

  parameter <- match.arg(parameter)

  if (parameter == "rho"){
    rho.num <- num
    rho.max <- par.max
    rho.seq <- seq(from = (-1) * rho.max, to = rho.max, length.out = rho.num)
    if (local){
      res <- parLapply(cl, rho.seq, function(t) D.prob(sigma_eta_h = param()$sigma_eta_h, rho = t, phi = param()$phi, Consumption = T))
    } else{
      res <- clusterMap(cl, function(t) D.prob(sigma_eta_h = param()$sigma_eta_h, rho = t, phi = param()$phi, Consumption = T), rho.seq)
    }
  }
  if (parameter == "phi"){
    phi.num <- num
    phi.max <- par.max
    phi.seq <- c(setdiff(seq(from =  1/phi.max, to = 1, length.out = floor(phi.num/2)), 1), seq(from = 1, to = phi.max, length.out = ceiling(phi.num/2)))
    if (local){
      res <- parLapply(cl, phi.seq, function(t) D.prob(sigma_eta_h = param()$sigma_eta_h, rho = param()$rho, phi = t, Consumption = T))
    } else{
      res <- clusterMap(cl, function(t) D.prob(sigma_eta_h = param()$sigma_eta_h, rho = param()$rho, phi = t, Consumption = T), phi.seq)
    }
  }

  return(res)
}


