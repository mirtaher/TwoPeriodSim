#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
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


