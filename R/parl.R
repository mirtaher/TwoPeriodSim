#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

parl <- function(x, parameter = c("Rho", "Phi"), other, sigma_eta_h, type = c("h", "w", "n", "c")){
  library(parallel)
  nc <- detectCores() - 1
  cl <- makeForkCluster(nc)
  parameter <- match.arg(parameter)
  if (parameter == "Rho"){
    res <- parLapply(cl, x, function(t) D.prob(sigma_eta_h = sigma_eta_h, Rho = t, Phi = other)[[type]])
  }
  if (parameter == "Phi"){
    res <- parLapply(cl, x, function(t) D.prob(sigma_eta_h = sigma_eta_h, Rho = other, Phi = t)[[type]])
  }

  return(res)
}
