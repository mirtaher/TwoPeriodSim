#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

test.U.lambda <-function(lambda, S){
  
  # just to test commit 
  
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
  
  
  # Cartesian product
  index <- expand.grid(1:N, 1:reps1)
  index2 <- expand.grid(1:N, 1:reps1, 1:reps2)
  
  A <- 1/beta * S + y2[,1,,] + y2[,2,,]
  
  c.h.S <- function(lambda, S){
    res <- (1/beta) / (((1 - lambda)/ lambda)^(1/chi) +1)
  }
  
  c.w.S <- function(lambda, S){
    res <- 1/beta - c.h.S(lambda, S)
  }
  
  c.h.lambda <- function(lambda, S){
    res <- A/(chi * lambda^2) * ((1 - lambda)/ lambda)^(1/chi - 1) * 1/(((1 - lambda)/ lambda)^(1/chi) +1)^2
  }
  
  c.w.lambda <- function(lambda, S){
    res <- (-1) * c.h.lambda(lambda, S)
  }
  
  
  U.2.m.vec <- mapply(function(h,w) U2.lam(lambda, h, w), c2.m.array.lam(lambda, S)[,1,,], c2.m.array.lam(lambda, S)[,2,,])
    
  
  # Taking expectation wrt information available at the first period
  U.2.m <- array(U.2.m.vec, dim = c(N,reps1, reps2))
  res <- apply(U.2.m, MARGIN = c(1,2), FUN = mean)
  return(res)
}