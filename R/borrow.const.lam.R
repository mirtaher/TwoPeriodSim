#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

Borrow.const.lam <- function(i, r){
  
  
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
  tau.d <- par$tau.d
  
  wt.h <- par$wt.h
  wt.w <- par$wt.w
  u <- par$u
  U <- par$U
  u.grad <- par$u.grad
  
  
  equal <- (y1[i,1,r] + y1[i,2,r])/2
  y2.h.min <- apply(y2[,1,,], MARGIN = c(1,2), FUN = min)
  y2.w.min <- apply(y2[,2,,], MARGIN = c(1,2), FUN = min)
  tol.cons <- 0.01
  
  b.mar <- (-beta) * (y2.h.min[i, r]  + y2.w.min[i, r])  * (1 - tol.cons)
  b.h <- (-beta)/delta * y2.h.min[i, r] * (1 - tol.cons)
  b.w <- (-beta)/(1 - delta) * y2.w.min[i, r] * (1 - tol.cons)
  
  res.min <- max(b.mar, b.h, b.w)
  
  s.max <- y1[i,1,r] + y1[i,2,r] * (1 - tol.cons)
 
  return(list("lower" = res.min, "upper" = s.max))
    
}





