#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

E.1.u.d.opt <- function(c2.d, type = c("u", "U"), spouse = c("h", "w")){
  
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
  
  if (!missing(spouse)){
    if (type == "u"){
      if (spouse == "h"){
        U.2.m.vec <- mapply(function(i, r1, r2) u(c2.d[i,1,r1,r2]), index2[,1], index2[,2], index2[,3])
      }
      if (spouse == "w"){
        U.2.m.vec <- mapply(function(i, r1, r2) u(c2.d[i,2,r1,r2]), index2[,1], index2[,2], index2[,3])
      }
      
    }
    
    if (type == "U"){
      U.2.m.vec <- mapply(function(i, r1, r2) U(c2.d[i,1,r1,r2], c2.d[i,2,r1,r2]), index2[,1], index2[,2], index2[,3])
      
    }
  } else {
    U.2.m.vec <- mapply(function(i, r1, r2) U(c2.d[i,1,r1,r2], c2.d[i,2,r1,r2]), index2[,1], index2[,2], index2[,3])
  }
  
  
  U.2.m <- array(U.2.m.vec, dim = c(N,reps1, reps2))
  res <- apply(U.2.m, MARGIN = c(1,2), FUN = mean)
  return(res)
}