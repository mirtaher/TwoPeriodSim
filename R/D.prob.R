#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

D.prob <- function(sigma_eta_h, Rho, Phi){
  
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
    c2.m <- aperm(array(mapply(function(i, r1, r2) period.2(c1.m[3,i,r1], y2[i,1, r1, r2], y2[i,2, r1, r2], Analytical = T),  index2[,1], index2[,2], index2[,3]), dim = c(2, N, reps1, reps2)), c(2,1,3,4))
    c2.d <- aperm(array(mapply(function(i, r1, r2) c(1/beta * delta * c1.m[3,i,r1] + y2[i,1, r1, r2], 1/beta * (1 - delta) * c1.m[3,i,r1] + y2[i,2, r1, r2] ) ,index2[,1], index2[,2], index2[,3]), dim = c(2, N, reps1, reps2)), c(2,1,3,4))
    
    # non-cooperative decision to divorce
    D.h <- ifelse(E.1.u.d.opt(c2.d, type = "u", spouse = "h") > E.1.u.m.opt(c2.m, type = "u", spouse = "h"), 1, 0)
    D.w <- ifelse(E.1.u.d.opt(c2.d, type = "u", spouse = "w") > E.1.u.m.opt(c2.m, type = "u", spouse = "w"), 1, 0)
    D.n <- ifelse(D.h | D.w, 1 , 0)
    
    # Cooperative decision to divorce 
    U.opt.m.per1 <- matrix(mapply(function(i, r) U(c1.m[1,i,r], c1.m[2,i,r]), index[,1], index[,2]), nrow = N, ncol = reps1)
    U.opt.m <- U.opt.m.per1  + beta * E.1.u.m.opt(c2.m, type = "U")
    U.opt.d <- U.opt.m.per1  + beta * E.1.u.d.opt(c2.d, type = "U")
    D.c <- ifelse(U.opt.d > U.opt.m, 1, 0)
    
    # Taking expectation
    res.h <- apply(D.h, MARGIN = 2, FUN = mean, na.rm = T)
    res.w <- apply(D.w, MARGIN = 2, FUN = mean, na.rm = T)
    res.n <- apply(D.n, MARGIN = 2, FUN = mean, na.rm = T)
    res.c <- apply(D.c, MARGIN = 2, FUN = mean, na.rm = T)
    res <- list("h" = res.h, "w" = res.w, "n" = res.n, "c" = res.c)
  }
  
  return(res)
}
















