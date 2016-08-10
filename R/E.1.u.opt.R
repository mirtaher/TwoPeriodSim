#' Expected Value Function
#'
#' This function calculates the expected second period value function given the information revealed in the first period.
#'  Essentially, it evaluates the objective function in the optimal second period consumption values and takes the expectation.
#' @param c2.opt The optimal consumption in the second
#' @param type Could be individual utility (u) or the joint utility of the couple (U) in a unitary framework
#' @param spouse If the type is individual utility (u), we should specify we mean husband or wife utility
#' @export

E.1.u.opt <- function(c2.opt, type = c("u", "U"), spouse = c("h", "w")){

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
        U.2.m.vec <- mapply(function(i, r1, r2) u(c2.opt[1,i,r1,r2]), index2[,1], index2[,2], index2[,3])
      }
      if (spouse == "w"){
        U.2.m.vec <- mapply(function(i, r1, r2) u(c2.opt[2,i,r1,r2]), index2[,1], index2[,2], index2[,3])
      }

    }

    if (type == "U"){
      U.2.m.vec <- mapply(function(i, r1, r2) U(c2.opt[1,i,r1,r2], c2.opt[2,i,r1,r2]), index2[,1], index2[,2], index2[,3])

    }
  } else {
    U.2.m.vec <- mapply(function(i, r1, r2) U(c2.opt[1,i,r1,r2], c2.opt[2,i,r1,r2]), index2[,1], index2[,2], index2[,3])
  }


  U.2.m <- array(U.2.m.vec, dim = c(N,reps1, reps2))
  res <- apply(U.2.m, MARGIN = c(1,2), FUN = mean)
  return(res)
}
