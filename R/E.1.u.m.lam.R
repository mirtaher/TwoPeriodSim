#' Expected Second Period Marriage Utility Given The First Period Info in the Collective Framework
#'
#' This function calculates the expected second period marriage utility given the information revealed at the first period for a given level of savings and a new sharing rule (for the second period) for an individual family.
#' Essentially, it calculates the second period optimal consumption for the given level of saving and updated sharing rule for different
#' relizations of the shocks in the second period and then substiute them in the individual utility functions (u) or the joint utility function(U).
#' Finally, it takes the average to return the expectations.
#' @param lambda Sharing rule
#' @param S Saving
#' @param type Could be individual utility (u) or the joint utility of the couple (U) in a unitary framework
#' @param spouse If the type is individual utility (u), we should specify we mean husband or wife utility
#' @param i The marriage index
#' @param r First period repetition. It is not necessary to be greater than one for the first period. It is needed for taking expectaions, which is required in the second period
#' @export


E.1.u.m.lam <-function(lambda, S, type = c("u", "U"), spouse = c("h", "w"), i, r){

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


  if (!missing(spouse)){
    if (type == "u"){
      if (spouse == "h"){
        U.2.m.vec <- mapply(function(h) u(h),
                            matrix(period.2.m.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,1])
      }
      if (spouse == "w"){
        U.2.m.vec <-  mapply(function(w) u(w),
                             matrix(period.2.m.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,2])
      }

    }

    if (type == "U"){
      U.2.m.vec <- mapply(function(h,w) U2.lam(lambda, h, w),
                          matrix(period.2.m.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,1],
                          matrix(period.2.m.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,2])

    }
  } else {
    U.2.m.vec <- mapply(function(h,w) U2.lam(lambda, h, w),
                        matrix(period.2.m.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,1],
                        matrix(period.2.m.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,2])
  }

  res <- mean(U.2.m.vec, na.rm = TRUE)
  return(res)
}
