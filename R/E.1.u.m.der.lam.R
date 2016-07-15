#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

E.1.u.m.der.lam <-function(lambda, S, var = c("S", "lambda"), type = c("u", "U"), spouse = c("h", "w"), i, r){

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

  A <- 1/beta * S + y2[i,1,r,] + y2[i,2,r,]

  c.h.S <- function(lambda, S){
    res <- (1/beta) / (((1 - lambda)/ lambda)^(1/chi) +1)
    return(res)
  }

  c.w.S <- function(lambda, S){
    res <- 1/beta - c.h.S(lambda, S)
    return(res)
  }

  c.h.lambda <- function(lambda, S){
    res <- A/(chi * lambda^2) * ((1 - lambda)/ lambda)^(1/chi - 1) * 1/(((1 - lambda)/ lambda)^(1/chi) +1)^2
    return(res)
  }

  c.w.lambda <- function(lambda, S){
    res <- (-1) * c.h.lambda(lambda, S)
    return(res)
  }


  if (var == "S"){
    if (!missing(spouse)){
      if (type == "u"){
        if (spouse == "h"){
          U.2.m.vec <- mapply(function(h, h.der) u.grad(h) * h.der ,
                              matrix(c2.m.array.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,1] ,
                              c.h.S(lambda, S))
        }
        if (spouse == "w"){
          U.2.m.vec <- mapply(function(w, w.der) u.grad(w) * w.der,
                              matrix(c2.m.array.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,2],
                              c.w.S(lambda, S))
        }

      }

      if (type == "U"){
        U.2.m.vec <- mapply(function(h, w, h.der, w.der) lambda * u.grad(h) * h.der + (1 - lambda) * u.grad(w) * w.der,
                            matrix(c2.m.array.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,1],
                            matrix(c2.m.array.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,2],
                            c.h.S(lambda, S),
                            c.w.S(lambda, S))

      }
    } else {
      U.2.m.vec <- mapply(function(h, w, h.der, w.der) lambda * u.grad(h) * h.der + (1 - lambda) * u.grad(w) * w.der,
                          matrix(c2.m.array.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,1],
                          matrix(c2.m.array.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,2],
                          c.h.S(lambda, S),
                          c.w.S(lambda, S))
    }

  }

  if (var == "lambda"){
    if (!missing(spouse)){
      if (type == "u"){
        if (spouse == "h"){
          U.2.m.vec <- mapply(function(h, h.der) u.grad(h) * h.der ,
                              matrix(c2.m.array.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,1] ,
                              c.h.lambda(lambda, S))
        }
        if (spouse == "w"){
          U.2.m.vec <- mapply(function(w, w.der) u.grad(w) * w.der,
                              matrix(c2.m.array.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,2],
                              c.w.lambda(lambda, S))
        }

      }

      if (type == "U"){
        U.2.m.vec <- mapply(function(h,w, h.der, w.der) u(h) - u(w) + lambda * (u.grad(h) * h.der - u.grad(w) * w.der) + u.grad(w) * w.der,
                            matrix(c2.m.array.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,1],
                            matrix(c2.m.array.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,2],
                            c.h.lambda(lambda, S),
                            c.w.lambda(lambda, S))

      }
    } else {
      U.2.m.vec <- mapply(function(h,w, h.der, w.der) u(h) - u(w) + lambda * (u.grad(h) * h.der - u.grad(w) * w.der) + u.grad(w) * w.der,
                          matrix(c2.m.array.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,1],
                          matrix(c2.m.array.lam(lambda, S, i, r, 1:reps2), ncol = 2)[,2],
                          c.h.lambda(lambda, S),
                          c.w.lambda(lambda, S))
    }
  }


  # Taking expectation wrt information available at the first period
  res <- mean(U.2.m.vec)
  return(res)
}
