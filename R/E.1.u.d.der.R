#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

E.1.u.d.der <-function(S, type = c("u", "U"), spouse = c("h", "w"), i, r){

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

  sigma_eta_h <- par$sigma_eta_h
  rho <- par$rho
  phi <- par$phi

  income <- TwoPeriodSim::incomeProcess(sigma_eta_h = sigma_eta_h, Rho = rho, Phi = phi)
  y1 <- income$y1
  y2 <- income$y2

  # Cartesian product
  index <- expand.grid(1:N, 1:reps1)
  index2 <- expand.grid(1:N, 1:reps1, 1:reps2)

  c2.d.array <- function(S){
    res <- 1/beta * delta * S + y2
    return(res)
  }

  if (!missing(spouse)){
    if (type == "u"){
      if (spouse == "h"){
        U.2.m.vec <- mapply(function(h) u.grad(h) * delta * 1/beta ,
                            matrix(period.2.d(S, i, r, 1:reps2), ncol = 2)[,1])
      }
      if (spouse == "w"){
        U.2.m.vec <- mapply(function(w) u.grad(w) * (1-delta) * 1/beta,
                            matrix(period.2.d(S, i, r, 1:reps2), ncol = 2)[,2])
      }
    }

    if (type == "U"){
      U.2.m.vec <- mapply(function(h,w) wt.h * u.grad(h) * delta * 1/beta + wt.w * u.grad(w) * (1-delta) * 1/beta,
                          matrix(period.2.d(S, i, r, 1:reps2), ncol = 2)[,1],
                          matrix(period.2.d(S, i, r, 1:reps2), ncol = 2)[,2])
    }
  } else {
    U.2.m.vec <- mapply(function(h,w) wt.h * u.grad(h) * delta * 1/beta + wt.w * u.grad(w) * (1-delta) * 1/beta,
                        matrix(period.2.d(S, i, r, 1:reps2), ncol = 2)[,1],
                        matrix(period.2.d(S, i, r, 1:reps2), ncol = 2)[,2])
  }


  res <- mean(U.2.m.vec)
  return(res)
}
