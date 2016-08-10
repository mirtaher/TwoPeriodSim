#' Initial parameters distribution
#'
#' This function simulate the income process for two spaouses
#' @export

E.1.u.m<-function(S, type = c("u", "U"), spouse = c("h", "w"), i, r, sigma_eta_h = param()$sigma_eta_h, rho = param()$rho, phi = param()$phi){

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
        U.2.m.vec <- mapply(function(h) u(h),  matrix(period.2.m(S, i, r, 1:reps2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), ncol = 2)[,1])
      }
      if (spouse == "w"){
        U.2.m.vec <-  mapply(function(w) u(w), matrix(period.2.m(S, i, r, 1:reps2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), ncol = 2)[,2])
      }

    }

    if (type == "U"){
      U.2.m.vec <- mapply(function(h,w) U(h, w),
                          matrix(period.2.m(S, i, r, 1:reps2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), ncol = 2)[,1],
                          matrix(period.2.m(S, i, r, 1:reps2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), ncol = 2)[,2])

    }
  } else {
    U.2.m.vec <- mapply(function(h,w) U(h, w),
                        matrix(period.2.m(S, i, r, 1:reps2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), ncol = 2)[,1],
                        matrix(period.2.m(S, i, r, 1:reps2, sigma_eta_h = sigma_eta_h, rho = rho, phi = phi), ncol = 2)[,2])
  }

  res <- mean(U.2.m.vec, na.rm = TRUE)
  return(res)
}
