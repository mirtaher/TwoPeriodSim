#' The Feasible Region Where Participation Constraints Satisfy
#'
#' This function identifies the pairs of saving and sharing rule where satisfy the participation constraints of both couples. Therefore,
#' under any pair of these, marriage can be a superior option over divorce. Given that we are assuming the sharing rule to be a half in the first period,
#' the optimal consumption in the first period is going to be equal. Therefore, given the value of saving, we can obtain the first period consumption levels.
#' These levels provide a feasible initial guess for the optimization routine. Then, thuru a simple two-dimensional grid search among all
#' feasible pairs we find the one leading to the smallest value of the objective function. Since, the length of serach area for saving is roughly 2 (between 1 and -1)
#' and 0.5 for the lambda (because we know either husband or wife should be compensated), with 200 and 50 for saving and sharing rule, respectively
#' we can pin down the optimal level of saving and lambda with two decimal point accuracy.
#' Given that the formal optimization procedure has to be global and therefore is so slow, this simple grid search is practically more desirable.
#' In the following I describe how the procedure works. First we obtain the max and min of saving to pin down the two ends of the saving grid. Then,
#' find the optimal savings under divorcea and marriage under the same sharing rule (unconstrained solution). If even with the same sharing rule the participation
#' contsraints of none of the spouses is binding, then the couple stays married with the old terms. Otherwise, we conduct a two-dimensional
#' search over saving and lambda that keeps both spouses better off in the marriage relative to divorce. We find the best of
#' these pairs as the rough approximation of the optimal saving and lambda. On the other hand, if the feasible set is empty the couple divorce
#' because there is no saving and sharing rule that dominates the expected utility of divorce for both spouses.
#' Notice if the income is missing the function return a vector of four NAs.
#'
#' @param i The marriage index
#' @param r First period repetition. It is not necessary to be greater than one for the first period. It is needed for taking expectaions, which is required in the second period
#' @param Extensive Returns extra ouput including the whole feasible region, and unconstrained (with no change of sharing rule) optimum
#' consumption and saving.
#' @param Optimize Returns only four elements vector of solution and no other output
#' @export


feasible.region <- function(i, r, Extensive = FALSE, Optimize = FALSE){
  library(nloptr)

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
  chi <- par$chi

  sigma_eta_h <- par$sigma_eta_h
  rho <- par$rho
  phi <- par$phi

  income <- TwoPeriodSim::incomeProcess(sigma_eta_h = sigma_eta_h, Rho = rho, Phi = phi)
  y1 <- income$y1
  y2 <- income$y2

  S0 <- borrow.const(i, r)
  S0.lower <- S0$lower
  S0.upper <- S0$upper

  miss <- is.na(y1[i,1,r]) | is.na(y1[i,2,r])
  if (miss){
    res <- rep(NA, 4)
  } else {
    uncon <- period.1.m(i, r)
    S.uncon <- uncon[3]

    sol.d <- period.1.d(i, r)
    s.d <- sol.d[3]

    S.length <- 200
    lam.length <- 50

    cond.h <- ifelse(E.1.u.d(S.uncon, type = "u", spouse = "h", i, r) > E.1.u.m(S.uncon, type = "u", spouse = "h", i, r), 1, 0)
    cond.w <- ifelse(E.1.u.d(S.uncon, type = "u", spouse = "w", i, r) > E.1.u.m(S.uncon, type = "u", spouse = "w", i, r), 1, 0)
    cond <- (cond.h | cond.w)

    if (is.na(cond)){
      res <- rep(NA, 4)
    } else {
      if (!cond) {
        res <- list("status" = "Stay Married with Old Terms", "c.h.uncon" = uncon[1], "c.w.uncon" = uncon[2], "s.uncon" = uncon[3])
        if (Optimize){
          res <- c(uncon[1], uncon[2], uncon[3], 0.5)
        }

      }
      if (cond) {
        if (cond.h){
          lam.seq <- seq(from = 0.99, to = 0.51, length.out = lam.length)
        }
        if(cond.w){
          lam.seq <- seq(from = 0.01, to = 0.49, length.out = lam.length)
        }
        s.seq <- seq(from = S0.lower, to = S0.upper, length.out = S.length )
        reg <- vector()
        temp <- c()
        for (is in 1:S.length){
          for(l in lam.seq){
            con.h <- ifelse(E.1.u.d(s.d, type = "u", spouse = "h", i, r) <= E.1.u.m.lam(l, s.seq[is], type = "u", spouse = "h", i, r), 1, 0)
            con.w <- ifelse(E.1.u.d(s.d, type = "u", spouse = "w", i, r) <= E.1.u.m.lam(l, s.seq[is], type = "u", spouse = "w", i, r), 1, 0)
            con <- con.h & con.w
            if (con) temp <- c( temp, c(s.seq[is], l) )
          }
        }
        if (  is.null(temp) ){
          res <- list("sol" = sol.d, "status" = "Divorce")

          if (Optimize){
            res <- c(sol.d, NA)
          }

        } else {
          sl <-matrix(temp, nrow = 2)
          c <- (y1[i,1,r] + y1[i,2,r] - sl[1, ])/2
          start.points <- cbind(c, c, sl[1, ], sl[2, ])
          values <- mapply(function(j) obj.lam( c.h = start.points[j, 1],
                                                c.w <- start.points[j, 2],
                                                S <- start.points[j, 3],
                                                lam <- start.points[j, 4],  i, r ), 1: nrow(start.points))

          sol <- start.points[which.min(values), ]
          res <- list("S0" = sol[3], "lam0" = sol[4], "c0" = sol[1], "status" = "Stay Married with New Terms")

          if (Optimize){
            res <- sol
          }

          if (Extensive){
            res <- list("S0" = sol[3], "lam0" = sol[4], "c0" = sol[1], "status" = "Stay Married with New Terms",
                        "region" = start.points, "c.h.uncon" = uncon[1], "c.w.uncon" = uncon[2], "s.uncon" = uncon[3])
          }

        }

      }
    }
  }
  return(res)
}







