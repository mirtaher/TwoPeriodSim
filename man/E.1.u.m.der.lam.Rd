% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/E.1.u.m.der.lam.R
\name{E.1.u.m.der.lam}
\alias{E.1.u.m.der.lam}
\title{Derivative of Expected Second Period Marriage Utility Given The First Period Info in the Collective Framework}
\usage{
E.1.u.m.der.lam(lambda, S, var = c("S", "lambda"), type = c("u", "U"),
  spouse = c("h", "w"), i, r)
}
\arguments{
\item{lambda}{Sharing rule}

\item{S}{Saving}

\item{var}{The variable with respect to which the derivative is taken}

\item{type}{Could be individual utility (u) or the joint utility of the couple (U) in a unitary framework}

\item{spouse}{If the type is individual utility (u), we should specify we mean husband or wife utility}

\item{i}{The marriage index}

\item{r}{First period repetition. It is not necessary to be greater than one for the first period. It is needed for taking expectaions, which is required in the second period}
}
\description{
This function calculates the derivative of the expected second period marriage utility with respect to saving and sharing rule (lambda) given the information revealed at the first period for a given level of savings and a new sharing rule (for the second period) for an individual family.
Essentially, it calculates the second period optimal consumption for the given level of saving and updated sharing rule for different
relizations of the shocks in the second period and then substiute them in the individual utility functions (u) or the joint utility function(U).
Finally, it takes the average to return the expectations.To take the derivative, we passed it thuru expectations and then used the second period analytical results to take the derivative.
}

