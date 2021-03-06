% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svr-functions.R
\name{gridSearch}
\alias{gridSearch}
\title{Function to perform grid search}
\usage{
gridSearch(FUN, init = 0, target, increment = 2, tol = 0.0025,
  tol.same = 1e-06, left = FALSE, iter.max = 25, ...)
}
\arguments{
\item{FUN}{Scalar function of scalar argument.}

\item{init}{Scalar, initial solution.}

\item{target}{Scalar, target p-value.}

\item{increment}{Scalar, determines the first increment of the
grid, and is refined over the iterations of the grid search.}

\item{tol}{Scalar, tolerance for determining when \code{FUN}
evaluated at the solution is sufficiently close to
\code{target}. Set to \code{2.5e-3} by default.}

\item{tol.same}{Scalar, tolerance for determining when solution can
no longer be improved. Set to \code{1e-06} by default.}

\item{left}{Boolean, set to \code{TRUE} if solution is left of the
initial solution declared in \code{init}.}

\item{iter.max}{Integer, maximum number of iterations to perform.}
}
\value{
Scalar, the value of the argument at which \code{FUN} is
    equal to \code{target}.
}
\description{
Perform a grid search to find the value at which argument
\code{FUN} is equal to argument \code{target}.
}
