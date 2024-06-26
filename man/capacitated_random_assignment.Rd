% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/immediate_cost.R
\name{capacitated_random_assignment}
\alias{capacitated_random_assignment}
\title{Capacitated Random Assignment}
\usage{
capacitated_random_assignment(model, obj_, rhs_, k, edx, ...)
}
\arguments{
\item{obj_}{Numeric vector, objective coefficients.}

\item{k}{Integer, number of bins.}

\item{n}{Integer, number of items.}

\item{x}{Numeric vector, bin capacities.}
}
\value{
List with 'x' (allocation) and 'objval' (objective value).
}
\description{
Randomly assigns items to bins while respecting capacity constraints.
}
\examples{
# Assign 10 items to 3 bins with capacities 5, 5, and 10 respectively
cap_assignment <- capacitated_random_assignment(10, 3, c(1, 2, 3), c(5, 5, 10))
cap_assignment$x
cap_assignment$objval
}
