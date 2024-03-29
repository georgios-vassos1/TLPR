% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/portf_eval.R
\name{eval_portfolio_t}
\alias{eval_portfolio_t}
\title{Single-stage portfolio contract evaluation with state-dependent constraints}
\usage{
eval_portfolio_t(env, A_, t, S.t, Q.t, ...)
}
\arguments{
\item{env}{The environment containing the parameters of the portfolio contract.}

\item{A_}{A vector containing the values of the assets.}

\item{t}{The time period for which the evaluation is performed.}

\item{S.t}{The state vector at time t.}

\item{Q.t}{The inflow constraint at time t.}

\item{...}{Additional arguments to be passed.}
}
\value{
A list containing the cost and allocation of the portfolio contract for the specified time period.
}
\description{
This function evaluates a single-stage portfolio contract considering state-dependent constraints.
}
\details{
This function evaluates a single-stage portfolio contract considering state-dependent constraints for the specified time period.
It computes the cost and allocation for each asset value provided in the vector \code{A_} based on the given state at time \code{t} and
the associated constraints.
}
