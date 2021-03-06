% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functionsa.R
\name{log_lkd_da}
\alias{log_lkd_da}
\title{Mixture likelihood function with data augmentation}
\usage{
log_lkd_da(y, w, mu, sigma, S)
}
\arguments{
\item{y}{vector of observations.}

\item{w}{vector of weights; they have to sum to one.}

\item{mu}{vector of location parameters for the Gaussian components.}

\item{sigma}{vector of standard deviation parameters for the Gaussian components.}

\item{S}{vector of length equal to the vector y; values in 1,...,K, the total number of components.}
}
\value{
a scalar representing the log-likelihood of the mixture of Gaussian components.
}
\description{
This function compute the log-likelihood of a mixture model for
a fixed number of Gaussian components with data augmentation, i.e.
including a latent variable representing the allocation.
}
\keyword{BayesMIC}
