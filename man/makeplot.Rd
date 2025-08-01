% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeplot.R
\name{makeplot}
\alias{makeplot}
\title{Plot facilities for fastCUB objects}
\usage{
makeplot(object)
}
\arguments{
\item{object}{An object of class "fastCUB"}
}
\description{
Plot facilities for objects of class "fastCUB".
}
\details{
Returns a plot comparing fitted
probabilities and observed relative frequencies for GEM models without covariates. If only one
explanatory dichotomous variable is included in the model for one or all components,
then the function returns a plot comparing the distributions of the responses conditioned to
the value of the covariate.
}
\keyword{device}
\keyword{models}
\keyword{package}
