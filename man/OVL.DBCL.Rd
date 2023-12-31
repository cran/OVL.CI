% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OVL.DBCL.R
\name{OVL.DBCL}
\alias{OVL.DBCL}
\title{OVL.DBCL}
\usage{
OVL.DBCL(x, y, alpha = 0.05, h_ini = -0.6)
}
\arguments{
\item{x}{controls}

\item{y}{cases}

\item{alpha}{confidence level}

\item{h_ini}{initial value in the optimization problem}
}
\value{
confidence interval
}
\description{
Parametric approach using the delta method after the Box-Cox transformation taking into account the variability of the estimated transformation parameter
}
\examples{
controls = rnorm(50,6,1)
cases = rnorm(100,6.5,0.5)
OVL.DBCL (controls,cases)
}
