% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OVL.LogitK.R
\name{OVL.LogitK}
\alias{OVL.LogitK}
\title{OVL.LogitK}
\usage{
OVL.LogitK(x, y, alpha = 0.05, B = 100, k = 1, h = 1)
}
\arguments{
\item{x}{controls}

\item{y}{cases}

\item{alpha}{confidence level}

\item{B}{bootstrap size}

\item{k}{kernel. When k=1 (default value) the kernel used in the estimation is the Gaussian kernel. Otherwise, the Epanechnikov kernel is used instead.}

\item{h}{bandwidth. When h=1 (default value) the cross-validation bandwidth is chosen. Otherwise, the bandwidth considered by Schmid and Schmidt (2006) is used instead.}
}
\value{
confidence interval
}
\description{
Kernel approach estimating the variance via bootstrap in the logit scale and back-transformed
}
\examples{
controls = rnorm(50,6,1)
cases = rnorm(100,6.5,0.5)
OVL.LogitK (controls,cases)
}
