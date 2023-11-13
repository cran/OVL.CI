#' @title Epanechnikov kernel density estimation
#' @description Estimates the density function using the Epanechnikov kernel
#' @param data vector of observations
#' @param points in which the function is evaluated
#' @param h bandwidth
#' @return density estimation
#' @export kernel.e.density
#' @examples
#' x = rnorm(100,1,2)
#' gridd = seq(-5,5,length.out=1000)
#' h = (4/3)^(1/5)*sd(x)*length(x)^(-1/5)
#' kernel.e.density (x,gridd,h)
kernel.e.density <- function(data,points,h){
ndata <- length(data)
npoints <- length(points)
matk <-  kernel.e((points%*%t(rep(1,ndata))-t(data%*%t(rep(1,npoints))))/h)
as.vector((matk%*%rep(1,ndata))/(ndata*h))
}
