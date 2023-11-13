#' @title Gaussian kernel
#' @description Evaluates the Gaussian kernel
#' @param u vector of observations
#' @return evaluation of the Gaussian kernel
#' @export kernel.g
#' @examples
#' x = rnorm(100,1,2)
#' kernel.g(x)
kernel.g <- function(u){
dnorm(u)
}
