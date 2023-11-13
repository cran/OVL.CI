#' @title Epanechnikov kernel
#' @description Evaluates the Epanechnikov kernel
#' @param u vector of observations
#' @return evaluation of the Epanechnikov kernel
#' @export kernel.e
#' @examples
#' x = rnorm(100,1,2)
#' kernel.e(x)
kernel.e <- function(u){
3/4*(1-u^2)*(u<1)*(u>-1)
}
