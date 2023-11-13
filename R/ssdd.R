#' @title Sample variance computation
#' @description Computes the sample variance of a vector of observations
#' @param x vector of observations
#' @return the sample variance
#' @export ssdd
#' @examples
#' x = rnorm(100,1,2)
#' ssdd(x)
ssdd <- function(x){
v <- sqrt(var(x)*(length(x)-1)/length(x))
return(v)
}


