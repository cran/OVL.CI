#' @title Auxiliary function
#' @description Evaluates an auxiliary function
#' @param mu1 sample mean of a vector x
#' @param mu2 sample mean of a vector y
#' @param sigma1 sample standard deviation of a vector x
#' @param sigma2 sample standard deviation of a vector y
#' @return evaluation of an auxiliary function
#' @export U
#' @examples
#' U(1,2,1,1)
U<-function(mu1,mu2,sigma1,sigma2){
vv<-sqrt((mu1-mu2)^2+(sigma2^2-sigma1^2)*log(sigma2^2/sigma1^2))
return(vv)
}
