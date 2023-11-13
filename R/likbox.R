#' @title Likelihood function of the BoxCox transformation
#' @description Computation of the likelihood function of the BoxCox transformation
#' @param h parameter of the Box-Cox transformation
#' @param data joint vector of controls (first) and cases
#' @param n length of the vector of controls
#' @return the likelihood function of the BoxCox transformation
#' @export likbox
#' @examples
#' h=-1.6
#' controls=rnorm(50,6,1)
#' cases=rnorm(100,6.5,0.5)
#' likbox(h,c(controls,cases),n=length(controls))
likbox<-function(h,data,n){
m <-length(data)-n
x<-data[1:n]
y<-data[(n+1):length(data)]
if (abs(h)<1e-5){
    xh<-log(x)
    yh<-log(y)
} else {
    xh<-((x^h)-1)/h
    yh<-((y^h)-1)/h
}
oout <- -n/2*log(sum((xh-sum(xh)/n)^2)/n)-m/2*log(sum((yh-sum(yh)/m)^2)/m)+(h-1)*(sum(log(x))+sum(log(y)))
return(-oout)
}

