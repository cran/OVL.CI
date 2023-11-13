#' @title OVL.K
#' @description Kernel approach estimating the variance via bootstrap
#' @param x controls
#' @param y cases
#' @param alpha confidence level
#' @param B bootstrap size
#' @param k kernel. When k=1 (default value) the kernel used in the estimation is the Gaussian kernel. Otherwise, the Epanechnikov kernel is used instead.
#' @param h bandwidth. When h=1 (default value) the cross-validation bandwidth is chosen. Otherwise, the bandwidth considered by Schmid and Schmidt (2006) is used instead.
#' @return confidence interval
#' @export OVL.K
#' @importFrom stats dnorm optim pnorm qnorm quantile sd var
#' @examples
#' controls = rnorm(50,6,1)
#' cases = rnorm(100,6.5,0.5)
#' OVL.K (controls,cases)
OVL.K<-function(x,y,alpha=0.05,B=100,k=1,h=1){
n1<-length(x)
n2<-length(y)
hx<-ifelse(h==1,ks::hscv(x),(4/3)^(1/5)*sd(x)*n1^-0.2)
hy<-ifelse(h==1,ks::hscv(y),(4/3)^(1/5)*sd(y)*n2^-0.2)
xo<-x
yo<-y
gridxy<-seq(min(c(x,y)),max(c(x,y)),length.out=min(5*(n1+n2),1000))
if(k==1){fkx<-kernel.g.density(x,gridxy,hx)}else{fkx<-kernel.e.density(x,gridxy,hx)}
if(k==1){fky<-kernel.g.density(y,gridxy,hy)}else{fky<-kernel.e.density(y,gridxy,hy)}
ff<-pmin(fkx,fky)
OVL<-(gridxy[2]-gridxy[1])*sum(ff)
OVL_ib<-numeric(B)
for (b in 1:B){
xn<-sample(xo,replace=TRUE)
yn<-sample(yo,replace=TRUE)
gridxy<-seq(min(c(xn,yn)),max(c(xn,yn)),length.out=min(5*(n1+n2),1000))
if(k==1){fkx<-kernel.g.density(x,gridxy,hx)}else{fkx<-kernel.e.density(x,gridxy,hx)}
if(k==1){fky<-kernel.g.density(y,gridxy,hy)}else{fky<-kernel.e.density(y,gridxy,hy)}
ff<-pmin(fkx,fky)
OVL_ib[b]<-(gridxy[2]-gridxy[1])*sum(ff)
}
var_OVL<-var(OVL_ib)
IC1<-OVL-qnorm(1-alpha/2)*sqrt(var_OVL)
IC1_aux<-IC1
if(IC1<0){IC1<-0}else{IC1<-IC1_aux}
IC2<-OVL+qnorm(1-alpha/2)*sqrt(var_OVL)
IC2_aux<-IC2
if(IC2>1){IC2<-1}else{IC2<-IC2_aux}
return(list(IC1,IC2))
}
