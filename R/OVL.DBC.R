#' @title OVL.DBC
#' @description Parametric approach using the delta method after the Box-Cox transformation
#' @param x controls
#' @param y cases
#' @param alpha confidence level
#' @param h_ini initial value in the optimization problem
#' @return confidence interval
#' @export OVL.DBC
#' @importFrom stats dnorm optim pnorm qnorm quantile sd var
#' @examples
#' controls = rnorm(50,6,1)
#' cases = rnorm(100,6.5,0.5)
#' OVL.DBC (controls,cases)
OVL.DBC<-function(x,y,alpha=0.05,h_ini=-0.6){
x_aux<-x
y_aux<-y
all_values<-c(x_aux,y_aux)
if (any(all_values<=0)){
  x<-x_aux+abs(min(all_values))+(max(all_values)-min(all_values))/2
  y<-y_aux+abs(min(all_values))+(max(all_values)-min(all_values))/2
} else {
  x<-x_aux
  y<-y_aux
}
xo<-x
yo<-y
hhat<-optim(h_ini,likbox,data=c(xo,yo),n=length(xo),method="BFGS")$par
if (abs(hhat)<1e-5){
    x<-log(xo)
    y<-log(yo)
} else {
    x<-((xo^hhat)-1)/hhat
    y<-((yo^hhat)-1)/hhat
}
if(ssdd(x)<ssdd(y)){
muestra1<-x
muestra2<-y
} else {
muestra1<-y
muestra2<-x
}
mu1_hat<-mean(muestra1)
mu2_hat<-mean(muestra2)
sigma1_hat<-ssdd(muestra1)
sigma2_hat<-ssdd(muestra2)
x1<-((mu1_hat*sigma2_hat^2-mu2_hat*sigma1_hat^2)-sigma1_hat*sigma2_hat*sqrt((mu1_hat-mu2_hat)^2+(sigma1_hat^2-sigma2_hat^2)*log(sigma1_hat^2/sigma2_hat^2)))/(sigma2_hat^2-sigma1_hat^2)
x2<-((mu1_hat*sigma2_hat^2-mu2_hat*sigma1_hat^2)+sigma1_hat*sigma2_hat*sqrt((mu1_hat-mu2_hat)^2+(sigma1_hat^2-sigma2_hat^2)*log(sigma1_hat^2/sigma2_hat^2)))/(sigma2_hat^2-sigma1_hat^2)
OVL<-1+pnorm((x1-mu1_hat)/sigma1_hat)-pnorm((x1-mu2_hat)/sigma2_hat)-pnorm((x2-mu1_hat)/sigma1_hat)+pnorm((x2-mu2_hat)/sigma2_hat)
n1<-length(muestra1)
n2<-length(muestra2)
z11<-(x1-mu1_hat)/sigma1_hat
z12<-(x1-mu2_hat)/sigma2_hat
z21<-(x2-mu1_hat)/sigma1_hat
z22<-(x2-mu2_hat)/sigma2_hat
deltax1mu1<-(sigma2_hat^2-sigma1_hat*sigma2_hat*(mu1_hat-mu2_hat)*U(mu1_hat,mu2_hat,sigma1_hat,sigma2_hat)^(-1))/(sigma2_hat^2-sigma1_hat^2)
deltax2mu1<-(sigma2_hat^2+sigma1_hat*sigma2_hat*(mu1_hat-mu2_hat)*U(mu1_hat,mu2_hat,sigma1_hat,sigma2_hat)^(-1))/(sigma2_hat^2-sigma1_hat^2)
deltax1mu2<-(-sigma1_hat^2+sigma1_hat*sigma2_hat*(mu1_hat-mu2_hat)*U(mu1_hat,mu2_hat,sigma1_hat,sigma2_hat)^(-1))/(sigma2_hat^2-sigma1_hat^2)
deltax2mu2<-(-sigma1_hat^2-sigma1_hat*sigma2_hat*(mu1_hat-mu2_hat)*U(mu1_hat,mu2_hat,sigma1_hat,sigma2_hat)^(-1))/(sigma2_hat^2-sigma1_hat^2)
deltax1sigma1<-1/(sigma2_hat^2-sigma1_hat^2)*(-2*sigma1_hat*mu2_hat-sigma2_hat*U(mu1_hat,mu2_hat,sigma1_hat,sigma2_hat)+sigma1_hat*sigma2_hat*((sigma2_hat^2-sigma1_hat^2)/sigma1_hat+sigma1_hat*log(sigma2_hat^2/sigma1_hat^2))*U(mu1_hat,mu2_hat,sigma1_hat,sigma2_hat)^(-1)+2*sigma1_hat*x1)
deltax2sigma1<-1/(sigma2_hat^2-sigma1_hat^2)*(-2*sigma1_hat*mu2_hat+sigma2_hat*U(mu1_hat,mu2_hat,sigma1_hat,sigma2_hat)-sigma1_hat*sigma2_hat*((sigma2_hat^2-sigma1_hat^2)/sigma1_hat+sigma1_hat*log(sigma2_hat^2/sigma1_hat^2))*U(mu1_hat,mu2_hat,sigma1_hat,sigma2_hat)^(-1)+2*sigma1_hat*x2)
deltax1sigma2<-1/(sigma2_hat^2-sigma1_hat^2)*(2*sigma2_hat*mu1_hat-sigma1_hat*U(mu1_hat,mu2_hat,sigma1_hat,sigma2_hat)-sigma1_hat*sigma2_hat*((sigma2_hat^2-sigma1_hat^2)/sigma2_hat+sigma2_hat*log(sigma2_hat^2/sigma1_hat^2))*U(mu1_hat,mu2_hat,sigma1_hat,sigma2_hat)^(-1)-2*sigma2_hat*x1)
deltax2sigma2<-1/(sigma2_hat^2-sigma1_hat^2)*(2*sigma2_hat*mu1_hat+sigma1_hat*U(mu1_hat,mu2_hat,sigma1_hat,sigma2_hat)+sigma1_hat*sigma2_hat*((sigma2_hat^2-sigma1_hat^2)/sigma2_hat+sigma2_hat*log(sigma2_hat^2/sigma1_hat^2))*U(mu1_hat,mu2_hat,sigma1_hat,sigma2_hat)^(-1)-2*sigma2_hat*x2)
var_OVL<-((dnorm(z11)/sigma1_hat-dnorm(z12)/sigma2_hat)*deltax1mu1+(dnorm(z22)/sigma2_hat-dnorm(z21)/sigma1_hat)*deltax2mu1+(dnorm(z21)-dnorm(z11))/sigma1_hat)^2*sigma1_hat^2/n1+((dnorm(z11)/sigma1_hat-dnorm(z12)/sigma2_hat)*deltax1mu2+(dnorm(z22)/sigma2_hat-dnorm(z21)/sigma1_hat)*deltax2mu2+(dnorm(z12)-dnorm(z22))/sigma2_hat)^2*sigma2_hat^2/n2+((dnorm(z11)/sigma1_hat-dnorm(z12)/sigma2_hat)*deltax1sigma1+(dnorm(z22)/sigma2_hat-dnorm(z21)/sigma1_hat)*deltax2sigma1+(dnorm(z21)*z21-dnorm(z11)*z11)/(sigma1_hat))^2*(n1-1)*sigma1_hat^2/(2*n1^2)+((dnorm(z11)/sigma1_hat-dnorm(z12)/sigma2_hat)*deltax1sigma2+(dnorm(z22)/sigma2_hat-dnorm(z21)/sigma1_hat)*deltax2sigma2+(dnorm(z12)*z12-dnorm(z22)*z22)/(sigma2_hat))^2*(n2-1)*sigma2_hat^2/(2*n2^2)
IC1<-OVL-qnorm(1-alpha/2)*sqrt(var_OVL)
IC1_aux<-IC1
if(IC1<0 & (is.na(IC1)==FALSE)){IC1<-0}else{IC1<-IC1_aux}
IC2<-OVL+qnorm(1-alpha/2)*sqrt(var_OVL)
IC2_aux<-IC2
if(IC2>1 & (is.na(IC2)==FALSE)){IC2<-1}else{IC2<-IC2_aux}
return(list(IC1,IC2))
}
