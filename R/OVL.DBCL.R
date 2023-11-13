#' @title OVL.DBCL
#' @description Parametric approach using the delta method after the Box-Cox transformation taking into account the variability of the estimated transformation parameter
#' @param x controls
#' @param y cases
#' @param alpha confidence level
#' @param h_ini initial value in the optimization problem
#' @return confidence interval
#' @export OVL.DBCL
#' @importFrom stats dnorm optim pnorm qnorm quantile sd var
#' @examples
#' controls = rnorm(50,6,1)
#' cases = rnorm(100,6.5,0.5)
#' OVL.DBCL (controls,cases)
OVL.DBCL<-function(x,y,alpha=0.05,h_ini=-0.6){
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
aux<-xo
xo<-yo
yo<-aux
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
s11<- -n1/sigma1_hat^2
s22<-n1/sigma1_hat^2-3/sigma1_hat^4*sum((muestra1-mu1_hat)^2)
s33<- -n2/sigma2_hat^2
s44<-n2/sigma2_hat^2-3/sigma2_hat^4*sum((muestra2-mu2_hat)^2)
Sigma_aux<-diag(c(s11,s22,s33,s44))
aux1<-1/sigma1_hat^2*sum(xo^hhat*log(xo)/hhat-(xo^hhat-1)/hhat^2)
aux2<-2/sigma1_hat^3*sum((muestra1-mean(muestra1))*(xo^hhat*log(xo)/hhat-(xo^hhat-1)/hhat^2))
aux3<-1/sigma2_hat^2*sum(yo^hhat*log(yo)/hhat-(yo^hhat-1)/hhat^2)
aux4<-2/sigma2_hat^3*sum((muestra2-mean(muestra2))*(yo^hhat*log(yo)/hhat-(yo^hhat-1)/hhat^2))
aux5<-mean(muestra1)/(hhat*sigma1_hat^2)*sum(xo^hhat*log(xo)^2)-1/(hhat^2*sigma1_hat^2)*sum(2*xo^(2*hhat)*log(xo)^2-xo^hhat*log(xo)^2+2*mean(muestra1)*xo^hhat*log(xo))+1/(hhat^3*sigma1_hat^2)*sum(4*xo^(2*hhat)*log(xo)-4*xo^hhat*log(xo)+2*mean(muestra1)*(xo^hhat-1))-3/(hhat^4*sigma1_hat^2)*sum((xo^hhat-1)^2)+mean(muestra2)/(hhat*sigma2_hat^2)*sum(yo^hhat*log(yo)^2)-1/(hhat^2*sigma2_hat^2)*sum(2*yo^(2*hhat)*log(yo)^2-yo^hhat*log(yo)^2+2*mean(muestra2)*yo^hhat*log(yo))+1/(hhat^3*sigma2_hat^2)*sum(4*yo^(2*hhat)*log(yo)-4*yo^hhat*log(yo)+2*mean(muestra2)*(yo^hhat-1))-3/(hhat^4*sigma2_hat^2)*sum((yo^hhat-1)^2)
Sigma_aux2<-rbind(Sigma_aux,c(aux1,aux2,aux3,aux4))
Sigma<-cbind(Sigma_aux2,c(aux1,aux2,aux3,aux4,aux5))
if (abs(det((-1)*Sigma)) >= 1e-25 & abs(det((-1)*Sigma)) <= 1e+25){
Sigma_lambda<-solve((-1)*Sigma)[1:4,1:4]
var_OVL<-t(c(((dnorm(z11)/sigma1_hat-dnorm(z12)/sigma2_hat)*deltax1mu1+(dnorm(z22)/sigma2_hat-dnorm(z21)/sigma1_hat)*deltax2mu1+(dnorm(z21)-dnorm(z11))/sigma1_hat),((dnorm(z11)/sigma1_hat-dnorm(z12)/sigma2_hat)*deltax1sigma1+(dnorm(z22)/sigma2_hat-dnorm(z21)/sigma1_hat)*deltax2sigma1+(dnorm(z21)*z21-dnorm(z11)*z11)/(sigma1_hat)),((dnorm(z11)/sigma1_hat-dnorm(z12)/sigma2_hat)*deltax1mu2+(dnorm(z22)/sigma2_hat-dnorm(z21)/sigma1_hat)*deltax2mu2+(dnorm(z12)-dnorm(z22))/sigma2_hat),((dnorm(z11)/sigma1_hat-dnorm(z12)/sigma2_hat)*deltax1sigma2+(dnorm(z22)/sigma2_hat-dnorm(z21)/sigma1_hat)*deltax2sigma2+(dnorm(z12)*z12-dnorm(z22)*z22)/(sigma2_hat))))%*%Sigma_lambda%*%c(((dnorm(z11)/sigma1_hat-dnorm(z12)/sigma2_hat)*deltax1mu1+(dnorm(z22)/sigma2_hat-dnorm(z21)/sigma1_hat)*deltax2mu1+(dnorm(z21)-dnorm(z11))/sigma1_hat),((dnorm(z11)/sigma1_hat-dnorm(z12)/sigma2_hat)*deltax1sigma1+(dnorm(z22)/sigma2_hat-dnorm(z21)/sigma1_hat)*deltax2sigma1+(dnorm(z21)*z21-dnorm(z11)*z11)/(sigma1_hat)),((dnorm(z11)/sigma1_hat-dnorm(z12)/sigma2_hat)*deltax1mu2+(dnorm(z22)/sigma2_hat-dnorm(z21)/sigma1_hat)*deltax2mu2+(dnorm(z12)-dnorm(z22))/sigma2_hat),((dnorm(z11)/sigma1_hat-dnorm(z12)/sigma2_hat)*deltax1sigma2+(dnorm(z22)/sigma2_hat-dnorm(z21)/sigma1_hat)*deltax2sigma2+(dnorm(z12)*z12-dnorm(z22)*z22)/(sigma2_hat)))
IC1<-OVL-qnorm(1-alpha/2)*sqrt(var_OVL)
IC1_aux<-IC1
if(IC1<0 & (is.na(IC1)==FALSE)){IC1<-0}else{IC1<-IC1_aux}
IC2<-OVL+qnorm(1-alpha/2)*sqrt(var_OVL)
IC2_aux<-IC2
if(IC2>1 & (is.na(IC2)==FALSE)){IC2<-1}else{IC2<-IC2_aux}
return(list(IC1,IC2))
}
else{
return(list(NA,NA))
}
}
