#' @title OVL.LogitBCAN
#' @description BCAN procedure carried out in the logit scale and back-transformed
#' @param x controls
#' @param y cases
#' @param alpha confidence level
#' @param B bootstrap size
#' @param h_ini initial value in the optimization problem
#' @return confidence interval
#' @export OVL.LogitBCAN
#' @importFrom stats dnorm optim pnorm qnorm quantile sd var
#' @examples
#' controls = rnorm(50,6,1)
#' cases = rnorm(100,6.5,0.5)
#' OVL.LogitBCAN (controls,cases)
OVL.LogitBCAN<-function(x,y,alpha=0.05,B=100,h_ini=-0.6){
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
OVL_ib<-numeric(B)
for (b in 1:B){
xo_ib<-sample(xo,replace=TRUE)
yo_ib<-sample(yo,replace=TRUE)
hhat<-optim(h_ini,likbox,data=c(xo_ib,yo_ib),n=length(xo_ib),method="BFGS")$par
if (abs(hhat)<1e-5){
    x<-log(xo_ib)
    y<-log(yo_ib)
} else {
    x<-((xo_ib^hhat)-1)/hhat
    y<-((yo_ib^hhat)-1)/hhat
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
OVL_ib[b]<-1+pnorm((x1-mu1_hat)/sigma1_hat)-pnorm((x1-mu2_hat)/sigma2_hat)-pnorm((x2-mu1_hat)/sigma1_hat)+pnorm((x2-mu2_hat)/sigma2_hat)
}
var_OVL<-var(OVL_ib,na.rm=TRUE)
var_OVL_log<-var_OVL/(OVL*(1-OVL))^2
OVL_log<-log(OVL/(1-OVL))
IC1o<-OVL_log-qnorm(1-alpha/2)*sqrt(var_OVL_log)
IC2o<-OVL_log+qnorm(1-alpha/2)*sqrt(var_OVL_log)
IC1<-exp(IC1o)/(1+exp(IC1o))
IC2<-exp(IC2o)/(1+exp(IC2o))
return(list(IC1,IC2))
}

