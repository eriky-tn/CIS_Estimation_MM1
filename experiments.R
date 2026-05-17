################################################################################
# Computes the numerical experiments presented in:
#
# Saroja Kumar Singh; Abhijit Datta Banik; Eriky S. Gomes & Frederico R. B. Cruz
# (2020). Estimating Traffic Intensity for Single Server Markovian Queueing 
# System, T.B.A (submited).
#
# Programmed by:
#
# Eriky S. Gomes & Frederico R. B. Cruz
# Universidade Federal de Minas Gerais
# E-mail: fcruz@est.ufmg.br
# ? 2020 Gomes & Cruz
# v.2021.10.18
install.packages("xlsx")
library("xlsx")
################################################################################
# number of customers distribution
################################################################################
Pm<-function(m,rho){
  prob<-1/(1+rho)*(rho/(1+rho))^m
  return(prob)}
Pm(m=0:5,rho=0.5)
################################################################################
# maximum likelihood estimate for rho (eq.4) 14/01/2021
################################################################################
MM1RhoMLEf<-function(samp) {
  n<-length(samp)
  y<-sum(samp)
  # cat("n=",n,"\t","y=",y,"\n")
  res<-y/n
  # return(ifelse(res<1,res,0.99))
  return(res)}
################################################################################
# maximum likelihood estimate for Lq (eq.7) 14/01/2021
################################################################################
MM1LqMLEf<-function(samp){
  # if(MM1RhoMLEf(samp)>=1) return(Inf)
  rhoMLE<-MM1RhoMLEf(samp)
  LqMLE<-(rhoMLE^2/(1-rhoMLE))
  return(LqMLE)}
################################################################################
# maximum likelihood estimate for Ls (eq.8) 14/01/2021
################################################################################
MM1LsMLEf<-function(samp){
  # if(MM1RhoMLEf(samp)>=1) return(Inf)
  rhoMLE<-MM1RhoMLEf(samp)
  LsMLE<-(rhoMLE/(1-rhoMLE))
  return(LsMLE)}
################################################################################
# Inverse functions for Lq and Ls
################################################################################
InvLq<-function(Lq){
  return((-Lq+sqrt(Lq^2+4*Lq))/2)}
InvLs<-function(Ls){
  return(Ls/(1+Ls))}
################################################################################
# testing MLE
################################################################################
set.seed(13579)
rho=0.22
vet<-numeric(1000)
for(i in 1:1000){
smm1<-rgeom(n=10,prob=1/(1+rho))
vet[i]<-MM1RhoMLEf(smm1)}
table(vet)
set.seed(13579)
rho=0.99
smm1<-rgeom(n=2000,prob=1/(1+rho))
MM1RhoMLEf(smm1)
set.seed(13579)
rho=0.1
Lq=2; rhoLq<-InvLq(Lq)
Ls=1; rhoLs<-InvLs(Ls)
smm1<-rgeom(n=2000,prob=1/(1+rho))
smm2<-rgeom(n=2000,prob=1/(1+rhoLq))
smm3<-rgeom(n=2000,prob=1/(1+rhoLs))
MM1RhoMLEf(smm1)
MM1LqMLEf(smm2)
MM1LsMLEf(smm3)
################################################################################
# Gauss hypergeometric function (eq.12+)
################################################################################
GaussHypf<-function(a,b,c,z) {
  #  cat("GaussHypf(",a,",",b,",",c,",",z,"):\n")
  GaussHypItg<-function(u,a,b,c,z) {
    U<-u^(b-1)*(1-u)^(c-b-1)*(1-u*z)^(-a)
    return(U)}
  GaussHyp1F2<-integrate(GaussHypItg,0.0,1.0,a,b,c,z)[[1]]/beta(b,c-b)
  #  cat("GaussHyp1F2=",GaussHyp1F2,"\n")
  return(GaussHyp1F2)}
################################################################################
# testing Gauss hypergeometric function (eq.12)
################################################################################
set.seed(13579)
rho=0.50
smm1<-rgeom(n=100,prob=1/(1+rho))
y<-sum(smm1)
n<-length(smm1)
alphap<-1.0
betap<-1.0
GaussHypf(n+y,y+alphap,y+alphap+betap,-1)
################################################################################
# beta prior - Eq. (11)
################################################################################
p1<-function(p,a,b){
  num<-function(u,a,b){
    U<-u^(a-1)*(1-u)^(b-1)
    return(U)}
  p1<-num(p,a,b)/beta(a,b)
  return(p1)}
################################################################################
# plotting Beta prior
################################################################################
setEPS()
postscript(paste("FiBe.eps",sep=""),width=10.5*0.75,height=8*0.75)
#dev.new(width=10.5*0.75,height=8*0.75)
par(mfrow=c(1,1))
rho<-seq(0.0,1.0,0.05)
plot(rho,p1(rho,a=2.5,b=1.5),type="n",xlab=expression(rho),
     ylab=expression(p[1](rho)))
lines(rho,p1(rho,a=2.5,b=1.5),type="b",lty=1,pch=1,col=1,lwd=2)
lines(rho,p1(rho,a=1.0,b=1.1),type="b",lty=2,pch=2,col=2,lwd=2)
lines(rho,p1(rho,a=2.0,b=2.1),type="b",lty=3,pch=3,col=3,lwd=2)
lines(rho,p1(rho,a=1.5,b=2.5),type="b",lty=4,pch=4,col=4,lwd=2)
legend("bottom",lty=c(1,2,3,4),pch=c(1,2,3,4),col=c(1,2,3,4),lwd=2,
       legend=c(
         expression(paste("Beta(",alpha,"=2.5; ",beta,"=1.5)")),
         expression(paste("Beta(",alpha,"=1.0; ",beta,"=1.1)")),
         expression(paste("Beta(",alpha,"=2.0; ",beta,"=2.1)")),
         expression(paste("Beta(",alpha,"=1.5; ",beta,"=2.5)"))))
graphics.off()
################################################################################
# beta posterior - Eq. (12) - version 2
################################################################################
q1<-function(p,x,a,b){
  n=length(x); y=sum(x)
  num<-function(u,a,b,n,y){
    U<-u^(y+a-1)*(1-u)^(b-1)*(1+u)^(-n-y)
    return(U)}
  #  a1<-num(p,a,b,n,y)
  #  a2<-p1(p,a,b)*p^y*(1+p)^(-n-y)
  #  cat("a1=",a1,"\ta2=",a2,"\n")
  #  a3<-beta(y+a,b)/GaussHypf(n+y,y+a,y+a+b,-1)
  #  a4<-integrate(num,0,1,a,b,n,y)[[1]]
  #  cat("a3=",a3,"\ta4=",a4,"\n")
  q1<-num(p,a,b,n,y)/beta(y+a,b)/GaussHypf(n+y,y+a,y+a+b,-1)
  return(q1)}
################################################################################
# testing beta posterior
################################################################################
set.seed(13579)
rho=0.50
smm1<-rgeom(n=100,prob=1/(1+rho))
y<-sum(smm1)
n<-length(smm1)
alphap<-1.0
betap<-1.0
q1(rho,smm1,alphap,betap)-
  p1(rho,alphap,betap)*rho^y*(1+rho)^(-n-y)/beta(y+alphap,betap)/
  GaussHypf(n+y,y+alphap,y+alphap+betap,-1)
################################################################################
# plotting posterior + prior
# notice that posterior concentrates around rho
################################################################################
par(mfrow=c(1,1))
rho<-seq(0.01,0.99,0.01)
plot(rho,q1(rho,x=smm1,a=alphap,b=betap),xlab=expression(rho),
     ylab=expression(q[1](rho)),type="l",lty=1,pch=1,col=1,lwd=2)
lines(rho,p1(rho,a=alphap,b=betap),type="l",lty=2,pch=2,col=2,lwd=2)
################################################################################
# beta: loss functions 
################################################################################
# SELF-------------------------------------------------------------------------
# estimators and functions
################################################################################
BSelfRhof<-function(x,a,b){
  aux<-function(p,x,a,b){
    U<-p*q1(p,x,a,b)
    return(U)}
  res<-integrate(aux,0,1,x,a,b)[[1]]
  return(res)} # rho estimate by definition
BSelfRhoEf<-function(x,a,b){
  y<-sum(x) 
  n<-length(x)
  res<-((y+a)*GaussHypf(n+y,y+a+1,y+a+b+1,-1))/
    ((y+a+b)*(GaussHypf(n+y,y+a,y+a+b,-1)))
  return(res)} # rho estimator by Eq. (13) 
BSelfLqf<-function(x,a,b){
  aux<-function(rho,x,a,b){
    U<-rho^2*(1-rho)^(-1)*q1(rho,x,a,b)
    return(U)}
  Lq=integrate(aux,0,0.99,x,a,b)[[1]]
  return(Lq)} # Lq estimate by definition  
BSelfLqEf<-function(x,a,b){ # only for b>1
  y=sum(x)
  n=length(x)
  Lq=beta(y+a+2,b-1)*GaussHypf(n+y,y+a+2,y+a+b+1,-1)/
    beta(y+a,b)/GaussHypf(n+y,y+a,y+a+b,-1)
  return(Lq)} # Lq estimator by Eq. (14)
BSelfLsf<-function(x,a,b){
  aux<-function(rho,x,a,b){
    U<-rho*(1-rho)^(-1)*q1(rho,x,a,b)
    return(U)}
  Lq=integrate(aux,0,1-1E-10,x,a,b)[[1]]
  return(Lq)} # Lq estimate by definition  
BSelfLsEf<-function(x,a,b){
  y=sum(x)
  n=length(x)
  Ls=beta(y+a+1,b-1)/beta(y+a,b)*GaussHypf(n+y,y+a+1,y+a+b,-1)/
    GaussHypf(n+y,y+a,y+a+b,-1)
  return(Ls)} # Ls estimator by Eq. (15)
################################################################################
# testing beta SELF estimators and risk
################################################################################
set.seed(13579)
rho=0.5
smm1<-rgeom(n=100,prob=1/(1+rho))
y<-sum(smm1)
n<-length(smm1)
a<-1.0  #alphap
b<-3.0 #betap
RhoE<-BSelfRhoEf(smm1,a,b)  # (eq.13)
RhoE-BSelfRhof(smm1,a,b)    # checking by definition
LqE<-BSelfLqEf(smm1,a,b)    # (eq.14)
LqE-BSelfLqf(smm1,a,b)      # checking by definition
LsE<-BSelfLsEf(smm1,a,b)    # (eq.15)
LsE-BSelfLsf(smm1,a,b)      # checking by definition
################################################################################
# PLF---------------------------------------------------------------------------
# estimators and functions 
################################################################################
BPLFRhof<-function(x,a,b){
  aux<-function(rho,x,a,b){
    U<-rho^2*q1(rho,x,a,b)
    return(U)}
  res<-sqrt(integrate(aux,0,1,x,a,b)[[1]])
  return(res)} # rho estimate by defintion
BPLFRhoEf<-function(x,a,b){
  y=sum(x)
  n=length(x)
  RhoE<-beta(y+a+2,b)/beta(y+a,b)*GaussHypf(n+y,y+a+2,y+a+b+2,-1)/
    GaussHypf(n+y,y+a,y+a+b,-1)
  return(sqrt(RhoE))}  # rho estimator by Eq. (16)
BPLFLqf<-function(x,a,b){
  aux<-function(rho,x,a,b){
    U<-rho^4*(1-rho)^(-2)*q1(rho,x,a,b)
    return(U)}
  res<-sqrt(integrate(aux,0,1-1E-5,x,a,b)[[1]])
  return(res)} # Lq estimate by definition
BPLFLqEf<-function(x,a,b){
  y=sum(x)
  n=length(x)
  res<-beta(y+a+4,b-2)*GaussHypf(n+y,y+a+4,y+a+b+2,-1)/
    beta(y+a,b)/GaussHypf(n+y,y+a,y+a+b,-1)
  res<-sqrt(res)
  return(res)} # Lq estimator by Eq. (17)
BPLFLsf<-function(x,a,b){
  aux<-function(rho,x,a,b){
    U<-rho^2*((1-rho)^(-2))*q1(rho,x,a,b)
    return(U)}
  res<-integrate(aux,0+1E-10,1-1e-10,x,a,b)[[1]]
  res<-sqrt(res)
  return(res)} # Ls estimate by definition
BPLFLsEf<-function(x,a,b){
  y=sum(x)
  n=length(x)
  res<-beta(y+a+2,b-2)*GaussHypf(n+y,y+a+2,y+a+b,-1)/
    beta(y+a,b)/GaussHypf(n+y,y+a,y+a+b,-1)
  res<-sqrt(res)
  return(res)} # Ls estimator by Eq. (18)
################################################################################
# testing Beta PLF estimators and risk
################################################################################
set.seed(13579)
rho=0.5
smm1<-rgeom(n=100,prob=1/(1+rho))
a<-1 # alphap
b<-3 # betap
RhoE<-BPLFRhoEf(smm1,a,b) # Eq. (16)
RhoE-BPLFRhof(smm1,a,b)   # checking by definition
LqE<-BPLFLqEf(smm1,a,b)   # Eq. (17)
LqE-BPLFLqf(smm1,a,b)     # checking by definition
LsE<-BPLFLsEf(smm1,a,b)   # Eq. (18)
LsE-BPLFLsf(smm1,a,b)     # checking by definition
################################################################################
#  inverted beta distribution (eq.19)
################################################################################
dIBf<-function(rho,a,b) {
  dIB<-rho^(a-1)*(1+rho)^(-b-a)/beta(a,b)
  #  cat("dIB=",dIB,"\n")
  return(dIB)}
################################################################################
# regularized incomplete beta function (eq.20+)
################################################################################
Ixi<-function(a,b,xi) {
  IepsItg<-function(u,a,b) {
    U<-u^(a-1)*(1-u)^(b-1)
    return(U)}
  Ixi<-integrate(IepsItg,0,xi,a,b)[[1]]/beta(a,b)
  #  cat("Ixi=",Ixi,"\n")
  return(Ixi)}
################################################################################
# testing inverted beta distribution and
# regularized incomplete beta function (eq.20)
################################################################################
(y + ?? + j + 2, n + ?? ??? j ??? 2)
set.seed(13579)
rho<-0.5
a<-2.0
b<-2.1
xi<-1
j=0
n<-500
smm1<-rgeom(n,prob=1/(1+rho))
y<-sum(smm1)
Ixi(y+a+j+2,n+b-j-2,xi/(1+xi))
integrate(dIBf,0,xi,1,n+b-j-2)[[1]]-Ixi(y+a+j+2,n+b-j-2,xi/(1+xi))

################################################################################
# incomplete inverted beta prior (eq.21)
################################################################################
p2<-function(rho,a,b){
  num<-rho^(a-1)*(1+rho)^(-b-a)
  p2<-num/beta(a,b)/Ixi(a,b,0.5)
  return(p2)}
################################################################################
# testing incomplete inverted beta prior
################################################################################
integrate(p2,0,1,a=1,b=1)
integrate(p2,0,1,a=1,b=2)
integrate(p2,0,1,a=1,b=4)
integrate(p2,0,1,a=2,b=1)
integrate(p2,0,1,a=2,b=2)
integrate(p2,0,1,a=2,b=4)
integrate(p2,0,1,a=4,b=1)
integrate(p2,0,1,a=4,b=2)
integrate(p2,0,1,a=4,b=4)
################################################################################
# plotting incomplete inverted beta prior
################################################################################
setEPS()
postscript(paste("FiII.eps",sep=""),width=10.5*0.75,height=8*0.75)
#dev.new(width=10.5*0.75,height=8*0.75)
par(mfrow=c(1,1))
rho<-seq(0.0,1.0,0.05)
plot(rho,p2(rho,a=1.1,b=3.0),type="n",xlab=expression(rho),
     ylab=expression(p[2](rho)),ylim=c(0,3))
lines(rho,p2(rho,a=1.0,b=1.1),type="b",lty=1,pch=1,col=1,lwd=2)
lines(rho,p2(rho,a=1.0,b=2.1),type="b",lty=1,pch=2,col=1,lwd=2)
lines(rho,p2(rho,a=1.0,b=3.0),type="b",lty=1,pch=3,col=1,lwd=2)
lines(rho,p2(rho,a=2.0,b=1.1),type="b",lty=2,pch=4,col=2,lwd=2)
lines(rho,p2(rho,a=2.0,b=2.1),type="b",lty=2,pch=5,col=2,lwd=2)
lines(rho,p2(rho,a=2.0,b=3.0),type="b",lty=2,pch=6,col=2,lwd=2)
lines(rho,p2(rho,a=3.0,b=1.1),type="b",lty=3,pch=7,col=3,lwd=2)
lines(rho,p2(rho,a=3.0,b=2.1),type="b",lty=3,pch=8,col=3,lwd=2)
lines(rho,p2(rho,a=3.0,b=3.0),type="b",lty=3,pch=9,col=3,lwd=2)
legend("top",lty=c(1,1,1,2,2,2,3,3,3),pch=1:9,
       col=c(1,1,1,2,2,2,3,3,3),lwd=2,
       legend=c(
         expression(paste("IIB(",alpha,"=1.0; ",beta,"=1.1)")),
         expression(paste("IIB(",alpha,"=1.0; ",beta,"=2.1)")),
         expression(paste("IIB(",alpha,"=1.0; ",beta,"=3.0)")),
         expression(paste("IIB(",alpha,"=2.0; ",beta,"=1.1)")),
         expression(paste("IIB(",alpha,"=2.0; ",beta,"=2.1)")),
         expression(paste("IIB(",alpha,"=2.0; ",beta,"=3.0)")),
         expression(paste("IIB(",alpha,"=3.0; ",beta,"=1.1)")),
         expression(paste("IIB(",alpha,"=3.0; ",beta,"=2.1)")),
         expression(paste("IIB(",alpha,"=3.0; ",beta,"=3.0)"))
       ))
graphics.off()
################################################################################
# inverted beta posterior (eq.22)
################################################################################
q2<-function(rho,x,a,b){
  y<-sum(x)
  n<-length(x)
  num<-rho^(y+a-1)*(1+rho)^(-n-y-b-a)
  den<-beta(y+a,n+b)*Ixi(y+a,n+b,0.5)
  q2<-num/den
  #  q22<-p2(rho,a+y,b+n)
  #  cat("q2=",q2,"q22=",q22,"\n")
  return(q2)}
################################################################################
# testing inverted beta posterior
################################################################################
set.seed(13579)
rho=0.50
smm1<-rgeom(n=100,prob=1/(1+rho))
a<-1.0
b<-1.5
q2(rho,smm1,a,b)
integrate(q2,0,1,x=smm1,a,b)[[1]]
################################################################################
# plotting posterior + prior
# notice that posterior concentrates around rho
################################################################################
par(mfrow=c(1,1))
rho<-seq(0.01,0.99,0.01)
plot(rho,q2(rho,x=smm1,a,b),xlab=expression(rho),
     ylab=expression(q[1](rho)),type="l",lty=1,pch=1,col=1,lwd=2)
lines(rho,p2(rho,a,b),type="l",lty=2,pch=2,col=2,lwd=2)
legend("topright",lty=c(1,2),col=c(1,2),lwd=2,
       legend=c("Posterior","Prior"))
################################################################################
# testing equivalence between
# integrate(itg) and B(y+a+j+1,b+n-j-1)*Ixi(a+y+j+1,b+n-j-1)
################################################################################
set.seed(13579)
rho=0.50
smm1<-rgeom(n=100,prob=1/(1+rho))
y<-sum(smm1)
n<-length(smm1)
a<-1.1
b<-1.1
j<-100
itg<-function(rho,a,b){
  U<-rho^(y+a+j+1-1)*(1+rho)^(-n-y-b-a)
  #  cat("U=",U,"\n")
  return(U)}
integrate(itg,0,1,a,b)[[1]]
beta(y+a+j+1,n+b-j-1)*Ixi(y+a+j+1,n+b-j-1,0.5)
################################################################################
# inverted beta: loss functions
################################################################################
# inverted: SELF ---------------------------------------------------------------
# estimators and functions
################################################################################
IBSelfRhof<-function(x,a,b){
  aux<-function(rho,x,a,b){
    U<-rho*q2(rho,x,a,b)
    return(U)}
  res<-integrate(aux,0,1,x,a,b)[[1]]
  return(res)} # rho estimate by definition
IBSelfRhoEf<-function(x,a,b){
  y<-sum(x)
  n<-length(x)
  RhoE<-beta(y+a+1,n+b-1)*Ixi(y+a+1,n+b-1,0.5)/
    beta(y+a,n+b)/Ixi(y+a,n+b,0.5)
  return(RhoE)} # rho estimator by Eq. (23)
IBSelfRhoER<-function(rhoE,x,a,b){
  y<-sum(x)
  n<-length(x)
  num<-beta(y+a+2,n+b-2)*Ixi(y+a+2,n+b-2,0.5)
  den<-beta(y+a,n+b)*Ixi(y+a,n+b,0.5)
  return(num/den-rhoE^2)} # risk
IBSelfLqf<-function(x,a,b){
  aux<-function(rho,x,a,b){
    U<-rho^2*(1-rho)^(-1)*q2(rho,x,a,b)
    return(U)}
  res<-integrate(aux,0,1-1e-10,x,a,b)[[1]]
  return(res)} # Lq estimate by definition
IBSelfLqEf<-function(x,a,b){
  y<-sum(x)
  n<-length(x)
  den<-beta(y+a,n+b)*Ixi(y+a,n+b,0.5)
  num<-0
  j<-0
  #changes: n=10
  while ((10+b-j-2)>0){
    term<-beta(y+a+j+2,n+b-j-2)*Ixi(y+a+j+2,n+b-j-2,0.5)
#   cat("term=",term,"\n")
    num<-num+term
    j<-j+1}
  return(num/den)} # Lq estimator by Eq. (24)
IBSelfLsf<-function(x,a,b){
  aux<-function(rho,x,a,b){
    U<-rho*(1-rho)^(-1)*q2(rho,x,a,b)
    return(U)}
  res<-integrate(aux,0,1-1e-10,x,a,b)[[1]]
  return(res)} # Ls estimate by definition
IBSelfLsEf<-function(x,a,b){
  y<-sum(x)
  n<-length(x)
  den<-beta(y+a,n+b)*Ixi(y+a,n+b,0.5)
  num<-0
  j<-0
  # changes 
  while ((10+b-j-1)>0){
    term<-beta(y+a+j+1,n+b-j-1)*Ixi(y+a+j+1,n+b-j-1,0.5)
    #    cat("term=",term,"\n")
    num<-num+term
    j<-j+1}
  return(num/den)} # Lq estimator by Eq. (25)
################################################################################
# testing inverted beta SELF estimators and risk
################################################################################
set.seed(13579)
rho=0.8
smm1<-rgeom(n=100,prob=1/(1+rho))
a<-1.1
b<-1.1
RhoE<-IBSelfRhoEf(smm1,a,b) # Eq. (23)
RhoE-IBSelfRhof(smm1,a,b)   # checking by definition
IBSelfRhoER(RhoE,smm1,a,b)  # risk
LqE<-IBSelfLqEf(smm1,a,b)   # Eq. (24)
LqE-IBSelfLqf(smm1,a,b)     # checking by definition
LsE<-IBSelfLsEf(smm1,a,b)   # Eq. (25)
LsE-IBSelfLsf(smm1,a,b)     # checking by definition

f1<-function(x){
  return(1/sqrt(x))
}

f2<-function(x){
  return(1/x)
}

integrate(f1,0,1)
integrate(f2,0,1)
################################################################################
# inverted: PLF ----------------------------------------------------------------
# estimators and functions 
################################################################################
# PLF --------------------------------------------------------------------------
# estimators and functions 
IBPlfRhof<-function(x,a,b){
  aux<-function(rho,x,a,b){
    U<-rho^2*q2(rho,x,a,b-a)
    return(U)}
  res<-integrate(aux,0,1,x,a,b)[[1]]
  return(sqrt(res))}  # rho estimate by definition
IBPlfRhoEf<-function(x,a,b){
  n<-length(x)
  y<-sum(x)
  rhoE<-beta(y+a+2,n+b-a-2)*Ixi(y+a+2,n+b-a-2,0.5)/
    beta(y+a,n+b-a)/Ixi(y+a,n+b-a,0.5)
  return(sqrt(rhoE))} # rho estimator by Eq. (26)
IBPlfLqf<-function(x,a,b){
  aux<-function(rho,x,a,b){
    U<-rho^4/(1-rho)^2*q2(rho,x,a,b)
    return(U)}
  res<-integrate(aux,0,0.99,x,a,b)[[1]]
  return(sqrt(res)) }  # Lq estimate by definition
IBPlfLqEf<-function(x,a,b){
  y<-sum(x)
  n<-length(x)
  num<-0
  den<-beta(y+a,n+b)*Ixi(y+a,n+b,0.5)
  j<-0
  #changes
  while(10+b-j-4>0){
    term<-(j+1)*beta(y+a+j+4,n+b-j-4)*Ixi(y+a+j+4,n+b-j-4,0.5)
    num<-num+term
    j<-j+1}
  return(sqrt(num/den))}    # Lq estimator by Eq. (27)
IBPlfLsf<-function(x,a,b){
  aux<-function(rho,x,a,b){
    U<-rho^2/(1-rho)^2*q2(rho,x,a,b)
    return(U)}
  res<-integrate(aux,0,0.99,x,a,b)[[1]]
  return(sqrt(res))}   # Ls estimate by definition
IBPlfLsEf<-function(x,a,b){
  y<-sum(x)
  n<-length(x)
  num<-0 
  j<-0
  den<-beta(y+a,n+b)*Ixi(y+a,n+b,0.5)
  #changes
  while(10+b-j-2>0){
    term<-(j+1)*beta(y+a+j+2,n+b-j-2)*Ixi(y+a+j+2,n+b-j-2,0.5)
    num<-num+term
    j<-j+1}
  return(sqrt(num/den))}  # Ls estimator by Eq. (28)
################################################################################
# testing inverted beta PLF estimators and risk
################################################################################
set.seed(13579)
rho=0.5
smm1<-rgeom(n=100,1/(1+rho))
a<-1.0
b<-1.0
RhoE<-IBPlfRhoEf(smm1,a,b)  # Eq. (26)
RhoE-IBPlfRhof(smm1,a,b)    # checking by definition
IBPlfRhoEf(smm1,a,b)-IBSelfRhoEf(smm1,a,b)  # risk
LqE<-IBPlfLqEf(smm1,a,b)      # Eq. (27)
LqE-IBPlfLqf(smm1,a,b)      # checking by definition
LsE<-IBPlfLsEf(smm1,a,b)    # Eq. (28)
LsE-IBPlfLsf(smm1,a,b)      # checking by definition
rm(RhoE,LqE,LsE)
################################################################################
# plotting Jeffreys prior
################################################################################
p3<-function(rho){
  res<-rho^(-1/2)*(1+rho)^(-1/2)
  return(res)}  # Eq. (31)
setEPS()
#postscript(paste("FiJF.eps",sep=""),width=10.5*0.75,height=8*0.75)
par(mfrow=c(1,1))
rho<-seq(0.0,1.0,0.05)
plot(rho,p3(rho),xlab=expression(rho),
     ylab=expression(p[3](rho)),type="b",lty=1,pch=1,col=1,lwd=2)
legend("top",lty=1,pch=1,col=1,lwd=2,legend=c("Jeffreys"))
graphics.off()
################################################################################
# Jeffreys posterior
################################################################################
q3<-function(rho,x){
  y<-sum(x)
  n<-length(x)
  res<-rho^(y-0.5)*(1+rho)^(-n-y-0.5)/beta(y+0.5,n)/Ixi(y+0.5,n,0.5)
  return(res)}  # Eq. (32)
################################################################################
# Jeffreys: loss functions
################################################################################
# Jeffreys: SELF ---------------------------------------------------------------
# estimators and functions
################################################################################
JSelfRhof<-function(x){
  y<-sum(x)
  n<-length(x)
  aux<-function(rho,x){
    U<-rho*q3(rho,x)
    return(U)}
  res<-integrate(aux,0,1,x)$value
  return(res)}   # rho estimate by definition 
JSelfRhoEf<-function(x){
  y<-sum(x)
  n<-length(x)
  res<-(y+0.5)/(n-1)*Ixi(y+1.5,n-1,0.5)/Ixi(y+0.5,n,0.5)
  return(res)}  # rho estimator Eq. (33)
JSelfRhoER<-function(x){
  y<-sum(x)
  n<-length(x)
  risk<-beta(y+2.5,n-2)/beta(y+0.5,n)*
    Ixi(y+2.5,n-2,0.5)/Ixi(y+0.5,n,0.5)-JSelfRhoEf(smm1)^2
  return(risk)}  # risk
JSelfLqf<-function(x){
  aux<-function(rho,x){
    U<-rho^2/(1-rho)*q3(rho,x)
    return(U)}
  res<-integrate(aux,0,0.99,x)[[1]]
  return(res)}    # Lq estimate by definition
JSelfLqEf<-function(x){
  y<-sum(x)
  n<-length(x)
  num<-0
  den<-beta(y+0.5,n)*Ixi(y+0.5,n,0.5)
  j<-0
  #changes
  while(10-j-2>0){
    term<-beta(y+j+2.5,n-j-2)*Ixi(y+j+2.5,n-j-2,0.5)
    num<-num+term
    j<-j+1}
  return(num/den)}   # Lq estimator by Eq. (34)
JSelfLsf<-function(x){
  aux<-function(rho,x){
    U<-rho/(1-rho)*q3(rho,x)
    return(U)}
  res<-integrate(aux,0,1-1e-10,x)$value
  return(res)}    # Ls estimate by definition
JSelfLsEf<-function(x){
  y<-sum(x)
  n<-length(x)
  num<-0
  den<-beta(y+0.5,n)*Ixi(y+0.5,n,0.5)
  j<-0
  # changes
  while(10-j-1>0){
    term<-beta(y+j+1.5,n-j-1)*Ixi(y+j+1.5,n-j-1,0.5)
    num<-num+term
    j<-j+1}
  return(num/den)}   # Ls estimator by Eq. (35)
################################################################################
# testing Jeffreys SELF estimators and risk
################################################################################
set.seed(13579)
rho<-0.5
smm1<-rgeom(100,1/(1+rho))
RhoE<-JSelfRhoEf(smm1)        # Eq. (33)
RhoE-JSelfRhof(smm1)          # checking by definition
JSelfRhoER(smm1)              # risk
LqE<-JSelfLqEf(smm1)          # Eq. (34)
LqE-JSelfLqf(smm1)            # checking by definition
LsE<-JSelfLsEf(smm1)          # Eq. (35)
LsE-JSelfLsf(smm1)            # checking by definition
################################################################################
# Jeffreys: PLF ----------------------------------------------------------------
# estimators and functions
################################################################################
JPlfRhof<-function(rho,x){
  aux<-function(rho,x){
    U<-rho^2*q3(rho,x)
    return(U)}
  res<-integrate(aux,0,1,x)[[1]]
  return(sqrt(res))}  # rho estimate by definition
JPlfRhoEf<-function(x){
  y<-sum(x)
  n<-length(x)
  rhoE<-(y+1.5)*(y+0.5)/(n-1)/(n-2)*
    Ixi(y+2.5,n-2,0.5)/Ixi(y+0.5,n,0.5)
  return(sqrt(rhoE))}     # rho estimator by Eq. (36)
JPlfRhoER<-function(x){
  risk<-2*(JPlfRhoEf(x)-JSelfRhoEf(x))
  return(risk)}     # risk
JPlfLqf<-function(x){
  aux<-function(rho,x) return(rho^4/(1-rho)^2*q3(rho,x))
  res<-integrate(aux,0,0.99,x)[[1]]
  return(sqrt(res))}       # Lq estimate by definition
JPlfLqEf<-function(x){
  y<-sum(x)
  n<-length(x)
  num<-0
  den<-beta(y+0.5,n)*Ixi(y+0.5,n,0.5)
  j<-0
  #changes
  while(10-j-4>0){
    term<-(j+1)*beta(y+j+4.5,n-j-4)*Ixi(y+j+4.5,n-j-4,0.5)
    num<-num+term
    j<-j+1}
  return(sqrt(num/den))}      # Lq estimator by Eq. (37)
JPlfLsf<-function(x){
  aux<-function(rho,x) return(rho^2/(1-rho)^2*q3(rho,x))
  res<-integrate(aux,0,0.99,x)$value
  return(sqrt(res))}       # Ls estimate by definition
JPlfLsEf<-function(x){
  y<-sum(x)
  n<-length(x)
  num<-0
  den<-beta(y+0.5,n)*Ixi(y+0.5,n,0.5)
  j<-0
  while(10-j-2>0){
    term<-(j+1)*beta(y+j+2.5,n-j-2)*Ixi(y+j+2.5,n-j-2,0.5)
    num<-num+term
    j<-j+1}
  return(sqrt(num/den))}      # Ls estimator by Eq. (38)
################################################################################
# testing Jeffreys PLF estimators and risk
################################################################################
set.seed(13579)
rho<-0.5
smm1<-rgeom(100,1/(1+rho))
RhoE<-JPlfRhoEf(smm1)     # Eq. (36)
RhoE-JPlfRhof(rho,smm1)   # checking by definition
JPlfRhoER(smm1)           # risk    
LqE<-JPlfLqEf(smm1)       # Eq. (37)
LqE-JPlfLqf(smm1)         # checking by definition 
LsE<-JPlfLsEf(smm1)       # Eq. (38)
LsE-JPlfLsf(smm1)         # checking by definition
################################################################################
# beta: predictive distribution of number of customers
################################################################################
P1f<-function(m,x,a,b){
  aux<-function(rho,m,x,a,b){
    U<-Pm(m,rho)*q1(rho,x,a,b)
    return(U)}
  res<-sapply(m,function(M){
    integrate(aux,0,1,M,x,a,b)[[1]]})
  return(res)}
P1Ef<-function(m,x,a,b){
  y<-sum(x)
  n<-length(x)
  res<-sapply(m,function(M){
    beta(y+a+M,b)/beta(y+a,b)*
      GaussHypf(n+y+M+1,y+a+M,y+a+b+M,-1)/
      GaussHypf(n+y,y+a,y+a+b,-1)})
  return(res)}  # predictive distribution (eq.39)
################################################################################
# testing predictive distribution
################################################################################
set.seed(13579)
rho<-0.5
smm1<-rgeom(100,1/(1+rho))
m<-c(0,1,2,3,4,5)
a<-1.1; b<-1.1
P1Ef(m,smm1,a,b)-P1f(m,smm1,a,b)  
################################################################################
# inverted beta: predictive distribution of number of customers 
######P#########################################################################
P2f<-function(m,x,a,b){
  aux<-function(rho,m,x,a,b){
    U<-Pm(m,rho)*q2(rho,x,a,b-a)
    return(U)}
  res<-sapply(m,function(M){
    integrate(aux,0,1,M,x,a,b)[[1]]})
  return(res)}
P2Ef<-function(m,x,a,b){
  y<-sum(x)
  n<-length(x)
  res<-sapply(m,function(M){
    beta(y+a+M,n+b+1-a)*Ixi(y+a+M,n+b+1-a,0.5)/
      beta(y+a,n+b-a)/Ixi(y+a,n+b-a,0.5)})
  return(res)}        # predictive distribution (eq.40)
################################################################################
# testing predictive distribution
################################################################################
set.seed(13579)
rho<-0.5
smm1<-rgeom(100,1/(1+rho))
smm1
m<-c(0,1,2,3,4,5)
a<-1.1; b<-1.1
P2Ef(m,smm1,a,b)-P2f(m,smm1,a,b)
################################################################################
# jeffreys: predictive distribution of number of customers
################################################################################
P3f<-function(m,x){
  aux<-function(rho,m,x){
    U<-Pm(m,rho)*q3(rho,x)
    return(U)}
  res<-sapply(m,function(M){
    integrate(aux,0,1,M,x)$value})
  return(res)}
P3Ef<-function(m,x){
  y<-sum(x)
  n<-length(x)
  res<-sapply(m,function(M){
    beta(y+M+0.5,n+1)/beta(y+0.5,n)*
      Ixi(y+M+0.5,n+1,0.5)/Ixi(y+0.5,n,0.5)})
  return(res)}      # predictive distribution (eq.41)
################################################################################
# testing predictive distribution
################################################################################
set.seed(13579)
rho<-0.5
smm1<-rgeom(100,1/(1+rho))
m<-c(0,1,2,3,4,5)
P3Ef(m,smm1)-P3f(m,smm1)
################################################################################
# bayes factor (eq.42)
################################################################################
# Being B= Beta, I=Inverted, J=Jeffreys
BFBI<-function(m,x,a,b) return(P1Ef(m,x,a,b)/P2Ef(m,x,a,b))
BFBJ<-function(m,x,a,b) return(P1Ef(m,x,a,b)/P3Ef(m,x))
BFIJ<-function(m,x,a,b) return(P2Ef(m,x,a,b)/P3Ef(m,x))
################################################################################
# testing Bayes factor
################################################################################
set.seed(13579)
rho<-0.10
smm1<-rgeom(100,1/(1+rho))
a<-1.1; b<-1.1
m<-c(0,1,2,3,4,5)
BFBI(m,smm1,a,b)
BFBJ(m,smm1,a,b)
BFIJ(m,smm1,a,b)
################################################################################
# generating predictive distribution and Bayes factor 14/01/2021
################################################################################
a<-1.0; b<-1.1
size=100
rho<-0.10
set.seed(13579)
smm1<-rgeom(size,1/(1+rho))
m<-c(0,1,2,3,4,5)
 {LinesPD<-cbind(P1Ef(m,smm1,a,b))
  LinesPD<-cbind(LinesPD,P2Ef(m,smm1,a,b))
  LinesPD<-cbind(LinesPD,P3Ef(m,smm1))
  LinesPD<-cbind(LinesPD,BFBI(m,smm1,a,b))
  LinesPD<-cbind(LinesPD,1/BFBI(m,smm1,a,b))
  LinesPD<-cbind(LinesPD,BFBJ(m,smm1,a,b))
  LinesPD<-cbind(LinesPD,1/BFBJ(m,smm1,a,b))
  LinesPD<-cbind(LinesPD,BFIJ(m,smm1,a,b))
  LinesPD<-cbind(LinesPD,1/BFIJ(m,smm1,a,b))
  TablePD<-rbind(LinesPD)}
################################################################################
# Monte Carlo function for point estimates 14/01/2021
################################################################################
MontCarlo<-function(size,rho,Fest,...){
  set.seed(13579)
  rep<-100
  sample<-numeric(size)
  est<-numeric(rep)
  for(i in 1:rep){
    sample<-rgeom(size,1/(1+rho))
    est[i]<- Fest(sample,...)}
  # organizing result
  output=list(avr=mean(est),var=var(est))
  return(output)}
################################################################################
# Testing estimation of rho and Monte Carlo function
################################################################################
size<-500
rho<-0.99
Lq=1
rholq<-InvLq(1)
MontCarlo(size,rho,MM1RhoMLEf)
a<-2.0;b<-2.1
MontCarlo(size,rholq,IBSelfLqEf,a,b)
################################################################################
# Monte Carlo table function for point estimates 14/01/2021
################################################################################
TabMontCarlo<-function(size,rho,Fest,...){
  tab<-matrix(nrow=length(rho),ncol=2*length(size))
  for(i in 1:length(rho)){
    for(j in 1:length(size)){
      result<-MontCarlo(size[j],rho[i],Fest,...)
      tab[i,(2*j-1):(2*j)]<-c(result$avr,result$var)}}
    # cat(i/length(rho)*100,"%. \n")}
  return(tab)}
################################################################################
# Monte Carlo Parallel Rho function
################################################################################
# MontCP<-function(rho,size,a,b,type){
#   set.seed(13579)
#   rep<-1000
#   est1<-numeric(rep)  # EMV
#   est2<-numeric(rep)  # BSelf
#   est3<-numeric(rep)  # BPLF
#   est4<-numeric(rep)  # ISelf
#   est5<-numeric(rep)  # IPLF
#   est6<-numeric(rep)  # JSelf
#   est7<-numeric(rep)  # JPLF
#   estTot<-rbind(est1,est2,est3,est4,est5,est6,est7)
#   j=1
#   if(type=="RHO"){
#     while(j<=rep){
#       sample<-rgeom(size,1/(1+rho))
#       estTot[1:7,j]<-c(MM1RhoMLEf(sample),BSelfRhoEf(sample,a,b),
#                        BPLFRhoEf(sample,a,b),IBSelfRhoEf(sample,a,b),
#                        IBPlfRhoEf(sample,a,b),JSelfRhoEf(sample),
#                        JPlfRhoEf(sample))
#       if(is.finite(mean(estTot[1:7,j]))) j=j+1}}
#   if(type=="LQ"){
#     while(j<=rep){
#       sample<-rgeom(size,1/(1+rho))
#       estTot[1:7,j]<-c(MM1LqMLEf(sample),BSelfLqEf(sample,a,b),
#                        BPLFLqEf(sample,a,b),IBSelfLqEf(sample,a,b),
#                        IBPlfLqEf(sample,a,b),JSelfLqEf(sample),
#                        JPlfLqEf(sample))
#       if(is.finite(mean(estTot[1:7,j]))) j=j+1}}
#   if(type=="LS"){
#     while(j<=rep){
#       sample<-rgeom(size,1/(1+rho))
#       estTot[1:7,j]<-c(MM1LsMLEf(sample),BSelfLsEf(sample,a,b),
#                        BPLFLsEf(sample,a,b),IBSelfLsEf(sample,a,b),
#                        IBPlfLsEf(sample,a,b),JSelfLsEf(sample),
#                        JPlfLsEf(sample))
#      if(is.finite(mean(estTot[1:7,j]))) j=j+1}}
#   res<-rbind(apply(estTot,1,mean),apply(estTot,1,var))
#   return(res)}
################################################################################
# testing Monte Carlo Parallel Rho function
################################################################################
# a<-2.0;   b<-2.1; rho<-0.10; size<-200
# test<-MontCP(rho,size,a,b,"RHO")
# test[,"est1"]
# MontCP(rho,size,a,b,"LQ")
# MontCP(rho,size,a,b,"LS")
# ################################################################################
# # Monte Carlo Parallel table function
# ################################################################################
# TabMontCP<-function(size,rho,a,b,type){
#   for(k in 1:7){ name<-paste0("tab",k)
#   assign(name,matrix(nrow=length(rho),ncol=2*length(size),dimnames=
#                        list(rho,rep(c("m","v"),length(size)))))}
#   term<-0
#   termT<-length(rho)*length(size)
#   for(i in 1:length(rho)){
#     for(j in 1:length(size)){
#       ests<-MontCP(rho[i],size[j],a,b,type)
#       tab1[i,(2*j-1):(2*j)]<-ests[,"est1"]
#       tab2[i,(2*j-1):(2*j)]<-ests[,"est2"]
#       tab3[i,(2*j-1):(2*j)]<-ests[,"est3"]
#       tab4[i,(2*j-1):(2*j)]<-ests[,"est4"]
#       tab5[i,(2*j-1):(2*j)]<-ests[,"est5"]
#       tab6[i,(2*j-1):(2*j)]<-ests[,"est6"]
#       tab7[i,(2*j-1):(2*j)]<-ests[,"est7"]
#       term=term+1
#       cat(term/termT*100,"% conclu?dos. \n")}}
#   return(list(EMV=tab1,BS=tab2,BP=tab3,IS=tab4,IP=tab5,JS=tab6,JP=tab7))}
# ################################################################################
# # testing Monte Carlo Parallel table function
# ################################################################################
# a<-1.0; b<-1.1
# rho<-c(0.01,0.1,0.2,0.5)
# size<-c(10,20,50,100)
# RhoE<-TabMontCP(size,rho,a,b,"RHO")
# RhoE
# rm(RhoE)
# ################################################################################
# # generating parallel tables
# ################################################################################
# # Rho
# a<-1.0; b<-1.1
# rho<-c(0.01,0.1,0.2,0.5,0.7,0.9,0.99)
# size<-c(10,20,50,100,200,500)
# r.RhoPE<-TabMontCP(size,rho,a,b,"RHO")
# r.RhoPE 
# # Lq 
# a<-2.0; b<-2.1
# Lq<-c(0.1,0.2,0.5,1.0,2.0)
# size<-c(10,20,50,100,200,500)
# rhoLq<-InvLq(Lq)
# r.LqPE<-TabMontCP(size,rhoLq,a,b,"LQ")
# r.LqPE
# # Ls
# a<-2.0; b<-2.1
# Ls<-c(0.1,0.2,0.5,1.0,2.0)
# size<-c(10,20,50,100,200,500)
# rhoLs<-InvLs(Ls)
# r.LsPE<-TabMontCP(size,rhoLs,a,b,"LS")
# r.LsPE
# save(r.RhoPE,r.LqPE,r.LsPE,file="TabMontCP.RData")
# load(file=choose.files())
################################################################################
# credible interval estimation
################################################################################
CredInt<-function(dpf,sample,...){
  delta<-1E-03
  sl<-0.05 # significance level
  x<-seq(0,1,by=delta)
  dens<-dpf(x,sample,...)
  #  cat("CredInt:dens=",dens,"\n")
  accum<-0.0
  i=1
  # inferior limit
  while ((accum<sl/2)&(i<length(dens))){
    i<-i+1
    eps<-dens[i-1]+dens[i]
    if (is.nan(eps)) {
      accum<-sl/2}
    else {
      accum<-accum+delta*eps/2}}
  n1<-x[i-1]
  # superior limit
  while ((accum<(1-sl/2))&(i<length(dens))){
    i=i+1
    eps<-dens[i-1]+dens[i]
    if (!is.nan(eps)) {
      accum<-accum+delta*eps/2}}
  n2<-x[i]
  return(c(n1,n2))}
################################################################################
# testing credible interval estimation
################################################################################
set.seed(13579)
rho=0.99
smm1<-rgeom(100,1/(1+rho))
a<-1.0; b<-1.1
CredInt(q1,smm1,a,b)
CredInt(q2,smm1,a,b)
CredInt(q3,smm1)
################################################################################
# Monte Carlo function for credible interval estimation
################################################################################
MontCI<-function(dpf,size,rho,...){
  rep<-1000
  length<-numeric(rep)
  res<-numeric(2)
  cicover<-0
  for(i in 1:rep){
    sample<-rgeom(size,1/(1+rho))
    res<-CredInt(dpf,sample,...)
    #    cat("MontCI:res=",res,"\n")
    length[i]<-res[2]-res[1]
    if (res[1]<=rho&res[2]>=rho) cicover = cicover+1}
  return(c(mean(length),cicover/rep))}
################################################################################
# Testing Monte Carlo for credible interval estimation
################################################################################
set.seed(13579)
rho<-0.99
size<-500
a<-1.0; b<-1.1
MontCI(q1,size,rho,a,b)
################################################################################
# Monte Carlo table function for credible interval estimation
################################################################################
TabCI<-function(dpf,rho,size,...){
  tab<-matrix(nrow=length(rho),ncol=2*length(size))
  for(i in 1:length(rho)){
    for(j in 1:length(size)){
      set.seed(13579)
      tab[i,(2*j-1):(2*j)]<-MontCI(dpf,size[j],rho[i],...)}}
  return(tab)}
################################################################################
# generating tables: rho
################################################################################
rho<-c(0.01,0.1,0.2,0.5,0.7,0.9,0.99)
size<-seq(0,400,10)
r.EMVRhoE<-TabMontCarlo(size,rho,MM1RhoMLEf)
# Beta
a<-1.0; b<-1.1  # b > a
r.BSelfRhoE<-TabMontCarlo(size,rho,BSelfRhoEf,a,b)
r.BPlfRhoE<-TabMontCarlo(size,rho,BPLFRhoEf,a,b)
# Inverted 
r.ISelfRhoE<-TabMontCarlo(size,rho,IBSelfRhoEf,a,b)
r.IPlfRhoE<-TabMontCarlo(size,rho,IBPlfRhoEf,a,b)
# Jeffreys
r.JSelfRhoE<-TabMontCarlo(size,rho,JSelfRhoEf)
r.JPlfRhoE<-TabMontCarlo(size,rho,JPlfRhoEf)
save(rho,size,r.EMVRhoE,
     r.BSelfRhoE,r.BPlfRhoE,
     r.ISelfRhoE,r.IPlfRhoE,
     r.JSelfRhoE,r.JPlfRhoE,
     file="TabRho.RData")
#load(file=choose.files())
r.EMVRhoE
r.BSelfRhoE
r.BPlfRhoE
r.ISelfRhoE
r.IPlfRhoE
r.JSelfRhoE
r.JPlfRhoE
################################################################################
# generating tables: Lq
################################################################################
Lq<-c(0.1,0.2,0.5,1.0,2.0)
size<-c(10,20,50,100,200,500)
rhoLq<-InvLq(Lq)
r.EMVLqE<-TabMontCarlo(size,rhoLq,MM1LqMLEf)
# Beta
a<-2.0; b<-2.1  # b > a
r.BSelfLqE<-TabMontCarlo(size,rhoLq,BSelfLqEf,a,b)  # b > 1.0
r.BPlfLqE<-TabMontCarlo(size,rhoLq,BPLFLqEf,a,b)    # b > 2.0
# Inverted 
r.ISelfLqE<-TabMontCarlo(size,rhoLq,IBSelfLqEf,a,b) 
r.IPlfLqE<-TabMontCarlo(size,rhoLq,IBPlfLqEf,a,b) 
# Jeffreys
r.JSelfLqE<-TabMontCarlo(size,rhoLq,JSelfLqEf)  
r.JPlfLqE<-TabMontCarlo(size,rhoLq,JPlfLqEf)  
save(rho,size,
     r.EMVLqE,
     r.BSelfLqE,r.BPlfLqE,
     r.ISelfLqE,r.IPlfLqE,
     r.JSelfLqE,r.JPlfLqE,
     file="TableLq.RData")
#load(file=choose.files())
r.EMVLqE
r.BSelfLqE
r.BPlfLqE
r.ISelfLqE
r.IPlfLqE
r.JSelfLqE
r.JPlfLqE
################################################################################
# generating tables: Ls
################################################################################
Ls<-c(0.1,0.2,0.5,1.0,2.0)
size<-c(10,20,50,100,200,500)
rhoLs<-InvLs(Ls)
r.EMVLsE<-TabMontCarlo(size,rhoLs,MM1LsMLEf)
# Beta
a<-2.0; b<-2.1  # b > a
r.BSelfLsE<-TabMontCarlo(size,rhoLs,BSelfLsEf,a,b)  # b > 1.0
r.BPlfLsE<-TabMontCarlo(size,rhoLs,BPLFLsEf,a,b)    # b > 2.0
# Inverted 
r.ISelfLsE<-TabMontCarlo(size,rhoLs,IBSelfLsEf,a,b) 
r.IPlfLsE<-TabMontCarlo(size,rhoLs,IBPlfLsEf,a,b) 
# Jeffreys
r.JSelfLsE<-TabMontCarlo(size,rhoLs,JSelfLsEf)  
r.JPlfLsE<-TabMontCarlo(size,rhoLs,JPlfLsEf)  
save(rho,size,
     r.EMVLsE,
     r.BSelfLsE,r.BPlfLsE,
     r.ISelfLsE,r.IPlfLsE,
     r.JSelfLsE,r.JPlfLsE,
     file="TableLs.RData")
#load(file=choose.files())
r.EMVLsE
r.BSelfLsE
r.BPlfLsE
r.ISelfLsE
r.IPlfLsE
r.JSelfLsE
r.JPlfLsE
################################################################################
# generating tables: credible interval
################################################################################
#rho<-c(0.01,0.1,0.2,0.5,0.7,0.9,0.99)
rho<-seq(0,1,0.1)
size<-c(10,20,50,100,200)
a=1.0; b=1.1
CredB<-TabCI(q1,rho,size,a,b)
CredI<-TabCI(q2,rho,size,a,b)
CredJ<-TabCI(q3,rho,size)
save(rho,size,CredB,CredI,CredJ,file="TabCI.RData")
#load(file=choose.files())
CredB
CredI
CredJ
################################################################################
THE END
################################################################################
