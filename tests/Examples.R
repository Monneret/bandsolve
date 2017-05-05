require(bandsolve)
require(bvpSolve)
require(microbenchmark)
require(pracma)
set.seed(100)
### ODE- Poisson equation in 1D ###
### Laplacien(y)=f(x) -> y''=f(x); centered second order discrete derivative
### Reformulate as standard ODE: y2'=f(x) && y1'=y2
dx=10^{-4}
gridx=seq(from=0,to=0.05,by=dx)
n=length(gridx)

f=(exp(-(gridx-0.025)^2)-1)*1000
c_1=-500
c_2=1.5
analsol=-886.227*gridx*erf(0.025-gridx)+22.1557*erf(0.025-gridx)+499.688*exp(gridx*(0.05-gridx))+c_2*gridx+c_1-500*gridx^2
plot(analsol)
secondorder=(2*analsol[2:500]-analsol[1:499]-analsol[3:501])*(1/dx^2)
plot(secondorder)
plot(f)
## Dirichlet boundary condition u(0)=u(end)=0

D0=rep(2,n-2)*(1/dx^2);
D1=rep(-1,n-3)*(1/dx^2);
A=cbind(D0,c(D1,0));

res=bandsolve(A,b=f[2:(n-1)],inplace=TRUE);

derivsyst<-function(x,State,Pars) {
  with(as.list(c(State, Pars)),{
  Dy1=y2
  Dy2=(exp(-(x-0.025)^2)-1)*1000
  res<-c(Dy1,Dy2)
  return(list(res))
  })
}

yini=c(y1=0,y2=NA)
yend=c(y1=0,y2=NA)
pars=c()
out<-bvptwp(yini,gridx,derivsyst,yend,pars)
plot(out)
plot(res$x)

### Smoother ###
n=1000
lambda=300
dx=1
gridx=seq(from=1,to=1000,by=dx)
x=0.2*sin(gridx/(40*pi))+10^(-6)*gridx^2-10^(-3)*gridx+1+rnorm(n = n,mean = 0,sd = 0.1)
xmem=x
names(xmem)=paste(1:n)
plot(gridx,xmem) ## noisy data

D=matrix(0,n-1,n-1)
diag(D)=-1
diag(D[c(-n),-1])=1
D=cbind(D,c(rep(0,n-2),1)) ## discrete derivative

y=smooth.spline(x)$y
points(gridx,y,type="p",pch=1, lwd = 0.1,col="blue",ylim=c(0.2,1.4))

A=mat2rot(diag(n)+lambda*t(D)%*%D)
res=bandsolve(A,x,inplace=TRUE)
lines(gridx,x,type="l",lwd = 5,col="red",ylim=c(0.2,1.4))

sum(abs(y-x))
x=xmem

smith<-function(){
  x=xmem
  smooth.spline(x)
}

band<-function(){
  x=xmem
  bandsolve(A,x,inplace=TRUE)
}
tmp=microbenchmark(RSpline=smith(),Bdsolve=band())
boxplot(tmp,ylab="time in ms",xlab="")
