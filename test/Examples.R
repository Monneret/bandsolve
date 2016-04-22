require(bandsolve)
require(rootsolve)
### ODE- Poisson equation in 1D ###
dx=10^{-6}
gridx=seq(from=0,to=0.05,by=dx)
n=length(gridx)

f=(exp(-(gridx-0.025)^2)/exp(-(0.025)^2)-1)*1000
model<-function(space,y,params){
  with(as.list(c(y,params)),{
    dX <-
    list(dX)
  })
  
}
## Dirichlet boundary condition u(0)=u(end)=0

D0=rep(2,n-2)*(1/dx^2);
D1=rep(-1,n-3)*(1/dx^2);
A=cbind(D0,c(D1,0));

bandsolve(A,u=1,l=1,b=f[2:(n-1)],sym=TRUE);

## comparer avec rootsolve

### Cubic spline ###
n=1000
lambda=300
dx=1
gridx=seq(from=1,to=1000,by=dx)
x=0.2*sin(gridx/(40*pi))+10^(-6)*gridx^2-10^(-3)*gridx+1+rnorm(n = n,mean = 0,sd = 0.1)
xmem=x
names(xmem)=paste(1:n)
plot(gridx,xmem)

D=matrix(0,n-1,n-1)
diag(D)=-1
diag(D[c(-n),-1])=1
D=cbind(D,c(rep(0,n-2),1))

A=as.rotated(diag(n)+lambda*t(D)%*%D,u=1,l=0)
res=bandsolve(A,x,u=1,l=0,sym=TRUE)
lines(gridx,x,type="l",lwd = 5,col="red",ylim=c(0.2,1.4))
## comparer avec smooth.spline, microbench avec les 2 et hop
smooth.spline(x)