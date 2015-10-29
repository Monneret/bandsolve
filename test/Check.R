## Check functions
refbandsolve4=function(D0,D1,D2,D3,D4,b) {
  A=diag(D0);
  A[-n,-1]=A[-n,-1]+diag(D1)
  A[-1,-n]=A[-1,-n]+diag(D1)
  A[-c(n,n-1),-c(1,2)]=A[-c(n,n-1),-c(1,2)]+diag(D2)
  A[-c(1:2),-c(n,n-1)]=A[-c(1:2),-c(n,n-1)]+diag(D2)
  A[-c(n,n-1,n-2),-c(1:3)]=A[-c(n,n-1,n-2),-c(1:3)]+diag(D3)
  A[-c(1:3),-c(n,n-1,n-2)]=A[-c(1:3),-c(n,n-1,n-2)]+diag(D3)
  A[-c(n,n-1,n-2,n-3),-c(1:4)]=A[-c(n,n-1,n-2,n-3),-c(1:4)]+diag(D4)
  A[-c(1:4),-c(n,n-1,n-2,n-3)]=A[-c(1:4),-c(n,n-1,n-2,n-3)]+diag(D4)
  x=solve(A,b)
  return(list(A=A,x=x))
}

n=2000;
D0=runif(n);
D1=-2*runif(n-1);
D2=3*runif(n-2);
D3=0.4*runif(n-3);
D4=-0.9*runif(n-4);
b=runif(n)
ref=refbandsolve4(D0,D1,D2,D3,D4,b);
x=bandsolve(D0,D1,D2,D3,D4,b)
plot(x[1:1000],ref$x[1:1000]); abline(0,1)

## BandSolveK
D=cbind(D0,c(D1,0),c(D2,0,0),c(D3,0,0,0),c(D4,0,0,0,0))
D=cbind(D0,c(D1,0))
x=bandsolveK(Ddata = D,bdata = b)
xslow=slowbandsolve(D,b)