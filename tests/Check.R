## Check functions
refbandsolve5=function(D0,D1,D2,D3,D4,D5,b) {
  A=diag(D0);
  A[-n,-1]=A[-n,-1]+diag(D1)
  A[-1,-n]=A[-1,-n]+diag(D1)
  A[-c(n,n-1),-c(1,2)]=A[-c(n,n-1),-c(1,2)]+diag(D2)
  A[-c(1:2),-c(n,n-1)]=A[-c(1:2),-c(n,n-1)]+diag(D2)
  A[-c(n,n-1,n-2),-c(1:3)]=A[-c(n,n-1,n-2),-c(1:3)]+diag(D3)
  A[-c(1:3),-c(n,n-1,n-2)]=A[-c(1:3),-c(n,n-1,n-2)]+diag(D3)
  A[-c(n,n-1,n-2,n-3),-c(1:4)]=A[-c(n,n-1,n-2,n-3),-c(1:4)]+diag(D4)
  A[-c(1:4),-c(n,n-1,n-2,n-3)]=A[-c(1:4),-c(n,n-1,n-2,n-3)]+diag(D4)
  A[-c(n,n-1,n-2,n-3,n-4),-c(1:5)]=A[-c(n,n-1,n-2,n-3,n-4),-c(1:5)]+diag(D5)
  A[-c(1:5),-c(n,n-1,n-2,n-3,n-4)]=A[-c(1:5),-c(n,n-1,n-2,n-3,n-4)]+diag(D5)
  x=solve(A,b)
  return(list(A=A,x=x))
}

n=10;
D0=runif(n);
D1=-2*runif(n-1);
D2=3*runif(n-2);
D3=0.4*runif(n-3);
D4=-0.9*runif(n-4);
D5=runif(n-5);
b=runif(n);
D=cbind(D0,c(D1,0),c(D2,0,0),c(D3,0,0,0),c(D4,0,0,0,0),c(D5,0,0,0,0,0))

ref5=refbandsolve5(D0,D1,D2,D3,D4,D5,b);
res5=bandsolve(A=D,b=b)
plot(res5$x,ref5$x,xlab="Bandsolve",ylab="Solve"); abline(0,1)
title(main="Consistency of solutions")



Ures=rot.as.mat(D,u=5,l=0);
Dres=diag(diag(Ures));
Lres=t(Ures);
diag(Ures)=1;
diag(Lres)=1;
Aret=Lres%*%Dres%*%Ures;
plot(diag(Aret)[1:1000],diag(ref5$A)[1:1000],xlab="Bandsolve",ylab="Solve"); abline(0,1);
title(main="Consistency of decompositions")