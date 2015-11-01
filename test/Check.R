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

n=1000;
D0=runif(n);
D1=-2*runif(n-1);
D2=3*runif(n-2);
D3=0.4*runif(n-3);
D4=-0.9*runif(n-4);
D5=runif(n-5);
b=runif(n);
D=cbind(D0,c(D1,0),c(D2,0,0),c(D3,0,0,0),c(D4,0,0,0,0),c(D5,0,0,0,0,0))


ref1=refbandsolve5(D0,D1,D2*0,D3*0,D4*0,D5*0,b);
ref2=refbandsolve5(D0,D1,D2,D3*0,D4*0,D5*0,b);
ref3=refbandsolve5(D0,D1,D2,D3,D4*0,D5*0,b);
ref4=refbandsolve5(D0,D1,D2,D3,D4,D5*0,b);
ref5=refbandsolve5(D0,D1,D2,D3,D4,D5,b);
ref=list(ref1,ref2,ref3,ref4,ref5);
res=list()

par(mfrow=c(2,3))
for (i in 1:5){
  res[[length(res)+1]]=bandsolve(D=D[,1:(i+1)],b=b)
  plot(res[[i]]$x,ref[[i]]$x); abline(0,1)
}
title(main="Consistency of solutions")

par(mfrow=c(2,3))
for (i in 1:4){
  L=matrix(0,n,n)
  diag(L)=1
  U=matrix(0,n,n)
  diag(U[-c(n:(n-i+1)),-c(1:i)])=D[-c(n:(n-i+1)),i+1]
  for (j in 1:i) {
    diag(L[-c(1:j),-c(n:(n-j+1))])=res[[i]][[j+1]];
    if (j>1) diag(U[-c(n:(n-j+2)),-c(1:(j-1))])=res[[i]][[j+i+1]] else diag(U)=res[[i]][[j+i+1]];
  }
  Aret=L%*%U
  plot(diag(Aret),diag(ref[[i]]$A)); abline(0,1)
}

L=matrix(0,n,n);
diag(L)=1;
U=matrix(0,n,n);
diag(U[-c(n:(n-4)),-c(1:5)])=D5;
for (j in 1:5) {
  diag(L[-c(1:j),-c(n:(n-j+1))])=res[[5]]$L[-c(n:(n-j+1)),j+1];
  if (j>1) diag(U[-c(n:(n-j+2)),-c(1:(j-1))])=res[[5]]$U[-c(n:(n-j+2)),j] else diag(U)=res[[5]]$U[,j];
}
Aret=L%*%U;
plot(diag(Aret),diag(ref[[5]]$A)); abline(0,1);
title(main="Consistency of decompositions")