require(bandsolve)
require(microbenchmark)
require(pracma)
require(spam)
require(limSolve)

set.seed(42)
n=1000;
D0=rep(1.25,n)
D1=rep(-0.5,n-1)
b=rnorm(n)

A=diag(D0);
A[-n,-1]=A[-n,-1]+diag(D1);
A[-1,-n]=A[-1,-n]+diag(D1);
det(A)
spam.A=as.spam(A)
rotated.A=rbind(c(0,D1),D0,c(D1,0))

#plot(solve(spam.A,b),bandsolve(D0,D1,b=b)); abline(0,1)
#plot(Solve.tridiag(D1,D0,D1,b)[,1],bandsolve(D0,D1,b=b)); abline(0,1)
#plot(Solve.banded(rotated.A,nup=1,nlow=1,b)[,1],bandsolve(D0,D1,b=b)); abline(0,1)


r=microbenchmark(SOLVE=solve(A,b),
                 SPAM.SOLVE=solve(spam.A,b),
                 PRA=trisolve(D0,D1,D1,b),
                 TRIDIAG=Solve.tridiag(D1,D0,D1,b)[,1],
                 BANDED=Solve.banded(rotated.A,nup=1,nlow=1,b)[,1],
                 FBS=bandsolve(t(rotated.A)[,-1],l=1,u=1,b=b,sym=TRUE)$x,times=100)
boxplot(r)

tmp=sapply(split(r$time,r$expr),mean)
tmp/tmp["FBS"]

#####

set.seed(421)
n=10000;
D0=rep(2.5,n)
U1=rep(0.2,n-1)
b=rnorm(n)

A=spam(0,nrow=n,ncol=n)
U1sp=spam(0,nrow=n+1,ncol=n+1)
diag(U1sp)=c(0,U1,0)
diag(A)=D0
A=A+U1sp[-1,-n]+U1sp[-n,-1]
rotated.A=rbind(c(0,U1),D0,c(U1,0))
rotated.Afbs=cbind(c(U1,0),D0,c(U1,0))


#plot(solve(A,b),bandsolve(D0,D1,b=b))

r=microbenchmark(SPAM.SOLVE=solve(A,b),
                 BANDED=Solve.banded(rotated.A,nup=1,nlow=1,b)[,1],
                 FBS=bandsolve(rotated.Afbs,l=1,u=1,b=b)$x,times=100) #FBS1=bandsolve1(D0data=D0, D1data=D1, bdata=b)$x

boxplot(r)

tmp=sapply(split(r$time,r$expr),mean)
tmp/tmp["FBS"] 

#mean(abs(x/y-1)) error relative ne depasse pas %eps

### bandsolve 1,2,3,4 vs bandsolveK

# n=2000;
# D0=runif(n);
# D1=-2*runif(n-1);
# D2=3*runif(n-2);
# D3=0.4*runif(n-3);
# D4=-0.9*runif(n-4);
# 
# b=runif(n);
# D=cbind(D0,c(D1,0),c(D2,0,0),c(D3,0,0,0),c(D4,0,0,0,0))
# 
# r1=microbenchmark(Band1234=bandsolve1(D0data = D0,D1data = D1,bdata = b),
#                  BandK=bandsolveK(Ddata=D[,1:2],bdata=b)$x,times=100)
# r2=microbenchmark(Band1234=bandsolve2(D0data = D0,D1data = D1,D2data = D2,bdata = b),
#                   BandK=bandsolveK(Ddata=D[,1:3],bdata=b)$x,times=100)
# r3=microbenchmark(Band1234=bandsolve3(D0data = D0,D1data = D1,D2data = D2,D3data = D3,bdata = b),
#                   BandK=bandsolveK(Ddata=D[,1:4],bdata=b)$x,times=100)
# r4=microbenchmark(Band1234=bandsolve4(D0data = D0,D1data = D1,D2data = D2,D3data = D3,D4data=D4,bdata = b),
#                   BandK=bandsolveK(Ddata=D,bdata=b)$x,times=100)
# #microbenchmark(Band1234=bandsolve4(D0data = D0,D1data = D1,D2data = D2,D3data = D3,D4data=D4,bdata = b),times=10)
# par(mfrow=c(2,2))
# boxplot(r1,main="bandsolve1")
# boxplot(r2,main="bandsolve2")
# boxplot(r3,main="bandsolve3")
# boxplot(r4,main="bandsolve4")
# par(mfrow=c(1,1))

### Bandwith: complexity
require(MASS)
time=vector(mode="numeric",length=51)
for (bd in (1:50)){
  A=mvrnorm(200,rep(0,bd+1),Sigma=diag(bd+1))
  A=cbind(A[,(bd+1):2],A)
  b=rnorm(200)
  time[bd+1]=mean(microbenchmark(FBS=bandsolve(A,l=bd,u=bd,b=b)$x,times=100)$time)/1000
}
plot(time)
lines((1:50)^2,col = "red",add=TRUE)

### Size of matrix: complexity
time=vector(mode="numeric",length=200)
for (sz in (11:210)){
  A=mvrnorm(sz,rep(0,3),Sigma=diag(3))
  A=cbind(A[,3:2],A)
  b=rnorm(sz)
  time[sz-10]=mean(microbenchmark(FBS=bandsolve(A,l=bd,u=bd,b=b)$x,times=100)$time)/1000
}
plot(time)
