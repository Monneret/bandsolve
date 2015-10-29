require(bandsolve)
require(microbenchmark)
require(pracma)
require(spam);
require(limSolve)


set.seed(421)
n=2000;
D0=rep(1,n)
D1=rep(-0.5,n-1)
b=rnorm(n)

A=diag(D0);
A[-n,-1]=A[-n,-1]+diag(D1);
A[-1,-n]=A[-1,-n]+diag(D1);

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
                 FBS=bandsolve(D0,D1,b=b),times=100)
boxplot(r)

tmp=sapply(split(r$time,r$expr),mean)
tmp/tmp["FBS"]

#####

set.seed(421)
n=20000;
D0=rep(1,n)
D1=rep(-0.5,n-1)
b=rnorm(n)

A=spam(0,nrow=n,ncol=n)
diag(A)=D0
U=spam(0,nrow=n+1,ncol=n+1)
diag(U)=c(0,D1,0)
A=A+U[-(n+1),-1]+U[-1,-(n+1)]

rotated.A=rbind(c(0,D1),D0,c(D1,0))


#plot(solve(A,b),bandsolve(D0,D1,b=b))

r=microbenchmark(SPAM.SOLVE=solve(A,b),
                 BANDED=Solve.banded(rotated.A,nup=1,nlow=1,b)[,1],
                 FBS=bandsolve(D0,D1,b=b),times=100)
boxplot(r)

tmp=sapply(split(r$time,r$expr),mean)
tmp/tmp["FBS"]

