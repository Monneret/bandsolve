require(microbenchmark)
require(pracma)
n=1000;
D0=runif(n);
D1=-0.2*runif(n-1);
b=runif(n)

A=diag(D0);
A[-n,-1]=A[-n,-1]+diag(D1);
A[-1,-n]=A[-1,-n]+diag(D1);

r=microbenchmark(SOLVE=solve(A,b),PRA=trisolve(D0,D1,D1,b),FBS=bandsolve(D0,D1,b=b),times=100)
boxplot(r)
