\name{bandsolve-package}
\alias{bandsolve-package}
\alias{bandsolve}
\docType{package}
\title{
\packageTitle{bandsolve}
}
\description{
\packageDescription{bandsolve}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{bandsolve}
\packageIndices{bandsolve}
A C++ implementation of a specific LU decomposition for band matrix. A a symmetric band matrix, b a vector, we either give the matrixes L and U such that A=LU, and also the solution x of the system Ax=b.
}
\author{
\packageAuthor{bandsolve}

Maintainer: \packageMaintainer{bandsolve}
}
\references{

}

\keyword{ package }
\seealso{

}
\examples{
n=2000;
D0=runif(n);
D1=-0.2*runif(n-1);
D2=0.3*runif(n-2);
D3=0.4*runif(n-3);
D4=-0.1*runif(n-4);
b=runif(n)
ref=bandsolve(D0,D1,D2,D3,D4,b)
}
