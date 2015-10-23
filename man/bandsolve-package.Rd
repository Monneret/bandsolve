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
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{bandsolve}

Maintainer: \packageMaintainer{bandsolve}
}
\references{

}

\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
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
