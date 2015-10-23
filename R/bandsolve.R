#' @useDynLib bandsolve
#' @import Rcpp
#' @export
#' @title Fast solver for linears systems involving symmetric band matrix.
#' @param D0 Diagonal as a vector of length n.
#' @param D1 First subdiagonal and superdiagonal as a vector of length n-1.
#' @param D2 Second subdiagonal and superdiagonal as a vector of length n-2.
#' @param D3 Third subdiagonal and superdiagonal as a vector of length n-3.
#' @param D4 Fourth subdiagonal and superdiagonal as a vector of length n-4.
#' @param b Right hand side of the linear system
#' @return Solution of the linear system.
#' @examples n=2000;
#' D0=runif(n);
#' D1=-0.2*runif(n-1);
#' D2=0.3*runif(n-2);
#' D3=0.4*runif(n-3);
#' D4=-0.1*runif(n-4);
#' b=runif(n)
#' ref=bandsolve(D0,D1,D2,D3,D4,b)

bandsolve  <- function(D0,D1,D2=NULL,D3=NULL,D4=NULL,b){
  if (is.null(D2)) {bandsolve1(D0data=D0, D1data=D1, bdata=b)} else
    if (is.null(D3)) {bandsolve2(D0data=D0, D1data=D1, D2data=D2, bdata=b)} else
      if (is.null(D4)) {bandsolve3(D0data=D0, D1data=D1, D2data=D2, D3data=D3, bdata=b)} else {
        bandsolve4(D0data=D0, D1data=D1, D2data=D2, D3data=D3, D4data=D4, bdata=b)}
}