#' @useDynLib bandsolve
#' @import Rcpp
#' @export
#' @title Fast solver for linears systems involving symmetric band matrix.
#' @param D Matrice rotationné
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

bandsolve  <- function(D,b){## Matrice rotationnée
  if (ncol(D)==2) {
    bandsolve1(D0data=D[,1], D1data=D[-nrow(D),2], bdata=b) 
  }else if (ncol(D)==3){ 
    bandsolve2(D0data=D[,1], D1data=D[-nrow(D),2], D2data=D[-c(nrow(D)-1,nrow(D)),3], bdata=b) 
  }else if (ncol(D)==4){ 
    bandsolve3(D0data=D[,1], D1data=D[-nrow(D),2], D2data=D[-c(nrow(D)-1,nrow(D)),3], D3data=D[-c(nrow(D)-2,nrow(D)-1,nrow(D)),4], bdata=b) 
  }else if (ncol(D)==5){ 
    bandsolve4(D0data=D[,1], D1data=D[-nrow(D),2], D2data=D[-c(nrow(D)-1,nrow(D)),3], D3data=D[-c(nrow(D)-2,nrow(D)-1,nrow(D)),4], D4data=D[-c(nrow(D)-3,nrow(D)-2,nrow(D)-1,nrow(D)),5],bdata=b)
  }else{
    bandsolveK(Ddata = D,bdata=b)
  }
}