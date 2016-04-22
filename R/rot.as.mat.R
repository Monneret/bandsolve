#' @useDynLib bandsolve
#' @export
#' 
#' @description Get back a square matrix based on his rotated row-wised version.
#' @title rot.as.mat
#' @param R Rotated matrix.
#' @param u upper bandwith.
#' @param l lower bandwith.
#' @return Band square matrix.
#' @examples n=2000;
#' D0=runif(n);
#' D1=-0.2*runif(n-1);
#' D2=0.3*runif(n-2);
#' ref=rot.as.mat(cbind(D0,c(D1,0),c(D2,0,0)),l=0,u=2)
#'

rot.as.mat <- function(R,l=0,u=0){
  N=nrow(R)
  M=matrix(0,N,N)
  if (l>0){
    for (i in 1:l){
      diag(M[-c(1:i),-c(N:(N-i+1))])=R[-(N-i+1:N),l-i+1]
    }
  }
  diag(M)=R[,l+1];
  if (u>0){
    for (j in 1:u){
      diag(M[-c(N:(N-j+1)),-c(1:j)])=R[-(N-j+1:N),l+1+j]
    }
  }
  return(M)
}