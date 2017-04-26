#' @useDynLib bandsolve
#' @export
#' 
#' @description Get back a symmetric square matrix based on his rotated row-wised version. 
#' The rotated form of the input is such each column correspond to a diagonal, where the first column is the main diagonal and next ones are the upper/lower-diagonal.
#' To match dimension, last element of these columns are discarded.
#' @title rot.mat
#' @title Get back a symmetric square matrix based on his rotated row-wised version.
#' @param R Rotated matrix.
#' @return Band square matrix.
#' @examples
#' require(bandsolve)
#' D0=1:5;
#' D1=c(0,1,0,0);
#' D2=rep(2,3);
#' 
#' ref=rot.mat(cbind(D0,c(D1,0),c(D2,0,0)))
#' ref
#' mat.rot(rot.mat(cbind(D0,c(D1,0),c(D2,0,0))))
#' 


rot.mat <- function(R){
  N=nrow(R)
  l=ncol(R)
  M=matrix(0,N,N)
  if (l>1){
  for (i in 1:(l-1)){
      diag(M[-c(1:i),-c(N:(N-i+1))])=R[-c(N:(N-i+1)),i+1]
      diag(M[-c(N:(N-i+1)),-c(1:i)])=R[-c(N:(N-i+1)),i+1]
  }
  diag(M)=R[,1];
  }
  return(M)
}