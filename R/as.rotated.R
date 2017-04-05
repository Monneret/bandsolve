#' @useDynLib bandsolve
#' @import Rcpp
#' @export
#' @description Rotate a band matrix to get the rotated row-wised matrix associated.
#' @title as.rotated
#' @title Rotate a band matrix to get the rotated row-wised matrix associated.
#' @param M Band square matrix.
#' @param u upper bandwith.
#' @param l lower bandwith.
#' @return Rotated matrix.
#' @examples n=2000;
#' A=diag(4)
#' A[2,3]=2
#' ref=as.rotated(A,l=0,u=1)
#' 

as.rotated <- function(M,l=0,u=0){
  N=ncol(M);
  R=matrix(0,N,l+u+1);
  if (l>0){
    for (i in 1:l){
      R[,l-i+1]=c(diag(M[-c(1:i),-c(N:(N-i+1))]),rep(0,i))
    }
  }
  R[,l+1]=diag(M);
  if (u>0){
    for (j in 1:u){
      R[,l+1+j]=c(diag(M[-c(N:(N-j+1)),-c(1:j)]),rep(0,j))
    }
  }
  return(R)
}