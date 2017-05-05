#' @useDynLib bandsolve
#' @import Rcpp
#' @export
#' @description Rotate a symmetric band matrix to get the rotated matrix associated.
#' Each column of the rotated matrix correspond to a diagonal. The first column is the main diagonal, the second one is the upper-diagonal and so on.
#' Artificial 0 are placed at the end of each column if necessary.
#' @title mat2rot
#' @title Rotate a band matrix to get the rotated row-wised matrix associated.
#' @param M Band square matrix or a list of diagonal.
#' @return Rotated matrix.
#' @examples
#' 
#' A=diag(4)
#' A[2,3]=2
#' A[3,2]=2
#' 
#' ## Original Matrix
#' A
#' ## Rotated version
#' R=mat2rot(A)
#' R
#' 
#' rot2mat(mat2rot(A))

mat2rot <- function(M){
  if (is.matrix(M)){
  N=ncol(M);
  l=0
  for (i in 1:N){
    lprime=which(M[i,]!=0)
    if (lprime[length(lprime)]-i>l) l=lprime[length(lprime)]-i;
  }
  R=matrix(0,N,l+1);
  R[,1]=diag(M);
  if (l>0){
    for (j in 1:l){
      if (j == N-1){
        R[,1+j]=c(M[-c(N:(N-j+1)),-c(1:j)],rep(0,j))
      } else {
        R[,1+j]=c(diag(M[-c(N:(N-j+1)),-c(1:j)]),rep(0,j))
      }
    }
  }
  return(R)
  } else if (is.list(M)){
    R=matrix(0,length(M[[1]]),length(M))
    for (i in 1:length(M)){
      R[,i]=c(M[[i]],rep(0,i-1))
    }
    return(R)
  }
}
