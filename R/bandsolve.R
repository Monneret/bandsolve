#' @useDynLib bandsolve
#' @import Rcpp
#' @export
#' @description Main function to solve efficiently and quickly a symmetric bandlinear system. Theses systems are solved much faster than standards system, dropping from complexity O(n³) to O(0.5*nk²), where k is the number of sub diagonal.
#' @title bandsolve
#' @param A Band square matrix in rotated form. The rotated form can be obtained with the function as.rotated: it's the visual rotation by 90 degrees of the matrix, where subdiagonal are discarded.
#' @param b right hand side of the equation. Can be either a vector or a matrix. If not supplied, the function return the inverse of A.
#' @param inplace Should results overwrite pre-existing data? Default set to false.
#' @return Solution of the linear problem.
#' @examples n=2000;
#' A=diag(4)
#' A[2,3]=2
#' ref=as.rotated(A,l=0,u=1)
#' bandsolve(ref)
#' \code{\link{LDL}},\code{\link{rot.as.mat}},\code{\link{as.rotated}}

bandsolve<-function(A,b=NULL,inplace=FALSE){
  if(nrow(A)==ncol(A)) stop("A should be a rotated matrix!");
  if (is.vector(b)){
    if(length(b)!=nrow(A)) stop("Dimension problem");
    if (inplace){
      return(fastbandsolve(A,b))
    } else {
      Amem=A
      bmem=b
      return(fastbandsolve(Amem,bmem))
    }
  } else if (is.matrix(b)){
    if(nrow(b)!=nrow(A)) stop("Dimension problem");
    if (inplace){
      return(fastbandsolveMatrix(A,b))
    } else {
      A=Amem
      Bmem=b
      return(fastbandsolveMatrix(Amem,Bmem))
    }
  } else if (is.null(b)){
    if (inplace){
      return(LDL(A))
    } else {
      Amem=A
      return(LDL(Amem))
    }
  } else {
    stop("b must either be a vector or a matrix")
  }
}