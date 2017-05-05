#' @useDynLib bandsolve
#' @import Rcpp
#' @export
#' @description Main function to solve efficiently and quickly a symmetric bandlinear system. Theses systems are solved much faster than standards system, dropping from complexity O(n³) to O(0.5*nk²), where k is the number of sub diagonal.
#' @title bandsolve
#' @param A Band square matrix in rotated form. The rotated form can be obtained with the function as.rotated: it's the visual rotation by 90 degrees of the matrix, where subdiagonal are discarded.
#' @param b right hand side of the equation. Can be either a vector or a matrix. If not supplied, the function return the inverse of A.
#' @param inplace Should results overwrite pre-existing data? Default set to false.
#' @return Solution of the linear problem.
#' @examples
#' 
#' A=diag(4)
#' A[2,3]=2
#' A[3,2]=2
#' R=mat2rot(A)
#' solve(A)
#' bandsolve(ref)
#' 
#' set.seed(100)
#' 
#' n=1000;
#' D0=rep(1.25,n)
#' D1=rep(-0.5,n-1)
#' b=rnorm(n)
#'
#' ## Comparison with solve
#' 
#' if(require(microbenchmark)){
#' A=diag(D0);
#' A[-n,-1]=A[-n,-1]+diag(D1);
#' A[-1,-n]=A[-1,-n]+diag(D1);
#' R=mat2rot(list(D0,D1))
#' r=microbenchmark(
#' SOLVE=solve(A,b),
#' BANDSOLVE=bandsolve(R,b=b,inplace=TRUE),times=100)
#'boxplot(r)
#'}

bandsolve<-function(A,b=NULL,inplace=FALSE){
  if((nrow(A)==ncol(A))&(A[nrow(A),ncol(A)]!=0)) stop("A should be a rotated matrix!");
  if(A[nrow(A),2]!=0) stop("A should be a rotated matrix!");
  if (is.vector(b)){
    if(length(b)!=nrow(A)) stop("Dimension problem");
    if (inplace){
      return(bandsolve_cpp(A,as.matrix(b))$x)
    } else {
      Amem=matrix(NA,nrow(A),ncol(A))
      Amem[]=A[]
      bmem=rep(NA,length(b))
      bmem[]=b[]
      return(bandsolve_cpp(Amem,as.matrix(bmem))$x)
    }
  } else if (is.matrix(b)){
    if(nrow(b)!=nrow(A)) stop("Dimension problem");
    if (inplace){
      return(bandsolve_cpp(A,b)$x)
    } else {
      Amem=matrix(NA,nrow(A),ncol(A))
      Amem=A[]
      Bmem=matrix(NA,nrow(b),ncol(b))
      Bmem[]=b[]
      return(bandsolve_cpp(Amem,Bmem)$x)
    }
  } else if (is.null(b)){
    B=diag(nrow(A))
    if (inplace){
      return(bandsolve_cpp(A,B)$x)
    } else {
      Amem=matrix(NA,nrow(A),ncol(A))
      Amem[]=A[]
      return(x=bandsolve_cpp(Amem,B)[[2]])
    }
  } else {
    stop("b must either be a vector or a matrix")
  }
}