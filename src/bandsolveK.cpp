#include <Rcpp.h>
using namespace Rcpp;

//' @title Fast solver for linears systems involving symmetric band matrix of length k.
//' @param Ddata Rotated row-wised matrix of dimensions n*k
//' @param bdata Right hand side of the linear system
//' @return Solution of the linear system.
//' @examples n=2000;
//' D0=runif(n);
//' D1=-0.2*runif(n-1);
//' D=cbind(D0,c(D1,0))
//' b=runif(n)
//' ref=bandsolveK(D,b)
//' @export
// [[Rcpp::export]]
NumericVector bandsolveK(NumericMatrix Ddata, NumericVector bdata) {
    // initialization
    int m=Ddata.ncol()-1;
    int N=Ddata.nrow();
    Rcpp::NumericMatrix L(N,m+1);
    Rcpp::NumericMatrix U(N,m+1);
    Rcpp::NumericVector x(N);
    Rcpp::NumericVector y(N);
    
    for (int i=0;i<N;i++){
      for (int j=0;j<m+1;j++){
        if (j==0){
          L(i,j)=1;
          U(i,j)=Ddata(i,j);
          } 
          else 
          {
          L(i,j)=Ddata(i,j);
          U(i,j)=Ddata(i,j);
          }
      }
    }
    
    // LU decomposition
    for (int i=1; i<N; i++){
      for (int k=1; k<m+1;k++) L(i-1,k)=L(i-1,k)/U(i-1,0);
      for (int j=0; j<m;j++) {
        if (i+j<N) {
          for (int k=0; k<m-j; k++) U(i+j,k)=U(i+j,k)-L(i-1,1+j)*U(i-1,k+j+1);
          if (j<m-1) for (int k=1; k<m-j; k++) L(i+j,k)=L(i+j,k)-L(i-1,k+j+1)*U(i-1,1+j);
        }
      }
    }
    
    // solve y=inv(L)b
    for (int k=0;k<N;k++) y[k]=bdata[k];
      
    for (int i=0; i<N; i++) {
      for (int j=0; j<m; j++)  {
        if (i>j) y[i]=y[i]-L(i-j-1,1+j)*y[i-j-1];
      }
    }
    
    // solve x=inv(U)y=inv(LU)b
    for (int k=0;k<N;k++) x[k]=y[k];
    
    for (int i=(N-1); i>=0; i--) {
      for (int j=0;j<m;j++) {
        if (i+j+1<N) 
          x[i]=x[i]-U(i,1+j)*x[i+j+1];
      }
      x[i]=x[i]/U(i,0);
    }

    return(x);
}