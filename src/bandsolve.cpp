#include <Rcpp.h>
#include "LDU.h"
#include "LDL.h"
using namespace Rcpp;

//' @description Fast solver for linears systems involving band matrix.
//' @title bandsolve
//' @param D Rotated row-wised matrix.
//' @name bandsolve
//' @param b Right hand side of the linear system
//' @return List with components with the solution of the linear system and the rotated matrixes L and U.
//' @examples n=2000;
//' D0=3+runif(n);
//' D1=-0.2*runif(n-1);
//' D=cbind(D0,c(D1,0))
//' b=runif(n)
//' sol=bandsolve(D,b,l=1,u=0)
//' @export

// [[Rcpp::export]]
List bandsolve(NumericMatrix D, int l, int u, NumericVector b,bool sym=true) {
      int n = D.nrow();
      
  if (sym==true) { //LDL case
    LDL(D);
    int K = D.ncol()-1;
    
    for (int i=2; i<=n; i++) { //solve b=inv(L)b
    int jmax=i-1;
    if (jmax>K) jmax=K;
    for (int j=1; j<=jmax; j++)
    b[i-1]-=D(i-j-1,j)*b[i-j-1];
  }
    // solve b=b/D
  for (int i=0; i<n; i++) b[i]/=D(i,0);
    // solve b=inv(t(L))b=inv(L*D*t(L))b
  for (int i=n-1; i>=1; i--) {
    int jmax=n-i;
    if (jmax>K) jmax=K;
    for (int j=1; j<=jmax; j++)
    b[i-1]-=D(i-1,j)*b[i+j-1];
  }
  }
  else
  { // LDU case
    LDU(D,l,u);
    
    for (int i=2; i<=n; i++) { //solve b=inv(L)b
    for (int k=std::max(1,i-l); k<=i-1; k++)
    b[i-1]-=D(k-1,l-i+k)*b[k-1]; //b_i-=sum_i-l^i-1 Lik bk
  }
    // solve b=b/D
  for (int i=0; i<n; i++) b[i]/=D(i,l);
    // solve b=inv(U)b=inv(L*D*U)b
  for (int i=n-1; i>=1; i--) {
    for (int k=std::min(n,i+u); k>=i+1; k--)
    b[i-1]-=D(i-1,k-i+l)*b[k-1]; //b_i-=sum_i+1^i+u Uik bk
    
  }
}
  
  return Rcpp::List::create(Rcpp::Named("D")=D,Rcpp::Named("x")=b);
}
