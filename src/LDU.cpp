#include <Rcpp.h>
using namespace Rcpp;

//' @description Fast inplace LDU decomposition of band matrix.
//' @title LDU
//' @param D Rotated row-wised matrix of dimensions n*(l+1+u), where the l first column corresponding to the sub-diagonal elements, the column l+1 corresponding to the diagonal and the u next to the super-diagonal elements.
//' @param l Lower bandwidth.
//' @param u Upper bandwidth.
//' @name LDU
//' @return List with D as solution of our LDU decomposition.
//' @examples n=2000;
//' D0=0.3*runif(n-1);
//' D1=3+runif(n);
//' D2=-0.2*runif(n-1);
//' D3=runif(n-2)
//' D=cbind(c(D0,0),D1,c(D2,0),c(D3,0,0))
//' res=LDU(D,l=1,u=2)

// [[Rcpp::export]]
List LDU(NumericMatrix D, int l,int u) {
  int n = D.nrow();
  
// LDU in-place 
  for (int i=1; i<=n; i++) {
    for (int t=std::max(std::max(i-u,i-l),1); t<i; t++) D(i-1,l)-=D(t-1,l-i+t)*D(t-1,l)*D(t-1,i-t+l);
    for (int k=1; k<=std::max(u,l); k++) {
      if ((k<=u)&&(i+k<=n)){
        for (int t=std::max(std::max(i-l,i-k-u),1); t<i; t++)  D(i-1,k+l)-=D(t-1,l-i+t)*D(t-1,l)*D(t-1,i+k-t+l);
        D(i-1,k+l)/=D(i-1,l);
      }
      if ((k<=l)&&(i+k<=n)){
        for (int t=std::max(std::max(i-u,i-k-l),1); t<i; t++)  D(i-1,l-k)-=D(t-1,l-i-k+t)*D(t-1,l)*D(t-1,i-t+l);
        D(i-1,l-k)/=D(i-1,l);
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("D")=D);
}
