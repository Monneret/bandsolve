#include <Rcpp.h>
using namespace Rcpp;

//' @title Fast solver for linears systems involving symmetric band matrix.
//'  @export
//'  @rdname bandsolve4
// [[Rcpp::export]]
Rcpp::List bandsolve1(NumericVector D0data, NumericVector D1data, NumericVector bdata) {
    if (D0data.size()!=(D1data.size()+1))
    Rcpp::stop("we must have length(D0)=length(D1)+1");
    if (D0data.size()!=bdata.size())
    Rcpp::stop("we must have length(D0)=length(b)");
    
    Rcpp::NumericVector D0(D0data);
    Rcpp::NumericVector D1(D1data);
    Rcpp::NumericVector b(bdata);
    int n = D0.size();
    Rcpp::NumericVector L1(n-1);
    Rcpp::NumericVector U0(n);
    Rcpp::NumericVector x(n);
    Rcpp::NumericVector y(n);
    // LU decomposition
    U0[0]=D0[0];
    
    for (int i=0; i<(n-1); i++) {
      L1[i]=D1[i]/U0[i];
      U0[i+1]=D0[i+1]-L1[i]*D1[i];
    }
    
    // solve y=inv(L)b
    
    y[0]=b[0];
    for (int i=1; i<n; i++) 
      y[i]=(b[i]-L1[i-1]*y[i-1]);
    // solve x=inv(U)y=inv(LU)b
    
    x[n-1]=y[n-1]/U0[n-1];
    for (int i=(n-1); i>0; i--)
    x[i-1]=(y[i-1]-D1[i-1]*x[i])/U0[i-1];
    return Rcpp::List::create(Rcpp::Named("x")=x,Rcpp::Named("L1")=L1,Rcpp::Named("U0")=U0);
}
