#include <Rcpp.h>
using namespace Rcpp;

//'  @export
//'  @rdname bandsolve4
// [[Rcpp::export]]
List bandsolve3(NumericVector D0data, NumericVector D1data, NumericVector D2data, NumericVector D3data, NumericVector bdata) {
    if (D0data.size()!=(D3data.size()+3))
    Rcpp::stop("we must have length(D0)=length(D3)+3"); 
    if (D0data.size()!=(D2data.size()+2))
    Rcpp::stop("we must have length(D0)=length(D2)+2");
    if (D0data.size()!=(D1data.size()+1))
    Rcpp::stop("we must have length(D0)=length(D1)+1");
    if (D0data.size()!=bdata.size())
    Rcpp::stop("we must have length(D0)=length(b)");
    
    Rcpp::NumericVector D0(D0data);
    Rcpp::NumericVector D1(D1data);
    Rcpp::NumericVector D2(D2data);
    Rcpp::NumericVector D3(D3data);
    Rcpp::NumericVector b(bdata);
    int n = D0.size();
    Rcpp::NumericVector L1(n-1);
    Rcpp::NumericVector L2(n-2);
    Rcpp::NumericVector L3(n-3);
    Rcpp::NumericVector U0(n);
    Rcpp::NumericVector U1(n-1);
    Rcpp::NumericVector U2(n-2);
    Rcpp::NumericVector x(n);
    Rcpp::NumericVector y(n);
    // initialization
    for (int i=0; i<(n-0); i++) U0[i]=D0[i];
    for (int i=0; i<(n-1); i++) U1[i]=D1[i];
    for (int i=0; i<(n-1); i++) L1[i]=D1[i];
    for (int i=0; i<(n-2); i++) U2[i]=D2[i];
    for (int i=0; i<(n-2); i++) L2[i]=D2[i];
    for (int i=0; i<(n-3); i++) L3[i]=D3[i];
    
    // LU decomposition
    for (int i=1; i<n; i++) {
      // process row i
      L1[i-1]=L1[i-1]/U0[i-1];
      U0[i]=U0[i]-L1[i-1]*U1[i-1];
      if (i<n-1) {
        // process row i+1
        L2[i-1]=L2[i-1]/U0[i-1];
        L1[i]=L1[i]-L2[i-1]*U1[i-1];
        U0[i+1]=U0[i+1]-L2[i-1]*U2[i-1];
        U1[i]=U1[i]-L1[i-1]*U2[i-1];
      }
      if (i<n-2) {
        L3[i-1]=L3[i-1]/U0[i-1];
        L2[i]=L2[i]-L3[i-1]*U1[i-1];
        L1[i+1]=L1[i+1]-L3[i-1]*U2[i-1];
        U0[i+2]=U0[i+2]-L3[i-1]*D3[i-1];
        U1[i+1]=U1[i+1]-L2[i-1]*D3[i-1];
        U2[i]=U2[i]-L1[i-1]*D3[i-1];
      }
    }
    
    // solve y=inv(L)b
    y[0]=b[0];
    y[1]=b[1]-L1[0]*y[0];
    y[2]=b[2]-L2[0]*y[0]-L1[1]*y[1];
    for (int i=3; i<n; i++) 
    y[i]=b[i]-L1[i-1]*y[i-1]-L2[i-2]*y[i-2]-L3[i-3]*y[i-3];
    // solve x=inv(U)y=inv(LU)b
    x[n-1]=y[n-1]/U0[n-1];
    x[n-2]=(y[n-2]-U1[n-2]*x[n-1])/U0[n-2];
    x[n-3]=(y[n-3]-U1[n-3]*x[n-2]-U2[n-3]*x[n-1])/U0[n-3];
    for (int i=(n-4); i>=0; i--)
    x[i]=(y[i]-U1[i]*x[i+1]-U2[i]*x[i+2]-D3[i]*x[i+3])/U0[i];
    return Rcpp::List::create(Rcpp::Named("x")=x,Rcpp::Named("L1")=L1,Rcpp::Named("L2")=L2,Rcpp::Named("L3")=L3,Rcpp::Named("U0")=U0,Rcpp::Named("U1")=U1,Rcpp::Named("U2")=U2);
}
