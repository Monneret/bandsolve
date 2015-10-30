#include <Rcpp.h>
using namespace Rcpp;

//' @title Fast solver for linears systems involving symmetric band matrix.
//' @param D0data Diagonal as a vector of length n.
//' @param D1Data First subdiagonal and superdiagonal as a vector of length n-1.
//' @param D2Data Second subdiagonal and superdiagonal as a vector of length n-2.
//' @param D3Data Third subdiagonal and superdiagonal as a vector of length n-3.
//' @param D4Data Fourth subdiagonal and superdiagonal as a vector of length n-4.
//' @param bdata Right hand side of the linear system
//' @return Solution of the linear system.
//' @examples n=2000;
//' D0=runif(n);
//' D1=-0.2*runif(n-1);
//' D2=0.4*runif(n-2);
//' D3=-0.1*runif(n-3);
//' D4=0.3*runif(n-4);
//' b=runif(n)
//' ref=bandsolve4(D0,D1,D2,D3,D4,b)
//'  @export
// [[Rcpp::export]]
NumericVector bandsolve4(NumericVector D0data, NumericVector D1data, NumericVector D2data, NumericVector D3data, NumericVector D4data, NumericVector bdata) {
    if (D0data.size()!=(D4data.size()+4))
    stop("we must have length(D0)=length(D3)+4"); 
    if (D0data.size()!=(D3data.size()+3))
    stop("we must have length(D0)=length(D3)+3"); 
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
    Rcpp::NumericVector D4(D4data);
    Rcpp::NumericVector b(bdata);
    int n = D0.size();
    Rcpp::NumericVector L1(n-1);
    Rcpp::NumericVector L2(n-2);
    Rcpp::NumericVector L3(n-3);
    Rcpp::NumericVector L4(n-4);
    Rcpp::NumericVector U0(n);
    Rcpp::NumericVector U1(n-1);
    Rcpp::NumericVector U2(n-2);
    Rcpp::NumericVector U3(n-3);
    Rcpp::NumericVector x(n);
    Rcpp::NumericVector y(n);
    // initialization
    for (int i=0; i<(n-0); i++) U0[i]=D0[i];
    for (int i=0; i<(n-1); i++) U1[i]=D1[i];
    for (int i=0; i<(n-1); i++) L1[i]=D1[i];
    for (int i=0; i<(n-2); i++) U2[i]=D2[i];
    for (int i=0; i<(n-2); i++) L2[i]=D2[i];
    for (int i=0; i<(n-3); i++) L3[i]=D3[i];
    for (int i=0; i<(n-3); i++) U3[i]=D3[i];
    for (int i=0; i<(n-3); i++) L4[i]=D4[i];
    // LU decomposition
    for (int i=1; i<n; i++) {
      // process row i
      L1[i-1]=L1[i-1]/U0[i-1];
      U0[i]=U0[i]-L1[i-1]*U1[i-1];
      
      if (i<n-1) {
        // process row i+1
        L2[i-1]=L2[i-1]/U0[i-1];
        L1[i]=L1[i]-L2[i-1]*U1[i-1];
        
        U1[i]=U1[i]-L1[i-1]*U2[i-1];
        U0[i+1]=U0[i+1]-L2[i-1]*U2[i-1];
      }
      if (i<n-2) {
        L3[i-1]=L3[i-1]/U0[i-1];
        L2[i]=L2[i]-L3[i-1]*U1[i-1];
        L1[i+1]=L1[i+1]-L3[i-1]*U2[i-1];
        
        U2[i]=U2[i]-L1[i-1]*U3[i-1];
        U1[i+1]=U1[i+1]-L2[i-1]*U3[i-1];
        U0[i+2]=U0[i+2]-L3[i-1]*U3[i-1];
      }
      if (i<n-3) {
        L4[i-1]=L4[i-1]/U0[i-1];
        L3[i]=L3[i]-L4[i-1]*U1[i-1];
        L2[i+1]=L2[i+1]-L4[i-1]*U2[i-1];
        L1[i+2]=L1[i+2]-L4[i-1]*U3[i-1];
        
        U3[i]=U3[i]-L1[i-1]*D4[i-1];
        U2[i+1]=U2[i+1]-L2[i-1]*D4[i-1];
        U1[i+2]=U1[i+2]-L3[i-1]*D4[i-1];
        U0[i+3]=U0[i+3]-L4[i-1]*D4[i-1];
    }
    }
    
    // solve y=inv(L)b
    y[0]=b[0];
    y[1]=b[1]-L1[0]*y[0];
    y[2]=b[2]-L2[0]*y[0]-L1[1]*y[1];
    y[3]=b[3]-L3[0]*y[0]-L2[1]*y[1]-L1[2]*y[2];
    for (int i=4; i<n; i++) 
    y[i]=b[i]-L1[i-1]*y[i-1]-L2[i-2]*y[i-2]-L3[i-3]*y[i-3]-L4[i-4]*y[i-4];
    // solve x=inv(U)y=inv(LU)b
    x[n-1]=y[n-1]/U0[n-1];
    x[n-2]=(y[n-2]-U1[n-2]*x[n-1])/U0[n-2];
    x[n-3]=(y[n-3]-U1[n-3]*x[n-2]-U2[n-3]*x[n-1])/U0[n-3];
    x[n-4]=(y[n-4]-U1[n-4]*x[n-3]-U2[n-4]*x[n-2]-U3[n-4]*x[n-1])/U0[n-4];
    for (int i=(n-5); i>=0; i--)
    x[i]=x[i]=(y[i]-U1[i]*x[i+1]-U2[i]*x[i+2]-U3[i]*x[i+3]-D4[i]*x[i+4])/U0[i];
    return(x);
}
