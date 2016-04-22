#ifndef LDU_H
#define LDU_H
#include <Rcpp.h>

Rcpp::List LDU(Rcpp::NumericMatrix D, int l, int u);

#endif