#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void print_me(IntegerVector x) {
  // several ways to print in Rcpp
  Rprintf("my vector is (%d, %d, %d)\n",
          x[0], x[1], x[2]);
  Rcpp::Rcout << "my vector is " << x << std::endl;
  Rf_PrintValue(x);
}
