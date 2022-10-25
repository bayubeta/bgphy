#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mcmcmove(NumericVector p, double sigma){
  int d = p.length();
  // draw d samples from normal
  NumericVector step = Rcpp::rnorm(d, 0.0, sigma);

  // add the differences into the original parameter values
  NumericVector prop = p + step;

  return prop;
}
