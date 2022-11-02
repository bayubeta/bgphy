#include <RcppDist.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// MCMC move with a fixed step size
// [[Rcpp::export]]
arma::mat fixedMCMC(arma::vec X, double sigma){

  // obtain d, number of dimensions
  double d = X.n_elem;

  // covariance of fixed move
  arma::mat Sfixed = ((pow(sigma, 2))/d) * arma::eye(d,d);

  // fixed move
  arma::mat Xfixed = rmvnorm(1, X, Sfixed);


  return Xfixed;
}



// adaptive MCMC move given a set of observations
// [[Rcpp::export]]
arma::mat adaptMCMC(arma::mat X, double beta, double sigma){
  // obtain n, number of observations
  int n = X.n_rows;

  // obtain d, number of dimensions
  double d = X.n_cols;

  // set the last value as the current
  arma::vec Xn = vectorise(X.row(n-1));

  // covariance of adaptive move
  arma::mat S = ((pow(2.38, 2))/d)*cov(X);

  // adaptive move
  arma::mat Xadapt = (1-beta)*rmvnorm(1, Xn, S);


  // fixed move
  arma::mat Xfixed = beta * fixedMCMC(Xn, sigma);


  return Xadapt + Xfixed;
}


/*** R
fixedMCMC(P1n, 0.01)
adaptMCMC(P1, 0.05, 0.01)
*/

