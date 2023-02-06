#include <RcppDist.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;





// Update previous covariance matrix with current data
// [[Rcpp::export]]
arma::mat covUpdate(arma::vec x, arma::mat S, arma::vec mean_n, int n){
  arma::vec delta = x - mean_n;
  return (n-1)*S/n + arma::kron(delta, delta.as_row())/(n+1);
}


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
arma::mat adaptMCMC(arma::vec X, arma::mat S, double beta, double sigma){

  // obtain d, number of dimensions
  double d = X.n_elem;

  // covariance of adaptive move
  arma::mat Sadapt = ((pow(2.38, 2))/d)*S;

  // adaptive move
  arma::mat Xadapt = (1-beta)*rmvnorm(1, X, Sadapt);

  // fixed move
  arma::mat Xfixed = beta * fixedMCMC(X, sigma);


  return Xadapt + Xfixed;
}

