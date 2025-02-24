#ifndef RFBASE_H
#define RFBASE_H

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
arma::field<arma::vec> fitvec(arma::vec x,
                              arma::mat Tree,
                              arma::uword maxit = 1000);
arma::field<arma::vec> poolvec(arma::vec x,
                              arma::mat Tree,
                              arma::uword maxit = 1000);
arma::field<arma::mat> fitmat(arma::mat X,
                              arma::mat Tree,
                              bool weight=false);
arma::field<arma::mat> fitfield(arma::mat X,
                                arma::field<arma::mat> Trees,
                                bool weight=false);
#endif

