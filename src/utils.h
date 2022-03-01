#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
arma::uvec select_rnd(arma::uword m, // number of elements
                      arma::uword n);
arma::field<arma::vec> QuickReg(arma::vec x, // no missing
                                arma::vec y, // no missing
                                double r);
#endif

