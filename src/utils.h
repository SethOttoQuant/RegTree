#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
arma::field<arma::uvec> select_rnd(arma::uword m, // number of elements
                                   arma::uword n);
arma::field<arma::vec> QuickReg(arma::vec x, // no missing
                                arma::vec y, // no missing
                                double r);
double welch_t(double m_1,  // mean 1
               double m_2, // mean 2
               double sig_1,  // squared error 1
               double sig_2,  // squared error 2
               double n_1, // n obs 1
               double n_2 // n obs 2
);
#endif

