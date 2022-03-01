// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace arma;

// [[Rcpp::export]]
arma::uvec select_rnd(arma::uword m, // number of elements
                      arma::uword n){ // number to select
  uvec idx = regspace<uvec>(0,m-1);
  uvec out = shuffle(idx);
  return(out.head(n));
}

// [[Rcpp::export]]
arma::field<arma::vec> QuickReg(arma::vec x, // no missing
                                arma::vec y, // no missing
                                double r){ // ridge parameter
  double n = x.n_elem+r;
  vec par(3);
  double sx = sum(x);
  double sy = sum(y);
  double xx = sum(square(x)) + r;
  double xy = as_scalar(trans(x)*y);
  double dt = n*xx - sx*sx;
  vec ab(2);
  ab(0) = (xx*sy - sx*xy)/dt;
  ab(1) = (-sx*sy + n*xy)/dt;
  vec e = y - ab(0) - ab(1)*x;
  double sig = sum(square(e));
  field<vec> out(2);
  par(span(0,1)) = ab; // ab.col(0);
  par(2) = sig;
  out(0) = par; // parameters: intercept, beta, sigma
  out(1) = e; // residuals
  return(out);
}