// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace arma;

//[[Rcpp::export]]
arma::uvec sort_em(arma::vec v
){
  uvec out = find(v);
  return(out);
}

//[[Rcpp::export]]
double welch_t(double m_1,  // mean 1
                  double m_2, // mean 2
                  double sig_1,  // squared error 1
                  double sig_2,  // squared error 2
                  double n_1, // n obs 1
                  double n_2 // n obs 2
){
  double s_1 = sig_1/(n_1-1);
  double s_2 = sig_2/(n_2-1);
  double s_delta = std::sqrt(s_1/n_1 + s_2/n_2);
  double t = abs(m_1 - m_2)/s_delta;
  double df = std::pow(s_delta,4)/(std::pow(s_1,2)/(std::pow(n_1,2)*(n_1-1)) + std::pow(s_2,2)/(std::pow(n_2,2)*(n_2-1)));
  return R::pt(t, df, 0, 0); // area above t-stat
  // return normcdf(t);
}

//[[Rcpp::export]]
double t_tst(double t,  // t stat
             double df // degree freedom
){
  double out = R::pt(t, df, 0, 0);  // area above t-stat
  return out;
}

// [[Rcpp::export]]
arma::field<arma::uvec> select_rnd(arma::uword m, // number of elements
                                   arma::uword n){ // number to select
  uvec idx = regspace<uvec>(0,m-1);
  uvec to_keep = shuffle(idx);
  field<uvec> out(2);
  out(0) = to_keep.head(n);
  out(1) = to_keep.tail(m-n);
  return(out);
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

// // [[Rcpp::export]]
// arma::mat var_distributions(arma::vec y,  // rhs data
//                             arma::uword min_obs,  // min_obs per node
//                             arma::uword max_obs,  // max_obs at first node
//                             arma::uword draws){  
//   
//   uvec obs = regspace<uvec>(min_obs, max_obs);
//   uword k = obs.n_elem;
//   mat out(k,1000,fill::zeros);
//   
// }

