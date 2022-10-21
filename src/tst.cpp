// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace arma;
using namespace Rcpp;
#include "utils.h"

// [[Rcpp::export]]
Rcpp::List tst_sort(arma::vec x){
  uvec idx = sort_index(x);
  uvec unsort_idx = sort_index(idx);
  vec x_sorted = x(idx);
  vec x_unsort = x_sorted(unsort_idx);
  Rcpp::List out(3);
  out[0] = x_sorted;
  out[1] = x_unsort;
  out[2] = idx;
  return out;
}

// [[Rcpp::export]]
arma::vec select_cut(arma::vec x, // predictor
                     arma::vec y, // response
                     double prior_shrink, // prior shrinkage par
                     double my=0,
                     double weight_pow=1){ // prior for y
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  double m_na=0, vnce_na = 0, m_leq=0, m_g=0;
  vec tmp_out(2, fill::zeros);
  if(not_obs.n_elem>0){
    m_na = (sum(y(not_obs)) + prior_shrink*my)/(not_obs.n_elem + prior_shrink);
    vnce_na = sum(square(y(not_obs) - m_na))/(not_obs.n_elem + prior_shrink);
  }
  if(obs.n_elem<(4)){  // hard coded min node size of 2
    return(tmp_out);  // never split on this guy; probability of selection is 0
  }
  x = x(obs);
  y = y(obs);
  uvec idx = sort_index(x); // ascending; may be duplicate values here (doesn't matter)
  // vec x_sorted = x(idx);
  vec yx = y(idx); // sort y by x values
  uword n = yx.n_elem; // number of obs
  uword idx_j;
  mat out(2, n-3, fill::zeros);
  for(uword j=0; j<=(n-4); j++){
    idx_j = j + 1;
    m_leq = (sum(yx(span(0,idx_j))) + prior_shrink*my)/(idx_j+1 + prior_shrink);
    m_g = (sum(yx(span(idx_j+1,n-1))) + prior_shrink*my)/(n-idx_j-1 + prior_shrink);
    out(0,j) = sum(square(yx(span(0,idx_j)) - m_leq))/(idx_j + 1 + prior_shrink) +
      sum(square(yx(span(idx_j+1,n-1)) - m_g))/(n-idx_j-1 + prior_shrink) + 
      vnce_na; // control prob of being selected. 
    out(1,j) = idx_j; // cut value index
  }
  // this bit is just to select the split with probability weight given by out.row(0)
  vec cdf = out.row(0).t();
  cdf /= mean(cdf);
  cdf = cumsum(pow(cdf, -weight_pow));
  double draw = randu()*cdf(n-4); // rand u var over full range of cdf
  uword j=0;
  while(draw>cdf(j)){
    j++;
  }
  return(out.col(j)); // mat with row 0 = variance, row 1 = cut value
}
// [[Rcpp::export]]
arma::field<arma::vec> cutpars(arma::vec x, // predictor
                               arma::vec y, // response
                               arma::uvec ind, // original index
                               double imin,
                               double my =0,
                               double prior_shrink=0){  // where to cut
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  double m_na=0, vnce_na = 0, m_leq=0, m_g=0;
  vec tmp_out(2, fill::zeros);
  if(not_obs.n_elem>0){
    m_na = (sum(y(not_obs)) + prior_shrink*my)/(not_obs.n_elem + prior_shrink);
    vnce_na = sum(square(y(not_obs) - m_na))/(not_obs.n_elem + prior_shrink);
  }
  field<vec> out(3);
  uword n = obs.n_elem; // number of obs
  x = x(obs);
  y = y(obs);
  ind = ind(obs);  // original index
  uvec idx = sort_index(x); // ascending; may be duplicate values here (doesn't matter)
  x = x(idx); // sorted by x values
  y = y(idx); // sort y by x values
  double cut = x(imin);
  m_leq = (sum(y(span(0,imin))) + prior_shrink*my)/(imin+1 + prior_shrink);
  m_g = (sum(y(span(imin+1,n-1))) + prior_shrink*my)/(n-imin-1 + prior_shrink);
  vec pars = {m_leq, m_g, 
              sum(square(y(span(0,imin)) - m_leq))/(imin + 1 + prior_shrink), 
              sum(square(y(span(imin+1,n-1)) - m_g))/(n-imin-1 + prior_shrink), 
              vnce_na, cut, imin+1, n-imin-1};
  out(0) = pars; // mean <=, mean >, sig <=, sig >, sig NA, cut, obs <=, obs > 
  out(1) = conv_to<vec>::from(ind(idx(span(0,imin)))); // index of original dataset <=
  out(2) = conv_to<vec>::from(ind(idx(span(imin+1,n-1)))); // index of original dataset >
  return(out); 
}
// find the best series in X to identify y
// [[Rcpp::export]]
arma::field<arma::vec> selectsplit(arma::mat X, // predictors
                                   arma::vec y, // response
                                   arma::uvec ind, // index of observations in original data
                                   double prior_shrink=1, // shrink estimate towards unconditional mean
                                   double my=0,   // unconditional mean
                                   double weight_pow=1){
  arma::mat candidates(2, X.n_cols);
  for(uword j=0; j<X.n_cols; j++){
    // Rcpp::Rcout << candidates(j) << endl;
    candidates.col(j) = select_cut(X.col(j), y, prior_shrink, my); 
  }
  Rcpp::Rcout << candidates << endl;
  // this bit is just to select the split with probability weight given by out.row(0)
  vec cdf = candidates.row(0).t();
  cdf /= mean(cdf);
  cdf = cumsum(pow(cdf, -weight_pow));
  double draw = randu()*cdf(cdf.n_elem-1); // rand u var over full range of cdf
  uword j=0;
  Rcpp::Rcout << cdf << endl;
  Rcpp::Rcout << draw << endl;
  while(draw>cdf(j)){
    Rcpp::Rcout << j << endl;
    j++;
  }
  field<vec> out = cutpars(X.col(j), y, ind, candidates(1,j), my, prior_shrink); // re-running cut here to keep things clean
  // clunky but effective
  vec par_out(9);
  par_out(0) = j;
  par_out(span(1,8)) = out(0);
  out(0) = par_out;
  return(out); // for find_cut() see function above
}

// // [[Rcpp::export]]
// arma::vec fast_cut(arma::vec x, // predictor
//                    arma::vec y, // response
//                    double prior_shrink, // prior shrinkage par
//                    double my=0,
//                    double weight_pow=1){ // prior for y
//   uvec obs = find_finite(x);
//   uvec not_obs = find_nonfinite(x);
//   double m_na=0, vnce_na = 0, m_leq=0, m_g=0;
//   vec tmp_out(2, fill::zeros);
//   if(not_obs.n_elem>0){
//     m_na = (sum(y(not_obs)) + prior_shrink*my)/(not_obs.n_elem + prior_shrink);
//     vnce_na = sum(square(y(not_obs) - m_na))/(not_obs.n_elem + prior_shrink);
//   }
//   if(obs.n_elem<(4)){  // hard coded min node size of 2
//     return(tmp_out);  // never split on this guy; probability of selection is 0
//   }
//   x = x(obs);
//   y = y(obs);
//   uvec idx = sort_index(x); // ascending; may be duplicate values here (doesn't matter)
//   vec x_sorted = x(idx);
//   vec yx = y(idx); // sort y by x values
//   uword n = yx.n_elem; // number of obs
//   uword idx_j;
//   mat out(2, n-3, fill::zeros);
//   for(uword j=0; j<=(n-4); j++){
//     idx_j = j + 1;
//     // Rcpp::Rcout << j << endl;
//     m_leq = (sum(yx(span(0,idx_j))) + prior_shrink*my)/(idx_j+1 + prior_shrink);
//     m_g = (sum(yx(span(idx_j+1,n-1))) + prior_shrink*my)/(n-idx_j-1 + prior_shrink);
//     out(0,j) = sum(square(yx(span(0,idx_j)) - m_leq))/(idx_j + 1 + prior_shrink) +
//       sum(square(yx(span(idx_j+1,n-1)) - m_g))/(n-idx_j-1 + prior_shrink) + 
//       vnce_na; // control prob of being selected. 
//     out(1,j) = x_sorted(idx_j); // cut value
//   }
//   // this bit is just to select the split with probability weight given by out.row(0)
//   vec cdf = out.row(0).t();
//   Rcpp::Rcout << cdf << "cdf" << endl;
//   cdf /= mean(cdf);
//   cdf = cumsum(pow(cdf, -weight_pow));
//   Rcpp::Rcout << cdf << endl;
//   double draw = randu()*cdf(n-4); // rand u var over full range of cdf
//   uword j=0;
//   Rcpp::Rcout << draw << endl;
//   while(draw>cdf(j)){
//     j++;
//   }
//   Rcpp::Rcout << out.col(j) << endl;
//   return(out.col(j)); // mat with row 0 = variance, row 1 = cut value
// }