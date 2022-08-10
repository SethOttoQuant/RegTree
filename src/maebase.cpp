// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace arma;
using namespace Rcpp;
#include "utils.h"
#include "rfbase.h"

// [[Rcpp::export]]
arma::vec fast_mae_cut(arma::vec x, // predictor
                   arma::vec y){ // response
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  double vnce_na = 0;  // actually MAE here... didn't change variable name
  if(not_obs.n_elem>0) vnce_na = sum(abs(y(not_obs) - median(y)));
  if(obs.n_elem<2){
    vec out = {0, 0, vnce_na, 0}; // no observations; not useful
    return(out);
  }
  x = x(obs);
  y = y(obs);
  uvec idx = sort_index(x); // ascending; may be duplicate values here (doesn't matter)
  vec x_sorted = x(idx);
  vec yx = y(idx); // sort y by x values
  uword n = yx.n_elem; // number of obs
  mat vnce(2,n-1); // store MAE resulting from split; again didn't change var name
  for(uword j=0; j<n-1; j++){
    vnce(0,j) = sum(abs(yx(span(0,j)) - median(yx(span(0,j))))); // <= variance
    vnce(1,j) = sum(abs(yx(span(j+1,n-1)) - median(yx(span(j+1,n-1))))); // > variance
  }
  double imin = index_min(sum(vnce,0));
  vec out = {vnce(0,imin), vnce(1,imin), vnce_na, imin};
  return(out); // less than cut mean, greater than cut mean, less variance, greater variance, NA variance, and cut. 
}
// [[Rcpp::export]]
arma::field<arma::vec> mae_cut(arma::vec x, // predictor
                       arma::vec y, // response
                       arma::uvec ind, // index
                       double imin){  // where to cut
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  uword n = obs.n_elem;
  uword nnot = not_obs.n_elem;
  field<vec> out(3); // where output goes
  double vnce_na = 0;  // actually MAE
  if(nnot>0) vnce_na = sum(abs(y(not_obs) - median(y))); // E(y) with no new data 
  if(n<5){
    vec par = {0, 0, 0 , 0, sum(abs(y - median(y))), 0, 0, 0}; // not enough observations; not useful
    out(0) = par;
    out(1) = conv_to<vec>::from(ind); // index of conditions unchanged
    out(2) = conv_to<vec>::from(ind); // index of conditions unchanged
    return(out);
  }
  x = x(obs);
  y = y(obs);
  ind = ind(obs);
  uvec idx = sort_index(x); // ascending; may be duplicate values here (doesn't matter)
  x = x(idx); // sorted by x values
  y = y(idx); // sort y by x values
  double cut = (x(imin) + x(imin+1))/2;
  vec pars = {median(y(span(0,imin))), median(y(span(imin+1, n-1))), 
              sum(abs(y(span(0,imin)) - median(y(span(0,imin))))), 
              sum(abs(y(span(imin+1,n-1)) - median(y(span(imin+1,n-1))))), 
              vnce_na, cut, imin+1, n-imin-1};
  out(0) = pars; // mean <=, mean >, sig <=, sig >, sig NA, cut, obs <=, obs > 
  out(1) = conv_to<vec>::from(ind(idx(span(0,imin)))); // index of original dataset <=
  out(2) = conv_to<vec>::from(ind(idx(span(imin+1,n-1)))); // index of original dataset >
  return(out); 
}
// find the best series in X to identify y
// [[Rcpp::export]]
arma::field<arma::vec> maesplit(arma::mat X, // predictors
                       arma::vec y, // response
                       arma::uvec ind, // index of observations in original data
                       arma::uword n){ // number of candidates to select
  // select candidates for splitting
  field<uvec> cnd = select_rnd(X.n_cols, n);
  uvec candidates = cnd(0); // randomly selected candidates
  mat splits(4,n); // where output goes
  vec tmp;
  for(uword j=0; j<n; j++){
    tmp = fast_mae_cut(X.col(candidates(j)), y); 
    splits.col(j) = tmp; // pars
  }
  vec tot_var = trans(sum(splits.rows(0,2),0)); // rows 2:4 contain MAE for high, low, and NA values
  uword min_idx = index_min(tot_var); // which series offers max reduction in MAE of y
  field<vec> out = mae_cut(X.col(candidates(min_idx)), y, ind, splits(3, min_idx)); // re-running cut here to keep things clean
  // clunky but effective
  vec par_out(9);
  par_out(0) = candidates(min_idx);
  par_out(span(1,8)) = out(0);
  out(0) = par_out;
  return(out); // for find_cut() see function above
}
// [[Rcpp::export]]
arma::mat maetree(arma::vec y, // response (no missing obs)
          arma::mat X, // predictors (missing obs OK)
          arma::uvec to_keep,
          arma::uword min_obs = 5, 
          arma::uword max_nodes= 1000){
  X = X.rows(to_keep); // this shuffles X and y but it shouldn't matter
  y = y(to_keep);
  double xnc = X.n_cols;
  uword n = ceil(xnc/3); // number of candidates to use at each split
  mat Tree(max_nodes, 9, fill::zeros);
  Tree(0,5) = median(y); // unconditional median
  Tree(0,6) = sum(abs(y-Tree(0,5))); // unconditional MAE
  Tree(0,7) = to_keep.n_elem;
  mat leaves; vec par; uvec fobs; mat tmp_tree;
  uvec leaf_idx; 
  field<uvec> I(max_nodes); // store index of values corresponding to each node
  I(0) = regspace<uvec>(0,X.n_rows-1); // index in original data (after bagging)
  field<vec> tmp; // temp output from maesplits()
  uword i = 0;
  uword j = 0;
  while(i+2<max_nodes){
    tmp = maesplit(X.rows(I(j)), y(I(j)), I(j), n);
    par = tmp(0);
    if(any(par(span(1,4)))){ //  enough observations to split
      Tree(j,0) = par(0); // which variable
      Tree(j,1) = par(6); // cut value
      Tree(j,2) = i+1; // if <, go to node i+1
      Tree(j,3) = i+2; // if >, go to node i+2
      Tree(j,8) = 0; // no longer terminal node
      Tree(i+1,4) = j; // node coming from
      Tree(i+2,4) = j; // node coming from
      Tree(i+1,5) = par(1); // mu if <
      Tree(i+2,5) = par(2); // mu if >
      Tree(i+1,6) = par(3); // sig if <
      Tree(i+2,6) = par(4); // sig if >
      Tree(i+1,7) = par(7); // nobs if <
      Tree(i+2,7) = par(8); // nobs if >
      Tree(i+1,8) = 1; // terminal node (leaf)
      Tree(i+2,8) = 1; // terminal node (leaf)
      I(i+1) = conv_to<uvec>::from(tmp(1)); // indexes <=
      I(i+2) = conv_to<uvec>::from(tmp(2)); // indexes >
      i += 2; // split result in 2 new rows
    }else{
      Tree(j,8) = 2; // terminal node with no more possible splits
    }
    tmp_tree = Tree.rows(0,i);
    leaf_idx = find(tmp_tree.col(8) == 1 && tmp_tree.col(7) > min_obs); // index of terminal nodes in Tree matrix with enough obs
    leaves = Tree.rows(leaf_idx); // terminal nodes (i.e. leaves). 
    if(leaf_idx.n_elem==0) break;
    j = leaf_idx(index_max(leaves.col(6))); // index of leaf with max volatility to work on next
  }
  return(Tree.rows(0,i));
}

// Draw 'draws' number of trees
// [[Rcpp::export]]
Rcpp::List maeforest(arma::vec y, // response (no missing obs)
                   arma::mat X, // predictors (missing obs OK)
                   arma::uword min_obs = 5,
                   arma::uword max_nodes = 1000,
                   arma::uword draws = 1000){
  field<mat> Trees(draws);
  double T = X.n_rows;
  vec oob(T); mat Tree;
  field<uvec> to_keep;
  mat OOB(T, draws);
  mat fc(T,X.n_cols);  // out of bag feature contributions
  cube FC(T, X.n_cols, draws);
  vec mse(draws);
  field<mat> tmp;
  for(uword j = 0; j<draws; j++){
    to_keep = select_rnd(T, ceil(0.632*T));
    Tree = maetree(y, X, to_keep(0), min_obs, max_nodes);
    tmp = fitmat(trans(X.rows(to_keep(1))), Tree);
    oob.fill(datum::nan);
    oob(to_keep(1)) = tmp(0);
    mse(j) = mean(square(y(to_keep(1)) - tmp(0)));
    fc.fill(datum::nan);
    fc.rows(to_keep(1)) = trans(tmp(1));
    OOB.col(j) = oob;
    FC.slice(j) = fc;
    Trees(j) = Tree;
  }
  Rcpp::List rtrn;
  rtrn["Trees"] = Trees;
  rtrn["OOB"] = OOB; // out of bag fit
  rtrn["FC"] = FC;
  rtrn["MSE"] = mse;
  return(rtrn);
}

