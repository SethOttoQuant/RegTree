// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace arma;
using namespace Rcpp;
// #include "utils.h"

// [[Rcpp::export]]
arma::uvec select_rnd(arma::uword m, // number of elements
                      arma::uword n){ // number to select
  uvec idx = regspace<uvec>(0,m-1);
  uvec out = shuffle(idx);
  return(out.head(n));
}

// [[Rcpp::export]]
arma::vec fast_cut(arma::vec x, // predictor
                   arma::vec y){ // response
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  double vnce_na = 0;
  if(not_obs.n_elem>0) vnce_na = sum(square(y(not_obs) - mean(y)));
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
  mat vnce(2,n-1); // store variance resulting from split
  for(uword j=0; j<n-1; j++){
    vnce(0,j) = sum(square(yx(span(0,j)) - mean(yx(span(0,j))))); // <= variance
    vnce(1,j) = sum(square(yx(span(j+1,n-1)) - mean(yx(span(j+1,n-1))))); // > variance
  }
  double imin = index_min(sum(vnce,0));
  vec out = {vnce(0,imin), vnce(1,imin), vnce_na, imin};
  return(out); // less than cut mean, greater than cut mean, less variance, greater variance, NA variance, and cut. 
}

// [[Rcpp::export]]
arma::field<arma::vec> cut(arma::vec x, // predictor
                       arma::vec y, // response
                       arma::uvec ind, // index
                       double imin){  // where to cut
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  uword n = obs.n_elem;
  uword nnot = not_obs.n_elem;
  field<vec> out(3);
  double vnce_na = 0;
  if(nnot>0) vnce_na = sum(square(y(not_obs) - mean(y))); // E(y) with no new data 
  if(n<5){
    vec par = {0, 0, 0 , 0, sum(square(y - mean(y))), 0, 0, 0}; // not enough observations; not useful
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
  // Rcpp::Rcout << imin << endl;// [[Rcpp::export]]
  vec pars = {mean(y(span(0,imin))), mean(y(span(imin+1, n-1))), 
              sum(square(y(span(0,imin)) - mean(y(span(0,imin))))), 
              sum(square(y(span(imin+1,n-1)) - mean(y(span(imin+1,n-1))))), 
              vnce_na, cut, imin+1, n-imin-1};
  out(0) = pars; // a<=, b<=, a>, b>, sig<=, sig>, sig NA, c
  out(1) = conv_to<vec>::from(ind(idx(span(0,imin)))); // index of original dataset <=
  out(2) = conv_to<vec>::from(ind(idx(span(imin+1,n-1)))); // index of original dataset >
  return(out); 
}

// find the best series in X to identify y
// [[Rcpp::export]]
arma::field<arma::vec> bestsplit(arma::mat X, // predictors
                       arma::vec y, // response
                       arma::uvec ind, // index of observations in original data
                       arma::uword n){ // number of candidates to select
  // select candidates for splitting
  uvec candidates = select_rnd(X.n_cols, n); // randomly selected candidates
  mat splits(4,n); 
  vec tmp;
  for(uword j=0; j<n; j++){
    tmp = fast_cut(X.col(candidates(j)), y); 
    splits.col(j) = tmp; // pars
  }
  vec tot_var = trans(sum(splits.rows(0,2),0)); // rows 2:4 contain variance for high, low, and NA values
  uword min_idx = index_min(tot_var); // which series offers max reduction in variance of y
  field<vec> out = cut(X.col(candidates(min_idx)), y, ind, splits(3, min_idx)); // re-running find_cut here to keep things clean
  // clunky but effective
  vec par_out(9);
  par_out(0) = candidates(min_idx);
  par_out(span(1,8)) = out(0);
  out(0) = par_out;
  return(out); // for find_cut() see function above
}


arma::uvec field_obs(arma::field<uvec> I, arma::uword n){
  uvec out(n);
  uvec tmp;
  for(uword j=0; j<n; j++){
    tmp = I(j);
    out(j) = tmp.n_elem;
  }
  return(out);
}

// Core function to call
// [[Rcpp::export]]
arma::mat RegTree(arma::vec y, // response (no missing obs)
          arma::mat X, // predictors (missing obs OK)
          arma::uword min_obs = 5, 
          arma::uword max_nodes= 1000){
          // double bag_rows = 0.632,
          // double bag_cols = 0.333){ 
  // Bag each tree by randomly selecting observations
  // bag_rows = std::min(1.0,std::abs(bag_rows)); // safety first
  // bag_rows = std::min(1.0,std::abs(bag_cols)); // safety first
  double T = X.n_rows;
  uvec to_keep = select_rnd(T, ceil(0.632*T));
  X = X.rows(to_keep); // this shuffles X and y but it shouldn't matter
  y = y(to_keep);
  double xnc = X.n_cols;
  uword n = ceil(xnc/3); // number of candidates to use at each split
  mat Tree(max_nodes, 9, fill::zeros);
  Tree(0,5) = mean(y); // unconditional mean
  Tree(0,6) = sum(square(y-Tree(0,5))); // unconditional var
  Tree(0,7) = to_keep.n_elem;
  mat leaves; vec par; uvec fobs; mat tmp_tree;
  uvec leaf_idx; 
  // double var_new = 0;
  field<uvec> I(max_nodes);
  I(0) = regspace<uvec>(0,X.n_rows-1); // index in original data (after bagging)
  field<vec> tmp; // temp output from best_splits()
  uword i = 0;
  uword j = 0;
  // double var_old = Tree(0,7); // volatility on previous node
  while(i+2<max_nodes){
    tmp = bestsplit(X.rows(I(j)), y(I(j)), I(j), n);
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
    // find the leaf with the highest variance for the next iteration
    // fobs = field_obs(I, i+3);
    tmp_tree = Tree.rows(0,i);
    leaf_idx = find(tmp_tree.col(8) == 1 && tmp_tree.col(7) > 5); // index of terminal nodes in Tree matrix with enough obs
    // Find the leaf with the maximum variance
    leaves = Tree.rows(leaf_idx); // terminal nodes (i.e. leaves). 
    if(leaf_idx.n_elem==0) break;
    j = leaf_idx(index_max(leaves.col(6))); // index of leaf with max volatility to work on next
  }
  // Rcpp::List rtrn;
  // mat tree = Tree.rows(0,i);
  // arma::field<uvec> Iout = I.rows(0,i);
  // rtrn["Tree"] = tree;
  // rtrn["I"] = Iout;
  // return(rtrn);
  return(Tree.rows(0,i));
}

// Draw 'draws' number of trees
// [[Rcpp::export]]
arma::field<arma::mat> RForest(arma::vec y, // response (no missing obs)
                                  arma::mat X, // predictors (missing obs OK)
                                  arma::uword min_obs = 5,
                                  arma::uword max_nodes = 1000,
                                  arma::uword draws = 1000){
  field<mat> Trees(draws);
  for(uword j = 0; j<draws; j++){
    Trees(j) = RegTree(y, X, min_obs, max_nodes);
  }
  return(Trees);
}

// Fit a single observation using the estimated tree
// [[Rcpp::export]]
double fitvec(arma::vec x,
               arma::mat Tree,
               arma::uword maxit = 1000){
  double j = 0; double it=0;
  while(Tree(j,8) != 1 && it<maxit){
    if(!std::isfinite(x(Tree(j,0)))){
      return(Tree(j,5));
    }else{
      if(x(Tree(j,0))>Tree(j,1)){
        j = Tree(j,3);
      }else{
        j = Tree(j,2);
      }
      it++;
    }
  }
  return(Tree(j,5));
}

// Fit a vector of observations using the estimated tree
// [[Rcpp::export]]
arma::vec fitmat(arma::mat X,
                  arma::mat Tree){
  vec Mu(X.n_cols);
  for(uword j=0; j<X.n_cols; j++){
    Mu(j) = fitvec(X.col(j), Tree);
  }
  return(Mu);
}


// Fit output from RegForest
// [[Rcpp::export]]
arma::vec fitfield(arma::mat X,
                        arma::field<arma::mat> Trees){
  mat Mu(X.n_rows, Trees.n_elem);
  X = trans(X); //transpose for FitMat
  for(uword j=0; j<Trees.n_elem; j++){
    Mu.col(j) = fitmat(X, Trees(j));
  }
  vec mu = mean(Mu,1); // take average response
  return(mu);
}

