// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace arma;
using namespace Rcpp;
#include "utils.h"
#include <random>

// [[Rcpp::export]]
double rand_geom(const double n){
  std::random_device gen;
  std::geometric_distribution<> d(n); 
  double draw = d(gen);
  return(draw);
}

// [[Rcpp::export]]
arma::vec fast_cut_alt(arma::vec x, // predictor
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
  return(out); // less variance, greater variance, NA variance, and cut. 
}
// [[Rcpp::export]]
arma::field<arma::vec> cut_alt(arma::vec x, // predictor
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
  vec pars = {mean(y(span(0,imin))), mean(y(span(imin+1, n-1))), 
              sum(square(y(span(0,imin)) - mean(y(span(0,imin))))), 
              sum(square(y(span(imin+1,n-1)) - mean(y(span(imin+1,n-1))))), 
              vnce_na, cut, imin+1, n-imin-1};
  out(0) = pars; // mean <=, mean >, sig <=, sig >, sig NA, cut, obs <=, obs > 
  out(1) = conv_to<vec>::from(ind(idx(span(0,imin)))); // index of original dataset <=
  out(2) = conv_to<vec>::from(ind(idx(span(imin+1,n-1)))); // index of original dataset >
  return(out); 
}
// find the best series in X to identify y
// [[Rcpp::export]]
arma::field<arma::vec> bestsplit_alt(arma::mat X, // predictors
                       arma::vec y, // response
                       arma::uvec ind, // index of observations in original data
                       double n){ // geometric dist par
  // select candidates for splitting
  // field<uvec> cnd = select_rnd(X.n_cols, n);
  uword k = X.n_cols;
  // uvec candidates = cnd(0); // randomly selected candidates
  mat splits(4,k); 
  vec tmp;
  for(uword j=0; j<k; j++){
    tmp = fast_cut_alt(X.col(j), y); 
    splits.col(j) = tmp; // pars
  }
  vec tot_var = trans(sum(splits.rows(0,2),0)); // rows 2:4 contain variance for high, low, and NA values
  uword c = rand_geom(n);
  if(c>=k){
    c=k-1;
  }
  uvec sort_idx = sort_index(tot_var); // expensive function... rethink
  uword min_idx = sort_idx(c); // selected splitting variable
  field<vec> out = cut_alt(X.col(min_idx), y, ind, splits(3, min_idx)); // re-running cut here to keep things clean
  // clunky but effective
  vec par_out(9);
  par_out(0) = min_idx; // index
  par_out(span(1,8)) = out(0);
  out(0) = par_out; // adding index to out(0) here
  return(out); 
}
// [[Rcpp::export]]
arma::mat regtree_alt(arma::vec y, // response (no missing obs)
          arma::mat X, // predictors (missing obs OK)
          arma::uvec to_keep,
          arma::uword min_obs = 5, 
          arma::uword max_nodes= 1000,
          double n = .5){
  X = X.rows(to_keep); // this shuffles X and y but it shouldn't matter
  y = y(to_keep);
  // double xnc = X.n_cols;
  // uword n = ceil(xnc/3); // number of candidates to use at each split
  mat Tree(max_nodes, 9, fill::zeros);
  Tree(0,5) = mean(y); // unconditional mean
  Tree(0,6) = sum(square(y-Tree(0,5))); // unconditional var
  Tree(0,7) = to_keep.n_elem;
  mat leaves; vec par; uvec fobs; mat tmp_tree;
  uvec leaf_idx; 
  field<uvec> I(max_nodes);
  I(0) = regspace<uvec>(0,X.n_rows-1); // index in original data (after bagging)
  field<vec> tmp; // temp output from best_splits()
  uword i = 0;
  uword j = 0;
  while(i+2<max_nodes){
    tmp = bestsplit_alt(X.rows(I(j)), y(I(j)), I(j), n);
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

// Fit a single observation using the estimated tree
// [[Rcpp::export]]
arma::field<arma::vec> fitvec_alt(arma::vec x,
                       arma::mat Tree,
                       arma::uword maxit = 1000){
  uword j=0; uword j_old=0; uword it=0; 
  double y = Tree(0,5); double y_old=Tree(0,5);
  field<vec> out(2); // out(0) is value of y, out(1) is feature contributions
  vec fc(x.n_elem, fill::zeros); // vector of feature contributions
  while(Tree(j,8) != 1 && it<maxit){
    if(!std::isfinite(x(Tree(j,0)))){
      break;
    }else{
      if(x(Tree(j,0))>Tree(j,1)){
        j = Tree(j,3);
        y = Tree(j,5);
        fc(Tree(j_old,0)) += y - y_old;
      }else{
        j = Tree(j,2);
        y = Tree(j,5);
        fc(Tree(j_old,0)) += y - y_old;
      }
      y_old = y;
      j_old = j;
      it++;
    }
  }
  out[0] = y;
  out[1] = fc;
  return(out);
}
// Fit a vector of observations using the estimated tree
// [[Rcpp::export]]
arma::field<arma::mat> fitmat_alt(arma::mat X,
                  arma::mat Tree){
  vec Mu(X.n_cols);
  field<vec> tmp;
  mat FC(X.n_rows, X.n_cols, fill::zeros);
  field<mat> out(2);
  for(uword j=0; j<X.n_cols; j++){
    tmp = fitvec_alt(X.col(j), Tree);
    Mu(j) = as_scalar(tmp(0));
    FC.col(j) = tmp(1);
  }
  out(0) = Mu;
  out(1) = FC;
  return(out);
}
// Fit output from RegForest
// [[Rcpp::export]]
arma::field<arma::mat> fitfield_alt(arma::mat X,
                  arma::field<arma::mat> Trees,
                  arma::vec weight){ // length must agree with slices of Tree, must have mean 1
  uword k = Trees.n_elem;
  vec Mu(X.n_rows, fill::zeros); // mean (ie prediction)
  mat FC(X.n_cols, X.n_rows, fill::zeros); // feature contribution
  field<mat> tmp;
  field<mat> out(2);
  X = trans(X); //transpose for FitMat
  for(uword j=0; j<Trees.n_elem; j++){
    tmp = fitmat_alt(X, Trees(j));
    Mu += tmp(0)*weight(j);
    FC += tmp(1)*weight(j);
  }
  out(0) = Mu/k; // take average response (div by num trees)
  out(1) = trans(FC)/k;
  return(out);
}

// Draw 'draws' number of trees
// [[Rcpp::export]]
Rcpp::List rforest_alt(arma::vec y, // response (no missing obs)
                   arma::mat X, // predictors (missing obs OK)
                   arma::uword min_obs = 5, // min node size
                   arma::uword max_nodes = 1000, // max number of nodes
                   arma::uword draws = 1000,
                   double geom_par = 0.5){ // number of trees
  field<mat> Trees(draws);
  double T = X.n_rows;  // total number of obs
  vec oob(T); mat Tree;
  field<uvec> to_keep;
  mat OOB(T, draws);
  mat fc(T,X.n_cols);  // out of bag feature contributions
  cube FC(T, X.n_cols, draws);  // store OOB feature contribution
  vec mse(draws); // store OOB MSE for each model
  field<mat> tmp;
  for(uword j = 0; j<draws; j++){  // new tree (model) for each iteration
    to_keep = select_rnd(T, ceil(0.632*T));  // idx for in and out of bag
    Tree = regtree_alt(y, X, to_keep(0), min_obs, max_nodes, geom_par); // to_keep(0) is in bag
    tmp = fitmat_alt(trans(X.rows(to_keep(1))), Tree);  // to_keep(1) is OOB
    mse(j) = mean(square(y(to_keep(1)) - tmp(0)));
    oob.fill(datum::nan);
    oob(to_keep(1)) = tmp(0);
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

