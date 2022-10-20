// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace arma;
using namespace Rcpp;
#include "utils.h"

// [[Rcpp::export]]
Rcpp::List tst_sort(arma::vec x){
  uvec idx = sort_index(x);
  x = x(idx);
  vec x_unsort = x(idx);
  Rcpp::List out;
  out[0] = x;
  out[1] = x_unsort;
  out[2] = idx;
  return out;
}

// [[Rcpp::export]]
arma::mat fast_cut(arma::vec x, // predictor
                   arma::vec y, // response
                   double prior_shrink, // prior shrinkage par
                   double my=0){ // prior for y
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  double m_na=0, vnce_na = 0, m_leq=0, m_g=0;
  mat tmp_out(1,2);
  tmp_out.fill(datum::nan);
  if(not_obs.n_elem>0){
    m_na = (sum(y(not_obs)) + prior_shrink*my)/(not_obs.n_elem + prior_shrink);
    vnce_na = sum(square(y(not_obs) - m_na))/(not_obs.n_elem + prior_shrink);
  }
  if(obs.n_elem<(4)){  // hard coded min node size of 2
    // Rcpp::Rcout << "bad split" << endl;
    return(tmp_out);  // never split this guy
  }
  x = x(obs);
  y = y(obs);
  uvec idx = sort_index(x); // ascending; may be duplicate values here (doesn't matter)
  vec x_sorted = x(idx);
  vec yx = y(idx); // sort y by x values
  uword n = yx.n_elem; // number of obs
  uword idx_j;
  mat out(2, n-3, fill::zeros);
  for(uword j=0; j<=(n-4); j++){
    idx_j = j + 1;
    // Rcpp::Rcout << j << endl;
    m_leq = (sum(yx(span(0,idx_j))) + prior_shrink*my)/(idx_j + prior_shrink);
    m_g = (sum(yx(span(idx_j+1,n-1))) + prior_shrink*my)/(n-idx_j + prior_shrink);
    out(0,j) = 1/std::pow(sum(square(yx(span(0,idx_j)) - m_leq))/(idx_j + prior_shrink) +
               sum(square(yx(span(idx_j+1,n-1)) - m_g))/(n-idx_j + prior_shrink) + 
               vnce_na, 10); // control prob of being selected. 
    out(1,j) = x_sorted(idx_j); // cut value
  }
  // this bit is just to select the split with probability weight given by out.row(0)
  vec cdf = cumsum(out.row(0));
  double winner = randu()*cdf(n-3); // rand u var over full range of cdf
  uword j=0;
  while(winner<cdf(j)){
    j++;
  }
  return(out.col(j)); // mat with row 0 = variance, row 1 = cut value
}
// [[Rcpp::export]]
arma::field<arma::vec> cut(arma::vec x, // predictor
                       arma::vec y, // response
                       arma::uvec ind, // index
                       double imin,
                       arma::uword min_obs){  // where to cut
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  //uword n = obs.n_elem;
  uword nnot = not_obs.n_elem;
  field<vec> out(3);
  double vnce_na = 0;
  if(nnot>0) vnce_na = sum(square(y(not_obs) - mean(y))); // E(y) with no new data
  uword n = obs.n_elem; // number of obs
  if(n<2*min_obs){ 
    // Rcpp::Rcout << "bad cut" << endl;
    vec par = {0, 0, 0 , 0, sum(square(y - mean(y))), 0, 0, 0}; // not enough observations; not useful
    out(0) = par;
    out(1) = conv_to<vec>::from(ind); // index of conditions unchanged
    out(2) = conv_to<vec>::from(ind); // index of conditions unchanged
    return(out);
  }
  x = x(obs);
  y = y(obs);
  ind = ind(obs);  // original index
  uvec idx = sort_index(x); // ascending; may be duplicate values here (doesn't matter)
  x = x(idx); // sorted by x values
  y = y(idx); // sort y by x values
  double cut = (x(imin) + x(imin+1))/2;  // cut value
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
arma::field<arma::vec> selectsplit(arma::mat X, // predictors
                       arma::vec y, // response
                       arma::uvec ind, // index of observations in original data
                       double prior_shrink, // shrink estimate towards unconditional mean
                       double my=0){ // unconditional mean
  arma::field<mat> candidates(X.n_cols);
  for(uword j=0; j<y.n_elem; j++){
    // Rcpp::Rcout << candidates(j) << endl;
    candidates(j) = fast_cut(X.col(j), y, prior_shrink, my); 
  }
  vec tot_var = trans(sum(splits.rows(0,2),0)); // rows 0:2 contain variance for high, low, and NA values
  uword min_idx = index_min(tot_var); // which series offers max reduction in variance of y
  field<vec> out = cut(X.col(candidates(min_idx)), y, ind, splits(3, min_idx), min_obs); // re-running cut here to keep things clean
  // clunky but effective
  vec par_out(9);
  par_out(0) = candidates(min_idx);
  par_out(span(1,8)) = out(0);
  out(0) = par_out;
  return(out); // for find_cut() see function above
}
// [[Rcpp::export]]
arma::mat regtree(arma::vec y, // response (no missing obs)
          arma::mat X, // predictors (missing obs OK), obs in rows series in cols
          arma::uword max_nodes= 10,
          double prior_shrink = 0){
  double xnc = X.n_cols;
  mat Tree(max_nodes, 10, fill::zeros);
  double my = mean(y);
  Tree(0,5) = my; // unconditional mean
  double vy = mean(square(y-my)); // unconditional var
  Tree(0,6) = vy;
  Tree(0,7) = y.n_elem;
  mat leaves; vec par; uvec fobs; mat tmp_tree;
  uvec leaf_idx; 
  field<uvec> I(max_nodes);
  I(0) = regspace<uvec>(0,X.n_rows-1); // index in original data 
  field<vec> tmp; // temp output from best_splits()
  uword i = 0;  // next nodes
  uword j = 0;  // this node
  while(i+2<max_nodes){
    tmp = bestsplit(X.rows(I(j)), y(I(j)), I(j), prior_shrink);
    par = tmp(0);
    if(any(par(span(1,4)))){ //  enough observations to split
      Tree(j,0) = par(0); // which variable
      Tree(j,1) = par(6); // cut value
      Tree(j,2) = i+1; // if <, go to node i+1
      Tree(j,3) = i+2; // if >, go to node i+2
      Tree(j,9) = 0; // no longer terminal node
      Tree(i+1,4) = j; // node coming from
      Tree(i+2,4) = j; // node coming from
      Tree(i+1,5) = par(1); // mu if <
      Tree(i+2,5) = par(2); // mu if >
      Tree(i+1,6) = par(3); // sig if <
      Tree(i+2,6) = par(4); // sig if >
      Tree(i+1,7) = par(7); // nobs if <
      Tree(i+2,7) = par(8); // nobs if >
      Tree(i+1,8) = 1 - welch_t(par(1), par(2), par(3), par(4), par(7), par(8)); // weight on this node (1 - for previous node)
      Tree(i+2,8) =  Tree(i+1,8);  // weight on this node (1 - for previous node)
      Tree(i+1,9) = 1; // terminal node (leaf)
      Tree(i+2,9) = 1; // terminal node (leaf)
      I(i+1) = conv_to<uvec>::from(tmp(1)); // indexes <=
      I(i+2) = conv_to<uvec>::from(tmp(2)); // indexes >
      i += 2; // split result in 2 new rows
    }else{
      Tree(j,9) = 2; // terminal node with no more possible splits
    }
    tmp_tree = Tree.rows(0,i);
    leaf_idx = find(tmp_tree.col(9) == 1 && tmp_tree.col(7) > max_obs); // index of terminal nodes in Tree matrix with enough obs
    leaves = Tree.rows(leaf_idx); // terminal nodes (i.e. leaves). 
    if(leaf_idx.n_elem==0) break;
    j = leaf_idx(index_max(leaves.col(6))); // index of leaf with max volatility to work on next
  }
  return(Tree.rows(0,i));
}

// Fit a single observation using the estimated tree
// [[Rcpp::export]]
arma::field<arma::vec> fitvec(arma::vec x,
                       arma::mat Tree,
                       arma::uword maxit = 1000){
  uword j=0; uword j_old=0; uword it=0;
  double y = Tree(0,5); double y_old=Tree(0,5);
  field<vec> out(2); // out(0) is value of y, out(1) is feature contributions
  vec fc(x.n_elem, fill::zeros); // vector of feature contributions
  while(Tree(j,9) == 0 && it<maxit){
    if(!std::isfinite(x(Tree(j,0)))){
      break;
    }else{
      if(x(Tree(j,0))>Tree(j,1)){
        j = Tree(j,3);
      }else{
        j = Tree(j,2);
      }
      y = Tree(j,5);
      fc(Tree(j_old,0)) += y - y_old;
      // Rcpp::Rcout << "j = " << j << "; y = " << y << "; y_old = " << y_old << endl;
      y_old = y;
      j_old = j;
      it++;
    }
  }
  out[0] = y;
  out[1] = fc;
  return(out);
}


// Fit a single observation using the estimated tree weighting by var
// [[Rcpp::export]]
arma::field<arma::vec> weightvec(arma::vec x,
                               arma::mat Tree,
                               arma::uword maxit = 1000){
  uword j=0; uword j_old=0; uword it=1;
  double y = Tree(0,5);  // unconditional mean
  double y_old = y;  // used for feature contributions
  field<vec> out(2); // out(0) is value of y, out(1) is feature contributions
  vec fc(x.n_elem, fill::zeros); // vector of feature contributions
  // ---------------- testing -----------------------
  // vec w_out(Tree.n_rows, fill::zeros);
  // vec m_out(Tree.n_rows, fill::zeros);
  // w_out(0) = 1;
  // m_out(0) = Tree(0,5);
  // ------------------------------------------------
  while(Tree(j,9) == 0 && it<maxit){
    // Rcpp::Rcout << it << endl;
    if(!std::isfinite(x(Tree(j,0)))){
      break;
    }else{
      if(x(Tree(j,0))>Tree(j,1)){ // if geq
        j = Tree(j,3); // new node
      }else{ // else if leq
        j = Tree(j,2);
      }
      y = Tree(j,5)*Tree(j,8) + (1-Tree(j,8))*y; // mean times weight
      fc(Tree(j_old,0)) += y - y_old; // fc times weight
      // w_out(it) = Tree(j,8);
      // m_out(it) = Tree(j,5);
      y_old = y;
      j_old = j;
      it++;
    }
  }
  out[0] = y;
  out[1] = fc;
  // out[2] = w_out;
  // out[3] = m_out;
  return(out);
}

// Fit a single observation using the estimated tree weighting by var
// [[Rcpp::export]]
arma::field<arma::vec> poolvec(arma::vec x,
                              arma::mat Tree,
                              arma::uword maxit = 1000){
  uword j=0; uword j_old=0; uword it=1;
  double w = (Tree(0,7)-1)/Tree(0,6); // weight for this node
  double y = Tree(0,5)*w;  // unconditional mean times weight
  double W = w; // this will be sum of all weights
  double y_old = y/W;  // used for feature contributions
  field<vec> out(4); // out(0) is value of y, out(1) is feature contributions
  vec fc(x.n_elem, fill::zeros); // vector of feature contributions
  // ---------------- testing -----------------------
  vec w_out(Tree.n_rows, fill::zeros);
  vec m_out(Tree.n_rows, fill::zeros);
  w_out(0) = w;
  m_out(0) = Tree(0,5);
  // ------------------------------------------------
  while(Tree(j,9) == 0 && it<maxit){
    if(!std::isfinite(x(Tree(j,0)))){
      break;
    }else{
      if(x(Tree(j,0))>Tree(j,1)){ // if geq
        j = Tree(j,3); // new node
      }else{ // else if leq
        j = Tree(j,2);
      }
      w = (Tree(j,7)-1)/Tree(j,6); // weight for this (new) node
      y += Tree(j,5)*w; // mean times weight
      W += w; // sum of all weights
      fc(Tree(j_old,0)) += w*(y/W - y_old); // fc times weight
      w_out(it) = w;
      m_out(it) = Tree(j,5);
      y_old = y/W;
      j_old = j;
      it++;
    }
  }
  out[0] = y/W;
  out[1] = fc/W;
  out[2] = w_out;
  out[3] = m_out;
  return(out);
}

// Fit a vector of observations using the estimated tree
// [[Rcpp::export]]
arma::field<arma::mat> fitmat(arma::mat X,
                              arma::mat Tree,
                              bool weight=false){ // bool T/F
  vec Mu(X.n_cols);
  field<vec> tmp;
  mat FC(X.n_rows, X.n_cols, fill::zeros);
  field<mat> out(2);
  for(uword j=0; j<X.n_cols; j++){
    if(weight){
      tmp = weightvec(X.col(j), Tree);
    }else{
      tmp = fitvec(X.col(j), Tree);
    }
    Mu(j) = as_scalar(tmp(0));
    FC.col(j) = tmp(1);
  }
  out(0) = Mu;
  out(1) = FC;
  return(out);
}
// Fit output from RegForest
// [[Rcpp::export]]
arma::field<arma::mat> fitfield(arma::mat X,
                                arma::field<arma::mat> Trees,
                                bool weight = false){ // bool T/F
  uword k = Trees.n_elem;
  vec Mu(X.n_rows, fill::zeros); // mean (ie prediction)
  mat FC(X.n_cols, X.n_rows, fill::zeros); // feature contribution
  field<mat> tmp;
  field<mat> out(2);
  X = trans(X); //transpose for FitMat
  for(uword j=0; j<Trees.n_elem; j++){
    tmp = fitmat(X, Trees(j), weight);
    Mu += tmp(0);
    FC += tmp(1);
  }
  out(0) = Mu/k; // take average response (div by num trees)
  out(1) = trans(FC)/k;
  return(out);
}

// Draw 'draws' number of trees
// [[Rcpp::export]]
Rcpp::List rforest(arma::vec y, // response (no missing obs)
                   arma::mat X, // predictors (missing obs OK)
                   arma::uword max_obs = 20,
                   arma::uword min_obs = 5,
                   arma::uword max_nodes = 1000,
                   arma::uword draws = 1000,
                   bool weight = false){ // T/F use pooled forecasts
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
    Tree = regtree(y, X, to_keep(0), max_obs, min_obs, max_nodes);
    tmp = fitmat(trans(X.rows(to_keep(1))), Tree, weight);  // oob results
    oob.fill(datum::nan);
    oob(to_keep(1)) = tmp(0); // fitted vals
    mse(j) = mean(square(y(to_keep(1)) - tmp(0))); // oob mse
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

