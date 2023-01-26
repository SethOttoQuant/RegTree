// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace arma;
using namespace Rcpp;
#include "utils.h"

// [[Rcpp::export]]
arma::vec fast_cut(arma::vec x, // predictor
                   arma::vec y, // response
                   arma::uword min_obs=1){ // min obs per node
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  double vnce_na = 0;
  if(not_obs.n_elem>0) vnce_na = sum(square(y(not_obs) - mean(y)));
  if(obs.n_elem<(2*min_obs)){
    // Rcpp::Rcout << "bad split" << endl;
    vec out = {0, 0, vnce_na, 0}; // not enough observations; not useful
    return(out);
  }
  x = x(obs);
  y = y(obs);
  uvec idx = sort_index(x); // ascending; may be duplicate values here (doesn't matter)
  vec x_sorted = x(idx);
  vec yx = y(idx); // sort y by x values
  uword n = yx.n_elem; // number of obs
  mat vnce(2,n-2*min_obs+1); // store variance resulting from split
  uword idx_j;
  for(uword j=0; j<=(n-2*min_obs); j++){
    idx_j = j + min_obs - 1;
    // Rcpp::Rcout << j << endl;
    vnce(0,j) = sum(square(yx(span(0,idx_j)) - mean(yx(span(0,idx_j))))); // <= variance
    vnce(1,j) = sum(square(yx(span(idx_j+1,n-1)) - mean(yx(span(idx_j+1,n-1))))); // > variance
  }
  double imin = index_min(sum(vnce,0));
  vec out = {vnce(0,imin), vnce(1,imin), vnce_na, imin + min_obs - 1};
  return(out); // less than cut mean, greater than cut mean, less variance, greater variance, NA variance, and cut. 
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
arma::field<arma::vec> bestsplit(arma::mat X, // predictors
                       arma::vec y, // response
                       arma::uvec ind, // index of observations in original data
                       arma::uword n,  // number of candidates to select
                       arma::uword min_obs){ // min number obs per node
  // select candidates for splitting
  field<uvec> cnd = select_rnd(X.n_cols, n);
  uvec candidates = cnd(0); // randomly selected candidates
  mat splits(4,n); 
  vec tmp;
  for(uword j=0; j<n; j++){
    // Rcpp::Rcout << candidates(j) << endl;
    tmp = fast_cut(X.col(candidates(j)), y, min_obs); 
    splits.col(j) = tmp; // pars
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
          arma::mat X, // predictors (missing obs OK)
          arma::uvec to_keep,  // rows to keep for this tree
          arma::uword max_obs=20,  // max obs in terminal node
          arma::uword min_obs = 5,  // min obs in terminal node
          arma::uword max_nodes= 1000){
  X = X.rows(to_keep); // this shuffles X and y but it shouldn't matter
  y = y(to_keep);
  double xnc = X.n_cols;
  uword n = ceil(xnc/3); // number of candidates to use at each split
  mat Tree(max_nodes, 10, fill::zeros);
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
    tmp = bestsplit(X.rows(I(j)), y(I(j)), I(j), n, min_obs);
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
      //Tree(i+1,8) = 1 - welch_t(par(1), par(2), par(3), par(4), par(7), par(8)); // weight on this node (1 - for previous node)
      //Tree(i+2,8) =  Tree(i+1,8);  // weight on this node (1 - for previous node)
      Tree(i+1,8) = 1; // terminal node (leaf)
      Tree(i+2,8) = 1; // terminal node (leaf)
      I(i+1) = conv_to<uvec>::from(tmp(1)); // indexes <=
      I(i+2) = conv_to<uvec>::from(tmp(2)); // indexes >
      i += 2; // split result in 2 new rows
    }else{
      Tree(j,8) = 2; // terminal node with no more possible splits
    }
    tmp_tree = Tree.rows(0,i);
    leaf_idx = find(tmp_tree.col(8) == 1 && tmp_tree.col(7) > max_obs); // index of terminal nodes in Tree matrix with enough obs
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
  while(Tree(j,8) == 0 && it<maxit){
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


// // Fit a single observation using the estimated tree weighting by var
// // [[Rcpp::export]]
// arma::field<arma::vec> weightvec(arma::vec x,
//                                arma::mat Tree,
//                                arma::uword maxit = 1000){
//   uword j=0; uword j_old=0; uword it=1;
//   double y = Tree(0,5);  // unconditional mean
//   double y_old = y;  // used for feature contributions
//   field<vec> out(2); // out(0) is value of y, out(1) is feature contributions
//   vec fc(x.n_elem, fill::zeros); // vector of feature contributions
//   // ---------------- testing -----------------------
//   // vec w_out(Tree.n_rows, fill::zeros);
//   // vec m_out(Tree.n_rows, fill::zeros);
//   // w_out(0) = 1;
//   // m_out(0) = Tree(0,5);
//   // ------------------------------------------------
//   while(Tree(j,9) == 0 && it<maxit){
//     // Rcpp::Rcout << it << endl;
//     if(!std::isfinite(x(Tree(j,0)))){
//       break;
//     }else{
//       if(x(Tree(j,0))>Tree(j,1)){ // if geq
//         j = Tree(j,3); // new node
//       }else{ // else if leq
//         j = Tree(j,2);
//       }
//       y = Tree(j,5)*Tree(j,8) + (1-Tree(j,8))*y; // mean times weight
//       fc(Tree(j_old,0)) += y - y_old; // fc times weight
//       // w_out(it) = Tree(j,8);
//       // m_out(it) = Tree(j,5);
//       y_old = y;
//       j_old = j;
//       it++;
//     }
//   }
//   out[0] = y;
//   out[1] = fc;
//   // out[2] = w_out;
//   // out[3] = m_out;
//   return(out);
// }

// // Fit a single observation using the estimated tree weighting by var
// // [[Rcpp::export]]
// arma::field<arma::vec> poolvec(arma::vec x,
//                               arma::mat Tree,
//                               arma::uword maxit = 1000){
//   uword j=0; uword j_old=0; uword it=1;
//   double w = (Tree(0,7)-1)/Tree(0,6); // weight for this node
//   double y = Tree(0,5)*w;  // unconditional mean times weight
//   double W = w; // this will be sum of all weights
//   double y_old = y/W;  // used for feature contributions
//   field<vec> out(4); // out(0) is value of y, out(1) is feature contributions
//   vec fc(x.n_elem, fill::zeros); // vector of feature contributions
//   // ---------------- testing -----------------------
//   vec w_out(Tree.n_rows, fill::zeros);
//   vec m_out(Tree.n_rows, fill::zeros);
//   w_out(0) = w;
//   m_out(0) = Tree(0,5);
//   // ------------------------------------------------
//   while(Tree(j,9) == 0 && it<maxit){
//     if(!std::isfinite(x(Tree(j,0)))){
//       break;
//     }else{
//       if(x(Tree(j,0))>Tree(j,1)){ // if geq
//         j = Tree(j,3); // new node
//       }else{ // else if leq
//         j = Tree(j,2);
//       }
//       w = (Tree(j,7)-1)/Tree(j,6); // weight for this (new) node
//       y += Tree(j,5)*w; // mean times weight
//       W += w; // sum of all weights
//       fc(Tree(j_old,0)) += w*(y/W - y_old); // fc times weight
//       w_out(it) = w;
//       m_out(it) = Tree(j,5);
//       y_old = y/W;
//       j_old = j;
//       it++;
//     }
//   }
//   out[0] = y/W;
//   out[1] = fc/W;
//   out[2] = w_out;
//   out[3] = m_out;
//   return(out);
// }

// Fit a vector of observations using the estimated tree
// [[Rcpp::export]]
arma::field<arma::mat> fitmat(arma::mat X,
                              arma::mat Tree){ // bool weight=false
  vec Mu(X.n_cols);
  field<vec> tmp;
  mat FC(X.n_rows, X.n_cols, fill::zeros);
  field<mat> out(2);
  for(uword j=0; j<X.n_cols; j++){
    // if(weight){
    //   tmp = weightvec(X.col(j), Tree);
    // }else{}
    tmp = fitvec(X.col(j), Tree);
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
                                arma::field<arma::mat> Trees){ // bool T/F
  uword k = Trees.n_elem;
  vec Mu(X.n_rows, fill::zeros); // mean (ie prediction)
  mat FC(X.n_cols, X.n_rows, fill::zeros); // feature contribution
  field<mat> tmp;
  field<mat> out(2);
  X = trans(X); //transpose for FitMat
  for(uword j=0; j<Trees.n_elem; j++){
    tmp = fitmat(X, Trees(j));
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
    Tree = regtree(y, X, to_keep(0), max_obs, min_obs, max_nodes);
    tmp = fitmat(trans(X.rows(to_keep(1))), Tree);  // oob results
    oob.fill(datum::nan);
    oob(to_keep(1)) = tmp(0); // fitted vals
    mse(j) = mean(square(y(to_keep(1)) - tmp(0)))/(to_keep(1).n_elem); // oob mse
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

