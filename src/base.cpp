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
arma::field<vec> quickreg(arma::vec x, // no missing
                          arma::vec y, // no missing
                          double r){ // ridge parameter
  mat X(x.n_elem, 2, fill::ones);
  X.col(1) = x;
  mat I = r*eye(2,2);
  vec par(3);
  mat ab = solve(trans(X)*X + I, trans(X)*y);
  vec e = y - X*ab;
  double sig = sum(square(e));
  field<vec> out(2);
  par(span(0,1)) = ab.col(0);
  par(2) = sig;
  out(0) = par; // parameters: intercept, beta, sigma
  out(1) = e; // residuals
  return(out);
}

// Output is:

// [[Rcpp::export]]
arma::field<vec> find_cut(arma::vec x, // predictor
                          arma::vec y, // response (residuals)
                          arma::uvec ind){ 
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  uword n = obs.n_elem;
  uword nnot = not_obs.n_elem;
  field<vec> out(5);
  double vnce_na = 0;
  if(nnot>0) vnce_na = sum(square(y(not_obs))); // E(y) with no new data 
  if(n<10){
    vec par = {0, 0, 0, 0, 0 , 0, vnce_na, 0}; // not enough observations; not useful
    out(0) = par;
    out(1) = y; // residuals unchanged
    out(2) = y; // residuals unchanged
    out(3) = conv_to<vec>::from(ind); // index of conditions unchanged
    out(4) = conv_to<vec>::from(ind); // index of conditions unchanged
    return(out);
  }
  x = x(obs);
  y = y(obs);
  ind = ind(obs);
  uvec idx = sort_index(x); // ascending; may be duplicate values here (doesn't matter)
  vec x_sorted = x(idx);
  vec yx = y(idx); // sort y by x values
  mat Pleq(3,n-9); // pars, <=. Order is a  b  sig
  mat Pg(3,n-9); // pars, >
  field<vec> Eleq(n-9); // residuals, <=
  field<vec> Eg(n-9);; // residuals, >
  field<vec> tmp; // output of quickreg
  vec vnce(n-9); // make sure indexes correspond with actual x and y data just to keep things easier
  for(uword j=4; j<n-5; j++){ // require at least 5 observations
    tmp = quickreg(x(span(0,j)), y(span(0,j)), nnot); // <=
    Pleq.col(j-4) = tmp(0); // pars for <=
    Eleq(j-4) = tmp(1); // residuals for <= (different length depending on split)
    tmp = quickreg(x(span(j+1,n-1)), y(span(j+1,n-1)), nnot); // >
    Pg.col(j-4) = tmp(0); // pars for >
    Eg(j-4) = tmp(1); // residuals for >
    vnce(j-4) = Pleq(2,j-4) + Pg(2,j-4); // total variance
  }
  uword imin = index_min(vnce);
  double cut = (x_sorted(imin+4) + x_sorted(imin+5))/2;
  // Rcpp::Rcout << imin << endl;// [[Rcpp::export]]
  vec pars = {Pleq(0,imin), Pleq(1,imin), Pg(0,imin), Pg(1,imin), Pleq(2,imin), Pg(2,imin), vnce_na, cut};
  out(0) = pars; // a<=, b<=, a>, b>, sig<=, sig>, sig NA, c
  out(1) = Eleq(imin); // residuals, <=
  out(2) = Eg(imin); // residuals, >
  out(3) = conv_to<vec>::from(ind(idx(span(0,imin+4)))); // index of original dataset <=
  out(4) = conv_to<vec>::from(ind(idx(span(imin+5,n-1)))); // index of original dataset >
  return(out); 
}

// find the best series in X to identify y
// [[Rcpp::export]]
arma::field<vec> best_split(arma::mat X, // predictors
                 arma::vec y, // response
                 arma::uvec ind, // index of observations in original data
                 arma::uword n){ // number of candidates to select
  // select candidates for splitting
  uvec candidates = select_rnd(X.n_cols, n); // randomly selected candidates
  mat splits(8,n); // results matrix (mean <, mean >, vol <, vol >, vol NA, cut)
  field<vec> tmp;
  for(uword j=0; j<n; j++){
    tmp = find_cut(X.col(candidates(j)), y, ind); // mildly inefficient since we're using only one output here
    splits.col(j) = tmp(0); // find best split for each series
  }
  vec tot_var = trans(sum(splits.rows(4,6),0)); // rows 2:4 contain variance for high, low, and NA values
  uword min_idx = index_min(tot_var); // which series offers max reduction in variance of y
  field<vec> out = find_cut(X.col(candidates(min_idx)), y, ind); // re-running find_cut here to keep things clean
  // clunky but effective
  vec par_out(9);
  par_out(0) = min_idx;
  par_out(span(1,8)) = out(0);
  out(0) = par_out;
  return(out); // for find_cut() see function above
}


// // for each line of tree matrix returns: variable, cut values, 0/1 less or greater than
// // [[Rcpp::export]]
// arma::mat node_conditions(arma::mat Tree,
//                           arma::uword j){ // must start at a terminal node
//   mat out(3,Tree.n_rows, fill::zeros);
//   uword it = 0; uword k = 0;
//   while(j>0){
//     k = Tree(j,4); // node from
//     out(0,it) = Tree(k,0); // variable
//     out(1,it) = Tree(k,1); // cut value
//     if(Tree(k,2) == j){
//       out(2,it) = 0; // less than
//     }else if(Tree(k,3) == j){
//       out(2,it) = 1; // greater than
//     }else{
//       Rcpp::stop("Bad tree");
//     }
//     j = k;
//     it++;
//   }
//   out = out.cols(0,it-1); // shed unneeded columns
//   return(out);
// }

// [[Rcpp::export]]
arma::uvec testfun(arma::vec x){
  uvec out = find(x != x || x<2);
  return(out);
}

arma::uvec field_obs(arma::field<vec> E, arma::uword n){
  uvec out(n);
  vec tmp;
  for(uword j=0; j<n; j++){
    tmp = E(j);
    out(j) = tmp.n_elem;
  }
  return(out);
}

// Core function to call
// [[Rcpp::export]]
Rcpp::List RegTree(arma::vec y, // response (no missing obs)
                   arma::mat X, // predictors (missing obs OK)
                   arma::uword max_nodes = 31,
                   double threshold = 0.01){ // required improvement in variance to continue
  // Bag each tree by randomly selecting observations
  double T = X.n_rows;
  // uvec to_keep = select_rnd(T, ceil(0.632*T));
  // X = X.rows(to_keep); // this shuffles X and y but it shouldn't matter
  // y = y(to_keep);
  double xnc = X.n_cols;
  uword n = xnc; // ceil(xnc/3); // number of candidates to use at each split
  mat Tree(max_nodes, 9, fill::zeros);
  Tree(0,5) = mean(y); // unconditional mean
  Tree(0,6) = 0; // beta is zero for the original node
  Tree(0,7) = sum(square(y-Tree(0,5))); // unconditional var
  mat leaves; vec par; uvec fobs; mat tmp_tree;
  uvec leaf_idx; double var_new = 0;
  field<vec> E(max_nodes);
  field<uvec> I(max_nodes);
  E(0) = y - Tree(0,5); // residuals are demeaned values of y
  I(0) = regspace<uvec>(0,X.n_rows-1); // index in original data (after bagging)
  field<vec> tmp; // temp output from best_splits()
  uword i = 0;
  uword j = 0;
  double var_old = Tree(0,7); // volatility on previous node
  while(i+2<max_nodes){
    tmp = best_split(X.rows(I(j)), E(j), I(j), n);
    par = tmp(0);
    Tree(j,0) = par(0); // which variable
    Tree(j,1) = par(8); // cut value
    Tree(j,2) = i+1; // if <, go to node i+1
    Tree(j,3) = i+2; // if >, go to node i+2
    Tree(j,8) = 0; // no longer terminal node
    Tree(i+1,4) = j; // node coming from
    Tree(i+2,4) = j; // node coming from
    Tree(i+1,5) = par(1); // a if <
    Tree(i+2,5) = par(3); // a if >
    Tree(i+1,6) = par(2); // b if <
    Tree(i+2,6) = par(4); // b if >
    Tree(i+1,7) = par(5); // sig if <
    Tree(i+2,7) = par(6); // sig if >
    Tree(i+1,8) = 1; // terminal node (leaf)
    Tree(i+2,8) = 1; // terminal node (leaf)
    E(i+1) = tmp(1); // residuals <=
    E(i+2) = tmp(2); // residuals >
    I(i+1) = conv_to<uvec>::from(tmp(3)); // indexes <=
    I(i+2) = conv_to<uvec>::from(tmp(4)); // indexes >
    // find the leaf with the highest variance for the next iteration
    fobs = field_obs(E, i+3);
    tmp_tree = Tree.rows(0,i+2);
    leaf_idx = find(tmp_tree.col(8) == 1 && fobs > 9); // index of terminal nodes in Tree matrix with enough obs
    // Find the leaf with the maximum variance
    leaves = Tree.rows(leaf_idx); // terminal nodes (i.e. leaves). 
    // Rcpp::Rcout << Tree.rows(0,10) << endl;
    var_new = sum(par(span(4,6))); // volatility at new node
    // Rcpp::Rcout << "old = " << var_old << "new = " << var_new << "next old = " << Tree(j,6) << endl;
    Rcpp::Rcout << (Tree(j,7) - var_new)/Tree(j,7) << endl;
    // if((var_old - var_new)/var_old < threshold) break; // if volatility doesn't improve, break loop
    i += 2; // split result in 2 new rows
    if(leaf_idx.n_elem==0) break;
    j = leaf_idx(index_max(leaves.col(7))); // index of leaf with max volatility to work on next
  }
  Rcpp::Rcout << "A" << endl;
  Rcpp::List rtrn;
  mat tree = Tree.rows(0,i);
  Rcpp::Rcout << "C" << endl;
  arma::field<vec> Eout = E.rows(0,i);
  Rcpp::Rcout << "D" << endl;
  arma::field<uvec> Iout = I.rows(0,i);
  rtrn["Tree"] = tree;
  rtrn["E"] = Eout;
  rtrn["I"] = Iout;
  return(rtrn);
}

// Fit a single observation using the estimated tree
// [[Rcpp::export]]
double FitVec(arma::vec x,
              arma::mat Tree,
              arma::uword maxit = 1000){
  double j = 0; double it=0; double y=0;
  while(Tree(j,8) != 1 && it<maxit){
    if(!std::isfinite(x(Tree(j,0)))){
      return(y);
    }else{
      y += Tree(j,5) + Tree(j,6)*x(Tree(j,0)); // need to use the split here too!!!!!
      if(x(Tree(j,0))>Tree(j,1)){
        j = Tree(j,3);
      }else{
        j = Tree(j,2);
      }
      it++;
    }
  }
  y += Tree(j,5) + Tree(j,6)*x(Tree(j,0)); // terminal node (leaf)
  return(y);
}

// Fit a vector of observations using the estimated tree
// [[Rcpp::export]]
arma::vec FitMat(arma::mat X,
                 arma::mat Tree){
  vec Mu(X.n_cols);
  for(uword j=0; j<X.n_cols; j++){
    Mu(j) = FitVec(X.col(j), Tree);
  }
  return(Mu);
}

// // Draw 'draws' number of trees
// // [[Rcpp::export]]
// arma::field<arma::mat> RegForest(arma::vec y, // response (no missing obs)
//                            arma::mat X, // predictors (missing obs OK)
//                            arma::uword max_nodes = 31,
//                            double threshold = 0.01,
//                            arma::uword draws = 500){
//   field<mat> Trees(draws);
//   for(uword j = 0; j<draws; j++){
//     Trees(j) = RegTree(y, X, max_nodes, threshold);
//   }
//   return(Trees);
// }

// Fit output from RegForest
// [[Rcpp::export]]
arma::vec FitField(arma::mat X,
                   arma::field<arma::mat> Trees){
  mat Mu(X.n_rows, Trees.n_elem);
  X = trans(X); //transpose for FitMat
  for(uword j=0; j<Trees.n_elem; j++){
    Mu.col(j) = FitMat(X, Trees(j));
  }
  vec mu = mean(Mu,1); // take average response
  return(mu);
}




