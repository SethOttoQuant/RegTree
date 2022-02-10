// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace arma;
using namespace Rcpp;
// #include "utils.h"


// -----  these functions should go into a new file with a header file utils.h... someday -----
arma::uvec selectrnd(arma::uword m, // number of elements
                      arma::uword n){ // number to select
  uvec idx = regspace<uvec>(0,m-1);
  uvec out = shuffle(idx);
  return(out.head(n));
}

// --------------------------------------------------------------------------------------

// [[Rcpp::export]]
arma::vec find_scut(arma::vec x, // predictor
                   arma::vec y){ // response
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  double vnce_na = 0;
  if(not_obs.n_elem>0) vnce_na = sum(square(y(not_obs) - mean(y)));
  if(obs.n_elem<2){
    vec out = {0, 0, 0 , 0, vnce_na, 0}; // no observations; not useful
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
  uword imin = index_min(sum(vnce,0));
  double cut = (x_sorted(imin) + x_sorted(imin+1))/2;
  vec out = {mean(yx(span(0,imin))), mean(yx(span(imin+1, n-1))), vnce(0,imin), vnce(1,imin), vnce_na, cut};
  return(out); // less than cut mean, greater than cut mean, less variance, greater variance, NA variance, and cut. 
  
}

// find the best series in X to identify y
// [[Rcpp::export]]
arma::vec best_ssplit(arma::mat X, // predictors
                     arma::vec y, // response
                     arma::uword n){ // number of candidates to select
  // select candidates for splitting
  uvec candidates = selectrnd(X.n_cols, n); // randomly selected candidates
  mat splits(6,n); // results matrix (mean <, mean >, vol <, vol >, vol NA, cut)
  for(uword j=0; j<n; j++){
    splits.col(j) = find_scut(X.col(candidates(j)), y); // find best split for each series
  }
  vec tot_var = trans(sum(splits.rows(2,4),0)); // rows 2:4 contain variance for high, low, and NA values
  uword min_idx = index_min(tot_var); // which series offers least variance in y
  vec out(7); // results vector 
  out(0) = candidates(min_idx); // index of the best one
  out(span(1,6)) = splits.col(min_idx); // other outputs (see below)
  return(out);
}


// for each line of tree matrix returns: variable, cut values, 0/1 less or greater than
// [[Rcpp::export]]
arma::mat node_conditions(arma::mat Tree,
                          arma::uword j){ // must start at a terminal node
  mat out(3,Tree.n_rows, fill::zeros); 
  uword it = 0; uword k = 0;
  while(j>0){
    k = Tree(j,4); // node from
    out(0,it) = Tree(k,0); // variable
    out(1,it) = Tree(k,1); // cut value
    if(Tree(k,2) == j){
      out(2,it) = 0; // less than
    }else if(Tree(k,3) == j){
      out(2,it) = 1; // greater than
    }else{
      Rcpp::stop("Bad tree");
    }
    j = k;
    it++;
  }
  out = out.cols(0,it-1); // shed unneeded columns
  return(out);
}

// unconstrained optimization problem
// [[Rcpp::export]]
arma::mat Std_Reg_Tree(arma::vec y, // response (no missing obs)
                  arma::mat X, // predictors (missing obs OK)
                  arma::vec depth_range){ // required improvement in variance to continue
  // Bag each tree by randomly selecting observations
  double T = X.n_rows;
  uvec to_keep = selectrnd(T, ceil(0.632*T));
  X = X.rows(to_keep); // this shuffles X and y but it shouldn't matter
  y = y(to_keep);
  double xnc = X.n_cols;
  uword max_nodes= ceil(as_scalar((depth_range(1)-depth_range(0))*randu<vec>(1) +
    depth_range(0))); // randomise depth of model
  uword n = ceil(xnc/3); // number of candidates to use at each split
  mat Tree(max_nodes, 8, fill::zeros);
  Tree(0,5) = mean(y); // unconditional mean
  Tree(0,6) = sum(square(y-Tree(0,4))); // unconditional var
  vec split; mat Xtmp = X; vec ytmp = y; mat leaves;
  mat cndtn; uvec idx; uvec leaf_idx; 
  uword i = 0;
  uword j = 0;
  while(i+2<max_nodes){
    split = best_ssplit(Xtmp, ytmp, n);
    Tree(j,0) = split(0); // which variable
    Tree(j,1) = split(6); // cut value
    Tree(j,2) = i+1; // if <, go to node i+1
    Tree(j,3) = i+2; // if >, go to node i+2
    Tree(j,7) = 0; // no longer terminal node
    Tree(i+1,4) = j; // node coming from
    Tree(i+2,4) = j; // node coming from
    Tree(i+1,5) = split(1); // mean value if <
    Tree(i+2,5) = split(2); // mean value if >
    Tree(i+1,6) = split(3); // vol if <
    Tree(i+2,6) = split(4); // vol if >
    Tree(i+1,7) = 1; // terminal node (leaf)
    Tree(i+2,7) = 1; // terminal node (leaf)
    // find the leaf with the highest variance for the next iteration
    leaf_idx = find(Tree.col(7)); // index of terminal nodes in Tree matrix
    leaves = Tree.rows(leaf_idx); // terminal nodes (i.e. leaves)
    // Rcpp::Rcout << Tree.rows(0,10) << endl;
    j = leaf_idx(index_max(leaves.col(6))); // index of leaf with max volatility
    i += 2; // split result in 2 new rows
    if(sum(split(span(3,4)))==0) break; // only 1 obs above and below
    // collect conditions for node j
    cndtn = node_conditions(Tree, j); // conditions, i.e. x1>c1, x5>c2, etc... Format is: variable, cut, < (0) or > (1)
    Xtmp = X;
    ytmp = y;
    for(uword l=0; l<cndtn.n_cols; l++){
      if(cndtn(2,l)==0){ // less than
        idx = find(Xtmp.col(cndtn(0,l)) <= cndtn(1,l));
      }else{ // greater than
        idx = find(Xtmp.col(cndtn(0,l)) > cndtn(1,l));
      }
      Xtmp = Xtmp.rows(idx); // keep only rows that meet condition
      ytmp = ytmp(idx); // keep only rows that meet condition
    }
  }
  Tree = Tree.rows(0,i);
  return(Tree);  
}

// Fit a single observation using the estimated tree
// [[Rcpp::export]]
double FitSVec(arma::vec x,
              arma::mat Tree,
              arma::uword maxit = 1000){
  double j = 0; double it=0;
  while(Tree(j,7) != 1 && it<maxit){
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
arma::vec FitSMat(arma::mat X,
                 arma::mat Tree){
  vec Mu(X.n_cols);
  for(uword j=0; j<X.n_cols; j++){
    Mu(j) = FitSVec(X.col(j), Tree);
  }
  return(Mu);
}

// Draw 'draws' number of trees
// [[Rcpp::export]]
arma::field<arma::mat> Rnd_Forest(arma::vec y, // response (no missing obs)
                                 arma::mat X, // predictors (missing obs OK)
                                 arma::vec depth_range,
                                 arma::uword draws = 1000){
  field<mat> Trees(draws);
  for(uword j = 0; j<draws; j++){
    Trees(j) = Std_Reg_Tree(y, X, depth_range);
  }
  return(Trees);
}

// Fit output from RegForest
// [[Rcpp::export]]
arma::vec Std_Fit_Field(arma::mat X,
                   arma::field<arma::mat> Trees){
  mat Mu(X.n_rows, Trees.n_elem);
  X = trans(X); //transpose for FitMat
  for(uword j=0; j<Trees.n_elem; j++){
    Mu.col(j) = FitSMat(X, Trees(j));
  }
  vec mu = mean(Mu,1); // take average response
  return(mu);
}


