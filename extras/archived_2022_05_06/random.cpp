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

// [[Rcpp::export]]
arma::uword rnd_idx(double to,
                    double from = 0){
  return ceil(as_scalar((to-from)*randu<vec>(1) + from));
}

// generate a totally random tree
// [[Rcpp::export]]
arma::mat rnd_tree(arma::mat X, // predictors (missing obs OK)
                   arma::vec depth_range){ // required improvement in variance to continue
  uword k = X.n_cols;
  uword max_nodes = rnd_idx(depth_range(1), depth_range(0));// randomize depth of model
  mat Tree(max_nodes, 6, fill::zeros);
  mat Xtmp = X; 
  mat cndtn; uvec idx; uvec leaf_idx; 
  uword i = 0;
  uword j = 0;
  uword count = 0;
  uword c_idx, r_idx;
  while(i+2<max_nodes){
    count = 0;
    c_idx = rnd_idx(k-1);
    r_idx = rnd_idx(Xtmp.n_rows-1);
    while(!std::isfinite(Xtmp(r_idx,c_idx))){ // if bad draw, draw again
      c_idx = rnd_idx(k-1); // random column
      r_idx = rnd_idx(Xtmp.n_rows-1); // random row
      count++;
      if(count>98){
        Tree(j,5) = 2; // do not operate on this node
        goto pick_node;
      }; // break out of nested loop
    }
    Tree(j,0) = c_idx; // which variable
    Tree(j,1) = Xtmp(r_idx,c_idx); // cut value
    Tree(j,2) = i+1; // if <, go to node i+1
    Tree(j,3) = i+2; // if >, go to node i+2
    Tree(j,5) = 0; // no longer terminal node
    Tree(i+1,4) = j; // node coming from
    Tree(i+2,4) = j; // node coming from
    Tree(i+1,5) = 1; // terminal node (leaf)
    Tree(i+2,5) = 1; // terminal node (leaf)
    i += 2; // split result in 2 new rows
    pick_node:
      // find terminal nodes to work on
      leaf_idx = find(Tree.col(5)==1); // index of terminal nodes in Tree matrix
      if(leaf_idx.n_elem==0) break; // all done
      j = leaf_idx(rnd_idx(leaf_idx.n_elem-1)); // index of new leaf to work on
      cndtn = node_conditions(Tree, j); // conditions, i.e. x1>c1, x5>c2, etc... Format is: variable, cut, < (0) or > (1)
      Xtmp = X;
      for(uword l=0; l<cndtn.n_cols; l++){
        if(cndtn(2,l)==0){ // less than
          idx = find(Xtmp.col(cndtn(0,l)) <= cndtn(1,l));
        }else{ // greater than
          idx = find(Xtmp.col(cndtn(0,l)) >= cndtn(1,l));
        }
        Xtmp = Xtmp.rows(idx); // keep only rows that meet condition
        // ytmp = ytmp(idx); // keep only rows that meet condition
      }
      // Tree(j, 5) = mean(ytmp);
      if(Xtmp.n_rows==1){
        Tree(j,5) = 2; // do not work on this node
        goto pick_node;
      }
  }
  Tree = Tree.rows(0,i);
  return(Tree);  
}

// Fit a single observation using the estimated tree
// [[Rcpp::export]]
double fitvec(arma::vec x, // out of sample obs to fit
               arma::mat X, // in sample to select
               arma::vec y, // result
               arma::mat Tree,
               arma::uword maxit = 1000){
  double j = 0; double it=0;
  uword idx; double cut=0; uvec ind;
  uvec I = regspace<uvec>(0, y.n_elem);
  while(Tree(j,5) == 0 && it<maxit){
    if(!std::isfinite(x(Tree(j,0)))){
      return(mean(y));
    }else{
      idx = Tree(j,0);
      cut = Tree(j,1);
      if(x(idx)>cut){
        j = Tree(j,3);
        ind = find(X.col(idx)>cut);
      }else{
        j = Tree(j,2);
        ind = find(X.col(idx)<=cut);
      }
      if(ind.n_elem==0) return(mean(y));
      X = X.rows(ind);
      y = y(ind);
      I = I(ind);
      it++;
    }
  }
  return(mean(y));
}

// Fit a vector of observations using the estimated tree
// [[Rcpp::export]]
arma::vec fitmat(arma::mat X_out, // out of sample values to fit
                 arma::mat X_in, // in sample values for estimates
                 arma::vec y,
                 arma::mat Tree){
  vec Mu(X_out.n_cols);
  for(uword j=0; j<X_out.n_cols; j++){
    Mu(j) = fitvec(X_out.col(j), X_in, y, Tree);
  }
  return(Mu);
}

// Draw 'draws' number of trees
// [[Rcpp::export]]
arma::field<arma::mat> rand_rand(arma::mat X, // predictors (missing obs OK)
                                 arma::vec depth_range,
                                 arma::uword draws = 1000){
  field<mat> Trees(draws);
  for(uword j=0; j<draws; j++){
    Trees(j) = rnd_tree(X, depth_range);
  }
  return(Trees);
}

// Fit output from RegForest
// [[Rcpp::export]]
arma::field<vec> rnd_fit(arma::mat X_out,
                  arma::mat X_in,
                  arma::vec y,
                  arma::field<arma::mat> Trees,
                  arma::uword keep = 10){
  uword k = Trees.n_elem; 
  mat Mu(X_out.n_rows, k);
  X_out = trans(X_out);
  for(uword j=0; j<Trees.n_elem; j++){
    Mu.col(j) = fitmat(X_out, X_in, y, Trees(j));
  }
  mat E = Mu.each_col() - y;
  vec mse = mean(square(E),1);
  uvec idx = sort_index(mse);
  idx = idx(span(0,keep-1));
  mat mu = Mu.cols(idx);
  // vec xx = sum(Mu%Mu,1);
  // mat XY = Mu.each_col() % y;
  // vec xy = sum(XY,1);
  // vec w = xy/xx;
  // uvec idx = find(w>0);
  // w = w(idx);
  // Mu = Mu.cols(idx);
  // vec mu = Mu*w;
  vec mu_out = mean(mu, 1);
  field<vec> out(2);
  out(0) = mu_out;
  out(1) = conv_to<vec>::from(idx);
  out(2) = mse;
  return(out);
}


