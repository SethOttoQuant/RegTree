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

//[[Rcpp::export]]
arma::vec join_vec(arma::vec a,
                   arma::vec b){
  uword na = a.n_elem;
  uword nb = b.n_elem;
  vec out(na+nb);
  if(na>0){
    out(span(0, na-1)) = a;
  }
  if(nb>0){
    out(span(na, na+nb-1)) = b;
  }
  return out;
}

//[[Rcpp::export]]
arma::uvec join_uvec(arma::uvec a,
                   arma::uvec b){
  uword na = a.n_elem;
  uword nb = b.n_elem;
  uvec out(na+nb);
  if(na>0){
    out(span(0, na-1)) = a;
  }
  if(nb>0){
    out(span(na, na+nb-1)) = b;
  }
  return out;
}

// Find the resulting variance and index for the best cut point for a variable x
// This is the fast/lightweight version to find the best cut, without the 
// other needed data to proceed to next split
// [[Rcpp::export]]
arma::vec find_split(arma::vec x, // rhs data (one series)
                        arma::vec y, // lhs data (one series)
                        arma::vec w){ // weights
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  uword n = obs.n_elem;
  uword nnot = not_obs.n_elem;
  vec out(2); // output is variance (0) and cut index (1)
  double vnce_na = 0;
  if(nnot>0) vnce_na = dot(square(y(not_obs) - mean(y)), w(not_obs));
  if(n<2){ // only one observation, cannot split
    out(0) = vnce_na; // total variance
    out(1) = 0; // cut value --- not applicable here
    return(out);
  }
  x = x(obs);
  y = y(obs);
  w = w(obs); // weights associated with indexes
  uvec idx = sort_index(x); // ascending; may be duplicate values here (doesn't matter)
  x = x(idx); // sorted by x values
  y = y(idx); // sort y by x values
  w = w(idx); // weights sorted by x values
  mat vnce(2,n-1); // store variance resulting from split
  for(uword j=0; j<n-1; j++){
    vnce(0,j) = dot(square(y(span(0,j)) - 
      dot(y(span(0,j)), w(span(0,j)))/sum(w(span(0,j))) ), w(span(0,j))); // <= variance
    vnce(1,j) = dot(square(y(span(j+1,n-1)) - 
      dot(y(span(j+1,n-1)), w(span(j+1,n-1)))/sum(w(span(j+1,n-1))) ), w(span(j+1,n-1))); // > variance
  }
  // Rcpp::Rcout << "turd" << endl;
  vec tot_vnce = trans(sum(vnce,0)); // total variance after cut (not including NA)
  uword imin = index_min(tot_vnce);
  out(0) = tot_vnce(imin) + vnce_na;
  out(1) = (x(imin) + x(imin+1))/2;
  return(out);
}

// [[Rcpp::export]]
arma::field<arma::vec> get_split_details(arma::vec x, // RHS data
                                         arma::vec y, // LHS data
                                         double c, // cut
                                         arma::vec ind, // row index of observation
                                         arma::vec w){ // corresponding weight
  uvec obs = find_finite(x);
  uvec not_obs = find_nonfinite(x);
  // uword n = obs.n_elem;
  uword nnot = not_obs.n_elem;
  vec x_obs = x(obs);
  vec y_obs = y(obs);
  vec w_obs = w(obs); 
  uvec lidx = find(x_obs<=c);
  uvec gidx = find(x_obs>c);
  double mu_less = dot(y_obs(lidx), w_obs(lidx))/sum(w_obs(lidx)); // weighted mean
  double mu_greater = dot(y_obs(gidx), w_obs(gidx))/sum(w_obs(gidx)); // weighted mean
  double sig_less = dot(square(y_obs(lidx) - mu_less), w_obs(lidx));
  double sig_greater = dot(square(y_obs(gidx) - mu_greater), w_obs(gidx));
  double sig_NA = 0;
  if(nnot>0){
    sig_NA = dot(square(y(not_obs) - mean(y)), w(not_obs));
  }
  uvec ind_less = unique(join_uvec(obs(lidx), not_obs));
  uvec ind_greater = unique(join_uvec(obs(gidx), not_obs));
  vec w_less = w;
  vec w_greater = w;
  double p_less; double p_greater; double p_sum; uword l;
  for(uword j=0; j<nnot; j++){
    l = not_obs(j);
    p_less = normpdf( (y(l)-mu_less)/(sig_less + datum::eps) );
    p_greater = normpdf( (y(l)-mu_greater)/(sig_greater + datum::eps) );
    p_sum = p_less + p_greater;
    w_less(l) = w_less(l)*(p_less/p_sum);
    w_greater(l) = w_greater(l)*(p_greater/p_sum);
  }
  field<vec> out(5);
  out(0) = {mu_less, mu_greater, sig_less, sig_greater, sig_NA}; // pars
  out(1) = ind(ind_less); // index <=
  out(2) = ind(ind_greater); // index >
  out(3) = w_less(ind_less); // weight <.
  out(4) = w_greater(ind_greater);  // weight >
  return out;
}

// find the best series in X to identify y
// [[Rcpp::export]]
arma::field<arma::vec> best_split(arma::mat X, // predictors
                       arma::vec y, // response
                       arma::vec ind, // index of observations in original data
                       arma::vec w, // weight corresponding to ind
                       arma::uword n){ // number of candidates to select
  // select candidates for splitting
  uvec candidates = select_rnd(X.n_cols, n); // randomly selected candidate series
  mat splits(2,n); // results matrix (mean <, mean >, vol <, vol >, vol NA, cut)
  vec tmp;
  for(uword j=0; j<n; j++){
    splits.col(j) = find_split(X.col(candidates(j)), y, w);
  }
  uword min_idx = index_min(splits.row(0)); // which series offers max reduction in variance of y
  field<vec>out = get_split_details(X.col(candidates(min_idx)), y,
                                     splits(1, min_idx), ind, w);
  // clunky but effective
  vec par_out(7);
  par_out(0) = candidates(min_idx); // index of which variable to split on
  par_out(1) = splits(1, min_idx);
  par_out(span(2,6)) = out(0);
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

// // [[Rcpp::export]]
// arma::uvec testfun(arma::vec x){
//   uvec out = find(x != x || x<2);
//   return(out);
// }

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
arma::mat RegTreeM(arma::vec y, // response (no missing obs)
          arma::mat X, // predictors (missing obs OK)
          arma::vec depth_range){
          // double bag_rows = 0.632,
          // double bag_cols = 0.333){
  // Bag each tree by randomly selecting observations
  // bag_rows = std::min(1.0,std::abs(bag_rows)); // safety first
  // bag_rows = std::min(1.0,std::abs(bag_cols)); // safety first
  double T = X.n_rows;
  uvec to_keep = select_rnd(T, ceil(0.632*T)); // select rows to keep
  uword max_nodes= ceil(as_scalar((depth_range(1)-depth_range(0))*randu<vec>(1) +
                              depth_range(0))); // randomise depth of model
  X = X.rows(to_keep); // this shuffles X and y but it shouldn't matter
  y = y(to_keep);
  double xnc = X.n_cols;
  uword n = ceil(xnc/3); // number of candidates (columns) to use at each split
  mat Tree(max_nodes, 8, fill::zeros);
  Tree(0,5) = mean(y); // unconditional mean
  Tree(0,6) = sum(square(y-Tree(0,5))); // unconditional var
  mat leaves; vec par; uvec fobs; mat tmp_tree;
  uvec leaf_idx;
  // double var_new = 0;
  field<vec> W(max_nodes);
  field<uvec> I(max_nodes);
  W(0) = ones<vec>(X.n_rows);
  I(0) = regspace<uvec>(0,X.n_rows-1); // index in original data (after bagging)
  field<vec> tmp; // temp output from best_splits()
  uword i = 0;
  uword j = 0; // leaf to work on this iteration
  while(i+2<max_nodes){
    tmp = best_split(X.rows(I(j)), y(I(j)), conv_to<vec>::from(I(j)), W(j), n);
    par = tmp(0);
    Tree(j,0) = par(0); // which variable
    Tree(j,1) = par(1); // cut value
    Tree(j,2) = i+1; // if <, go to node i+1
    Tree(j,3) = i+2; // if >, go to node i+2
    Tree(j,7) = 0; // no longer terminal node
    Tree(i+1,4) = j; // node coming from
    Tree(i+2,4) = j; // node coming from
    Tree(i+1,5) = par(2); // mu if <=
    Tree(i+2,5) = par(3); // mu if >
    Tree(i+1,6) = par(4); // sig if <=
    Tree(i+2,6) = par(5); // sig if >
    Tree(i+1,7) = 1; // terminal node (leaf)
    Tree(i+2,7) = 1; // terminal node (leaf)
    W(i+1) = tmp(3); // residuals <=
    W(i+2) = tmp(4); // residuals >
    I(i+1) = conv_to<uvec>::from(tmp(1)); // indexes <=
    I(i+2) = conv_to<uvec>::from(tmp(2)); // indexes >
    // find the leaf with the highest variance for the next iteration
    fobs = field_obs(W, i+3);
    tmp_tree = Tree.rows(0,i+2);
    leaf_idx = find(tmp_tree.col(7) == 1 && fobs > 1); // index of terminal nodes in Tree matrix with enough obs
    // Find the leaf with the maximum variance
    leaves = Tree.rows(leaf_idx); // terminal nodes (i.e. leaves).
    // Rcpp::Rcout << Tree.row(j) << endl;
    // var_new = sum(par(span(4,6))); // volatility at new node
    // Rcpp::Rcout << "old = " << var_old << "new = " << var_new << "next old = " << Tree(j,6) << endl;
    // Rcpp::Rcout << (Tree(j,7) - var_new)/Tree(j,7) << endl;
    // if((var_old - var_new)/var_old < threshold) break; // if volatility doesn't improve, break loop
    i += 2; // split result in 2 new rows
    if(leaf_idx.n_elem==0) break;
    j = leaf_idx(index_max(leaves.col(6))); // index of leaf with max volatility to work on next
  }
  // Rcpp::List rtrn;
  // mat tree = Tree.rows(0,i);
  // arma::field<vec> Eout = E.rows(0,i);
  // arma::field<uvec> Iout = I.rows(0,i);
  // rtrn["Tree"] = tree;
  // rtrn["E"] = Eout;
  // rtrn["I"] = Iout;
  return(Tree.rows(0,i));
}

// Fit a single observation using the estimated tree
// [[Rcpp::export]]
double FitVecM(arma::vec x,
              arma::mat Tree,
              arma::uword maxit = 1000){
  uword j = 0; uword it=0; double y=0; uword i = 0;
  while(Tree(j,7) != 1 && it<maxit){
    if(!std::isfinite(x(Tree(j,0)))){
      return(y);
    }else{
      i = Tree(j,0);
      if(x(Tree(j,0))>Tree(j,1)){
        j = Tree(j,3);
      }else{
        j = Tree(j,2);
      }
      it++;
    }
  }
  y = Tree(j,5); // terminal node (leaf)
  return(y);
}

// Fit many observations using the estimated tree
// [[Rcpp::export]]
arma::vec FitMatM(arma::mat X,
                  arma::mat Tree){
  vec Mu(X.n_cols);
  for(uword j=0; j<X.n_cols; j++){
    Mu(j) = FitVecM(X.col(j), Tree);
  }
  return(Mu);
}

// Draw 'draws' number of trees
// [[Rcpp::export]]
arma::field<arma::mat> RegForestM(arma::vec y, // response (no missing obs)
                           arma::mat X, // predictors (missing obs OK)
                           arma::vec depth_range, // try 15 too
                           arma::uword draws = 1000){
  field<mat> Trees(draws);
  for(uword j = 0; j<draws; j++){
    Trees(j) = RegTreeM(y, X, depth_range);
  }
  return(Trees);
}

// Fit output from RegForest
// [[Rcpp::export]]
arma::vec FitFieldM(arma::mat X,
                    arma::field<arma::mat> Trees){
  mat Mu(X.n_rows, Trees.n_elem);
  X = trans(X); //transpose for FitMat
  for(uword j=0; j<Trees.n_elem; j++){
    Mu.col(j) = FitMatM(X, Trees(j));
  }
  vec mu = mean(Mu,1); // take average response
  return(mu);
}




