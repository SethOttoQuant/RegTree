// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace arma;
using namespace Rcpp;
#include "utils.h"
#include "rfbase.h"


//[[Rcpp::export]]
arma::uword find_unique(arma::uvec v,
                        double k
){
  uword j;
  for(j=0; j<v.n_elem; j++){
    if(v(j) == k) break;
  }
  return(j);
}

// This function randomly draws a cut point where cut points with a higher
// reduction in MSE have a higher probability of being selected according to
// the weight_pow parameter
// [[Rcpp::export]]
arma::vec sim_cut(arma::vec x, // predictor
                     arma::vec y, // response
                     double prior_shrink, // prior shrinkage par
                     double my,
                     double weight_pow){ // determines probability of cut points
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
  // Rcpp::Rcout<< cdf << endl;
  double draw = randu()*cdf(n-4); // rand u var over full range of cdf
  uword j=0;
  while(draw>cdf(j)){
    j++;
  }
  return(out.col(j)); // mat with row 0 = variance, row 1 = cut value
}
// [[Rcpp::export]]
arma::field<arma::vec> sim_pars(arma::vec x, // predictor
                               arma::vec y, // response
                               arma::uvec ind, // original index
                               double imin,
                               double my,
                               double prior_shrink){  // where to cut
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
// draw a series in X to identify y; probability of drawing "best" is determined
// by the weight_pow parameter
// [[Rcpp::export]]
arma::field<arma::vec> sim_split(arma::mat X, // predictors
                                   arma::vec y, // response
                                   arma::uvec ind, // index of observations in original data
                                   double prior_shrink, // shrink estimate towards unconditional mean
                                   double my,   // unconditional mean
                                   double weight_pow){
  arma::mat candidates(2, X.n_cols);
  for(uword j=0; j<X.n_cols; j++){
    // Rcpp::Rcout << candidates(j) << endl;
    candidates.col(j) = sim_cut(X.col(j), y, prior_shrink, my, weight_pow); 
  }
  // this bit is just to select the split with probability weight given by out.row(0)
  vec cdf = candidates.row(0).t();
  cdf /= mean(cdf);
  cdf = cumsum(pow(cdf, -weight_pow));
  double draw = randu()*cdf(cdf.n_elem-1); // rand u var over full range of cdf
  uword j=0;
  while(draw>cdf(j)){
    j++;
  }
  field<vec> out = sim_pars(X.col(j), y, ind, candidates(1,j), my, prior_shrink); // re-running cut here to keep things clean
  // clunky but effective
  vec par_out(9);
  par_out(0) = j;
  par_out(span(1,8)) = out(0);
  out(0) = par_out;
  return(out); // for find_cut() see function above
}

struct TreeReturn
{
  arma::mat Tree;
  arma::field<uvec> I;
};

TreeReturn simtree(arma::vec y, // response (no missing obs)
          arma::mat X, // predictors (missing obs OK), obs in rows series in cols
          arma::mat Tree_in, // mat Tree(max_nodes, 10, fill::zeros);
          arma::field<uvec> I_in, // index corresponding to each node of Tree: 
          double prior_shrink = 0,
          double weight_pow=2, // increase to increase probability of selecting best split variable
          arma::uword max_nodes=31,
          arma::uword j=0){  // prune from this node onward
  double my=mean(y); double vy = mean(square(y-my));
  TreeReturn rtrn;
  mat leaves; vec par; mat tmp_tree;
  uvec leaf_idx; 
  field<vec> tmp; // temp output from best_splits()
  // uword i = j;  // next nodes  
  mat Tree(max_nodes, 10, fill::zeros);
  field<uvec> I(max_nodes);
  uword i = 0;
  if(j==0){
    Tree(0,5) = my; // unconditional mean
    Tree(0,6) = vy;
    Tree(0,7) = y.n_elem;
    I(0) = regspace<uvec>(0,X.n_rows-1);
  }else if(Tree_in(j,8)==1){ // selected a terminal node; no change to model
    rtrn.Tree = Tree_in;
    rtrn.I = I_in;
    return(rtrn);
  }else{
    // Rcpp::Rcout << Tree_in << endl;
    // --------------------------------------
    // Pruning
    // Need to wipe out any nodes that come after the start node and reorganize the tree 
    // --------------------------------------
    uvec to_wipe(Tree_in.n_rows, fill::zeros);
    to_wipe(0) = Tree_in(j, 2);  // go to <
    to_wipe(1) = Tree_in(j, 3);  // go to > 
    Tree_in(j,8) = 1; // now a terminal node
    uword g=0;
    uword h=2;
    // Rcpp::Rcout << "to_wipe(g): " << to_wipe(g) << endl;
    while(g<Tree_in.n_rows){  // get rid of all nodes that follow the start node
      if(Tree_in(to_wipe(g), 8) == 0){
        to_wipe[h] = Tree_in(to_wipe[g], 2);  
        to_wipe[h+1] = Tree_in(to_wipe[g], 3);
        h+=2;
      }
      g++;
      if(g==h) break;  // all lower nodes selected
    }
    // Rcpp::Rcout << to_wipe << endl;
    vec tmp_vec = Tree_in.col(9);
    tmp_vec(to_wipe.rows(0,h-1)).zeros(); // Tree_in.rows(to_wipe.rows(0, h-1)) = zeros<mat>(h, 10); // fill rows with zeros
    uvec keep_idx = find(tmp_vec); // get index of rows to keep; this should preserve order
    // Rcpp::Rcout << keep_idx << endl;
    i = keep_idx.n_elem-1; // where to start adding new nodes
    for(uword rw=0; rw<i+1; rw++){
      Tree.row(rw) = Tree_in.row(keep_idx(rw));
      I(rw) = I_in(keep_idx(rw));
      Tree(rw, 2) = find_unique(keep_idx, Tree(rw, 2)); // update go-to row
      Tree(rw, 3) = find_unique(keep_idx, Tree(rw, 3)); // update go-to row
      Tree(rw, 4) = find_unique(keep_idx, Tree(rw, 4)); // update from row
    }
    // select terminal node to start on
    tmp_tree = Tree.rows(0,i);
    leaf_idx = find(tmp_tree.col(8) == 1 && tmp_tree.col(7) > 3); // index of terminal nodes in Tree matrix with enough obs
    leaves = Tree.rows(leaf_idx); // terminal nodes (i.e. leaves). 
    j = leaf_idx(floor(randu()*leaf_idx.n_elem)); // index_max(leaves.col(6)) totally randomized for now, in the future might want to be more selective.
  }
  // --------------------------------------
  // Calculate new Tree
  // --------------------------------------
  while(i+2<max_nodes){
    tmp = sim_split(X.rows(I(j)), y(I(j)), I(j), prior_shrink, my, weight_pow);
    par = tmp(0);
    if(any(par(span(1,4)))){ //  enough observations to split
      Tree(j,0) = par(0); // which variable
      Tree(j,1) = par(6); // cut value
      Tree(j,2) = i+1; // if <, go to node i+1
      Tree(j,3) = i+2; // if >, go to node i+2
      Tree(j,8) = 0; // no longer terminal node: 1=terminal node
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
    // draw which terminal node to work on next
    tmp_tree = Tree.rows(0,i);
    leaf_idx = find(tmp_tree.col(8) == 1 && tmp_tree.col(7) > 3); // index of terminal nodes in Tree matrix with enough obs
    leaves = Tree.rows(leaf_idx); // terminal nodes (i.e. leaves). 
    if(leaf_idx.n_elem==0) break;
    j = leaf_idx(floor(randu()*leaf_idx.n_elem)); // index_max(leaves.col(6)) totally randomized for now, in the future might want to be more selective.
  }
  for(uword h=0; h<max_nodes; h++){
    Tree(h,9) = max_nodes - h; // descending index for sorting (because blank is zero)
  }
  rtrn.Tree = Tree.rows(0,i);
  rtrn.I = I.rows(0,i);
  return(rtrn);
}

// [[Rcpp::export]]
Rcpp::List Rsimtree(arma::vec y, // response (no missing obs)
                   arma::mat X, // predictors (missing obs OK), obs in rows series in cols
                   double prior_shrink = 0,
                   double weight_pow=2, // increase to increase probability of selecting best split variable
                   arma::uword max_nodes=21,
                   arma::uword j=0){
  
  mat Tree_in(0, 10, fill::zeros);
  arma::field<uvec> I_in(0);
  TreeReturn out = simtree(y, // response (no missing obs)
                           X, // predictors (missing obs OK), obs in rows series in cols
                           Tree_in, // mat Tree(max_nodes, 10, fill::zeros);
                           I_in, // index corresponding to each node of Tree: 
                           prior_shrink,
                           weight_pow,
                           max_nodes,
                           j);
  
  Rcpp::Rcout << out.Tree << endl;
  
  out = simtree(y, X, out.Tree, out.I, 0, 2, 21, 3);
  
  Rcpp::List rtrn(2);
  rtrn(0) = out.Tree;
  rtrn(1) = out.I;
  return(rtrn);
}

// [[Rcpp::export]]
Rcpp::List simforest(arma::vec y, // response (no missing obs)
                    arma::mat X, // predictors (missing obs OK), obs in rows series in cols
                    arma::uword max_nodes = 21,
                    double prior_shrink = 0,
                    double weight_pow=2,
                    arma::uword burn = 100,
                    arma::uword reps = 1000){
  // uword max_nodes = 2*splits + 1;
  uword T = y.n_elem;
  uvec terminal, internal;
  vec tmp_variance;
  double var=0;
  double a=0;
  mat Tree_in(0, 10, fill::zeros), Tree_candidate, Tree_terminal;
  field<uvec> I_in(0);
  uword j=0; 
  TreeReturn out;
  double sig = sum(square(y - mean(y)))/(2*T); // whatever.... 
  double k = -log(sqrt(sig*2*datum::pi)); double k2 = 1/(2*sig);
  double ll_new=0; 
  double ll_old=0;
  // Burn in
  for(uword it=0; it<burn; it++){
    out = simtree(y, // response (no missing obs)
                   X, // predictors (missing obs OK), obs in rows series in cols
                   Tree_in, // mat Tree(max_nodes, 10, fill::zeros);
                   I_in, // index corresponding to each node of Tree: 
                   prior_shrink,
                   weight_pow,
                   max_nodes,
                   j);
    Tree_candidate= out.Tree;
    terminal = find(Tree_candidate.col(8)>0); // terminal nodes
    tmp_variance = Tree_candidate.col(6); // variance at each node
    var = sum(tmp_variance(terminal)); // total variance of the model
    ll_new = T*k - k2*var;
    if(it==0 or ll_new>ll_old){
      a=1;
    }else{
      a = exp(ll_new - ll_old);
    }
    // Rcpp::Rcout << "\r" << a ;
    if(a > randu()){
      Tree_in = out.Tree;
      I_in = out.I;
      ll_old = ll_new;
    }
    internal = find(Tree_in.col(8)==0); //non-terminal nodes
    j = internal(floor(randu()*internal.n_elem)); // randomly select a new starting node
  }
  field<mat> Trees(reps);
  field<mat> tmp;
  mat fit(T,reps, fill::zeros); 
  cube FC(T, X.n_cols, reps, fill::zeros);
  vec mse(reps, fill::zeros);
  for(uword it=0; it<reps; it++){
    out = simtree(y, // response (no missing obs)
                  X, // predictors (missing obs OK), obs in rows series in cols
                  Tree_in, // mat Tree(max_nodes, 10, fill::zeros);
                  I_in, // index corresponding to each node of Tree: 
                  prior_shrink,
                  weight_pow,
                  max_nodes,
                  j);
    Tree_candidate= out.Tree;
    terminal = find(Tree_candidate.col(8)>0); // terminal nodes
    tmp_variance = Tree_candidate.col(6); // variance at each node
    var = sum(tmp_variance(terminal)); // total variance of the model
    ll_new = T*k - k2*var;
    if(it==0 or ll_new>ll_old){
      a=1;
    }else{
      a = exp(ll_new - ll_old);
    }
    // Rcpp::Rcout << "\r a: " << a << " it: " << it ;
    if(a > randu()){
      Tree_in = out.Tree;
      I_in = out.I;
      ll_old = ll_new;
    }
    internal = find(Tree_in.col(8)==0); //non-terminal nodes
    j = internal(floor(randu()*internal.n_elem)); // randomly select a new starting node
    Trees(it) = Tree_in;
    tmp = fitmat(trans(X), Tree_in);
    fit.col(it) = tmp(0);
    mse(it) = sum(square(y - tmp(0)));
    FC.slice(j) = trans(tmp(1));
  }
  Rcpp::List rtrn;
  rtrn["Trees"] = Trees;
  rtrn["OOB"] = fit; // not out of bag, but keeping naming the same
  rtrn["FC"] = FC;
  rtrn["MSE"] = mse;
  return(rtrn);
}


