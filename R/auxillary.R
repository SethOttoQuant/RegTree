

# Calls for C++ functions
RegForest <- function(y,X, min_obs=15, max_nodes = 31, draws = 1000) Reg_Forest(y,X,min_obs,max_nodes,draws)
RndForest <- function(y, X, min_obs=5, max_nodes=1000, draws = 1000) rforest(y, X, min_obs, max_nodes, draws)
MAEForest <- function(y, X, min_obs=5, max_nodes=1000, draws = 1000) maeforest(y, X, min_obs, max_nodes, draws)
AltForest <- function(y, X, min_obs=5, max_nodes=1000, draws = 1000, geom_par = .5) rforest_alt(y, X, min_obs, max_nodes, draws, geom_par)
FitField <- function(X,Trees) Fit_Field(X,Trees)
FitFieldWeight <- function(X,Trees,weight){
  out <- Fit_Field_Weight(X,Trees,weight)
  colnames(out[[2]]) <- colnames(X)
  return(out)
} 
StdFitField <- function(X,Trees){
  out <- fitfield_alt(X,Trees)
  colnames(out[[2]]) <- colnames(X)
  return(out)
}
StdFitFieldWeight <- function(X,Trees,weight){ 
  out <- fitfield_alt(X,Trees,weight)
  colnames(out[[2]]) <- colnames(X)
  return(out)
}
RegTree <- function(y,X,max_nodes = 31) Reg_Tree(y,X,max_nodes)
StdRegTree <- function(y,X,min_obs=5,max_nodes=1000) regtree(y,X,min_obs,max_nodes)
quickreg <- function(x,y,r=0) QuickReg(x,y,r)

fake_pca <- function(W){
  W[is.na(W)] <- 0
  S = t(W)%*%W/NROW(W)
  eg <- eigen(S, symmetric = TRUE)
  C <- W%*%eg$vectors
  return(C)
}


# backfill <- function(x){
#   x_obs <- which(!is.na(x))
#   x_na <- which(is.na(x))
#   x_na <- x_na[x_na<max(x_obs)] # don't fill NA values in tail
#   k <- length(x_na)
#   if(k==0) return(x)
#   for(j in seq(k)){
#     idx <- min(x_obs[x_obs>x_na[j]])
#     x[x_na[j]] <- x[idx] 
#   }
#   return(x)
# }

backfill <- function(x) {   
  ind = which(!is.na(x))      # get positions of nonmissing values
  if(is.na(x[length(x)]))             # if it begins with a missing, add the 
    ind = c(ind, length(x))        # first position to the indices
  rep(x[ind], times = diff(c(0, ind)))   # repeat the values at these indices
  # diffing the indices + length yields how often 
}  

# must include ref_date and value... other columns optional
distribute_weekly <- function(df){
  tmp <- merge(data.table("ref_date"=seq.Date(min(df$ref_date)-6, max(df$ref_date), by = "day")), df, all = TRUE)
  tmp[ , value := backfill(value)/7]
  if("pub_date"%in%names(tmp)) tmp[ , pub_date := backfill(pub_date)]
  if("series_name"%in%names(tmp)) tmp[ , series_name := backfill(series_name)]
  if("pub_lag"%in%names(tmp)) tmp[ , pub_lag := as.numeric(pub_date - ref_date)]
  return(tmp)
}

pretty_plot <- function(X, x = NULL, lwd = 2, xlab = "Date", ylab = "", legend_pos = "bottomleft", title = "", ...){
  
  if("ref_date"%in%colnames(X)){
    X <- data.frame(X)
    x <- X[ , "ref_date"]
    X <- X[ , !colnames(X)%in%"ref_date"]
  }else if(is.null(x)) x <- 1:NROW(X)
  k <- NCOL(X)
  color <- rep(c("black", "deeppink", "deepskyblue3", "green3", "tan3", "turquoise", "darkgrey", "steelblue", "violetred", "yellow3"),5)
  color <- color[1:k]
  matplot(x, X, type = 'l', lty = 1, col = color, lwd = lwd, xlab = xlab, ylab = ylab, ...)
  grid(col = "lightslategrey")
  if(!is.null(colnames(X))){
    legend(legend_pos, colnames(X),
           col = color, lwd = lwd) #  bty = "n"
  }
  title(title)
}

first_split <- function(Trees, X_names = NULL){
  idx <- sapply(Trees, function(j) j[1,1]) + 1
  tab <- table(idx)
  if(!is.null(X_names)) names(tab) <- X_names[as.numeric(names(tab))]
  return(sort(tab, decreasing = TRUE))
}

get_all_splits <- function(Tree){
  Tree <- Tree[Tree[,9] != 1, ,drop=FALSE]
  return(Tree[,1])
}

all_splits <- function(Trees, X_names = NULL, regression = TRUE){
  idx <- sapply(Trees, get_all_splits)
  if(is.list(idx)) idx <- do.call("c", idx)
  idx <- idx + 1
  tab <- table(idx)
  if(!is.null(X_names)) names(tab) <- X_names[as.numeric(names(tab))]
  return(sort(tab, decreasing = TRUE))
}

reg_forest <- function(y, X, min_obs="auto", max_nodes = "auto", draws = 1000, 
                       steps = 1, type = c("standard", "regression", "alt", "mae"), return_trees = TRUE, 
                       orthogonal = FALSE, weight_by_mse = TRUE, geom_par=.5, weight_pow = 2){
  type <- match.arg(type)
  y <- c(y)
  X <- as.matrix(X)
  if(length(y) != NROW(X)) stop("Length of 'y' and rows of 'X' to not agree")
  y_finite <- is.finite(y) # fit model on periods in which y is finite
  if(!is.null(steps)){
    last_period <- max(which(y_finite)) + steps # the period we want to predict
    k <- min(NROW(X), last_period)
    y <- y[seq(k)]
    X <- X[seq(k),]
    y_finite <- y_finite[seq(k)]
    X <- as.matrix(X[ ,is.finite(X[k, ])]) # if X not observed in the period we wish to predict, drop it
  }else{ # if NULL ignore time series dimension of the problem and predict using everything
    k <- length(y)
  }
  if(orthogonal){
    X <- fake_pca(X)
  }
  if(max_nodes == "auto"){
    if(type == "regression"){
      max_nodes <- 40
    }else{
      max_nodes <- 1000
    }
  }
  if(min_obs == "auto"){
    if(type == "regression"){
      min_obs <- 15
    }else{
      min_obs <- 5
    }
  }
  if(type == "regression"){
    rf_out <- RegForest(y[y_finite], X[y_finite, ], min_obs, max_nodes, draws) # estimate model
  }else if(type=="standard"){
    rf_out <- RndForest(y[y_finite], X[y_finite, ], min_obs, max_nodes, draws) # estimate model
  }else if(type=="alt"){
    rf_out <- AltForest(y[y_finite], X[y_finite, ], min_obs, max_nodes, draws, geom_par) # estimate model
  }else if(type=="mae"){
    rf_out <- MAEForest(y[y_finite], X[y_finite, ], min_obs, max_nodes, draws) # estimate model
  }
  
  if(weight_by_mse){
    weight = (rf_out$MSE/mean(rf_out$MSE))^-weight_pow # standardize first to avoid very large numbers
    weight = weight/mean(weight) # standardize again
    oob <- rowMeans(t(apply(rf_out$OOB, 1, function(j) j*weight)), na.rm = TRUE)
    for(j in seq(draws)) rf_out$FC[,,j] <- rf_out$FC[,,j]*weight[j]
    fc_oob <- apply(rf_out$FC, c(1,2), mean, na.rm=TRUE) # contributions for each period
  }else{
    oob <- rowMeans(rf_out$OOB, na.rm = TRUE)
    fc_oob <- apply(rf_out$FC, c(1,2), mean, na.rm=TRUE) # contributions for each period
  }
  mse <- mean((y[y_finite]-oob)^2)
  
  cnames <- colnames(X)
  if(!orthogonal){
    fstsplt <- first_split(rf_out$Trees, cnames)
    allsplt <- all_splits(rf_out$Trees, cnames, regression)
  }else{
    fstsplt <- NULL
    allsplt <- NULL
  }
  
  if(type=="regression"){
    if(weight_by_mse){
      InSamp <- FitFieldWeight(X, rf_out$Trees, weight) # in sample fit
      fit <- InSamp[[1]]
      fc_fit <- InSamp[[2]]
    }else{
      InSamp <- FitField(X, rf_out$Trees) # in sample fit
      fit <- InSamp[[1]]
      fc_fit <- InSamp[[2]]
    }
  }else{
    if(weight_by_mse){
      InSamp <- StdFitFieldWeight(X, rf_out$Trees, weight) # in sample fit
      fit <- InSamp[[1]]
      fc_fit <- InSamp[[2]]
    }else{
      InSamp <- StdFitField(X, rf_out$Trees) # in sample fit
      fit <- InSamp[[1]]
      fc_fit <- InSamp[[2]]
    }
  }

  out <- list(in_samp_fit = fit, true_vals = y, mean_abs_feature_contribution = colMeans(abs(fc_fit)))
  if(return_trees) out$Trees <- rf_out$Trees
  if(weight_by_mse) out$weight <- weight
  out$MSE <- rf_out$MSE
  out_of_sample <- rep(NA, k)
  out_of_sample[y_finite] <- oob
  out_of_sample[!y_finite] <- fit[!y_finite]
  
  FC <- matrix(NA, NROW(X), NCOL(X))
  FC[y_finite, ] <- fc_oob
  FC[!y_finite, ] <- fc_fit[!y_finite, ]
  colnames(FC) <- colnames(X)
    
  out$out_of_sample <- out_of_sample
  out$mse <- mse
  out$oob_feature_contribution <- FC
  out$first_split <- fstsplt
  out$all_splits <- allsplt
  out$X_names <- cnames
  out$idx <- seq(k)
  return(out)
}

skew_loss <- function(x, y){
  y[y>0] <- y[y>0]^x[1]
  y[y>0] <- y[y>0]*abs(x[2])
  loss <- (100*mean(y, na.rm=TRUE))^2 + (sum(y^3, na.rm=TRUE))^2
  return(loss)
}

make_symetric <- function(y){
  y_finite <- is.finite(y)
  y_in <- y[y_finite]
  y_med <- median(y_in)
  y_in <- y_in - y_med
  out <- optim(par=c(1,1), fn=skew_loss, y=y_in)
  x <- out$par
  y_in[y_in>0] <- y_in[y_in>0]^x[1]
  y_in[y_in>0] <- y_in[y_in>0]*abs(x[2])
  y[y_finite] <- y_in
  return(list(
    y=y,
    x=x,
    med=y_med
  ))
}


# x <- X[359, ]

# FitVec(x, Tree)
# pretty_plot(cbind(fit, y[seq(last_period)]))

