

# Calls for C++ functions
RegForest <- function(y,X,max_nodes = 31, draws = 1000) Reg_Forest(y,X,max_nodes,draws)
RndForest <- function(y, X, min_obs=5, max_nodes=1000, draws = 1000) rforest(y, X, min_obs, max_nodes, draws)
FitField <- function(X,Trees) Fit_Field(X,Trees)
StdFitField <- function(X,Trees) fitfield(X,Trees)
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

pretty_plot <- function(X, x = NULL, lwd = 2, xlab = "Date", ylab = "", legend_pos = "bottomleft", title = ""){
  
  if("ref_date"%in%colnames(X)){
    X <- data.frame(X)
    x <- X[ , "ref_date"]
    X <- X[ , !colnames(X)%in%"ref_date"]
  }else if(is.null(x)) x <- 1:NROW(X)
  k <- NCOL(X)
  color <- rep(c("black", "deeppink", "deepskyblue3", "green3", "tan3", "turquoise", "darkgrey", "steelblue", "violetred", "yellow3"),5)
  color <- color[1:k]
  matplot(x, X, type = 'l', lty = 1, col = color, lwd = lwd, xlab = xlab, ylab = ylab)
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

reg_forest <- function(y, X, min_obs=5, max_nodes = "auto", draws = 1000, 
                       steps = 1, regression = FALSE, return_trees = TRUE, 
                       orthogonal = FALSE){
  y <- c(y)
  X <- as.matrix(X)
  if(length(y) != NROW(X)) stop("Length of 'y' and rows of 'X' to not agree")
  y_finite <- is.finite(y) # fit model on periods in which y is finite
  last_period <- max(which(y_finite)) + steps # the period we want to predict
  k <- min(NROW(X), last_period)
  X <- X[ ,is.finite(X[k, ])] # if X not observed in the period we wish to predict, drop it
  if(orthogonal){
    X <- fake_pca(X)
  }
  if(max_nodes == "auto"){
    if(regression){
      max_nodes = 30
    }else{
      max_nodes = 1000
    }
  }
  if(regression){
    Trees <- RegForest(y[y_finite], X[y_finite, ], max_nodes, draws) # estimate model
  }else{
    Trees <- RndForest(y[y_finite], X[y_finite, ], min_obs, max_nodes, draws = draws) # estimate model
  }
  
  cnames <- colnames(X)
  if(!orthogonal){
    fstsplt <- first_split(Trees, cnames)
    allsplt <- all_splits(Trees, cnames, regression)
  }else{
    fstsplt <- NULL
    allsplt <- NULL
  }
  
  if(regression){
    fit <- FitField(X[seq(k), ], Trees) # in sample fit
  }else{
    fit <- StdFitField(X[seq(k), ], Trees) # in sample fit
  }
  
  fit_out <- rep(NA, )
  
  out <- list(fit = fit, true_vals = y[seq(k)])
  if(return_trees) out$Trees <- Trees
  out$first_split <- fstsplt
  out$all_splits <- allsplt
  out$X_names <- cnames
  out$idx <- seq(k)
  return(out)
}


# x <- X[359, ]

# FitVec(x, Tree)
# pretty_plot(cbind(fit, y[seq(last_period)]))

