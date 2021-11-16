

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

all_splits <- function(Trees, X_names = NULL){
  idx <- sapply(Trees, get_all_splits)
  if(is.list(idx)) idx <- do.call("c", idx)
  idx <- idx + 1
  tab <- table(idx)
  if(!is.null(X_names)) names(tab) <- X_names[as.numeric(names(tab))]
  return(sort(tab, decreasing = TRUE))
}

reg_forest <- function(y, X, max_nodes = 64, draws = 1000, steps = 1, regression = TRUE, return_trees = TRUE){
  y <- c(y)
  X <- as.matrix(X)
  y_finite <- is.finite(y) # fit model on periods in which y is finite
  last_period <- max(which(y_finite)) + steps # the period we want to predict
  k <- min(NROW(X), last_period)
  X <- X[ ,is.finite(X[k, ])] # if X not observed in the period we wish to predict, drop it
  
  if(regression){
    Trees <- RegForest(y[y_finite], X[y_finite, ], max_nodes, draws) # estimate model
  }else{
    Trees <- RndForest(y[y_finite], X[y_finite, ], max_nodes, draws) # estimate model
  }
  
  
  cnames <- colnames(X)
  fstsplt <- first_split(Trees, cnames)
  allsplt <- all_splits(Trees, cnames)
  
  if(regression){
    fit <- FitField(X[seq(k), ], Trees) # in sample fit
  }else{
    fit <- StdFitField(X[seq(k), ], Trees) # in sample fit
  }
  
  out <- list(fit = fit, true_vals = y[seq(k)])
  if(return_trees) out$Trees <- Trees
  out$first_split <- fstsplt
  out$all_splits <- allsplt
  out$X_names <- cnames
  return(out)
}


# x <- X[359, ]

# FitVec(x, Tree)
# pretty_plot(cbind(fit, y[seq(last_period)]))

