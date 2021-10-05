



pretty_plot <- function(X, x = NULL, lwd = 2, xlab = "Date", ylab = "", legend_pos = "bottomleft", title = ""){
  
  if("ref_date"%in%colnames(X)){
    X <- data.frame(X)
    x <- X[ , "ref_date"]
    X <- X[ , !colnames(X)%in%"ref_date"]
  }else if(is.null(x)) x <- 1:NROW(X)
  k <- ncol(X)
  color <- rep(c("black", "deeppink", "deepskyblue3", "green3", "tan3", "turquoise"),5)
  color <- color[1:k]
  matplot(x, X, type = 'l', lty = 1, col = color, lwd = lwd, xlab = xlab, ylab = ylab)
  grid(col = "lightslategrey")
  legend(legend_pos, colnames(X),
         col = color, lwd = lwd) #  bty = "n"
  title(title)
}

reg_forest <- function(y, X, max_nodes = 64, threshold = 0.02, draws = 500, steps = 1, treehugger = TRUE){
  y <- c(y)
  X <- as.matrix(X)
  y_finite <- is.finite(y) # fit model on periods in which y is finite
  last_period <- max(which(y_finite)) + steps # the period we want to predict
  X <- X[ ,is.finite(X[last_period, ])] # if X not observed in the period we wish to predict, drop it
  
  Trees <- RegForest(y[y_finite], X[y_finite, ], max_nodes, threshold, draws) # estimate model
  
  fit <- FitField(X[seq(last_period), ], Trees) # fit model
  out <- list(fit = fit, true_vals = y[seq(last_period)])
  if(treehugger) out$Trees <- Trees
  return(out)
}


# x <- X[359, ]

# FitVec(x, Tree)
# pretty_plot(cbind(fit, y[seq(last_period)]))

