



pretty_plot <- function(x = NULL, X, lwd = 2, xlab = "Date", ylab = "", legend_pos = "bottomleft", title = ""){
  
  if(is.null(x)) x <- 1:NROW(X)
  k <- ncol(X)
  color <- rep(c("black", "deeppink", "deepskyblue3", "green3", "tan3", "turquoise"),5)
  color <- color[1:k]
  matplot(x, X, type = 'l', lty = 1, col = color, lwd = lwd, xlab = xlab, ylab = ylab)
  grid(col = "lightslategrey")
  legend(legend_pos, colnames(X),
         col = color, lwd = lwd, bty = "n")
  title(title)
}
                            
