# An example of random forest with regression at nodes

library(RegTree)

# (1) Simulate some data
x <- 1:40
y <- rep(0,40)
y[1:20] <- x[1:20] + 5*rnorm(20) 
y[21:40] <- -x[21:40] + 5*rnorm(20) 
y[1:10] <- y[1:10] - 10
y[11:20] <- y[11:20] + 10
y[21:30] <- y[21:30] + 10

X <- cbind(rnorm(40), rnorm(40), x, rnorm(40), rnorm(40))



tree <- RegTree(y, X, bag_rows = 1, bag_cols = 1, max_nodes = 7)

# splits
y1 <- tree[1,6] + tree[2,6] + (1:20)*tree[2,7]
y2 <- tree[1,6] + tree[3,6] + (21:40)*tree[3,7]
y3 <- tree[1,6] + tree[2,6] + (1:10)*tree[2,7] + tree[4,6] + (1:10)*tree[4,7]
y4 <- tree[1,6] + tree[2,6] + (11:20)*tree[2,7] + tree[5,6] + (11:20)*tree[5,7]
y5 <- tree[1,6] + tree[3,6] + (21:30)*tree[3,7] + tree[6,6] + (21:30)*tree[6,7]
y6 <- tree[1,6] + tree[3,6] + (31:40)*tree[3,7] + tree[7,6] + (31:40)*tree[7,7]

r1 <- rep(mean(y[1:20]), 20)
r2 <- rep(mean(y[21:40]), 20)
r3 <- rep(mean(y[1:10]), 10)
r4 <- rep(mean(y[11:20]), 10)
r5 <- rep(mean(y[21:30]), 10)
r6 <- rep(mean(y[31:40]), 10)

plot(x,y)
lines(1:20, y1, col = "deeppink", lwd = 2)
lines(21:40, y2, col = "deeppink", lwd = 2)
lines(1:10, y3, col = "deepskyblue3", lwd = 2)
lines(11:20, y4, col = "deepskyblue3", lwd = 2)
lines(21:30, y5, col = "green3", lwd = 2)
lines(31:40, y6, col = "green3", lwd = 2)
grid()
legend("bottomleft", c("First Split", "Second Split", "Third Split"), lty = 1, lwd = 2,
       col = c("deeppink", "deepskyblue3", "green3"))
title("Regression Leaves")
dev.print(png, filename = "C:/Users/seton/Dropbox/system2/inflation/regression_leaf.png", width = 1000, height = 900, res = 150)



plot(x,y)
lines(1:20, r1, col = "deeppink", lwd = 2)
lines(21:40, r2, col = "deeppink", lwd = 2)
lines(1:10, r3, col = "deepskyblue3", lwd = 2)
lines(11:20, r4, col = "deepskyblue3", lwd = 2)
lines(21:30, r5, col = "green3", lwd = 2)
lines(31:40, r6, col = "green3", lwd = 2)
grid()
legend("bottomleft", c("First Split", "Second Split", "Third Split"), lty = 1, lwd = 2,
       col = c("deeppink", "deepskyblue3", "green3"))
title("Standard Leaves")
dev.print(png, filename = "C:/Users/seton/Dropbox/system2/inflation/standard_leaf.png", width = 1000, height = 900, res = 150)



