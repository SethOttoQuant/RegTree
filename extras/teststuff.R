x <- matrix(rnorm(50), 10, 5)
y <- 2*x[,1] - x[,2] + x[,3] + rnorm(10)
ind <- c(3,4,5,7,8,9,10,15,16,17)
oot <- selectsplit(x, y, ind, 1, mean(y), weight_pow = 2)


q <- c(5.3872,   2.9550,   7.2107 ,  6.6894,   6.9126)
q <- q/mean(q)
cumsum(q^-2)


xsort <- sort.int(x, index.return = TRUE)
xsort$ix


yx <- y[xsort$ix]


sum((yx[1:5]-mean(yx[1:5]))^2)/5

xsort$x

x <- rnorm(10)
y <- 2*x + rnorm(10)

sue <- fast_cut(x, y, 0, my=mean(y), 10)




bob <- tst_sort(x)

bob[[1]]
bob[[2]]
x
bob[[3]]