x <- rnorm(10)
y <- 2*x + rnorm(10)

xsort <- sort.int(x, index.return = TRUE)
xsort$ix


yx <- y[xsort$ix]


sum((yx[1:5]-mean(yx[1:5]))^2)/5

xsort$x

sue <- fast_cut(x, y, 0, my=mean(y), 10)




bob <- tst_sort(x)

bob[[1]]
bob[[2]]
x
bob[[3]]