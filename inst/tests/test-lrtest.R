n <- 100
w.x <- rep(1, n)
w.y <- rep(1, n)
x <- rnorm(n, 10)
y <- rnorm(n, 8)
lrtest(w.x, w.y, x, y)

w.x <- rep(1, n)
w.y <- rep(0, n)
x <- rnorm(n, 10)
y <- integer(0)
lrtest(w.x, w.y, x, y)
