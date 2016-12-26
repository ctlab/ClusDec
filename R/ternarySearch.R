ternarySearch <- function(fun, range, eps=1e-9, maxiter=100, ...) {
    iter <- 0
    l <- range[1]
    r <- range[2]

    while (r - l > eps && iter < maxiter) {
        diff3 <- (r - l) / 3
        m1 <- l + diff3
        m2 <- l + 2 * diff3

        f1 <- fun(m1, ...)
        f2 <- fun(m2, ...)
        iter <- iter + 1

        if (f1 < f2) {
            r <- m2
            next
        }
        if (f1 > f2) {
            l <- m1
            next
        }
        r <- m2
    }
    return(l)
}
