##' Create a simple linear model
##' 
##' This creates a simple linear model.
##' 
##' @param y numeric.  A vector of the responses
##' @param X matrix.  The model matrix of predictors.
##' @param V weights.  Option weights to use to fit the model
##' @export
##' @examples
##' X <- matrix(c(rep(1, 20), 1:20), ncol = 2)
##' y <- X %*% c(1, 0.7) + rnorm(20)
##' o <- linmod(y, X)
linmod <- function(y, X, V = diag(length(y)) ){
    n <- length(y)
    rank <- qr(X)$rank
    df <- n - rank
    Vinv <- solve(V)
    Xt <- t(X)
    B <- solve(Xt %*% Vinv %*% X) %*% Xt %*% Vinv %*%y
    yhat <- X %*% B
    residuals <- y - yhat
    S2 <- drop(t(residuals) %*% Vinv %*% residuals / df)
    betas.cov <- S2 * solve(Xt %*% Vinv %*% X)
    results <- list(X = X,
                    y = y,
                    yhat = yhat,
                    residuals = residuals,
                    betas = B,
                    MSE = S2,
                    df = df,
                    betas.cov = betas.cov)
    class(results) <- "linmod"
    return(results)
}