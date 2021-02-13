library(MASS)
set.seed(10L)
n <- 500
mu <- c(1, 2)
Sigma <- matrix(c(1, -2, 3, 2), nrow = 2)
Sigma <- Sigma %*% t(Sigma)
sdMat <- diag(1/sqrt(diag(Sigma)))
corrMat <- sdMat %*% Sigma %*% sdMat
beta <- c(-2, 1, 0.5)
X <- as.matrix(mvrnorm(n = n, mu = mu, Sigma = Sigma))
x <- cbind(1, X)
u <- rnorm(n, 0, 12)
y <- x[, 1:3] %*% beta + u
simdata <- data.frame(cbind(y, x[, -1]))
colnames(simdata) <- c('y', 'x1', 'x2')
library("devtools")
use_data(simdata, overwrite = TRUE)
