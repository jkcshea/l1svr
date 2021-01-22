rm(list = ls())
library(MASS)
source('svr-functions.R')

## Generate random data
set.seed(10L)
n <- 1000
mu <- c(1, 2)
Sigma <- matrix(c(1, -2, 3, 2), nrow = 2)
Sigma <- Sigma %*% t(Sigma)
sdMat <- diag(1/sqrt(diag(Sigma)))
corrMat <- sdMat %*% Sigma %*% sdMat
beta <- c(-2, 1, 0.5)
X <- as.matrix(mvrnorm(n = n, mu = mu, Sigma = Sigma))
x <- cbind(1, X,
           sample(x = c(1, 2, 3),
                  prob = c(0.2, 0.3, 0.5),
                  size = n,
                  replace = TRUE))
## Construct hollow residuals


u <- rnorm(n, 0, 12)
y <- x[, 1:3] %*% beta + u

## Generate the data set a user would have
simdata <- data.frame(cbind(y, x[, -1]))
colnames(simdata) <- c('y', 'x1', 'x2')

## Test that regression works
fit1 <- lm(formula = y ~ 1 + x1 + x2, data = simdata)
fit1
summary(fit1)


## Set the SVM parameters
tmpEpsilon <- 7
tmpLambda <- 20

## Now try your wrapper function, point estimates only
res1 <- l1svr(formula = y ~ 1 +  x1 + x2, data = simdata,
      epsilon = tmpEpsilon, lambda = tmpLambda, inference = FALSE)
res1
summary(res1)

## Homoskedastic inference
res2 <- l1svr(formula = y ~ 1 +  x1 + x2, data = simdata,
      epsilon = tmpEpsilon, lambda = tmpLambda, inference = TRUE)
res2
summary(res2)

## Heteroskedastic inference
res3 <- l1svr(formula = y ~ 1 +  x1 + x2, data = simdata,
              epsilon = tmpEpsilon, lambda = tmpLambda, inference = TRUE,
              heteroskedastic = TRUE, h = 2, kappa = 1.75)
summary(res3)


table(res3$dualPositive)
table(res3$dualNegative)



test <- function(x) {
    x^2
}


gridSearch(test, init = 0, target = 4, increment = -1, left = TRUE)
gridSearch(test, init = 0, target = 4, increment = -1, left = FALSE)


## -----------------------------------------------------------------------------
## Check confidence interval
## -----------------------------------------------------------------------------

source('svr-functions.R')
## Homoskedastic inference
res4 <- l1svr(formula = y ~ 0 +  x1 + x2, data = simdata,
              epsilon = tmpEpsilon, lambda = tmpLambda,
              inference = TRUE,
              confidence = FALSE)

source('svr-functions.R')
summary(res4)
fullRes <- l1svr(formula = y ~ 0 +  x1 + x2, data = simdata,
                 epsilon = tmpEpsilon, lambda = tmpLambda,
                 inference = TRUE,
                 confidence = FALSE,
                 confidence.iter = 20,
                 confidence.level = 0.9)
summary(fullRes)


source('svr-functions.R')
fit3 <- l1svr(formula = y ~ 1 + x1 + x2,
              data = simdata,
              epsilon = 7,
              lambda = 20,
              confidence = TRUE,
              confidence.level = 0.95,
              confidence.iter = 20)
summary(fit3)

fit3$ci


## You want to vary the amount you concentrate out, that is all you
## have to do. So for each variable, you do the grid search over the
## values that are concenterated out.

