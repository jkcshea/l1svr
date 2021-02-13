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
x <- cbind(1, as.matrix(mvrnorm(n = n, mu = mu, Sigma = Sigma)))

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
source('svr-functions.R')
res1 <- l1svr(formula = y ~ 1 +  x1 + x2, data = simdata,
              epsilon = tmpEpsilon, lambda = tmpLambda, inference = FALSE,
              solver = 'lpsolveapi')

res1
summary(res1)


source('svr-functions.R')
res1 <- l1svr(formula = y ~ 1 +  x1 + x2, data = simdata,
              epsilon = tmpEpsilon, lambda = tmpLambda, inference = FALSE,
              solver = 'lpsolveapi', h = 2, kappa = 1)

## Homoskedastic inference
source('svr-functions.R')
res2 <- l1svr(formula = y ~ 1 +  x1 + x2, data = simdata,
              epsilon = tmpEpsilon, lambda = tmpLambda, inference = TRUE,
              solver = 'gurobi')
summary(res2)


source('svr-functions.R')
res2 <- l1svr(formula = y ~ 1 +  x1 + x2, data = simdata,
              epsilon = tmpEpsilon, lambda = tmpLambda, inference = FALSE,
              solver = 'gurobi',  confidence.level = 0.95, confidence.iter = 5,
              h = 2, kappa = 1.75)
summary(res2)
res2$ci


resid <- y - x %*% res2$coef
dualDiff <- res2$dualNegative - res2$dualPositive
plot(x = resid, y = dualDiff)

par(mfrow = c(1, 2))
plot(x = resid, y = +res2$primalNegative)
plot(x = resid, y = -res2$primalPositive)

ttt <- res2$primalPositive - res2$primalNegative
cbind(ttt, resid)
hist(ttt -  resid)


head(resid)
tail(res2$primalPositive)
tail(res2$primalNegative)


load('cplexDual.Rdata')
load('lpDual.Rdata')

length(cplexDual)
length(lpDual)

max(abs(cplexDual - lpDual))


## Heteroskedastic inference
res3 <- l1svr(formula = y ~ 1 +  x1 + x2, data = simdata,
              epsilon = tmpEpsilon, lambda = tmpLambda, inference = TRUE,
              h = 2, kappa = 1.75)
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
              confidence.level = 0.95)

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
              confidence.iter = 10)
summary(fit3)

fit3$ci


## You want to vary the amount you concentrate out, that is all you
## have to do. So for each variable, you do the grid search over the
## values that are concenterated out.

