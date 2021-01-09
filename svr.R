#' SVR-l1 regression function
#'
#' This function carries out a SVR regression under l1 penalty.
#'
#' @param Y (vector) Dependent variable.
#' @param X (2d array) Independent variables, excluding column of
#'     constants for the intercept.
#' @param epsilon (real) Parameter defining soft-thresholding rule.
#' @param lambda (real) Tuning parameter to scale l1 penalty.
#' @return (list) A list including the coefficient estimates, and
#'     solutions to the primal and dual problem that can be used to
#'     identify the support vectors and recover the coefficient
#'     estimates.
svmRegress <- function(Y, X, epsilon, lambda) {
    X <- as.matrix(X)
    Y <- as.vector(Y)
    ## Primal problem:
    ## Declare all required maticex
    cvecPrimal <- c(0,
              rep(lambda, ncol(X)),
              rep(lambda, ncol(X)),
              rep(1, nrow(X)),
              rep(1, nrow(X)))
    AmatPrimal <- cbind(
        c(rep(1, nrow(X)), rep(-1, nrow(X))),
        rbind(-X, X),
        rbind(X, -X),
        rbind(-diag(nrow(X)), matrix(0, nrow = nrow(X), ncol = nrow(X))),
        rbind(matrix(0, nrow = nrow(X), ncol = nrow(X)), -diag(nrow(X))))
    bvecPrimal <- c(rep(epsilon, nrow(X)) + Y,
                  rep(epsilon, nrow(X)) - Y)
    ubPrimal <- c(Inf, rep(Inf, ncol(X)), rep(Inf, ncol(X)),
                  rep(Inf, nrow(X)), rep(Inf, nrow(X)))
    lbPrimal <- c(-Inf, rep(0, ncol(X)), rep(0, ncol(X)),
                  rep(0, nrow(X)), rep(0, nrow(X)))
    sensePrimal <- rep("<=", 2 * nrow(X))
    ## Pass matrices into Gurobi
    modelPrimal = list()
    modelPrimal$modelsense <- "min"
    modelPrimal$obj   <- cvecPrimal
    modelPrimal$A     <- AmatPrimal
    modelPrimal$rhs   <- bvecPrimal
    modelPrimal$sense <- sensePrimal
    modelPrimal$ub    <- ubPrimal
    modelPrimal$lb    <- lbPrimal
    resultPrimal <- gurobi::gurobi(modelPrimal, list(outputflag = 0))
    ## Prepare solutions
    solutions <- resultPrimal$x
    intercept <- solutions[1]
    betaMinus  <- solutions[2:(1 + ncol(X))]
    betaPlus <- solutions[(2 + ncol(X)):(1 + 2 * ncol(X))]
    auxiliaryPrimal <- solutions[(2 * ncol(X) + 2):length(solutions)]
    auxiliaryPrimal1 <- auxiliaryPrimal[1:nrow(X)] ## xi minus
    auxiliaryPrimal2 <- auxiliaryPrimal[(nrow(X) + 1):length(auxiliaryPrimal)] ## xi plus
    beta = c(intercept, betaPlus - betaMinus)
    ## Dual problem:
    ## Declare all required matrix
    cvecDual <- bvecPrimal
    AmatDual <- t(AmatPrimal)
    bvecDual <- cvecPrimal
    ubDual <- rep(0, 2 * nrow(X))
    lbDual <- rep(-Inf, 2 * nrow(X))
    senseDual <- c("=",
                   rep("<=", ncol(X)),
                   rep("<=", ncol(X)),
                   rep("<=", nrow(X)),
                   rep("<=", nrow(X)))
    ## Pass matrices into Gurobi
    modelDual = list()
    modelDual$modelsense <- "max"
    modelDual$obj   <- cvecDual
    modelDual$A     <- AmatDual
    modelDual$rhs   <- bvecDual
    modelDual$sense <- senseDual
    modelDual$ub    <- ubDual
    modelDual$lb    <- lbDual
    resultDual <- gurobi::gurobi(modelDual, list(outputflag = 0))
    ## Prepare solutions
    auxiliaryDual1 <- resultDual$x[1:nrow(X)] ## alpha minus
    auxiliaryDual2 <- resultDual$x[(nrow(X) + 1):(2 * nrow(X))] ## alpha plus
    ## Return output
    return(list(beta = beta,
                auxiliaryPrimal1 = auxiliaryPrimal1,
                auxiliaryPrimal2 = auxiliaryPrimal2,
                auxiliaryDual1  = auxiliaryDual1,
                auxiliaryDual2  = auxiliaryDual2))

}

#' Function to perform inference
#'
#' This function performs the SVR regression and calculates the p-value.
#'
#' @param Xr Matrix of covariates, excluding the covariate of interest
#'     for which we are conducting inference.
#' @param Zr Vector of the covariate of interest for which we are
#'     conducting inference.
#' @param Yr Vector of the outcome, where we have concentrated out
#'     \code{Zr}.
#' @param Y Vector of the outcome.
#' @param epsilon Real scalar, bandwidth for SVR.
#' @param lambda Real scalar, tuning parmaeter for SVR that controls
#'     the penalty.
#' @param nullGamma Real scalar, the coefficient on \code{Zr} under
#'     the null.
#' @param heteroskedastic Boolean, indicate whether data has
#'     heteroskedastic errors.
#' @param hc Real scalar, tuning parameter for density estimation.
#' @param kappa Real scalar, tuning parameter for density estimation.
#' @param tau Real scalar, tuning parameter for density estimation.
#' @param wald Boolean, set to \code{TRUE} if Wald confidence
#'     intervals are desired.
#' @return pvalue.
svmInference <- function(Xr, Zr, Yr, Y, epsilon, lambda = 0, nullGamma,
                         heteroskedastic = FALSE, hc = 0.5, kappa = 1,
                         tau = 0.5, wald = FALSE) {
    if (epsilon == 0 && wald) {
        stop("Wald inference curently not supported for QR.")
    }
    n <- length(Y)
    Xrc <- cbind(1, Xr)
    if (!(heteroskedastic & epsilon == 0)) {
        results <- svmRegress(Y = Yr, X = Xr, epsilon = epsilon,
                              lambda = lambda)
        alphaPlus <- results$auxiliaryDual1 - results$auxiliaryDual2
        alphaMinus <- results$auxiliaryDual1 + results$auxiliaryDual2
    }
    if (epsilon != 0) {
        if (!heteroskedastic && !wald) {
        rho <- mean(Yr <= Xrc %*% results$beta +
                    nullGamma * Zr - epsilon)
        } else {
            if (!wald) {
                testStatNumer <- t(Zr) %*% alphaPlus / sqrt(n)
                svmFit <- svmRegress(Yr, cbind(Xr, Zr), epsilon, lambda)
                uHat <- Y - (cbind(Xrc, Zr) %*% svmFit$beta - epsilon)
                hn <- hc * n ^ (-1 / 3)
                cn <- kappa * (qnorm(tau + hn) - qnorm(tau - hn))
                fX <- sweep(x = Xrc, MARGIN = 1,
                            STATS = as.vector(abs(uHat) < cn),
                            FUN = "*")
                fZX <- t(Zr) %*% fX / n
                fXX <- t(Xrc) %*% fX / n
                A <- fZX %*% solve(fXX)
                Ztilde <- Zr - Xrc %*% t(A)
                testStatDenom <- sqrt(2 * t(Ztilde * as.vector(uHat < 0)) %*%
                                      Ztilde / n)
                testStat <- testStatNumer / testStatDenom
                pvalue <- as.numeric(pnorm(-abs(testStat)) * 2)
                return(pvalue)
            } else {
                svmFit <- svmRegress(Y, cbind(Xr, Zr), epsilon, lambda)
                uHat <- Y - (cbind(Xrc, Zr) %*% svmFit$beta - epsilon)
                hn <- hc * n ^ (-1 / 3)
                cn <- kappa * (qnorm(tau + hn) - qnorm(tau - hn))
                fX <- sweep(x = cbind(Xrc, Zr), MARGIN = 1,
                            STATS = as.vector(abs(uHat) < cn),
                            FUN = "*")
                fXX <- t(cbind(Xrc, Zr)) %*% fX / n
                FXX <- t(sweep(cbind(Xrc, Zr), MARGIN = 1,
                               STATS = as.vector(uHat < 0), FUN = "*")) %*%
                    cbind(Xrc, Zr) / n
                avar <- 0.5 * solve(fXX) %*% FXX %*% solve(fXX)
                avar <- avar[ncol(avar), ncol(avar)]
                ci <- c(svmFit$beta[length(svmFit$beta)] - 1.96 *
                        sqrt(avar / n),
                        svmFit$beta[length(svmFit$beta)] + 1.96 *
                        sqrt(avar / n))
                pval <- 1 - pnorm(abs((sqrt(n) *
                                       svmFit$beta[length(svmFit$beta)]) /
                                      sqrt(avar)))
                return(list(ci = ci,
                            pval = pval))
            }
        }
    }
    if (epsilon == 0) {
        if (!heteroskedastic) {
            rho <- 0.5
        } else {
            n <- length(Yr)
            rqFit <- rq(Yr ~ Xr)
            bhat <- rqFit$dual - (1 - tau)
            testStatNumer <- t(bhat) %*% Zr / sqrt(n)
            rqFit <- rq(Y ~ Xr + Zr, tau = 0.5)
            uHat <- Y - cbind(Xrc, Zr) %*% rqFit$coefficients
            hn <- hc * n ^ (-1 / 3)
            cn <- kappa * (qnorm(tau + hn) - qnorm(tau - hn))
            fX <- sweep(x = Xrc, MARGIN = 1, STATS = as.vector(abs(uHat) < cn),
                     FUN = "*")
            fZX <- t(Zr) %*% fX / n
            fXX <- t(Xrc) %*% fX / n
            A <- solve(fXX) %*% t(fZX)
            testStatDenom <- sqrt(tau * (1 - tau) * t(Zr - Xrc %*% A) %*%
                                  (Zr - Xrc %*% A) / n)
            testStat <- testStatNumer / testStatDenom
            pvalue <- as.numeric(pnorm(-abs(testStat)) * 2)
            return(pvalue)
        }
    }
    M <- diag(rep(1, nrow(Xrc))) - Xrc %*% solve(t(Xrc) %*% Xrc) %*% t(Xrc)
    T <- (1 / sqrt(nrow(Xrc))) * t(Zr) %*% alphaPlus /
        sqrt((1 / nrow(Xrc)) * t(Zr) %*% M %*% Zr * rho * 2)
    pvalue <- as.numeric(pnorm(-abs(T)) * 2)
    return(pvalue)
}
