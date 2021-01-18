l1svr <- function(formula, data, epsilon, lambda, inference = TRUE,
                  confidence = FALSE,
                  heteroskedastic = FALSE,
                  h = NULL, kappa = NULL) {
    origCall <- match.call()
    ## Check that a response variable is included in the formula
    tmpTerms <- terms(formula)
    if (attr(tmpTerms, 'response') == 0) {
        stop('Dependent variable not declared in formula.',
             call. = FALSE)
    }
    ## Check whether there is intercept
    tmpInt <- as.logical(attr(tmpTerms, 'intercept'))
    ## Check that tuning parameters are provided if heteroskedasticity
    ## is declared.
    if (heteroskedastic) {
        if (is.null(h) | is.null(kappa)) {
            stop(gsub('\\s+', ' ',
                      "If 'heteroskedastic = TRUE', then arguments 'h'
                       and 'kappa' must be provided."),
                 call. = FALSE)
        }
        if (!is.numeric(h) | !(is.numeric(kappa) & kappa > 0)) {
            stop(gsub('\\s+', ' ',
                      "Arguments 'h' and 'kappa'
                       must be strictly positive scalars."),
                 call. = FALSE)
        }
    } else {
        if (!is.null(h) | !is.null(kappa)) {
            warning(gsub('\\s+', ' ',
                         "Arguments 'h' and 'kappa' are only used when
                          'heteroskedastic = TRUE'."),
                    call. = FALSE)
        }
    }
    if (!inference && (heteroskedastic == TRUE |
                       !is.null(h) | !is.null(kappa))) {
        warning(gsub('\\s+', ' ',
                     "Arguments 'heteroskedastic', 'h', and 'kappa' are
                      only used when 'inference = TRUE'."),
                call. = FALSE)
    }
    ## Construct design matrix
    tmpMf <- model.frame(formula = formula, data = data)
    tmpMt <- terms(x = formula, data = data, rhs = 1)
    fullX <- model.matrix(tmpMt, tmpMf)
    n <- nrow(fullX)
    Y <- as.vector(tmpMf[, all.vars(formula)[1]])
    origNames <- colnames(fullX)
    ## Remove vector of constants/intercept from design matrix
    if (tmpInt) {
        X <- matrix(fullX[, -1], ncol = ncol(fullX) - 1)
        colnames(X) <- origNames[-1]
    } else {
        X <- fullX
    }
    ## Estimate coefficients
    coefEstimates <- svmRegress(Y, X, epsilon, lambda, tmpInt)
    coefEstimates$N <- n
    coefEstimates$call <- origCall
    class(coefEstimates) <- 'l1svr'
    ## Generate output without inference
    if (!inference) {
        return(coefEstimates)
    } else {
        ## Conduct inference
        ## Obtain residuals
        U <- Y - fullX %*% coefEstimates$coef
        ## Concentrate out each variable one by one
        tmpXB <- sweep(x = fullX, MARGIN = 2,
                       STATS = rep(0, times = ncol(fullX)), FUN = '*')
        concY <- sweep(x = -tmpXB, MARGIN = 1, STATS = Y, FUN = '+')
        ## Conduct inference for each term
        pVec <- NULL
        for (i in 1:ncol(fullX)) {
            if (tmpInt) {
                if (i == 1) {
                    tmpP <- svmInference(Xr = X, Zr = rep(1, n),
                                         Yr = concY[, 1], Y = Y, U = U,
                                         epsilon = epsilon, lambda = lambda,
                                         intercept = FALSE, nullGamma = 0,
                                         heteroskedastic = heteroskedastic,
                                         hc = h, kappa = kappa)
                    pVec <- c(pVec, tmpP)
                } else {
                    Xr <- matrix(X[, -(i - 1)], ncol = ncol(X) - 1)
                    colnames(Xr) <- colnames(X)[-(i - 1)]
                    tmpP <- svmInference(Xr = Xr, Zr = X[, (i - 1)],
                                         Yr = concY[, i], Y = Y, U = U,
                                         epsilon = epsilon, lambda = lambda,
                                         intercept = TRUE, nullGamma = 0,
                                         heteroskedastic = heteroskedastic,
                                         hc = h, kappa = kappa)
                    pVec <- c(pVec, tmpP)
                }
            } else {
                Xr <- matrix(X[, -i], ncol = ncol(X) - 1)
                colnames(Xr) <- colnames(X)[-i]
                tmpP <- svmInference(Xr = Xr, Zr = X[, i],
                                     Yr = concY[, i], Y = Y, U = U,
                                     epsilon = epsilon, lambda = lambda,
                                     intercept = FALSE, nullGamma = 0,
                                     heteroskedastic = heteroskedastic,
                                     hc = h, kappa = kappa)
                pVec <- c(pVec, tmpP)
            }
        }
        names(pVec) <- origNames
        coefEstimates$pvalues <- pVec
        ## ---------------------------------------
        ## TESTING
        ## confidence <- TRUE
        alpha <- 0.05
        if (confidence) {
            confidenceMat <- NULL
            for (i in 1:ncol(fullX)) {
                args <- list(FUN = svmInference,
                             init = coefEstimates$coef[i],
                             target = alpha,
                             tol = 5e-3,
                             Y = Y, U = U,
                             epsilon = epsilon, lambda = lambda,
                             intercept = FALSE,
                             heteroskedastic = heteroskedastic,
                             hc = h, kappa = kappa)
                if (tmpInt) {
                    if (i == 1) {
                        args$intercept <- FALSE
                        args$Xr <- X
                        args$Zr <- rep(1, n)
                        args$Yr <- concY[, 1]
                    } else {
                        Xr <- matrix(X[, -(i - 1)], ncol = ncol(X) - 1)
                        colnames(Xr) <- colnames(X)[-(i - 1)]
                        args$intercept <- TRUE
                        args$Xr <- Xr
                        args$Zr = X[, (i - 1)]
                        args$Yr = concY[, i]
                    }
                } else {
                    Xr <- matrix(X[, -i], ncol = ncol(X) - 1)
                    colnames(Xr) <- colnames(X)[-i]
                    args$intercept <- FALSE
                    args$Xr <- Xr
                    args$Zr <- X[, i]
                    args$Yr <- concY[, i]
                }
                ## Find lower bound
                print('coefEstimate$coef[i]')
                print(coefEstimates$coef)
                print(coefEstimates$coef[i])
                args$left <- TRUE
                args$increment <- -abs(coefEstimates$coef[i])
                save(args, file = 'test-args.Rdata')
                stop('end of test')
                lb <- do.call(gridSearch, args)
                print('found lb')
                print(lb)
                ## Find upper bound
                args$left <- FALSE
                args$increment <- abs(coefEstimates$coef[i])
                ub <- do.call(gridSearch, args)
                print('found ub')
                print(ub)
                ## Save results
                confidenceMat <- rbind(confidenceMat, c(lb, ub))
            }
            print('this is the confidence mat')
            print(confidenceMat)
        }
        ## ---------------------------------------
        return(coefEstimates)
    }
}


#' SVR-l1 regression function
#'
#' This function carries out a SVR regression under l1 penalty.
#'
#' @param Y (vector) Dependent variable.
#' @param X (2d array) Independent variables, excluding column of
#'     constants for the intercept.
#' @param epsilon (real) Parameter defining soft-thresholding rule.
#' @param lambda (real) Tuning parameter to scale l1 penalty.
#' @param intercept (boolean) Set to \code{TRUE} if an intercept
#'     should be included in the regression.
#' @return (list) A list including the coefficient estimates, and
#'     solutions to the primal and dual problem that can be used to
#'     identify the support vectors and recover the coefficient
#'     estimates.
svmRegress <- function(Y, X, epsilon, lambda, intercept = TRUE) {
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
    if (!intercept) {
        AmatPrimal <- AmatPrimal[, -1]
        cvecPrimal <- cvecPrimal[-1]
        ubPrimal <- ubPrimal[-1]
        lbPrimal <- lbPrimal[-1]
    }
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
    tmpShift <- 0
    if (intercept) {
        tmpShift <- 1
        betaIntercept <- solutions[1]
    }
    betaMinus  <- solutions[(1 + tmpShift):(ncol(X) + tmpShift)]
    betaPlus <- solutions[(1 + tmpShift + ncol(X)):(tmpShift + 2 * ncol(X))]
    auxiliaryPrimal <- solutions[(2 * ncol(X) + 1 + tmpShift):length(solutions)]
    primalNegative <- auxiliaryPrimal[1:nrow(X)] ## xi minus
    primalPositive <- auxiliaryPrimal[(nrow(X) + 1):length(auxiliaryPrimal)] ## xi plus
    if (intercept) {
        beta <- c(betaIntercept, betaPlus - betaMinus)
        names(beta) <- c('(Intercept)', colnames(X))
    } else {
        beta <- betaPlus - betaMinus
        names(beta) <- colnames(X)
    }
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
    if (!intercept) {
        senseDual <- senseDual[-1]

    }
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
    dualNegative <- resultDual$x[1:nrow(X)] ## alpha minus
    dualPositive <- resultDual$x[(nrow(X) + 1):(2 * nrow(X))] ## alpha plus
    ## Return output
    return(list(coefficients = beta,
                primalNegative = primalNegative,
                primalPositive = primalPositive,
                dualNegative  = dualNegative,
                dualPositive  = dualPositive))
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
#' @param U Vector of the residuals.
#' @param epsilon Real scalar, bandwidth for SVR.
#' @param lambda Real scalar, tuning parmaeter for SVR that controls
#'     the penalty.
#' @param intercept (boolean) Set to \code{TRUE} if an intercept
#'     should be included in the regression.
#' @param nullGamma Real scalar, the coefficient on \code{Zr} under
#'     the null.
#' @param heteroskedastic Boolean, indicate whether data has
#'     heteroskedastic errors.
#' @param hc Real scalar, tuning parameter for density estimation.
#' @param kappa Real scalar, tuning parameter for density estimation.
#' @return pvalue.
svmInference <- function(Xr, Zr, Yr, Y, U, epsilon, lambda = 0,
                         intercept = TRUE, nullGamma = 0,
                         heteroskedastic = FALSE, hc = 0.5, kappa = 1) {
    tau <- 0.5
    n <- length(Yr)
    uHat <- U + epsilon
    if (intercept)  Xrc <- cbind(1, Xr)
    if (!intercept) Xrc <- Xr
    if (!(heteroskedastic & epsilon == 0)) { ## i.e., if not hetero.
                                             ## median regression
        concResults <- svmRegress(Y = Yr, X = Xr, epsilon = epsilon,
                                  lambda = lambda, intercept = intercept)
        alphaPlus <- concResults$dualNegative - concResults$dualPositive
    }
    if (epsilon != 0) {
        if (!heteroskedastic) {
            rho1 <- mean(Yr <= Xrc %*% concResults$coefficients +
                         nullGamma * Zr - epsilon)
            rho2 <- mean(Yr >= Xrc %*% concResults$coefficients +
                         nullGamma * Zr + epsilon)
            rho <- rho1 + rho2
        } else {
            testStatNumer <- t(Zr) %*% alphaPlus / sqrt(n)
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
        }
    }
    if (epsilon == 0) {
        if (!heteroskedastic) {
            rho <- 0.5
        } else {
            if (intercept)  rqFit <- quantreg::rq(Yr ~ 1 + Xr)
            if (!intercept) rqFit <- quantreg::rq(Yr ~ 0 + Xr)
            bhat <- rqFit$dual - (1 - tau)
            testStatNumer <- t(bhat) %*% Zr / sqrt(n)
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
        sqrt((1 / nrow(Xrc)) * t(Zr) %*% M %*% Zr * rho)
    pvalue <- as.numeric(pnorm(-abs(T)) * 2)
    return(pvalue)
}


#' Format result for display
#'
#' This function simply takes a number and formats it for being
#' displayed. Numbers less than 1 in absolute value are rounded to 6
#' significant figure. Numbers larger than
#'
#' @param x The scalar to be formated
#' @return A scalar.
fmtResult <- function(x) {
    if (abs(x) < 1) {
        fx <- signif(x, digits = 7)
    } else if (abs(x) >= 1 & abs(x) < 1e+7) {
        fx <- signif(round(x, digits = 4), digits = 7)
    } else {
        fx <- formatC(x, format = "e", digits = 7)
    }
    if (is.numeric(fx)) fx <- as.character(fx)
    return(fx)
}

#' Print results
#'
#' This function uses the print method on the l1svr return list.
#'
#' @param x an object returned from '\code{l1svr}'.
#' @param ... additional arguments.
#' @return basic set of results.
#' @export
print.l1svr <- function(x, ...) {
    stopifnot(inherits(x, "l1svr"))
    if (!is.null(x$coefficients)) {
        cat('\nCall:\n')
        print(x$call)
        cat('\nCoefficient estimates:\n')
        print(x$coefficients)
        cat('\n')
    }
}

#' Summarize results
#'
#' This function uses the summary method on the l1svr return list.
#'
#' @param object an object returned from '\code{l1svr}'.
#' @param ... additional arguments.
#' @return summarized results.
#' @export
summary.l1svr <- function(object, ...) {
    stopifnot(inherits(object, "l1svr"))
    if (!is.null(object$coefficients) && is.null(object$pvalues)) {
        cat('\nCall:\n')
        print(object$call)
        cat('\nCoefficient estimates:\n')
        print(object$coefficients)
        ## cat('\nObs: ', object$N, '\n')
    }
    if (!is.null(object$coefficients) && !is.null(object$pvalues)) {
        cat('\nCall:\n')
        print(object$call)
        cat('\nCoefficient estimates:\n')
        resultsMat <- rbind(object$coefficients, object$pvalues)
        resultsMat <- data.frame(t(resultsMat))
        signifVec <- rep('', times = length(object$pvalues))
        signifVec[which(object$pvalues < 1e-1)] <- '.  '
        signifVec[which(object$pvalues < 5e-2)] <- '*  '
        signifVec[which(object$pvalues < 1e-2)] <- '** '
        signifVec[which(object$pvalues < 1e-3)] <- '***'
        resultsMat$.tmp.sig.vec <- signifVec
        colnames(resultsMat) <- c('Estimate', 'Pr(>|t|)', '')
        print(resultsMat)
        cat('---\n')
        cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
        ## cat('Obs: ', object$N, '\n')
    }
}

#' Function to solve equations via grid search
#'
#' This function solves an equation using a grid search.
#'
#' @param FUN A function with a scalar argument.
#' @param init Scalar, the starting point for the grid search.
#' @param target Scalar, the target value for the function.
#' @param increment Scalar, starting increment for the grid
#'     search. This increment is progessively halved.
#' @param tol Scalar, tolerance to determine when the target value is
#'     reached.
#' @param left Boolean, indicate the direction of search. For example,
#'     \code{x^2 = 4} has two solutions, one being \code{-2} and the
#'     other being \code{2}. If \code{init = 0} and \code{left =
#'     TRUE}, then this function would return \code{-2}. If instead
#'     \code{left = FALSE}, then this function would return \code{2}.
#' @return Scalar, the argument value for which \code{FUN} is equal to
#'     \code{target}.
#'
gridSearch <- function(FUN, init = 0, target, increment = 2,
                       tol = 5e-3, left = FALSE, iter.max = 10,
                       ...) {
    exponent <- 0
    initVal <- FUN(init, ...)
    diff <- target - initVal
    if (left) diff <- -diff
    if (diff > 0) {
            tooSmall <- TRUE
            direction <- 1
    }
    if (diff < 0) {
            tooSmall <- FALSE
            direction <- -1
    }
    iter.count <- 1
    while (abs(diff) > tol & iter.count <= iter.max) {
        print('ORIGINAL init')
        print(init)
        print('stuff added')
        print(direction * sign(increment) * abs(increment) ^ exponent)
        init <- init + direction * sign(increment) * abs(increment) ^ exponent
        value <- FUN(init, ...)
        diff <- target - value
        if (left) diff <- -diff        
        if (diff > 0) {
            if (tooSmall == FALSE) {
                direction <- direction * -1
                exponent <- exponent - 1
            }
            tooSmall <- TRUE
        }
        if (diff < 0) {
            if (tooSmall == TRUE) {
                direction <- direction * -1
                exponent <- exponent - 1
            }
            tooSmall <- FALSE
        }
        print('init value')
        print(init)
        print('pvalue')
        print(value)
        iter.count <- iter.count + 1
    }
    if (abs(diff) > tol) {
        warning(gsub('\\s+', ' ',
                     'Calculation of confidence intervals terminated,
                      iteration maximum reached.'),
                call. = FALSE)
    }
    return(init)
}
