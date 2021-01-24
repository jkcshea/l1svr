#' Function to perform SVM regression under l1 regularization
#'
#' @param formula Regression formula.
#' @param Y (vector) Dependent variable.
#' @param X (2d array) Independent variables, excluding column of
#'     constants for the intercept.
#' @param epsilon (real) Parameter defining soft-thresholding rule.
#' @param lambda (real) Tuning parameter to scale l1 penalty.
#' @param intercept (boolean) Set to \code{TRUE} if an intercept
#'     should be included in the regression.
#' @param solver character, name of the linear programming package in
#'     R used to obtain the bounds on the treatment effect. The
#'     function supports \code{'gurobi'}, \code{'cplexapi'},
#'     \code{'lpsolveapi'}. The name of the solver should be provided
#'     with quotation marks.
l1svr <- function(formula, data, epsilon, lambda,
                  h = NULL, kappa = NULL,
                  solver = 'gurobi',
                  inference = TRUE,
                  confidence.level = NULL,
                  confidence.iter = 25,
                  confidence.tol = 1e-3,
                  confidence.same = 1e-06) {
    origCall <- match.call()
    ## Check that the solver is valid
    solver <- tolower(solver)
    if (! solver %in% c("gurobi",
                        "cplexapi",
                        "lpsolveapi")) {
        stop(gsub("\\s+", " ",
                  paste0("Estimator is incompatible with linear
                          programming package '", solver,"'. Please
                          install one of the
                          following linear programming packages instead:
                          gurobi (version 7.5-1 or later);
                          cplexAPI (version 1.3.3 or later);
                          lpSolveAPI (version 5.5.2.0 or later).")),
             call. = FALSE)
    }
    ## Check that the user actually has the solver
    missingSolver <- FALSE
    if (solver == 'gurobi') {
        if (!requireNamespace("gurobi", quietly = TRUE)) {
            missingSolver <- TRUE
        }
    }
    if (solver == 'lpsolveapi') {
        if (!requireNamespace("lpSolveAPI", quietly = TRUE)) {
            missingSolver <- TRUE
        }
    }
    if (solver == 'cplexapi') {
        if (!requireNamespace("cplexAPI", quietly = TRUE)) {
            missingSolver <- TRUE
        }
    }
    if (missingSolver) {
        stop(gsub("\\s+", " ",
                  paste0("The solver '", solver, "' is not installed.
                  Please install one of the following packages required
                  for estimation:
                  gurobi (version 7.5-1 or later);
                  cplexAPI (version 1.3.3 or later);
                  lpSolveAPI (version 5.5.2.0 or later).")),
             call. = FALSE)
    }
    ## Check that a response variable is included in the formula
    tmpTerms <- terms(formula)
    if (attr(tmpTerms, 'response') == 0) {
        stop('Dependent variable not declared in formula.',
             call. = FALSE)
    }
    ## Check whether there is intercept
    tmpInt <- as.logical(attr(tmpTerms, 'intercept'))
    ## Determine whether confidence intervals should be constructed
    if (!is.null(confidence.level)) {
        confidence <- TRUE
        if (!is.numeric(confidence.level)) {
            stop(gsub('\\s+', ' ',
                      "The argument 'confidence.level' must be numeric and
                       lie between 0 and 1."),
                 call. = FALSE)
        }
        if (confidence.level <= 0 | confidence.level >= 1) {
            stop(gsub('\\s+', ' ',
                      "The argument 'confidence.level' must
                       lie between 0 and 1."),
                 call. = FALSE)
        }
    } else {
        confidence <- FALSE
    }
    if (!confidence) {
        if (hasArg(confidence.iter) |
            hasArg(confidence.tol) |
            hasArg(confidence.same)) {
            warning(gsub('\\s+', ' ',
                         "To enable estimation of confidence intervals,
                          set 'confidence.level' to a value between 0 and 1.
                          Otherwise, the arguments 'confidence.iter',
                          'confidence.tol', and 'confidence.same' are
                          ignored."),
                    call. = FALSE)
        }
    }
    ## Determine heteroskedasticity. Heteroskedasticity is assumed as
    ## long as the user inputs an argument into either h or kappa.
    if (!is.null(h) | !is.null(kappa)) {
        heteroskedastic <- TRUE
    } else {
        heteroskedastic <- FALSE
    }
    if (heteroskedastic) {
        if (is.null(h) | is.null(kappa)) {
            stop(gsub('\\s+', ' ',
                      "If the errors are homoskedastic, then both arguments
                       'h' and 'kappa' must not be provided.
                       If the errors are heteroskedastic, then
                       both arguments 'h' and 'kappa' must be provided."),
                 call. = FALSE)
        }
        if (!(is.numeric(h) && h > 0) | !(is.numeric(kappa) && kappa > 0)) {
            stop(gsub('\\s+', ' ',
                      "Arguments 'h' and 'kappa'
                       must be strictly positive scalars."),
                 call. = FALSE)
        }
    }
    if (!inference && (heteroskedastic == TRUE |
                       confidence == TRUE)) {
        warning(gsub('\\s+', ' ',
                     "Arguments 'h', 'kappa', and those beginning with
                      'confidence.' are all ignored if
                      'inference = FALSE'."),
                call. = FALSE)
        heteroskedastic <- FALSE
        confidence <- FALSE
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
    coefEstimates <- svmRegress(Y, X, epsilon, lambda, tmpInt, solver)
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
        ## Conduct inference for each term
        pVec <- NULL
        for (i in 1:ncol(fullX)) {
            if (tmpInt) {
                if (i == 1) {
                    tmpP <- svmInference(Xr = X, Zr = rep(1, n), Y = Y, U = U,
                                         epsilon = epsilon, lambda = lambda,
                                         solver = solver,
                                         intercept = FALSE,
                                         nullGamma = 0,
                                         heteroskedastic = heteroskedastic,
                                         hc = h, kappa = kappa)
                    pVec <- c(pVec, tmpP)
                } else {
                    Xr <- matrix(X[, -(i - 1)], ncol = ncol(X) - 1)
                    colnames(Xr) <- colnames(X)[-(i - 1)]
                    tmpP <- svmInference(Xr = Xr, Zr = X[, (i - 1)],
                                         Y = Y, U = U,
                                         epsilon = epsilon, lambda = lambda,
                                         solver = solver,
                                         intercept = TRUE, nullGamma = 0,
                                         heteroskedastic = heteroskedastic,
                                         hc = h, kappa = kappa)
                    pVec <- c(pVec, tmpP)
                }
            } else {
                Xr <- matrix(X[, -i], ncol = ncol(X) - 1)
                colnames(Xr) <- colnames(X)[-i]
                tmpP <- svmInference(Xr = Xr, Zr = X[, i],
                                     Y = Y, U = U,
                                     epsilon = epsilon, lambda = lambda,
                                     solver = solver,
                                     intercept = FALSE, nullGamma = 0,
                                     heteroskedastic = heteroskedastic,
                                     hc = h, kappa = kappa)
                pVec <- c(pVec, tmpP)
            }
        }
        names(pVec) <- origNames
        coefEstimates$pvalues <- pVec
        if (confidence) {
            confidenceMatLb <- NULL
            confidenceMatUb <- NULL
            ## ## Below is code for an alternative method for
            ## ## determining the confidence intervals.
            ## if (optimAlt) {
            ## args <- list(f = ciIterOptim,
            ##                  target = 1 - confidence.level,
            ##                  intercept = tmpInt, Y = Y,
            ##                  X = X, U = U, epsilon = epsilon,
            ##                  lambda = lambda,
            ##                  heteroskedastic = heteroskedastic, h = h,
            ##              kappa = kappa)
            ## }
            args <- list(FUN = ciIter,
                         target = 1 - confidence.level,
                         iter.max = confidence.iter,
                         tol = confidence.tol,
                         tol.same = confidence.same,
                         intercept = tmpInt, Y = Y,
                         X = X, U = U, epsilon = epsilon,
                         lambda = lambda,
                         heteroskedastic = heteroskedastic, h = h,
                         kappa = kappa)
            cat('\n')
            for (i in 1:ncol(fullX)) {
                cat(paste0('Estimating confidence interval for ',
                    colnames(fullX)[i], '...\n'))
                ## ## Below is code for an alternative method for
                ## ## determining the confidence intervals.
                ## if (optimAlt) {
                ##     args$index <- i
                ##     ## Lower bound
                ##     t0 <- Sys.time()
                ##     args$interval <-
                ##         c(coefEstimates$coef[i] - 5 * abs(coefEstimates$coef[i]),
                ##           coefEstimates$coef[i])
                ##     lb <- do.call(optimize, args)
                ##     ## Upper bound
                ##     t1 <- Sys.time()
                ##     args$interval <-
                ##         c(coefEstimates$coef[i],
                ##           coefEstimates$coef[i] + 5 * abs(coefEstimates$coef[i]))
                ##     ub <- do.call(optimize, args)
                ## }
                ## Find lower bound
                args$index <- i
                args$init <- coefEstimates$coef[i]
                args$left <- TRUE
                t0 <- Sys.time()
                args$increment <- abs(coefEstimates$coef[i])
                lb <- do.call(gridSearch, args)
                ## Find the upper bound
                t1 <- Sys.time()
                args$left <- FALSE
                ub <- do.call(gridSearch, args)
                ## Store results
                confidenceMatLb <- rbind(confidenceMatLb, lb)
                confidenceMatUb <- rbind(confidenceMatUb, ub)
            }
            cat('\n')
            confidenceMatLb <- data.frame(confidenceMatLb)
            confidenceMatUb <- data.frame(confidenceMatUb)
            confidenceMatLb$status <- factor(confidenceMatLb$status,
                                              levels = c(1, 2, 3),
                                              labels = c('Optimal',
                                                         'Iter. limit',
                                                         'Prec. limit'))
            confidenceMatUb$status <- factor(confidenceMatUb$status,
                                              levels = c(1, 2, 3),
                                              labels = c('Optimal',
                                                         'Iter. limit',
                                                         'Prec. limit'))
            rownames(confidenceMatLb) <- colnames(fullX)
            rownames(confidenceMatUb) <- colnames(fullX)
            coefEstimates$ci <- list(lower = confidenceMatLb,
                                     upper = confidenceMatUb,
                                     level = confidence.level,
                                     tol = confidence.tol,
                                     heteroskedastic = heteroskedastic)
            if ('Iter. limit' %in% confidenceMatLb$status |
                'Iter. limit' %in% confidenceMatUb$status) {
                warning(gsub('\\s+', ' ',
                             'Estimation for some confidence intervals was
                              terminated,
                              iteration maximum reached.'),
                        call. = FALSE)
            }
            if ('Prec. limit' %in% confidenceMatLb$status |
                'Prec. limit' %in% confidenceMatUb$status) {
                warning(gsub('\\s+', ' ',
                             'Estimation for some confidence intervals was
                              terminated,
                              precision limit reached.'),
                        call. = FALSE)
            }
        }
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
#' @param solver character, name of the linear programming package in
#'     R used to obtain the bounds on the treatment effect. The
#'     function supports \code{'gurobi'}, \code{'cplexapi'},
#'     \code{'lpsolveapi'}. The name of the solver should be provided
#'     with quotation marks.
#' @return (list) A list including the coefficient estimates, and
#'     solutions to the primal and dual problem that can be used to
#'     identify the support vectors and recover the coefficient
#'     estimates.
svmRegress <- function(Y, X, epsilon, lambda, intercept = TRUE, solver) {
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
    if (solver == 'gurobi') {
        resultPrimal <- runGurobi(modelPrimal)
    } else {
        modelPrimal <- lpSetupSolver(modelPrimal, solver)
        if (solver == 'lpsolveapi') {
            resultPrimal <- runLpSolveAPI(modelPrimal, 'min')
        }
        if (solver == 'cplexapi') {
            resultPrimal <- runCplexAPI(modelPrimal, cplexAPI::CPX_MIN)
        }
    }
    ## Prepare primal solutions
    solutionPrimal <- resultPrimal$optx
    tmpShift <- 0
    if (intercept) {
        tmpShift <- 1
        betaIntercept <- solutionPrimal[1]
    }
    betaMinus  <- solutionPrimal[(1 + tmpShift):(ncol(X) + tmpShift)]
    betaPlus <- solutionPrimal[(1 + tmpShift + ncol(X)):(tmpShift + 2 * ncol(X))]
    auxiliaryPrimal <- solutionPrimal[(2 * ncol(X) + 1 + tmpShift):length(solutionPrimal)]
    primalNegative <- auxiliaryPrimal[1:nrow(X)] ## xi minus
    primalPositive <- auxiliaryPrimal[(nrow(X) + 1):length(auxiliaryPrimal)] ## xi plus
    if (intercept) {
        beta <- c(betaIntercept, betaPlus - betaMinus)
        names(beta) <- c('(Intercept)', colnames(X))
    } else {
        beta <- betaPlus - betaMinus
        names(beta) <- colnames(X)
    }
    ## Prepare dual solutions
    if (solver != 'lpsolveapi') {
        solutionDual <- resultPrimal$opty
        dualNegative <- solutionDual[1:nrow(X)] ## alpha minus
        dualPositive <- solutionDual[(nrow(X) + 1):(2 * nrow(X))] ## alpha plus
    } else {
        ## lpsolveAPI has rather odd dual solutions, and the
        ## documentation is not good enough to determine why. To avoid
        ## inconsistency, manually declare and solve the dual problem.
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
        resultDual <- runLpSolveAPI(modelDual, 'min')
        ## Prepare solutions
        dualPositive <- resultDual$optx[1:nrow(X)] ## alpha plus
        dualNegative <- resultDual$optx[(nrow(X) + 1):(2 * nrow(X))] ## alpha minus
    }
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
svmInference <- function(Xr, Zr, Y, U, epsilon, lambda = 0,
                         solver = 'gurobi',
                         intercept = TRUE, nullGamma = 0,
                         heteroskedastic = FALSE, hc = 0.5, kappa = 1) {
    Yr <- Y - Zr * nullGamma
    tau <- 0.5
    n <- length(Yr)
    uHat <- U + epsilon
    if (intercept)  Xrc <- cbind(1, Xr)
    if (!intercept) Xrc <- Xr
    if (!(heteroskedastic & epsilon == 0)) { ## i.e., if not hetero.
                                             ## median regression
        concResults <- svmRegress(Y = Yr, X = Xr, epsilon = epsilon,
                                  lambda = lambda,
                                  solver = solver,
                                  intercept = intercept)
        alphaPlus <- concResults$dualNegative - concResults$dualPositive
    }
    if (epsilon != 0) {
        if (!heteroskedastic) {
            rho1 <- mean(Yr <= (Xrc %*% concResults$coefficients +
                         nullGamma * Zr - epsilon))
            rho2 <- mean(Yr >= (Xrc %*% concResults$coefficients +
                         nullGamma * Zr + epsilon))
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
            rho <- 1
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
    if (!is.null(object$ci)) {
        ciTable <- data.frame(cbind(object$ci$lower$bound,
                                    object$ci$upper$bound))
        statusLower <- rep('  ', nrow(object$ci$lower))
        statusLower[which(object$ci$lower$status == 'Iter. limit')] <- ' *'
        statusLower[which(object$ci$lower$status == 'Prec. limit')] <- '**'
        statusUpper <- rep('  ', nrow(object$ci$upper))
        statusUpper[which(object$ci$upper$status == 'Iter. limit')] <- '* '
        statusUpper[which(object$ci$upper$status == 'Prec. limit')] <- '**'
        ciTable <- cbind(statusLower, ciTable, statusUpper)
        colnames(ciTable) <- c('', 'Lower', 'Upper', '')
        rownames(ciTable) <- rownames(object$ci$lower)
        cat(paste0(round(object$ci$level * 100, digits = 2),
                   '% confidence intervals:\n'))
        print(ciTable)
        cat('---\n')
        cat("Conf. interval codes:\n")
        cat("'*'  Iteration limit reached\n")
        cat("'**' Precision limit reached\n")
    }
}

ciIter <- function(gamma0, index, intercept, Y, X, U, epsilon, lambda,
               heteroskedastic, h, kappa) {
    args <- list(nullGamma = gamma0, Y = Y, U = U,
                 epsilon = epsilon, lambda = lambda,
                 heteroskedastic = heteroskedastic,
                 hc = h, kappa = kappa)
    if (intercept) {
        if (index == 1) {
            args$intercept <- FALSE
            args$Xr <- X
            args$Zr <- rep(1, n)
        } else {
            Xr <- matrix(X[, -(index - 1)], ncol = ncol(X) - 1)
            colnames(Xr) <- colnames(X)[-(index - 1)]
            args$intercept <- TRUE
            args$Xr <- Xr
            args$Zr = X[, (index - 1)]
        }
    } else {
        Xr <- matrix(X[, -index], ncol = ncol(X) - 1)
        colnames(Xr) <- colnames(X)[-index]
        args$intercept <- FALSE
        args$Xr <- Xr
        args$Zr <- X[, index]
    }
    do.call(svmInference, args)
}


gridSearch <- function(FUN, init = 0, target, increment = 2,
                       tol = 2.5e-3, tol.same = 1e-06,
                       left = FALSE, iter.max = 25,
                       ...) {
    exponent <- 0
    if (left) {
        farDirection <- -1
    }
    if (!left) {
        farDirection <- +1
    }
    x <- init + farDirection * abs(increment)
    value <- FUN(x, ...)
    diff <- target - value
    if (diff < 0) { ## i.e. haven't gone far enough
        tooFar <- FALSE
        dirMod <- 1
    }
    if (diff > 0) { ## i.e. gone too far
        tooFar <- TRUE
        dirMod <- -1
        exponent <- exponent - 1
    }
    iter.count <- 1
    sameX <- 0
    sameX.lim <- 3
    minDiff <- Inf
    while (abs(diff) > tol & iter.count <= iter.max & sameX < sameX.lim) {
        newX <- x + dirMod * farDirection * abs(increment) * 2.01 ^ exponent
        if (abs(newX - x) < tol.same) sameX <- sameX + 1
        x <- newX
        if (left) {
            if (x > init) x <- init
        } else {
            if (x < init) x <- init
        }
        value <- FUN(x, ...)
        diff <- target - value
        if (diff <= 0) { ## i.e. haven't gone far enough from initial point
            if (tooFar) {
                dirMod <- 1
                exponent <- exponent - 1
            }
            tooFar <- FALSE
        }
        if (diff >= 0) { ## i.e. gone too far from initial point
            if (!tooFar) {
                dirMod <- -1
                exponent <- exponent - 1
            }
            tooFar <- TRUE
        }
        if (abs(diff) < minDiff) {
            minDiff <- abs(diff)
            bestX <- x
            bestPvalue <- value
        }
        iter.count <- iter.count + 1
    }
    iter.count <- iter.count - 1
    status <- 1
    if (iter.count == iter.max && abs(diff) > tol) {
        ## warning(gsub('\\s+', ' ',
        ##              'Calculation of confidence intervals terminated,
        ##               iteration maximum reached.'),
        ##         call. = FALSE)
        status <- 2
    }
    if (sameX == sameX.lim && abs(diff) > tol) {
        ## warning(gsub('\\s+', ' ',
        ##              'Calculation of confidence intervals terminated,
        ##               precision limit reached.'),
        ##         call. = FALSE)
        status <- 3
    }
    return(c(bound = unname(bestX),
             pvalue = bestPvalue,
             iters = iter.count,
             status = status))
}

#' Alternative function for determining confidence intervals
#'
#' This function should be used with the R function \code{optimize} to
#' find the bounds of the confidence interval.
#'
#' @param target Target, equal to the size of the test.
#' @param gamma0 Scalar, a conjectured value of the coefficient for
#'     the covariate of interest, i.e. the covariate for which the
#'     confidence interval is being constructed.
#' @param index Integer, indexes which covariate is of interest.
#' @param intercept Boolean, indicates whether an intercept term is
#'     included in the regression.
#' @param Y Vector of dependent variable.
#' @param X Matrix of covariates.
#' @param U Vector of residuals from regressing \code{Y} on \code{X}
#'     using l1svr.
#'
#' @inheritParams ivmteEstimate
ciIterOptim <- function(target, gamma0, index, intercept, Y, X, U,
                        epsilon, lambda, heteroskedastic, h, kappa) {
    args <- list(nullGamma = gamma0, Y = Y, U = U,
                 epsilon = epsilon, lambda = lambda,
                 heteroskedastic = heteroskedastic,
                 hc = h, kappa = kappa)
    if (intercept) {
        if (index == 1) {
            args$intercept <- FALSE
            args$Xr <- X
            args$Zr <- rep(1, n)
        } else {
            Xr <- matrix(X[, -(index - 1)], ncol = ncol(X) - 1)
            colnames(Xr) <- colnames(X)[-(index - 1)]
            args$intercept <- TRUE
            args$Xr <- Xr
            args$Zr = X[, (index - 1)]
        }
    } else {
        Xr <- matrix(X[, -index], ncol = ncol(X) - 1)
        colnames(Xr) <- colnames(X)[-index]
        args$intercept <- FALSE
        args$Xr <- Xr
        args$Zr <- X[, index]
    }
    pvalue <- do.call(svmInference, args)
    return(abs(pvalue - target))
}


#' Running Gurobi LP solver
#'
#' This function solves the LP problem using the Gurobi package. The
#' object generated by \code{\link{lpSetup}} is compatible with the
#' \code{gurobi} function. See \code{\link{runCplexAPI}} for
#' additional error code labels.
#' @param lpobj list of matrices and vectors defining the linear
#'     programming problem.
#' @param solver.options list, each item of the list should
#'     correspond to an option specific to the LP solver selected.
#' @return a list of the output from Gurobi. This includes the
#'     objective value, the solution vector, and the optimization
#'     status (status of \code{1} indicates successful optimization) .
runGurobi <- function(lpobj, solver.options = list(outputflag = 0)) {
    result <- gurobi::gurobi(lpobj, solver.options)
    status <- 0
    if (result$status == "OPTIMAL") status <- 1
    if (result$status == "INFEASIBLE") status <- 2
    if (result$status == "INF_OR_UNBD") status <- 3
    if (result$status == "UNBOUNDED") status <- 4
    if (result$status == "NUMERIC") status <- 5
    if (result$status == "SUBOPTIMAL") status <- 6
    optx <- result$x
    return(list(objval = result$objval,
                optx = result$x,
                opty = result$pi,
                status = status))
}


#' Running cplexAPI LP solver
#'
#' This function solves the LP problem using the cplexAPI package. The
#' object generated by \code{\link{lpSetup}} is not compatible with
#' the \code{cplexAPI} functions. This function adapts the object to
#' solve the LP problem. See \code{\link{runGurobi}} for additional
#' error code labels.
#' @param lpobj list of matrices and vectors defining the linear
#'     programming problem.
#' @param lpdir input either CPX_MAX or CPX_MIN, which sets the LP
#'     problem as a maximization or minimization problem.
#' @param solver.options list, each item of the list should
#'     correspond to an option specific to the LP solver selected.
#' @return a list of the output from CPLEX. This includes the
#'     objective value, the solution vector, and the optimization
#'     status (status of \code{1} indicates successful optimization).
runCplexAPI <- function(lpobj, lpdir, solver.options = NULL) {
    ## Declare environment and set options
    env  <- cplexAPI::openEnvCPLEX()
    prob <- cplexAPI::initProbCPLEX(env)
    cplexAPI::chgProbNameCPLEX(env, prob, "sample")
    if (!is.null(solver.options)) {
        for(i in seq(length(solver.options))) {
            eval(parse(text = solver.options[[i]]))
        }
    }
    ## Declare LP prblem
    sense <- lpobj$sense
    cnt <- apply(lpobj$A, MARGIN = 2, function(x) length(which(x != 0)))
    beg <- rep(0, ncol(lpobj$A))
    beg[-1] <- cumsum(cnt[-length(cnt)])
    ind <- unlist(apply(lpobj$A, MARGIN = 2, function(x) which(x != 0) - 1))
    val <- c(lpobj$A)
    val <- val[val != 0]
    cplexAPI::copyLpwNamesCPLEX(env = env,
                                lp = prob,
                                nCols = ncol(lpobj$A),
                                nRows = nrow(lpobj$A),
                                lpdir = lpdir,
                                objf = lpobj$obj,
                                rhs = lpobj$rhs,
                                sense = sense,
                                matbeg = beg,
                                matcnt = cnt,
                                matind = ind,
                                matval = val,
                                lb = lpobj$lb,
                                ub = lpobj$ub)
    cplexAPI::lpoptCPLEX(env, prob)
    solution <- cplexAPI::solutionCPLEX(env, prob)
    cplexAPI::delProbCPLEX(env, prob)
    cplexAPI::closeEnvCPLEX(env)
    status <- 0
    if (typeof(solution) == "S4") {
        if (attr(solution, "class") == "cplexError") {
            status <- 5
            solution <- list()
            solution$objval <- NA
            solution$x <- NA
        }
    }  else {
        if (solution$lpstat == 1) status <- 1
        if (solution$lpstat == 2) status <- 4
        if (solution$lpstat == 3) status <- 2
        if (solution$lpstat == 4) status <- 3
        if (solution$lpstat == 5) status <- 7
        if (solution$lpstat == 6) status <- 6
    }
    cplexDual <- solution$pi
    save(cplexDual, file = 'cplexDual.Rdata')
    return(list(objval = solution$objval,
                optx   = solution$x,
                opty   = solution$pi,
                status = status))
}

#' Running lpSolveAPI
#'
#' This function solves the LP problem using the \code{lpSolveAPI}
#' package. The object generated by \code{\link{lpSetup}} is not
#' compatible with the \code{lpSolveAPI} functions. This function
#' adapts the object to solve the LP problem. See
#' \code{\link{runGurobi}} and \code{\link{runCplexAPI}} for
#' additional error code labels.
#' @param lpobj list of matrices and vectors defining the linear
#'     programming problem.
#' @param modelsense input either 'max' or 'min' which sets the LP
#'     problem as a maximization or minimization problem.
#' @param solver.options list, each item of the list should
#'     correspond to an option specific to the LP solver selected.
#' @return a list of the output from \code{lpSolveAPI}. This includes
#'     the objective value, the solution vector, and the optimization
#'     status (status of \code{1} indicates successful optimization).
runLpSolveAPI <- function(lpobj, modelsense, solver.options = NULL) {
    lpmodel <- lpSolveAPI::make.lp(nrow(lpobj$A), ncol(lpobj$A))
    for (j in 1:ncol(lpobj$A)) {
        lpSolveAPI::set.column(lprec = lpmodel,
                               column = j,
                               x = lpobj$A[, j])
    }
    lpSolveAPI::set.constr.value(lprec = lpmodel,
                                 rhs = lpobj$rhs)
    sense <- lpobj$sense
    sense[sense == "<"]  <- "<="
    sense[sense == ">"]  <- ">="
    sense[sense == "=="] <- "="
    lpSolveAPI::set.constr.type(lprec = lpmodel,
                                types = sense)
    lpSolveAPI::set.objfn(lprec = lpmodel,
                          obj = lpobj$obj)
    lpSolveAPI::lp.control(lprec = lpmodel,
                           sense = modelsense)
    if (!is.null(solver.options)) {
        eval(solver.options)
    }
    lpSolveAPI::set.bounds(lprec = lpmodel,
                           lower = lpobj$lb,
                           upper = lpobj$ub)
    solved <- lpSolveAPI::solve.lpExtPtr(lpmodel)
    status <- 0
    if (solved == 0) status <- 1
    if (solved == 1) status <- 6
    if (solved == 2) status <- 2
    if (solved == 3) status <- 4
    if (solved == 5) status <- 5
    ## Remove extraneous dual output from lpSolveAPI
    lpDual <- lpSolveAPI::get.dual.solution(lpmodel)
    lpDual <- lpDual[2:(length(lpSolveAPI::get.variables(lpmodel)) + 1)]
    save(lpDual, file = 'lpDual.Rdata')
    return(list(objval = lpSolveAPI::get.objective(lpmodel),
                optx   = lpSolveAPI::get.variables(lpmodel),
                opty   = lpDual,
                status = status))
}

#' Configure LP environment to be compatible with solvers
#'
#' This alters the LP object so the model will be compatible with
#' specific solvers.
#' @param env List, the LP object.
#' @param solver Character, the LP solver.
#' @return Nothing, as this modifies an environment variable to save
#'     memory.
#' @export
lpSetupSolver <- function(model, solver) {
    if (solver == "cplexapi") {
        model$sense[model$sense == "<"]  <- "L"
        model$sense[model$sense == "<="] <- "L"
        model$sense[model$sense == ">"]  <- "G"
        model$sense[model$sense == ">="] <- "G"
        model$sense[model$sense == "="]  <- "E"
        model$sense[model$sense == "=="] <- "E"
        model$ub[model$ub == Inf] <- cplexAPI::CPX_INFBOUND
        model$lb[model$lb == -Inf] <- -cplexAPI::CPX_INFBOUND
    }
    if (solver == "lpsolveapi") {
        model$sense[model$sense == "<"]  <- "<="
        model$sense[model$sense == ">"]  <- ">="
        model$sense[model$sense == "=="] <- "="
    }
    return(model)
}
