## Sum of powered score (SPU) test (perm, version 1, vector used in permutation, coded in C)
##
## It gives the p-values of the SPU test and aSPU test based on based on the permutation of residuals.  (This is version 1, matrix version is faster but if it doesn't work, we should use version 1, vector version, coded in C)
##
## @param Y phenotype data. It can be disease lables; =0 for controls, =1 for cases.
##     or It can be any quantitative traits. Vector with length n (number of observations)
##
## @param X genotype data; each row for a subject, and each column
##     for an SNP. The value of each element is the # of the copies
##     for an allele. Matrix with dimension n by k (n : number of observation, k : number of genotype data)
##
## @param cov covariates. Matrix with dimension n by p (n :number of observation, p : number of covariates)
##
## @param model Use "gaussian" for quantitative trait (Default)
##    , and Use "binomial" for binary trait.
##
## @param pow power used in SPU test. Vector of g number of power.
##
## @param n.perm number of permutation
##
## @export
## @return Test Statistics and p-values for SPU tests and aSPU test.
##
## @examples
##
## data(exdat)
## out <- aSPUperm2C(exdat$Y, exdat$X, cov = NULL, model = "binomial",
##                   pow = c(1:8, Inf), n.perm = 1000)
## out
##


aSPUperm2C <- function(Y, X, cov = NULL, model=c("gaussian","binomial"), pow=c(1:8, Inf), n.perm=1000){

    model = match.arg(model)
    n <- length(Y)
    if (is.null(X) && length(X) > 0) 
        X = as.matrix(X, ncol = 1)
    k <- ncol(X)
    if (is.null(cov)) {
        XUs <- Xg <- X
        r <- Y - mean(Y)
        U <- as.vector(t(Xg) %*% r)
    } else {
        tdat1 <- data.frame(trait = Y, cov)
        if (is.null(colnames(cov))) {
            colnames(tdat1) = c("trait", paste("cov", 1:dim(cov)[2], 
                sep = ""))
        } else {
            colnames(tdat1) = c("trait", colnames(cov))
        }
        fit1 <- glm(trait ~ ., family = model, data = tdat1)
        pis <- fitted.values(fit1)
        XUs <- matrix(0, nrow = n, ncol = k)
        Xmus = X
        for (i in 1:k) {
            tdat2 <- data.frame(X1 = X[, i], cov)
            fit2 <- glm(X1 ~ ., data = tdat2)
            Xmus[, i] <- fitted.values(fit2)
            XUs[, i] <- (X[, i] - Xmus[, i])
        }
        r <- Y - pis
        U <- t(XUs) %*% r
    }
    Ts = rep(NA, length(pow))
    for (j in 1:length(pow)) {
        if (pow[j] < Inf) {
            Ts[j] = sum(U^pow[j])
        } else Ts[j] = max(abs(U))
    }
    
    if(max(pow) == Inf) {
        pow[which(pow ==Inf)] = -1
    }

    s <- sample(1:10^5, 1)
    Results = aSPUpermEngine(t(XUs), r, pow, n.perm, Ts, s)

    minp0 <- Results$minp0
    pPerm0 <- Results$pPerm0
    P0s <- Results$P0

    Paspu <- (sum(minp0 <= min(pPerm0)) + 1)/(n.perm + 1)
    pvs <- c(pPerm0, Paspu)
    Ts <- c(Ts, min(pPerm0))

    if(min(pow) == -1) {
        pow[which(pow == -1 )] = Inf
    }

    names(Ts) <- c(paste("SPU", pow, sep = ""), "aSPU")
    names(pvs) = names(Ts)
    list(Ts = Ts, pvs = pvs)
}





