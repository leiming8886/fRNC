fitBumModel <- function (x, starts = 10)
{
    if (is.null(names(x))) {
        warning("Please name the p-values with the gene names!")
        names(x) = as.character(1:length(x))
    }
    fit <- bumOptim(x = x, starts)
    return(fit)
}

bumOptim <- function (x, starts = 1, labels = NULL)
{
    if (is.null(names(x)) && is.null(labels)) {
        warning("Please name the p-values with the gene names or give labels!")
        names(x) <- as.character(1:length(x))
    }
    if (!is.null(labels)) {
        names(x) <- labels
    }
    a <- runif(starts, 0.3, 0.7)
    lambda <- runif(starts, 0.3, 0.7)
    value <- Inf
    best <- list()
    for (i in 1:starts) {
        test.optim <- try(opt <- optim(c(lambda[i], a[i]), fn = .fbumnLL,
            gr = .fpLL, x = x, lower = rep(1e-05, 3), method = "L-BFGS-B",
            upper = rep(1 - 1e-05, 3)))
        if ((!class(test.optim) == "try-error") && all(opt$par >=
            1e-05) && all(opt$par <= 1 - 1e-05)) {
            value <- opt$value
            best <- opt
        }
		if(value > 0){
		break
		}
    }
    if (length(best) == 0) {
        return(warning("BUM model could not be fitted to data"))
    }
    else {
        if (any(opt$par == 1e-05) || any(opt$par == 1 - 1e-05)) {
            #warning("One or both parameters are on the limit of the defined parameter space")
        }
        ret <- list(lambda = best$par[1], a = best$par[2], negLL = best$value,
            pvalues = x)
        class(ret) <- "bum"
        return(ret)
    }
}
scoreFunction <-function (fb, fdr = 0.001){


    return((fb$a - 1) * (log(fb$pvalues) - log(fdrThreshold(fdr,
                                                            fb))))
}
fdrThreshold <- function (fdr, fb){
    pihat <- fb$lambda + (1 - fb$lambda) * fb$a
    return(((pihat - (fdr * fb$lambda))/(fdr * (1 - fb$lambda)))^(1/(fb$a -
                                                                         1)))
}

.fbumnLL <- function (parms, x){
    sum(log(fbum(x, parms[1], parms[2])))
}

fbum <- function (x, lambda, a)
{
    lambda + (1 - lambda) * a * x^(a - 1)
}
.fpLL <- function (parms, x) {
    l <- parms[1]
    a <- parms[2]
    dl <- -sum((1 - a * x^(a - 1))/(a * (1 - l) * x^(a - 1) +
                                        l))
    da <- -sum((a * (1 - l) * x^(a - 1) * log(x) + (1 - l) *
                    x^(a - 1))/(a * (1 - l) * x^(a - 1) + l))
    return(c(dl, da))
}




