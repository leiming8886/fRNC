#' @title Integrate multiple p-values into joint p-values
#'
#' @description The function integrate multiple p-values into the joint p-value of p-values based on the order statistics of p-values.
#' An joint p-value #is given by the kth order statistic.
#' @param pvalmatrix Numeric matrix of p-values, columns represent different sets of p-values
#' @return matrix The matrix of two columns: gene, weight(p-value)
#' @export
integPvals <- function (pvalmatrix)
{
    old.pvals <- pvalmatrix
	order <- ncol(old.pvals)
    if (!is.matrix(pvalmatrix)) {
        return("Input is not matrix")
    }
    if (order > dim(pvalmatrix)[2]) {
        return("order is larger than array dimensions")
    }
    for (j in 1:dim(pvalmatrix)[1]) {
        pvalmatrix[j, ] <- sort(as.numeric(pvalmatrix[j, ]),
            na.last = TRUE)
    }
    x.vec <- as.numeric(pvalmatrix[, order])
    n <- dim(pvalmatrix)[2]
    concat.pvals <- pbeta(x.vec, order, n - order + 1)
    names(concat.pvals) <- rownames(pvalmatrix)
    gene2weight <- data.frame(gene = rownames(pvalmatrix),weigth = concat.pvals)
    colnames(gene2weight) <- c("gene","weight")
    return(gene2weight)
}
