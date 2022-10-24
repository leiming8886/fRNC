#' @title Run module search function
#'
#' @description runmodule constructs a node-weighted ncRNA network, performs module searching, generates simulation data from random networks,
#' normalizes module scores using simulation data, removes un-qualified modules, and orders resultant modules according to their significance.
#'
#' @param network A data frame containing a symbolic edge list of the ncRNA network in which the columns must contain "node_gene_ID", "type", "target_gene_Name"
#'
#' @param gene2weight  A weigth data frame containing three columns:"type","gene", "weight" 
#' the first "type" the type of the gene identifier; lncRNA, miRNA, circRNA and RBP
#' the second gene is unique, gene identifier (should be coordinate with the node symbol used in ncRNA network);
#'  the third weight is gene-based  p-value or corrected p-value derived from differentially gene analysis or survival analysis
#' @param method	a character string indicating which the search method is to be computed
#'  . One of "global" (default, refer to Heinz method), "local ( refer to GS method)": can be abbreviated
#' @param expr1	the expression matrix of the case sample
#' @param expr2	the expression matrix of the control sample
#' @param  min.size An integer: the min numbel of size of the module for user settings in the method of "global", default 5.
#' @param d An integer used to define the order of neighbour genes to be
#'   searched  in  the method of  the method "local" . This parameter is default set up as 2
#' @param r A float indicating the cut-off for increment during module expanding
#'   process in  the method of  the method "local". Greater r will generate smaller module. Default is 0.1.
#' @param  maxsize An integer: the max numbel of size of the module for user settings in the method of "global", default 15.
#' @param FDR Numeric value, from the false discovery rate a p-value threshold is calculated. P-values below this threshold are considered to be significant
#'  The FDR can be used to control the size of the maximum scoring module
#' @param seletN a vector: gene identifier IDs, or a gene identifier ID, for example "MIMAT0000461",
#' or c("MIMAT0000461", "ENSG00000250742")
#' @param issymbol Boolean value, whether to set the node attribute "symbol"(gene symbol) in the network.

#' @return \code{runmodule} returns a list containing relevant data and results,
#'   including:
#'  \tabular{ll}{
#'    \code{GNCW} \tab the node-weighted network used for searching \cr
#'    \code{module} \tab list of genes comprising each module, named for the seed gene if the method is "local" 
#'    or the igraph class of the module if the method is "global" \cr
#'    \code{module.score.matrix} \tab contains Zm, Zn \cr
#'  }
#'
#' @references Hongbo Shi, Jiayao Li, Qiong Song et al. (2019) Systematic identification and analysis of dysregulated miRNA and transcription factor feed-forward loops in hypertrophic cardiomyopathy
#' @references Peilin Jia, Siyuan Zheng, Jirong Rong, Wei Zheng, Zhongming Zhao. (2011) Bioformatics. dmGWAS: dense module searching for genome-wide association studies in protein-protein interaction networks.
#' @references Daniela Beisser, Gunnar W. Klau, Thomas Dandekar  et al. (2019)  BioNet: an R-Package for the functional analysis of biological networks
#' @examples
#' \dontrun{
#' data("dataN")
#' gene2weight <- combinp(dataN[,c("type","logFC","PValue")])
#' interac <- interStringency(type = "transcription", spec ="hg",
#' 				stringency = "strict")
#' interac <- interac[,c("node_gene_ID","type","target_gene_ID")]
#'  res.list_global <- runmodule(network = interac, gene2weight, 
#'  								method = "global",FDR = 1e-14)
#'  res.list_local <- runmodule(network = interac, gene2weight, 
#'  		method = "local",maxsize=15, seletN = "MIMAT0000461")
#'  
#' }
#'
#' @export

runmodule <- function (network, gene2weight, method = c("global","local"), expr1 = NULL, expr2 = NULL, d = 2, r = 0.1, seletN = NULL, FDR = 1e-14, lambda = 0.5, min.size = 5, maxsize = 15, issymbol = TRUE) {
	########test###########
	###################
    library(igraph)
    if (is.null(expr1) | is.null(expr2)) {
		if(method == "local")
			res_no_edge <- runmodule_1(network = network, gene2weight= gene2weight, maxsize = maxsize, method = "local", d = d, r = r, seletN = seletN)
		if(method == "global")
			res_no_edge <- runmodule_1(network = network, gene2weight= gene2weight, method = method, FDR = FDR, issymbol = TRUE)
        return(res_no_edge)
    }
	nodes_user <- seletN
	if (is.null(nodes_user)){
		cat("please input seed Nodes "," ...\n", sep = "")
		return("ERROR")
	}
    GWPI <- generate_network(expr1= expr1, expr2 = expr2, network=network, geneweight = gene2weight)
    rm(expr1)
    rm(expr2)
    rm(gene2weight)
    rm(network)
    gc()
    cat("start searching at ", format(Sys.time(), "%H:%M, %b %d %Y"), 
        " ...\n", sep = "")
	#############node as seed to generate 
	sublist = list()
	nodes_user <- seletN
    for (node in nodes_user) {
        ng <- seedQueryJ_edge(GWPI, node, search_r= 1, r=r, lambda = lambda, maxsize = maxsize)
        if (vcount(ng) >= min.size)
            sublist[[node]] <- ng
    }
	dm.result <- sublist

	#############
    if (length(dm.result) == 0) 
        stop("No dense modules identified!\n")
    cat("extracting modules...\n", sep = "")
    genesets <- list()
    for (k in 1:length(dm.result)) {
        node = names(dm.result[k])
        g = dm.result[[k]]
        genesets[[node]] <- V(g)$name
    }
    seed.genes <- names(genesets)
    cat("removing identical modules...\n", sep = "")
    identical.idx <- list()
    for (k in 1:length(genesets)) {
        tmp.idx <- c(k)
        for (kt in 1:length(genesets)) {
            if (kt == k) 
                (next)()
            genesk <- genesets[[k]]
            genest <- genesets[[kt]]
            if (length(genesk) != length(genest)) 
                (next)()
            overlap = intersect(genesk, genest)
            if (length(overlap) == length(genest)) {
                tmp.idx <- c(tmp.idx, kt)
            }
        }
        if (length(tmp.idx) > 1) {
            tmp.idx <- sort(tmp.idx)
            identical.idx[[seed.genes[k]]] <- tmp.idx
        }
    }
    toremove.idx <- c()
    for (k in 1:length(identical.idx)) {
		if(length(identical.idx) == 0)
			break
        tmp.idx <- identical.idx[[k]]
        toremove.idx <- c(toremove.idx, tmp.idx[-1])
    }
    toremove.idx <- unique(toremove.idx)
	if(is.null(toremove.idx)){
		genesets.clear <- genesets
	}else {
		genesets.clear <- genesets[-toremove.idx]
		dm.result <- dm.result[-toremove.idx]
	}
    cat("permutation on random network...\n", sep = "")
    genesets.length <- unique(as.numeric(lapply(dm.result, vcount)))
    genesets.length.null.dis <- list()
    cat("module size: ")
    for (k in min(genesets.length):max(genesets.length)) {
        cat(k, ".", sep = "")
        genesets.length.null.dis[[as.character(k)]] <- random_network(size = k, 
            G = GWPI, lambda)
    }
    genesets.length.null.stat <- list()
    for (k in min(genesets.length):max(genesets.length)) {
        l.zperm <- genesets.length.null.dis[[as.character(k)]]
        k.mean <- mean(l.zperm)
        k.sd <- sd(l.zperm)
        genesets.length.null.stat[[as.character(k)]] = c(k.mean, 
            k.sd)
    }
    calculate_score <- function(g) {
        if (ecount(g) > 0) 
            subsum <- (1 - lambda) * sum(V(g)$weight)/sqrt(vcount(g)) + 
                lambda * sum(as.numeric(E(g)$weight))/sqrt(ecount(g))
        else subsum <- (1 - lambda) * sum(V(g)$weight)/sqrt(vcount(g))
        subsum
    }
    ms <- data.frame(gene = names(genesets.clear), Sm = -9, Sn = -9)
    for (k in 1:length(genesets.clear)) {
        ms[k, 2] <- calculate_score(dm.result[[k]])
        tmp <- genesets.length.null.stat[[as.character(vcount(dm.result[[k]]))]]
        ms[k, 3] = (ms[k, 2] - tmp[1])/tmp[2]
    }
    ms_ordered = ms[order(ms[, 3], decreasing = T), ]
    res.list <- list()
    res.list[["GNCW"]] = GWPI
    res.list[["module"]] = genesets.clear
    #res.list[["genesets.length.null.dis"]] = genesets.length.null.dis
    #res.list[["genesets.length.null.stat"]] = genesets.length.null.stat
    res.list[["module.score.matrix"]] = ms
    #res.list[["ordered.module.score.matrix"]] = ms_ordered
	save(res.list, file = paste("Lambda_", lambda, "_estimated_by_default_result.RData", sep = ""))
    cat("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"), " ...\n", sep = "")
    return(res.list)
}
