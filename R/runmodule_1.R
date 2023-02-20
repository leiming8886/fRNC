runmodule_1 <- function (network, gene2weight, maxsize = 15, method = c("global","local"), d = 2, r = 0.1, seletN = NULL, FDR = 1e-14, issymbol = TRUE)
{
	#library(igraph)
	#library(stats)
	maxs <- maxsize
	if (min(gene2weight[, c("weight")]) <= 0 | max(gene2weight[, c("weight")]) >=
			1) {
			stop("P values out of range 0 < p < 1")
	}
	if(!missing(method) & length(method)>1) stop("Only one 'method' allowed.")
	method <- match.arg(method)
	cat("Number of genes with node weight: ", length(gene2weight[, c("gene")]), "\n", sep = "")
	rawG <- graph.data.frame(network[,c("node_gene_ID","target_gene_ID")], directed = F)
	rawG <- set_edge_attr(rawG, "types", value= as.character(network$type))
	rawG <- simplify(rawG,edge.attr.comb=toString)


	if(method == "local"){
		g.weight <- sapply(as.numeric(gene2weight[, c("weight")]), function(x) qnorm(1 -
		x))
	}
	if(method == "global"){
		pvals <- gene2weight$weight
		names(pvals) <- rownames(gene2weight)
		fb <- fitBumModel(pvals)
		g.weight <- scoreFunction(fb, fdr = FDR)
	}
	intG <- integGM(rawG, as.character(gene2weight[, c("gene")]), g.weight, as.character(gene2weight[, c("type")]), issymbol = issymbol )
	GNCW <- simplify(intG, edge.attr.comb=toString)
	cat("start searching at ", format(Sys.time(), "%H:%M, %b %d %Y"),
		" ...\n", sep = "")
	if( method == "global"){
		scores <- V(GNCW)$weight
		names(scores)<- V(GNCW)$name
		module <- FastHeinz(GNCW, scores)
		genes.idx <- seq(1, length(V(GNCW)$name))
		graph.g.weight = data.frame(GNCWene = V(GNCW)$name, gain.weight = V(GNCW)$weight)
		genes <- V(module)$name
		idx <- match(genes, graph.g.weight[, 1])
		idx <- idx[!is.na(idx)]
		ZM <- sum(graph.g.weight[idx, 2])/sqrt(length(idx))
		l.zperm <- c()
		for (j in 1:1e+05) {
            idx.pseudo = sample(genes.idx, size = length(V(module)))
            l.zperm <- c(l.zperm, sum(graph.g.weight[idx.pseudo,
                2])/sqrt(length(idx.pseudo)))
        }
		k.mean <- mean(l.zperm)
		k.sd <- sd(l.zperm)
		ZN = (ZM - k.mean)/k.sd
		res.list <- list()
		res.list[["GNCW"]] = GNCW
		res.list[["module"]] = module
		res.list[["module.score.matrix"]] = data.frame(Zm=ZM,Zn=ZN)
		cat("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"),
		" ...\n", sep = "")
		return(res.list)

	}
	if(method == "local"){
		sublist = list()
		nodes_user <- seletN
		if (is.null(nodes_user)){
				cat("please input Nodes "," ...\n", sep = "")
		}
		for (node in nodes_user) {#V(GNCW)$name
			if(!node %in% V(GNCW)$name){
			cat("Error: the Node  ",node," have no expression value or the interaction\n", sep = "")
			next
			}
			ng <- seedQuery_GS(GNCW, node, maxsize = maxs, search_r = d, r = r)
			if (vcount(ng) >= 5)
			sublist[[node]] <- ng
		}

		dm.result <- sublist
		cat("extracting modules...\n", sep = "")
		genesets <- list()
		for (k in 1:length(dm.result)) {
			node = names(dm.result[k])
			g = dm.result[[k]]
			genesets[[node]] <- V(g)$name
		}
		seed.genes <- names(genesets)
		#cat("removing identical modules...\n", sep = "")
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
		}
		cat("permutation on random network...\n", sep = "")
		genesets.length <- c()
		for (k in 1:length(genesets.clear)) {
			genes <- genesets.clear[[k]]
			genesets.length <- c(genesets.length, length(genes))
		}

		genesets.length <- unique(genesets.length)
		genes.idx <- seq(1, length(V(GNCW)$name))
		graph.g.weight = data.frame(GNCWene = V(GNCW)$name, gain.weight = V(GNCW)$weight)
		genesets.length.null.dis <- list()
		length.max = max(genesets.length)
		length.min = min(genesets.length)
		for (k in length.min:length.max) {
			l.zperm <- c()
			for (j in 1:1e+05) {
				idx.pseudo = sample(genes.idx, size = k)
				l.zperm <- c(l.zperm, sum(graph.g.weight[idx.pseudo,
					2])/sqrt(length(idx.pseudo)))
			}
			genesets.length.null.dis[[as.character(k)]] = l.zperm
			cat(k, ".", sep = "")
		}
		genesets.length.null.stat <- list()
		for (k in length.min:length.max) {
			l.zperm <- genesets.length.null.dis[[as.character(k)]]
			k.mean <- mean(l.zperm)
			k.sd <- sd(l.zperm)
			genesets.length.null.stat[[as.character(k)]] = c(k.mean,
				k.sd)
		}
		zim <- data.frame(gene = names(genesets.clear), Zm = -9,
			Zn = -9, zcount = -9)
		for (k in 1:length(genesets.clear)) {
			genes <- genesets.clear[[k]]
			idx <- match(genes, graph.g.weight[, 1])
			idx <- idx[!is.na(idx)]
			zim[k, 2] <- sum(graph.g.weight[idx, 2])/sqrt(length(idx))
			tmp <- genesets.length.null.stat[[as.character(length(idx))]]
			zim[k, 3] = (zim[k, 2] - tmp[1])/tmp[2]
			zim[k, 4] = sum(genesets.length.null.dis[[as.character(length(idx))]] >=
				zim[k, 2])
		}
		zom = zim[order(zim[, 3], decreasing = T), ]
		res.list <- list()
		res.list[["GNCW"]] = GNCW
		#res.list[["graph.g.weight"]] = graph.g.weight
		res.list[["module"]] = genesets.clear
		#res.list[["genesets.length.null.dis"]] = genesets.length.null.dis
		#res.list[["genesets.length.null.stat"]] = genesets.length.null.stat
		#res.list[["zi.matrix"]] = zim
		res.list[["module.score.matrix"]] = zom[,c("gene", "Zm", "Zn")]
		save(res.list, file = "RESULT.list.RData")
		cat("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"),
			" ...\n", sep = "")
		return(res.list)
	}
}
