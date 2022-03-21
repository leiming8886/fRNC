generate_network <- function (expr1, expr2 = NULL, network, geneweight) 
{
    require(igraph)
    if (!require(igraph)) 
        stop("igraph must be pre-installed!\n")
    if (sum(is.na(geneweight$weight)) > 0) 
        stop("Have missing P values!\n")
    if (min(geneweight$weight) <= 0 | max(geneweight$weight) >= 1) 
        stop("P values out of range 0<p<1\n")
    edgeweight <- generate_edge_weight(expr1, expr2, network, 
        geneweight)
    gene1 <- as.character(unique(geneweight$gene))
    gene2 <- as.character(rownames(edgeweight))
    gene <- intersect(gene1, gene2)
    geneweight <- geneweight[is.element(geneweight$gene, gene), 
        ]
    index <- is.element(rownames(edgeweight), gene)
    edgeweight <- edgeweight[index, index]
	rawG <- graph.data.frame(network[,c("node_gene_ID","target_gene_ID")], directed = F)
	rawG <- set_edge_attr(rawG, "types", value= as.character(network$type))
	rawG <- simplify(rawG,edge.attr.comb=toString)
	##################################
	gene_type_all <- geneweight$type
    names(gene_type_all) <- geneweight$gene
	genes <- intersect(geneweight$gene, V(rawG)$name)
    rawG <- induced.subgraph(rawG, genes)
    V(rawG)$type <- as.character(gene_type_all[V(rawG)$name])
    symbol <- fRNC::IDsymbol$gene_Name
    names(symbol) <- fRNC::IDsymbol$genes
    V(rawG)$symbol <- as.character(symbol[V(rawG)$name])
	
	##################################
	###
    gene <- as.character(geneweight$gene)
    g.weight <- sapply(as.numeric(geneweight$weight), function(x) qnorm(1 - 
        x))
    names(g.weight) <- gene
    gene <- intersect(gene, V(rawG)$name)
    subG <- induced.subgraph(rawG, gene)
    if (is.element("weight", list.vertex.attributes(subG))) 
        cat("Warning: previous G node weight replaced!\n")
    V(subG)$weight <- g.weight[V(subG)$name]
    subG <- simplify(subG, edge.attr.comb=toString)
    if (is.element("weight", list.edge.attributes(subG))) 
        cat("Warning: previous edge weight replaced!\n")
    index1 <- match(get.edgelist(subG)[, 1], rownames(edgeweight))
    index2 <- match(get.edgelist(subG)[, 2], rownames(edgeweight))
    index <- cbind(index1, index2)
    subG <- set.edge.attribute(subG, "weight", index = E(subG), 
        value = edgeweight[index])
    subG <- simplify(subG, edge.attr.comb=toString)
    rm(edgeweight)
    rm(rawG)
    gc()
	del_e <- c()
	if (sum(is.na(E(subG)$weight)) > 0) {
		del_e <- E(subG)
		del_e <- del_e[is.na(del_e$weight)]
		subG <- delete_edges(subG, del_e)
		subG <- simplify(subG, edge.attr.comb=toString)
		}
    edge2weight <- as.numeric(E(subG)$weight)
    index <- E(subG)$weight == Inf
    if (any(index)) 
        E(subG)$weight[index] <- max(edge2weight[is.finite(edge2weight)])
    index <- E(subG)$weight == (-Inf)
    if (any(index)) 
        E(subG)$weight[index] <- min(as.numeric(edge2weight[is.finite(edge2weight)]))
    cat("The final background network has ", vcount(subG), " nodes and ", 
        ecount(subG), " edges.\n", sep = "")
    return(subG)
}












