integGM <- function (G, genes, weights, type, issymbol = TRUE, simplify = TRUE)
{
    #if (!require(igraph)) {
    #    stop("igraph must be pre-installed!\n")
    #}
    if (is.element("weight", list.vertex.attributes(G))) {
        cat("Warning: previous G node weight replaced!\n")
    }

    names(weights) <- genes
    names(type) <- genes
    genes <- intersect(genes, V(G)$name)
    subG <- induced.subgraph(G, genes)
    V(subG)$weight <- weights[V(subG)$name]
    V(subG)$type <- type[V(subG)$name]
    if (issymbol){
        symbol <- fRNC::IDsymbol$gene_Name
        names(symbol) <- fRNC::IDsymbol$genes
        V(subG)$symbol <- as.character(symbol[V(subG)$name])
    }
    #if(simplify)
    #    subG = simplify(subG, edge.attr.comb=toString)

    return(subG)
}
