largestComp <- function (network)
{
    #if (!require(igraph)) {
    #    stop("igraph must be pre-installed!\n")
    #}
    if (is(network, "graphNEL")) {
        cc <- RBGL::connectedComp(network)
        idx <- which.max(Biobase::listLen(cc))
        return(graph::subGraph(cc[[idx]], network))
    }
    else if (is(network, "igraph")) {
        clust <- clusters(network)
        cid <- which.max(clust$csize)
        lg <- induced.subgraph(network, V(network)[clust$membership ==
            cid])
        return(lg)
    }
}
