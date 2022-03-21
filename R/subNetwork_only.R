#' @title Create a subGraph
#'
#' @description The function creates a subgraph with the nodes given in the nodeList
#'
#' @param nodeList Character vector of nodes, contained in the subgraph.
#' @param network Graph that is used for subgraph extraction
#'
#' @return A graph object.
#'
#' @examples
#' library(igraph)
#' edgel <- cbind(c("a1", "a2", "a3", "a4", "a5", "a6", "a7"), c("b", "c", "d", "e", "f", "a", "b"))
#' g <- graph.edgelist(edgel, directed=TRUE)
#' node.list <- c("a1", "b", "c", "a1")
#' graph <- subNetwork_only(nodeList=node.list, network=g)
#'
#'
#' @export

subNetwork_only <- function (nodeList, network)
{
    #library(stats)
    #library(methods)
    if (is(network, "igraph")) {
        mapping <- seq(1, (length(V(network))))
        if (is.null(V(network)$name)) {
            V(network)$name <- as.character(V(network))
        }
        names(mapping) <- V(network)$name
        nodeList = mapping[nodeList]
        if (any(is.na(nodeList))) {
            nodeList = na.omit(nodeList)
            warning("Not all nodes found in network")
        }
        subgr <- induced.subgraph(network, vids = nodeList)
    }
    else {
        subgr <- graph::subGraph(graph::nodes(network)[graph::nodes(network) %in%
            nodeList], network)
    }
    return(subgr)
}
