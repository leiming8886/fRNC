#' @title save of the subnetwork
#'
#' @description The function plots a network from graphNEL or igraph format. It is used to visualize the modules. For further plotting options use the plot.igraph
#' function of the igraph package. The shapes of the nodes can be changed according to the scores argument, then negative scores appear squared
#' The color of the nodes can be changed according to the diff.expr argument. Negative(positive) values lead to green(red) nodes.
#'
#' @param network A graph in igraph or graphNEL format.
#' @param name Name of the network, only needed for the XGMML format.
#' @param file File name to save.
#' @param type Type in which graph shall be saved.
#' @examples
#' library(igraph)
#' edgel <- cbind(c("1", "2", "3", "4", "5", "6", "7"), c("b", "c", "d", "e", "f", "a", "b"))
#' g <- graph.edgelist(edgel, directed=TRUE)
#' V(g)$type <- c(rep("lncRNA",4),rep("miRNA",4),rep("circRNA",5))
#' saveNetwork(g,file ="g", type = "XGMML")
#' @references Daniela Beisser, Gunnar W. Klau, Thomas Dandekar et al. (2010) BioNet: an R-Package for the functional analysis of biological networks
#'
#' @export

saveNetwork <- function (network, name = "network", file, type = c("table",
    "XGMML", "sif", "tab"))
{
    file <- .cleanFile(file)
    type <- match.arg(type)
    if (is.null(V(network)$ID)) {
        V(network)$ID <- as.character(V(network)$name)
    }
    if (is.null(V(network)$name)) {
        #V(network)$name <- as.character(V(network))
        V(network)$name <- as.character(V(network)$symbol)
    }

    if ("symbol" %in% list.vertex.attributes(network)){
        V(network)$name <- as.character(V(network)$symbol)
    }
    if (type == "XGMML") {
        requireNamespace("XML")
        addNode <- XML::addNode
        if (is(network, "graphNEL")) {
            network <- igraph.from.graphNEL(network)
        }
        top <- .XGMML.destription(name = name)
        print("...adding nodes")
        nodes <- .XGMML.nodes(network = network)
        top <- XML::append.xmlNode(top, nodes)
        print("...adding edges")
        edges <- .XGMML.edges(network = network)
        top <- XML::append.xmlNode(top, edges)
        print("...writing to file")
        XML::saveXML(top, file = paste(file, ".xgmml",
            sep = ""), encoding = "UTF-8")
        if ("package:XML" %in% search()) {
            detach("package:XML")
        }
        addNode <- graph::addNode
    }
    if (type == "table") {
        if (is(network, "graphNEL")) {
            network <- igraph.from.graphNEL(network)
        }
        .graph.table(network = network, file = file)
    }
    if (type == "sif") {
        if (is(network, "graphNEL")) {
            network <- igraph.from.graphNEL(network)
        }
        edges <- .graph.sif(network = network, file = file)
        if (length(list.edge.attributes(network)) != 0) {
            .graph.eda(network = network, file = file, edgelist.names = edges)
        }
        if (length(list.vertex.attributes(network)) != 0) {
            .graph.noa(network = network, file = file)
        }
    }
    if (type == "tab") {
        if (is(network, "igraph")) {
            nE <- ecount(network)
            network <- simplify(network, remove.multiple = TRUE)
            if (nE != ecount(network)) {
                warning("Multiple edges are not allowed for the graphNEL format, they had to be removed")
            }
            network <- igraph.to.graphNEL(network)
            .saveGraph.tab(graph = network, filename = file)
        }
    }

}

