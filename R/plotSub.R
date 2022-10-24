#' @title Plot of the subnetwork
#'
#' @description The function plots a network from graphNEL or igraph format. It is used to visualize the modules. For further plotting options use the plot.igraph
#' function of the igraph package. The shapes of the nodes can be changed according to the scores argument, then negative scores appear squared
#' The color of the nodes can be changed according to the diff.expr argument. Negative(positive) values lead to green(red) nodes.
#'
#' @param network A graph in igraph or graphNEL format.
#' @param layout Layout algorithm, e.g. layout.fruchterman.reingold or layout.kamada.kawai.
#' @param labels Labels for the nodes of the network
#' @param diff.expr Named numerical vector of log2FC of the nodes in the network for coloring of the nodes.
#' @param scores Named numerical vector of scores of the nodes for the shape of the node in the network.
#' @param main Main title of the plot.
#' @param vertex.size Numerical value or verctor for the size of the vertices.

#' @examples
#' library(igraph)
#' edgel <- cbind(c("1", "2", "3", "4", "5", "6", "7"),
#' c("b", "c", "d", "e", "f", "a", "b"))
#' g <- graph.edgelist(edgel, directed=TRUE)
#' V(g)$type <- c(rep("lncRNA",4),rep("miRNA",4),rep("circRNA",5))
#' plotSub(g)
#'



#' @references Daniela Beisser, Gunnar W. Klau, Thomas Dandekar et al. (2010) BioNet: an R-Package for the functional analysis of biological networks
#'
#' @export

plotSub <- function (network, layout = layout.fruchterman.reingold, labels = NULL,
                            diff.expr = NULL, scores = NULL, main = NULL, vertex.size = NULL)
{
  if (is(network, "graphNEL")) {
    network <- igraph.from.graphNEL(network)
  }
  if (is.null(V(network)$name)) {
    #V(network)$name <- as.character(V(network))
    V(network)$name <- as.character(V(network)$symbol)
  }
  if ("symbol" %in% list.vertex.attributes(network)){
  V(network)$name <- as.character(V(network)$symbol)
  }
  if (is.null(labels)) {
    if ("geneSymbol" %in% list.vertex.attributes(network)) {
      labels <- V(network)$geneSymbol
    }
    else {
      labels <- V(network)$name
    }
  }
  shapes <- rep("circle", length(V(network)))
  names(shapes) <- V(network)$name
#"circle"     "square"     "csquare"    "rectangle"
  #if (!is.null(type) && !is.null(names(type))) {
    #shapes[intersect(names(which(scores < 0)), V(network)$name)] <- "csquare"





  #}

  mytriangle <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.size <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }

    symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
            stars=cbind(vertex.size, vertex.size, vertex.size),
            add=TRUE, inches=FALSE)
  }
  # clips as a circle
  add.vertex.shape("triangle", clip=vertex.shapes("circle")$clip,
                   plot= mytriangle)

  if ( "type" %in% list.vertex.attributes(network)) {
    #scores <- V(network)$score
    #names(scores) <- V(network)$name
    #shapes[names(which(scores < 0))] <- "csquare"
    shapes <- sapply(V(network)$type, function(x) {
      color <- NULL
      if ( x == "lncRNA")
        color = "circle"
      else if (x == "miRNA")
        color = "square"
      else if (x == "circRNA")
        color = "triangle"
      else if (x == "rbp")
        color = "triangle"
      else
        color = "circle"
      return(color)
    }
    )
  }

  if (!is.null(diff.expr) && !is.null(names(diff.expr))) {
    coloring <- node.color(network, diff.expr)
  }else {
    coloring <- "SkyBlue2"
  }
  if (is.null(diff.expr) && "diff.expr" %in% list.vertex.attributes(network)) {
    diff.exprs = V(network)$diff.expr
    names(diff.exprs) <- V(network)$name
    coloring <- node.color(network, diff.exprs)
  }
  max.labels <- max(nchar(labels))
  network.size = length(V(network))
  vertex.size2 <- 8
  cex = 0.6
  if (network.size < 50) {
    if (max.labels > 2) {
      labels.dist <- 0.5
    }
    else {
      vertex.size2 <- 15
      labels.dist <- 0
    }
  }
  if (network.size < 100 && network.size >= 50) {
    if (max.labels > 2) {
      labels.dist <- 0.5
    }
    else {
      labels.dist <- 0
    }
  }
  if (network.size >= 100) {
    if (max.labels > 3) {
      labels.dist <- 0.5
    }
    else {
      labels.dist <- 0
    }
  }
  if (!is.null(vertex.size)) {
    vertex.size2 <- vertex.size
    labels.dist <- vertex.size/15
  }
  plot(network, layout = layout, vertex.size = vertex.size2,
       vertex.label = labels, vertex.label.cex = cex, vertex.label.dist = labels.dist,
       vertex.color = coloring, vertex.label.family = "sans",
       vertex.shape = shapes, main = main)
}
