PathScore <- function (path, graph1, graph2, node.score) {
    sum(c(node.score[V(graph1)[path]$name], node.score[V(graph2)[unique(unlist(V(graph1)[path]$clusters))]$name]))
}
