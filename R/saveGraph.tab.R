.saveGraph.tab <- function (graph, filename)
{
    en <- graph::edgeNames(graph)
    x <- gsub("\\~", "\t", en)
    cat(x, file = paste(filename, ".txt", sep = ""),
        sep = "\n")
}
