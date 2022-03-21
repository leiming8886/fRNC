.graph.table <- function (network, file)
{
    attrib <- list.vertex.attributes(network)
    if (length(attrib) != 0) {
        node.attribs <- matrix(data = NA, ncol = length(attrib),
            nrow = length(V(network)))
        for (i in 1:length(attrib)) {
            node.attribs[, i] <- get.vertex.attribute(network,
                attrib[i])
        }
        node.attribs <- cbind(V(network)$name, node.attribs)
        colnames(node.attribs) <- c("id", attrib)
        write.table(node.attribs, file = paste(file, "_n.txt",
            sep = ""), row.names = FALSE, sep = "\t")
    }
    attrib <- list.edge.attributes(network)
    if (length(attrib) != 0) {
        edge.attribs <- matrix(data = NA, ncol = length(attrib),
            nrow = length(E(network)))
        for (i in 1:length(attrib)) {
            edge.attribs[, i] <- get.edge.attribute(network,
                attrib[i])
        }
        edgelist.names <- get.edgelist(network, names = TRUE)
        edge.attribs <- cbind(edgelist.names, edge.attribs)
        colnames(edge.attribs) <- c("nodeA", "nodeB",
            attrib)
        write.table(edge.attribs, file = paste(file, "_e.txt",
            sep = ""), row.names = FALSE, sep = "\t")
    }
    if (length(attrib) == 0) {
        edgelist.names <- get.edgelist(network, names = TRUE)
        colnames(edgelist.names) <- c("node", "target")
        write.table(edgelist.names, file = paste(file, "_e.txt",
                                               sep = ""), row.names = FALSE, sep = "\t")
    }


}
