.graph.sif <-function (network, file)
{
    edgelist.names <- get.edgelist(network, names = TRUE)
    edgelist.names <- cbind(edgelist.names[, 1], rep("g2g",
                                                     length(E(network))), edgelist.names[, 2])
    write.table(edgelist.names, row.names = FALSE, col.names = FALSE,
                file = paste(file, ".sif", sep = ""), sep = "\t",
                quote = FALSE)
    return(edgelist.names)
}
