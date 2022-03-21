.graph.eda <- function (network, file, edgelist.names)
{
    attrib <- list.edge.attributes(network)
    for (i in 1:length(attrib)) {
        if (is(get.edge.attribute(network, attrib[i]))[1] ==
            "character") {
            type <- "String"
        }
        if (is(get.edge.attribute(network, attrib[i]))[1] ==
            "integer") {
            type <- "Integer"
        }
        if (is(get.edge.attribute(network, attrib[i]))[1] ==
            "numeric") {
            type <- "Double"
        }
        eda <- cbind(cbind(edgelist.names[, 1], edgelist.names[, 3]), rep("=",
            length(E(network))), get.edge.attribute(network,
            attrib[i]))
        first.line <- paste(attrib[i], " ",
            type, ")", sep = "")
        write(first.line, file = paste(file, "_", attrib[i],
            ".EA", sep = ""), ncolumns = 1, append = FALSE,
            sep = " ")
        write.table(eda, row.names = FALSE, col.names = FALSE,
            file = paste(file, "_", attrib[i], ".EA",
                sep = ""), sep = " ", append = TRUE,
            quote = FALSE)
    }





}

