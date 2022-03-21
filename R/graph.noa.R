.graph.noa <- function (network, file)
{
    attrib <- list.vertex.attributes(network)
    for (i in 1:length(attrib)) {
        if (is(get.vertex.attribute(network, attrib[i]))[1] ==
            "character") {
            type <- "String"
        }
        if (is(get.vertex.attribute(network, attrib[i]))[1] ==
            "integer") {
            type <- "Integer"
        }
        if (is(get.vertex.attribute(network, attrib[i]))[1] ==
            "numeric") {
            type <- "Double"
        }
        noa <- cbind(V(network)$name, rep("=", length(V(network))),
            get.vertex.attribute(network, attrib[i]))
        first.line <- paste(attrib[i], " ",
            type, ")", sep = "")
        write(first.line, file = paste(file, "_", attrib[i],
            ".NA", sep = ""), ncolumns = 1, append = FALSE,
            sep = " ")
        write.table(noa, row.names = FALSE, col.names = FALSE,
            file = paste(file, "_", attrib[i], ".NA",
                sep = ""), sep = " ", append = TRUE,
            quote = FALSE)
    }
}
