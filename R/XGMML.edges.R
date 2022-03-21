.XGMML.edges <- function (network)
{
    requireNamespace("XML")
    c.edge <- rep("edge", length(E(network)))
    edges <- lapply(c.edge, XML::xmlNode)
    edgelist.names <- get.edgelist(network, names = TRUE)
    edgelist.names <- paste(edgelist.names[, 1], edgelist.names[,
        2], sep = " (g2g) ")
    edgelist.ids <- get.edgelist(network, names = FALSE)
    attrib <- list.edge.attributes(network)
    edge.attribs <- matrix(data = NA, nrow = length(attrib),
        ncol = length(E(network)))
    for (i in 1:length(attrib)) {
        #i
        if (length(attrib) == 0){
            break
        }
        if (is(get.edge.attribute(network, attrib[i]))[1] ==
            "character") {
            type <- "string"
        }
        if (is(get.edge.attribute(network, attrib[i]))[1] ==
            "integer") {
            type <- "integer"
        }
        if (is(get.edge.attribute(network, attrib[i]))[1] ==
            "numeric") {
            type <- "real"
        }
        edge.attribs[i, ] <- paste("att type=", "\"",
            type, "\"", " name=", "\"", attrib[i],
            "\"", " value=", "\"", get.edge.attribute(network,
                attrib[i]), "\"", sep = "")
    }
    edge.attribs <- matrix(lapply(edge.attribs, XML::xmlNode),
        nrow = length(attrib), ncol = length(E(network)))
    for (i in 1:length(E(network))) {
        edges[[i]] <- XML::addAttributes(edges[[i]], label = edgelist.names[i],
            source = edgelist.ids[i, 1], target = edgelist.ids[i,
                2])
        if(length(attrib) != 0){
        edges[[i]] <- XML::append.xmlNode(edges[[i]], c(edge.attribs[,
            i]))
        }
    }
    return(edges)
}
