.XGMML.nodes <- function (network) 
{
    requireNamespace("XML")
    c.node <- rep("node", length(V(network)))
    nodes <- lapply(c.node, XML::xmlNode)
    attrib <- list.vertex.attributes(network)
    node.attribs <- matrix(data = NA, nrow = length(attrib), 
        ncol = length(V(network)))
    for (i in 1:length(attrib)) {
        if (is(get.vertex.attribute(network, attrib[i]))[1] == 
            "character") {
            type <- "string"
        }
        if (is(get.vertex.attribute(network, attrib[i]))[1] == 
            "integer") {
            type <- "integer"
        }
        if (is(get.vertex.attribute(network, attrib[i]))[1] == 
            "numeric") {
            type <- "real"
        }
        node.attribs[i, ] = paste("att type=", "\"", 
            type, "\"", " name=", "\"", attrib[i], 
            "\"", " value=", "\"", get.vertex.attribute(network, 
                attrib[i]), "\"", sep = "")
    }
    node.attribs <- matrix(lapply(node.attribs, XML::xmlNode), 
        nrow = length(attrib), ncol = length(V(network)))
    if (is.null(V(network)$name)) {
        V(network)$name <- as.character(V(network))
    }
    node.label <- V(network)$name
    node.id <- as.vector(V(network))
    for (i in 1:length(V(network))) {
        nodes[[i]] <- XML::addAttributes(nodes[[i]], label = node.label[i], 
            id = node.id[i])
        nodes[[i]] <- XML::append.xmlNode(nodes[[i]], node.attribs[, 
            i])
    }
    return(nodes)
}
