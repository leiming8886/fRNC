.XGMML.destription <- function (name)
{
    requireNamespace("XML")
    top <- XML::xmlNode("graph", attrs = c(label = name,
        `xmlns:dc` = "http://purl.org/dc/elements/1.1/",
        `xmlns:xlink` = "http://www.w3.org/1999/xlink",
        `xmlns:rdf` = "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        `xmlns:cy` = "http://www.cytoscape.org",
        xmlns = "http://www.cs.rpi.edu/XGMML"))
    top <- XML::append.xmlNode(top, XML::xmlNode("att",
        attrs = c(name = "documentVersion", value = "1.1")))
    d <- XML::xmlNode("rdf:Description", attrs = c(`rdf:about` = "http://www.cytoscape.org/"))
    d <- XML::append.xmlNode(d, XML::xmlNode("dc:type",
        "RPB ncRNA Interaction"))
    d <- XML::append.xmlNode(d, XML::xmlNode("dc:description",
        "N/A"))
    d <- XML::append.xmlNode(d, XML::xmlNode("dc:identifier",
        "N/A"))
    d <- XML::append.xmlNode(d, XML::xmlNode("dc:date",
        Sys.time()))
    d <- XML::append.xmlNode(d, XML::xmlNode("dc:title",
        name))
    d <- XML::append.xmlNode(d, XML::xmlNode("dc:format",
        "fRNC-Cytoscape-XGMML"))
    c <- XML::xmlNode("att", attrs = c(name = "networkMetadata"),
        XML::xmlNode("rdf:RDF", d))
    top <- XML::append.xmlNode(top, c)
    return(top)
}
