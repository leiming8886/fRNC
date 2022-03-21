FastHeinz <- function (network, scores)
{
    #library(methods)
    #library(stats)
    #if (!requireNamespace("igraph", quietly = TRUE)) {
    #    stop("Package \"igraph\" needed for this function to work. Please install it.",
    #         call. = FALSE)
    #}
    net.flag <- FALSE
    if (is.null(names(scores))) {
        warning("scores were empty")
    }
    if (any(is.na(names(scores)))) {
        warning("scores has NA values, and will be removed")
        scores <- scores[!is.na(names(scores))]
    }
    if (is(network, "graphNEL")) {
        #convert grahpNEL format of the network to igraph
        network <- igraph.from.graphNEL(network)
        net.flag <- TRUE
    }
    if (is.null(V(network)$name)) {
        #convert mode of name of vertex to character
        V(network)$name <- as.character(V(network))
    }
    V(network)$score <- scores[V(network)$name]
    #extract pos nodes from network
    pos.nodes <- names(scores[which(scores > 0)])
    if (length(pos.nodes) == 0) {
        warning("the network have no positive nodes")
        module <- graph.empty(n = 0, directed = FALSE)
        return(module)
    }
    if (length(pos.nodes) == 1) {
        module <- subNetwork_only(pos.nodes, network)
        if (net.flag) {
            nE <- ecount(module)
            module <- simplify(module, remove.multiple = TRUE)
            if (nE != ecount(module)) {
                warning("Multiple edges between two nodes had to be removed")
            }
            module <- igraph.to.graphNEL(module)
        }
        return(module)
    }
    pos.subgraph <- subNetwork_only(pos.nodes, network)
    conn.comp.graph <- decompose.graph(pos.subgraph)
    score.comp <- unlist(lapply(lapply(conn.comp.graph, get.vertex.attribute,
        "score"), sum))
    conn.comp.graph <- conn.comp.graph[order(score.comp, decreasing = TRUE)]
    score.comp <- sort(score.comp, TRUE)
    for (i in 1:length(conn.comp.graph)) {
        conn.comp.graph[[i]]$score <- score.comp[i]
    }
    v.id <- seq(1, vcount(network))
    names(v.id) <- V(network)$name
    edgelist <- get.edgelist(network, FALSE)
    edgelist1 <- edgelist[, 1]
    edgelist2 <- edgelist[, 2]
    for (i in 1:length(conn.comp.graph)) {
        new.id <- length(V(network)) + i
        for (j in as.character(v.id[V(conn.comp.graph[[i]])$name])) {
            edgelist1[which(edgelist1 == j)] <- new.id
            edgelist2[which(edgelist2 == j)] <- new.id
        }
    }
    new.ids <- seq(length(V(network)) + 1, length(V(network)) +
        length(conn.comp.graph))
    new.names <- paste("cluster", seq(1:length(conn.comp.graph)),
        sep = "")
    names(new.ids) <- new.names
    v.id <- c(v.id, new.ids)
    v.name <- names(v.id)
    names(v.name) <- v.id
    new.edgelist <- cbind(v.name[as.character(edgelist1)], v.name[as.character(edgelist2)])
    interactome2 <- graph.edgelist(new.edgelist, FALSE)
    E(interactome2)$weight <- rep(0, length(E(interactome2)))
    interactome2 <- simplify(interactome2, remove.loops = TRUE,
        remove.multiple = TRUE)
    score1 <- scores[V(interactome2)$name]
    names(score1) <- V(interactome2)$name
    score1[which(is.na(score1))] <- 0
    score.degree <- score1/(degree(interactome2) + 1)
    V(interactome2)$score.degree <- score.degree
    E(interactome2)$weight <- -(V(interactome2)[get.edgelist(interactome2,
        FALSE)[, 1]]$score.degree + V(interactome2)[get.edgelist(interactome2,
        FALSE)[, 2]]$score.degree)
    node.score <- scores[V(interactome2)$name]
    names(node.score) <- V(interactome2)$name
    node.score.cluster <- sapply(conn.comp.graph, get.graph.attribute,
        "score")
    names(node.score.cluster) <- new.names
    node.score[grep("cluster", names(node.score))] <- node.score.cluster[names(node.score[grep("cluster",
        names(node.score))])]
    decomp.graphs <- decompose.graph(interactome2)
    sum.pos <- lapply(decomp.graphs, function(x) {
        sum(node.score[names(which(node.score[V(x)$name] > 0))])
    })
    interactome2 <- decomp.graphs[[which.max(sum.pos)]]
    rm(decomp.graphs)
    mst <- minimum.spanning.tree(interactome2, weights = E(interactome2)$weight)
    mst.cluster.id <- grep("cluster", V(mst)$name)
    names(mst.cluster.id) <- V(mst)[mst.cluster.id]$name
    mst.cluster.id <- mst.cluster.id[order(as.numeric(matrix(unlist(strsplit(names(mst.cluster.id),
        "cluster")), nrow = 2)[2, ]))]
    all.ids <- c()
    if (length(mst.cluster.id) == 1) {
        neg.node.ids.2 = c()
    }
    else {
        for (j in 1:(length(mst.cluster.id) - 1)) {
            path <- unlist(get.all.shortest.paths(mst, from = mst.cluster.id[j],
                to = mst.cluster.id[(j + 1):length(mst.cluster.id)]))
            all.ids <- c(all.ids, path)
        }
        all.ids <- unique(all.ids)
        sub.interactome2 <- subNetwork_only(V(mst)[all.ids]$name,
            interactome2)
        neg.node.ids <- which(node.score[V(sub.interactome2)$name] <
            0)
        for (i in neg.node.ids) {
            V(sub.interactome2)[i]$clusters <- list(neighbors(sub.interactome2,
                v = i)[grep("cluster", V(sub.interactome2)[neighbors(sub.interactome2,
                v = i)]$name)])
        }
        score.neg.nodes <- c()
        for (i in neg.node.ids) {
            if (length(V(sub.interactome2)[i]$clusters[[1]]) >
                0) {
                score.neg.nodes <- c(score.neg.nodes, sum(node.score[V(sub.interactome2)[c(i,
                  V(sub.interactome2)[i]$clusters[[1]])]$name]))
            }
            else {
                score.neg.nodes <- c(score.neg.nodes, node.score[V(sub.interactome2)[i]$name])
            }
        }
        neg.node.ids.2 <- neg.node.ids[score.neg.nodes > 0]
    }
    if (length(neg.node.ids.2) == 0) {
        module <- unlist(lapply(conn.comp.graph, get.vertex.attribute,
            "name")[as.numeric(matrix(unlist(strsplit(names(node.score.cluster)[which.max(node.score.cluster)],
            "cluster")), nrow = 2)[2, ])])
        module <- subNetwork_only(module, network)
        if (net.flag) {
            nE <- ecount(module)
            module <- simplify(module, remove.multiple = TRUE)
            if (nE != ecount(module)) {
                warning("Multiple edges between two nodes had to be removed")
            }
            module <- igraph.to.graphNEL(module)
        }
        return(module)
    }
    subg <- largestComp(induced.subgraph(sub.interactome2, neg.node.ids.2))
    mst.subg <- minimum.spanning.tree(subg, E(subg)$weight)
    max.score <- 0
    best.path <- c()
    for (i in 1:(length(V(mst.subg)))) {
        path <- get.all.shortest.paths(mst.subg, from = V(mst.subg)[i])
        path.score <- unlist(lapply(path$res, PathScore,
            graph1 = mst.subg, graph2 = sub.interactome2, node.score = node.score))
        best.pos <- which.max(path.score)
        if (path.score[[best.pos]] > max.score) {
            best.path <- path$res[[best.pos]]
            max.score <- path.score[[best.pos]]
        }
    }
    if (length(best.path) != 1) {
        cluster.list <- V(mst.subg)[best.path]$clusters
        names.list <- as.character(1:length(cluster.list))
        names(cluster.list) <- names.list
        names(best.path) <- names.list
        for (i in names.list) {
            res <- lapply(cluster.list, intersect, cluster.list[[i]])
            if (length(intersect(unlist(cluster.list[as.character(which(as.numeric(names.list) <
                as.numeric(i)))]), unlist(cluster.list[as.character(which(as.numeric(names.list) >
                as.numeric(i)))]))) > 0) {
                if (length(setdiff(res[[i]], unique(unlist(res[names(res) !=
                  i])))) == 0) {
                  cluster.list <- cluster.list[names(cluster.list) !=
                    i]
                  names.list <- names.list[names.list != i]
                }
            }
        }
        best.path <- best.path[names.list]
    }
    module <- V(mst.subg)[best.path]$name
    pos.cluster <- V(sub.interactome2)[unique(unlist(V(mst.subg)[best.path]$clusters))]$name
    module <- c(module, unlist(lapply(conn.comp.graph, get.vertex.attribute,
        "name")[as.numeric(matrix(unlist(strsplit(pos.cluster,
        "cluster")), nrow = 2)[2, ])]))
    module <- subNetwork_only(module, network)
    if (net.flag) {
        nE <- ecount(module)
        module <- simplify(module, remove.multiple = TRUE)
        if (nE != ecount(module)) {
            warning("Multiple edges between two nodes had to be removed")
        }
        module <- igraph.to.graphNEL(module)
    }
    return(module)
}
