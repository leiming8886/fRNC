generate_edge_weight <- function (expr1, expr2, network, geneweight) 
{
	

    if (!require(igraph)) 
        stop("igraph must be pre-installed!\n")

	cat("Number of genes with node weight: ", length(geneweight$gene), "\n", sep = "")
	cat("Number of genes with expression value in case samples: ", 
		nrow(expr1), "\n", sep = "")
	cat("Number of genes with expression value in control sampels: ", 
		nrow(expr2), "\n", sep = "")
	gene1 <- as.character(unique(geneweight$gene))
	gene2 <- as.character(unique(rownames(expr1)))
	gene3 <- as.character(unique(rownames(expr2)))
	gene4 <- unique(unlist(network[,c("node_gene_ID","target_gene_ID")]))
	common_gene <- intersect(intersect(gene1, gene4), intersect(gene2, 
		gene3))
	expr1 <- expr1[is.element(rownames(expr1), common_gene), ]
	expr2 <- expr2[is.element(rownames(expr2), common_gene), ]
	expr1 <- expr1[order(rownames(expr1)), ]
	expr2 <- expr2[order(rownames(expr2)), ]
	expr1_temp <- t(expr1)
	colnames(expr1_temp) <- rownames(expr1)
	expr1_temp <- cor(expr1_temp, use = "pairwise.complete.obs",method = "spearman")
	diag(expr1_temp) <- 0
	expr2_temp <- t(expr2)
	colnames(expr2_temp) <- rownames(expr2)
	expr2_temp <- cor(expr2_temp, use = "pairwise.complete.obs",method = "spearman")
	diag(expr2_temp) <- 0
	ncase <- ncol(expr1)
	ncontrol <- ncol(expr2)
	adjust_factor <- 1/(1.06/(ncase - 3) + 1.06/(ncontrol - 3))^0.5
	expr1_temp <- 0.5 * log((1 + expr1_temp)/(1 - expr1_temp))
	expr2_temp <- 0.5 * log((1 + expr2_temp)/(1 - expr2_temp))
	rewire <- adjust_factor * (expr1_temp - expr2_temp)
	rm(expr1_temp)
	rm(expr2_temp)
	gc()
	diag(rewire) <- 0
	temp <- rewire[lower.tri(rewire)]
	temp <- temp[is.finite(temp)]
	mu <- mean(temp)
	sigma <- sd(temp)
	rewire <- (rewire - mu)/sigma
	rewire <- abs(rewire)
	rewire <- qnorm(1 - (pnorm(rewire, lower.tail = F) * 
		2))
	return(rewire)

}
