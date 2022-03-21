#' @title Generate corrected p-value based on p-value and logFC in the expression matrix
#'
#' @description based on the formula: corrected p-value = 2*( 1-pnorm( (-log10( p-value )) * abs( log2FC )) ), corrected p-value was generated
#'
#' @param node_attr A data frame containing three columns: type, logFC and p value, and the row name is the gene identifier.
#' @param islog Boolean value, whether to use the logFC, if FALSE, the weight is the p-value, or "TRUE", the corrected p-value is used
#' @return A data matrix containing three columns: type, gene, weight(corrected p-value), the row name is the gene identifier.
#'
#' @examples
#' data("dataN")
#' result <- combinp(dataN[,c("type","logFC","PValue")])
#'
#' @references Hongbo Shi, Jiayao Li, Qiong Song et al. (2019) Systematic identification and analysis of dysregulated miRNA and transcription factor feed-forward loops in hypertrophic cardiomyopathy
#' @export
combinp <- function(node_attr = NULL,islog = T){
	#
	#node_attr <- read.table(dir,row.names=1, sep= "\t")
	#gene2weight <- seq(1, length(node_attr[,1]))
	if(islog){
	  node_p2 <-apply(node_attr[,c("logFC", "PValue")],1,function(x) 2*( 1-pnorm( (-log10(unlist(x[2]))) * abs(unlist(x[1]))) ))
	  gene2weight <- data.frame(type = node_attr[,c("type")], gene=names(node_p2), weigth = node_p2)
	  rownames(gene2weight)<- rownames(node_attr)
	  }
	else{
	  gene2weight <- data.frame(type = node_attr[,c("type")], gene=rownames(node_attr), weigth=node_attr[,c("PValue")])
	  rownames(gene2weight)<- rownames(node_attr)
	  }

	colnames(gene2weight) <- c("type","gene","weight")
	min_p <- min(gene2weight[which(gene2weight$weight != 0),"weight"])
	gene2weight[which(gene2weight$weight == 0),"weight"] <- min_p
	gene2weight <-gene2weight[which(gene2weight$weight != 1),]
	return(gene2weight)
}

