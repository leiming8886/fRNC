#' @title extract lncRNA-circRNA ceRNA pair
#'
#' @description extract lncRNA-circRNA ceRNA pair
#'
#' @param interac the interaction data matrix, in which the column name contain "node_gene_ID", "type" and  "target_gene_ID"
#' @param N_mi  a numeric value, shared miRNAs by a lncRNA and circRNA
#'
#' @return interac_temp the interaction data matrix, in which the column name contain "node_gene_ID", "type" and  "target_gene_ID"
#' @examples
#' \dontrun{
#' interac <- interStringency(type = "ncRNA", stringency = "strict")
#' interac_strict_p_50 <- cepair(interac, N_mi= 50)
#' }
#' @export

cepair <- function (interac,N_mi= 50){


#  cat("start searching at ", format(Sys.time(), "%H:%M, %b %d %Y"),
#      " ...\n", sep = "")
m2c <- interac[which(interac$type== "M2C"),]
m2l <- interac[which(interac$type== "M2L"),]
miRNA_1 <- unique(m2c$node_gene_ID)
miRNA_2 <- unique(m2l$node_gene_ID)
circRNA <- unique(m2c$target_gene_ID)
lncRNA <- unique(m2l$target_gene_ID)
miRNA <- unique(c(miRNA_1,miRNA_2))
N = length(miRNA)
#length(circRNA)
#length(lncRNA)
#phyper(q-1, n1, N-n1, n2, lower.tail=F)
#the number and hypergeometric probability of shared miRNAs by a lncRNA-mRNA pair
output <- NULL
for (n1 in lncRNA){
  #everyline <- c()
  #n1="ENSG00000163597"
  #n2="NM_000022"
  mi_n1 <- m2l[which(m2l$target_gene_ID== n1),c("node_gene_ID")]
  mi_n1 <- unique(mi_n1)
  for ( n2 in circRNA){
      mi_n2 <- m2c[which(m2c$target_gene_ID== n2),c("node_gene_ID")]
      mi_n2 <- unique(mi_n2)
      #continue
      q = intersect(mi_n1,mi_n2)
      if(is.null(q)){
        next
      }
      p_value <- phyper(length(q)-1, length(mi_n1), N-length(mi_n1), length(mi_n2), lower.tail=F)
      if (length(q)> N_mi && p_value < 0.05){# length(q)>50 && p_value < 0.05
        output <-  c(output,c(n1,n2))
        output <- unique(output)
      }
  }
}
location <- interac$target_gene_ID %in% output
interac_temp <- interac[location,]

#cat("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"),
#    " ...\n", sep = "")
return(interac_temp)
}
