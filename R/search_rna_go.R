search_rna_go <- function(RNA){
  #circRNA <- match.arg(RNA)
  load("data/rna2gene_GO.rda")
  if(! RNA %in% rna2gene_GO$RNA_ID){
    cat("The RNA has no corresponding GO function information ")
    return("NA")
  }
  infer_go <- rna2gene_GO[which(rna2gene_GO$RNA_ID == RNA),]
  return(infer_go)
}

