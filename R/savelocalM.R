#' @title save and plot module with local method
#'
#' @description save and plot module for the method local result in the function runmodule
#'
#' @param res.list_local the methd local result

#' @return the plot and the format XGMML of the each module, filenames is the seed node
#' @examples
#' \dontrun{
#' data("dataN")
#' gene2weight <- combinp(dataN[,c("type","logFC","PValue")])
#' interac <- interStringency(type = "transcription", spec ="hg",
#' 				stringency = "strict")
#' interac <- interac[,c("node_gene_ID","type","target_gene_ID")]
#' res.list_local <- runmodule(network = interac, gene2weight,
#'  method = "local",maxsize=15, seletN = "MIMAT0000461")
#' savelocalM(res.list_local)
#' }
#'
#' @export

savelocalM <- function(res.list_local){
    list_local <- res.list_local$module
   Map(function(x, i) {
    list_n <- as.character(x)
    node_plot <- subNetwork_only(list_n, res.list_local$GNCW)
    saveNetwork(node_plot, file= i, type = "XGMML")
    png(filename = paste0(i, ".jpg"), width = 2400, height = 1800, res = 200)
    try(plotSub(node_plot))
    dev.off()
  }, list_local , names(list_local))
}
