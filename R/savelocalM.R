#' @title save and plot module
#'
#' @description save and plot module for the method local result in the function runmodule
#'
#' @param res.list_local the methd local result

#' @return the plot and the format XGMML of the each module, filenames is the seed node

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
