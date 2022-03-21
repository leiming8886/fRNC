
node.color <- function (network, colors)
{
  colors <- colors[V(network)$name]
  colors2 <- colors
  if (max(abs(colors)) < 5) {
    colors <- colors * 5
  }
  if (any(colors > 0)) {
    max.red <- max(ceiling(abs(colors[which(colors > 0)])))
    reds <- colorRampPalette(colors = c("white", "red"))
    red.vec <- reds(max.red)
    colors2[which(colors > 0)] <- red.vec[ceiling(abs(colors[which(colors >
                                                                     0)]))]
  }
  if (any(colors < 0)) {
    max.green <- max(ceiling(abs(colors[which(colors < 0)])))
    greens <- colorRampPalette(colors = c("white",
                                          "green"))
    green.vec <- greens(max.green)
    colors2[which(colors < 0)] <- green.vec[ceiling(abs(colors[which(colors <
                                                                       0)]))]
  }
  return(colors2)
}



