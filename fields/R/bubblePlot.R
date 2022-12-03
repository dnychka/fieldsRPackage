# bubblePlot fix

bubblePlot<- function (x, y, z, col = viridis::viridis(256), zlim = NULL, 
                            horizontal = FALSE, legend.cex = 1, legend.lab = NULL, legend.line = 2, 
                            legend.shrink = 0.9, legend.width = 1.2, legend.mar = ifelse(horizontal, 
                                                                                         3.1, 5.1), axis.args = NULL, legend.args = NULL, size = 1, 
                            add = FALSE, legendLayout = NULL, highlight = TRUE, highlight.color = "grey30", 
                            ...) 
{
  x <- as.matrix(x)
  if (dim(x)[2] == 2) {
    z <- y
    y <- x[, 2]
    x <- x[, 1]
  }
  ctab = color.scale(z, col, zlim = zlim)
  if (!add & is.null(legendLayout)) {
    legendLayout <- setupLegend(legend.shrink = legend.shrink, 
                                legend.width = legend.width, legend.mar = legend.mar, 
                                horizontal = horizontal)
  }
  if (!add) {
    plot(x, y, col = ctab, cex = size, pch = 16, ...)
  }
  else {
    points(x, y, cex = size, col = ctab, pch = 16)
  }
  if (highlight) {
    points(x, y, cex = size, col = highlight.color)
  }
  big.par <- par(no.readonly = TRUE)
  mfg.save <- par()$mfg
  if (!add | !is.null(legendLayout)) {
    levelsZ <- attr(ctab, "levelsZ")
    if ((is.null(axis.args)) & (!is.null(levelsZ))) {
      axis.args = list(at = 1:length(levelsZ), labels = levelsZ)
    }
    addLegend(legendLayout, col = attr(ctab, "col"), zlim = attr(ctab, 
                                                                 "zlim"), axis.args = axis.args, legend.args = legend.args, 
              legend.cex = legend.cex, legend.lab = legend.lab, 
              legend.line = legend.line)
  }
  par(plt = big.par$plt, xpd = FALSE)
  par(mfg = mfg.save, new = FALSE)
}
