# bubblePlot fix

bubblePlot<- function (x, y, z, col = viridisLite::viridis(256), zlim = NULL, 
                            horizontal = FALSE, legend.cex = 1, legend.lab = NULL, legend.line = 2, 
                            legend.shrink = 0.9, legend.width = 1.2, legend.mar = ifelse(horizontal, 
                                                                                         3.1, 5.1), axis.args = NULL, legend.args = NULL, size = 1, 
                            add = FALSE, legendLayout = NULL, highlight = FALSE, highlight.color = "grey30", 
                            bubbleType="circle",
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
  
  # only square and circle supported!
  if( bubbleType=="square"){
  pchSolid<- 15
  pchOutline<- 0
  }
  else{
    pchSolid<- 16
      pchOutline<- 1
  }
    
    
  if (!add) {
    plot(x, y, col = ctab, cex = size, pch = pchSolid, ...)
  }
  else {
    points(x, y, cex = size, col = ctab, pch = pchSolid)
  }
  if (highlight) {
    points(x, y, cex = size,
           pch = pchOutline, col = highlight.color)
  }
  # save the graphical parameter so they can be reset
  big.par <- par(no.readonly = TRUE)
  mfg.save <- par()$mfg
  # figure out some legend parameters if they are not already
  # specified
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
