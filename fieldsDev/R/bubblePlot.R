#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2026 Colorado School of Mines
# 1500 Illinois St., Golden, CO 80401
# Contact: Douglas Nychka,  douglasnychka@gmail.com,
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2
##END HEADER
# bubblePlot fix

bubblePlot<- function (x, y, z, col = viridisLite::viridis(256),
                       zlim = NULL, noLegend=FALSE,
                            horizontal = FALSE, 
                         legend.cex = 1, legend.lab = NULL,
                       legend.line = 2, 
                            legend.shrink = 0.8, legend.width = 1.2, 
                       legend.mar = ifelse(horizontal, 3.1, 5.1),
                       axis.args = NULL, legend.args = NULL, size = 1, 
                            add = FALSE,  highlight = FALSE,
                            highlight.color = "grey30", 
                            bubbleType="circle", 
                            ...) 
{
  x <- as.matrix(x)
  if (dim(x)[2] == 2) {
    z <- y
    y <- x[, 2]
    x <- x[, 1]
  }
  if( is.null( zlim)){
    # coerce to numeric in case 
    # this is from a table (see help example)
    zlim<- range( as.numeric(z), na.rm =TRUE)
  }
  ctab = color.scale(z, col, zlim = zlim)
  
  # only square and circle supported!
  if( bubbleType=="square"){
  pchSolid<- 15
  pchOutline<- 0
  }
  else{
    pchSolid<- 16
      pchOutline<- 1
  }
    
  # only change margin if not adding and need a legend ...
  if (!add & !noLegend) {
    oldMar<- par()$mar
    # this reserves some room in outer margin for the color bar.
    LEGENDSPACE<- legend.width + legend.mar # in units of characters
 if( horizontal){
      newMar<-c(LEGENDSPACE,0,0,0 ) + oldMar
 }
else{
  newMar<-c(0,0,0,LEGENDSPACE ) + oldMar
}
    par(mar=newMar)
  }
  
  if( !add){ 
    # yep this simple!
    plot(x, y, col = ctab, cex = size, pch = pchSolid, ...)
  }
  else {
    points(x, y, cex = size, col = ctab, pch = pchSolid)
  }
  
  if (highlight) {
    points(x, y, cex = size,
           pch = pchOutline, col = highlight.color)
  }
  
  if ( noLegend==FALSE) {
    imagePlot( legend.only=TRUE,
      col=col, zlim = zlim,horizontal=horizontal,
      axis.args = axis.args, legend.args = legend.args, 
      legend.cex = legend.cex, legend.lab = legend.lab, 
      legend.line = legend.line, legend.width=legend.width)
  # 
  # imagePlot when it exits automatically resets the graphical parameters
  # to the original plot.
  }
  
}
