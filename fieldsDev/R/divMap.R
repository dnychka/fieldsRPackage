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
divMap <- function(x,y,z,
                   zForce=NULL, zCap=NULL, col= NULL,
                   worldCol="magenta1",lty=1, lwd=1, 
                   brewerPal = 'Spectral',
                   rev=TRUE, 
                   wrapVal=c(0,360),...){
  #
  # Nathan's preferred plot of geophysical fields 
  #
	if(!is.null(zCap)){
		z[which(z>zCap)] <- zCap
		z[which(z < -zCap)] <- -zCap
	}
	zMax <- max(abs(z),zForce, zCap, na.rm=T)
	zr <- c(-zMax,zMax)

	if(is.null(col)){
            pal <- designer.colors(256,
                                   RColorBrewer::brewer.pal(11,brewerPal)
                                   )
	}

	if(rev) pal <- rev(pal)
	
	imagePlot(x,y, matrix(z,length(x),length(y)),
	          xlab='',ylab='', zlim=zr, col=pal, ...)
	world(add=T,wrap=wrapVal, col=worldCol, lty=lty, lwd=lwd)
}


