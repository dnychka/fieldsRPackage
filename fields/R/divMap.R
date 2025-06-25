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


