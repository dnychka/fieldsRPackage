divMap <- function(x,y,z,zForce=NULL, zCap=NULL, pal= NULL,
                   brewerPal = 'Spectral', rev=TRUE,wrapVal=c(0,360),...){
	require(fields)
	require(RColorBrewer)

	if(!is.null(zCap)){
		z[which(z>zCap)] <- zCap
		z[which(z < -zCap)] <- -zCap
	}


	zMax <- max(abs(z),zForce,zCap,na.rm=T)
	zr <- c(-zMax,zMax)

	if(is.null(pal)){
            pal <- designer.colors(256,
                                   RColorBrewer::brewer.pal(11,brewerPal)
                                   )
	}

	if(rev) pal <- rev(pal)

	worldMap(x,y,z,zlim=zr,col=pal,wrapVal=wrapVal,...)
}

worldMap <- function(x,y,z, wrapVal=c(0,360),...){
	require(fields)

	imagePlot(x,y,matrix(z,length(x),length(y)),
		xlab='',ylab='',...)
	world(add=T,wrap=wrapVal)
}
