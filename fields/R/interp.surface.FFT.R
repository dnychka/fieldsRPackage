
interp.surface.FFT<-function( obj, M){
#a
  
z<- obj$z
m1<- nrow( z)
m2<- ncol( z)
#Require that  m1 and m2 are odd
if( (m1%%2 != 1) | (m2%%2 != 1)){
  cat( m1, m2, fill=TRUE)
  stop("Need an odd number of grid points for FFT interp method")
}

zf<- fft( z)
# refined grid sizes
N1<- M*m1
N2<- M*m2
temp<- matrix( 0, N1,N2)
# this stuffing is touchy and that is why the odd
# number of grid values is needed. 
n1<- (m1-1)/2
n2<- (m2-1)/2

temp[1:(n1+1),       1:(n2+1)]<- zf[ 1:(n1+1),  1:(n2+1)]
temp[1:(n1+1),    N2-(n2-1):0]<- zf[ 1:(n1+1), (n2+2):m2]
temp[N1-(n1-1):0,    1:(n2+1)]<- zf[(n1+2):m1,  1:(n2+1)]
temp[N1-(n1-1):0, N2-(n2-1):0]<- zf[(n1+2):m1, (n2+2):m2]

look<- Re(fft( temp, inverse=TRUE))/ (m1*m2)
# trim interpolated array to be within bounds of source image
look<- look[1: (N1 - M +1), 1: (N2 - M +1) ]

xOut<- seq(obj$x[1], obj$x[m1], length.out=(N1 - M +1) ) 
yOut<- seq(obj$y[1], obj$y[m2], length.out=(N2 - M +1) )

out<- list( x= xOut, y=yOut, z=look )

return(out)

}
  
