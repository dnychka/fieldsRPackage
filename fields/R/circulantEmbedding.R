#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2024 Colorado School of Mines
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
"circulantEmbedding" <- function(obj) {
    m<- obj$m   
    M<- obj$M
    prodM<- prod( M)
    Z <- fft( array(rnorm(prodM), M) )
    
    if( any( Re(obj$wght) < 0)){
      stop("Weight function has negtive values")
    }
    out<-  Re(
               fft(sqrt(obj$wght) * Z, inverse = TRUE)
             )/sqrt(prodM)
    # this arcania 
    # handles general fields with more than 2 dimensions
    out<- array( out[makeMultiIndex(m)], m)
  
    return( out)
}
