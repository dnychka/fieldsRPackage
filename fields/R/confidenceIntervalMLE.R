#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2022 Colorado School of Mines
# 1500 Illinois St., Golden, CO 80401
# Contact: Douglas Nychka,  douglasnychka@gmail.edu,
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
confidenceIntervalMLE<- function( obj, CILevel, verbose=FALSE){
   if( is.na( CILevel)){
     cat("CI not found")
     return(NULL)
   }
  MLEInfo<- obj$MLEInfo
  transformedMLEPars<- MLEInfo$optimResults$par
  sampleFisherInfo<- solve(-1*MLEInfo$optimResults$hessian)
  if( verbose){
    print( sampleFisherInfo)
  }
  testD<- diag(sampleFisherInfo)
  
  if( any( testD<=0)){
    SE<- rep( NA, length( testD) )
  }
  else{
  SE<-  sqrt( diag(sampleFisherInfo ))
  }
  
  CITmp1<- transformedMLEPars - qnorm(CILevel)*SE
  CITmp2<- transformedMLEPars + qnorm(CILevel)*SE
  tableCI<- 
    cbind(
      MLEInfo$parTransform( CITmp1, inv=TRUE),
      MLEInfo$parTransform( CITmp2, inv=TRUE)
    )
  percentCI<- paste0(100* CILevel,"%")
   colnames(tableCI ) <- c(
    paste0("lower",percentCI), 
    paste0("upper",percentCI)
  )
  return( tableCI)
}