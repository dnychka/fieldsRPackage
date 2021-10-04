#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2021 Colorado School of Mines
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
profileMLE<- function (obj, parName, parGrid=NULL, gridN=15,
                       cov.params.start=NULL, GCV=FALSE, REML=FALSE,
                       verbose=FALSE ){
  if( class(obj)[1]!= "spatialProcess"){
    stop("only implemented tfor spatialProcess objects")
  }
  if( obj$CASE==0){
    stop("object needs to be have MLE parameters")
  }
  if(is.null(parGrid) ){
  parGrid<- data.frame(obj$summary[parName] *
                         seq( .5,2,length.out=gridN)
                       )
  names(parGrid)<- parName
  }
  
  # remove profile parameter as having a starting value
  
  if( is.null(cov.params.start)){
  cov.params.startTmp<-  as.list(obj$MLEInfo$pars.MLE)
  #cov.params.startTmp<- obj$cov.params.start
  cov.params.startTmp[parName]<- NULL
  }
  else{
    cov.params.startTmp<-cov.params.start
  }
  
# just evaluate other parameters at MLEs and don't optimze. 
# E.g. just use  MLE lambda hat for all choice in in parGrid. 
# if( fast){
#    cov.params.startTmp<-NULL
#    
#  }
#  cov.argsTemp<- obj$cov.argsFull
#  cov.argsTemp[parName]<-  NULL
  
  profileInfo<-  mKrigMLEGrid(obj$x, obj$y,  
               weights = obj$weights,
               Z = obj$Z, 
               mKrig.args = obj$mKrig.args,
               cov.function = obj$cov.function, 
               cov.args  = obj$cov.args,
               par.grid = parGrid, 
               na.rm = TRUE,
               verbose = verbose,
               REML = REML,
               GCV  =  GCV,
               cov.params.start = cov.params.startTmp)
  return( profileInfo)
}