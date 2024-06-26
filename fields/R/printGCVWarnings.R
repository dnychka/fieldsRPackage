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
printGCVWarnings<- function( Table, method="all"){
 ind<- Table$Warning
 if( method=="user"){
   #print nothing, expect that the user knows what he/she is doing!
  return()
 }
 if( method == "all"){
 	kIndex<- 1:6
 }
 else{
    kIndex<- match( method,c("GCV",
                          "GCV.model",
                          "GCV.one",
                          "RMSE",
                          "pure error",
                         "REML")
                         
                   )
      }
 methodList<- c(
 "(GCV) Generalized Cross-Validation ",
 "(GCV.model) Generalized Cross-Validation on replicate means ",
 "(GCV.one) Generalized Cross-Validation on individual observations ",
 "(RMSE) Matching estimate of tau to supplied rmse ",
 "Matching estimate of tau to that from replicated observations",
 "(REML) Restricted maximum likelihood ",
 "user supplied lambda"
  )
 if( any( ind[kIndex])){
 	cat("Warning: ", fill=TRUE)
 	cat("Grid searches over lambda (nugget and sill variances) with  minima at the endpoints: ", fill=TRUE) }
 for( k in kIndex){
 if( ind[k]){
    whichEnd<- ifelse(Table[k,2],"left","right")
 	cat( " ", methodList[k], fill =TRUE)
 	cat( "   minimum at ", whichEnd, "endpoint",
 	                " lambda  = ", Table[k,6] ,
    "(eff. df=", Table[k,7] , ")", fill = TRUE )
 	     }
 }	     
 
}
