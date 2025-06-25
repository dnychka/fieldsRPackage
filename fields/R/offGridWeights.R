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
offGridWeights<-function(s, gridList, NNSize=2,
                         mKrigObject=NULL, 
                         Covariance=NULL, 
                         covArgs=NULL,
                         aRange=NULL, sigma2=NULL, 
                         giveWarnings=TRUE,
                         debug=FALSE,
                         verbose=FALSE
                   )
  #
  #  This is a wrapper function to distinguish between the 
  # 1D and 2D cases. The 1D is useful to follow the algorithm 
  # details but typically not needed in practice. The 2D is 
  # the workhorse for fast approximation simulation of a stationary
  # spatial process on a regular grid combined with irregular 
  # observations
  # 
  {
  

   if( length( gridList)==1){
     out<- offGridWeights1D( s=s, gridList=gridList, NNSize=NNSize,
                      mKrigObject=mKrigObject, 
                      Covariance=Covariance,
                      covArgs=covArgs,
                      aRange=aRange, 
                      sigma2=sigma2, 
                      giveWarnings=giveWarnings,
                      debug=debug,
                      verbose=verbose)
     out$gridList<- gridList
  return( out)
    }
  # 
  if( length( gridList)==2){
     marginInfo<- addMarginsGridList(s, gridList, NNSize)
    
    if(marginInfo$gridExpand){
     
    warning("gridList expanded to include enough 
            nearest neighbors")
     if( verbose){
       cat("NNSize", NNSize, fill=TRUE)
       cat(  "x min/max",
             min(s[,1]), max( s[,2]), fill=TRUE)
       cat("Old grid", gridList$x,fill=TRUE)
       cat("New grid", marginInfo$gridListNew$x,fill=TRUE)
       cat(  "y min/max",
             min(s[,2]), max( s[,2]), fill=TRUE)
       cat("Old grid", gridList$y,fill=TRUE)
       cat("New grid", marginInfo$gridListNew$y,fill=TRUE)
     }
      
      gridList<- marginInfo$gridListNew
      }
  
    out<- offGridWeights2D( s=s, gridList=gridList, NNSize=NNSize,
                     mKrigObject=mKrigObject, 
                     Covariance=Covariance,
                     covArgs=covArgs,
                     aRange=aRange, 
                     sigma2=sigma2, 
                     giveWarnings=giveWarnings,
                     debug=debug,
                     verbose=verbose)
    
    return( out )
  }
  
   stop("GridList should be either 1 or 2 components")
  
  }
