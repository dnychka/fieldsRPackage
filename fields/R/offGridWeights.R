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
##END HEADER\
offGridWeights<-function(s, gridList, np=2,
                         mKrigObject=NULL, 
                         Covariance=NULL, covArgs=NULL,
                         aRange=NULL, sigma2=NULL, 
                         giveWarnings=TRUE,
                         debug=FALSE
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
  callList<- list( s=s, gridList=gridList, np=np,
                   mKrigObject=mKrigObject, 
                   Covariance=Covariance,
                   covArgs=covArgs,
                   aRange=aRange, 
                   sigma2=sigma2, 
                   giveWarnings=giveWarnings,
                   debug=debug )
  
  if( length( gridList)==1){
    # always find prediction variances  for 1D.
     return( do.call("offGridWeights1D", callList) )
   }
  
  if( length( gridList)==2){
    return( do.call("offGridWeights2D", callList) )
  }
  
   stop("GridList should be either 1 or 2 components")
  
  }
