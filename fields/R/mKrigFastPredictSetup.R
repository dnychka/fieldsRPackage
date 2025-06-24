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

mKrigFastPredictSetup <- function(mKrigObject,
                          gridList, 
                          NNSize=4,
                          giveWarnings=TRUE,
                          verbose=FALSE
                                  ) {
  # the main reason to pass new args to the covariance is to increase
  # the temp space size for sparse multiplications
  # other optional arguments that typically describe the covariance function 
  # from mKrig are passed along in the list object$args
   xObs<- mKrigObject$x
   if(ncol( xObs) !=2){
     cat( " dim s", ncol( xObs), fill=TRUE )
     stop("fast predict only implemented for 2 D")
   }
  
  # adjust grid if needed to include a margin of np+1 grid points beyond xObs
   
    marginInfo<- addMarginsGridList(xObs, gridList, NNSize)

  # these are the slightly larger grids by adding margins.
    gridListNew<-  marginInfo$gridListNew
    approxGridObject<- approximateCovariance2D( xObs,
                                    gridListNew,
                                mKrigObject = mKrigObject,
                                         np = NNSize,
                               giveWarnings = giveWarnings,
                               verbose=verbose)
                                   
    cov.obj<- stationary.image.cov( setup=TRUE, 
                                    grid=gridListNew, 
                                    cov.function=mKrigObject$cov.function.name,
                                    cov.args= mKrigObject$args)
    return( 
       list(offGridObject =approxGridObject,
                  cov.obj = cov.obj,
               marginInfo = marginInfo
            )
          )
  }
