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

mKrigFastPredictSetup <- function(mKrigObject,
                          gridList, 
                          NNSize=4,
                          giveWarnings=TRUE
                                  ) {
  # the main reason to pass new args to the covariance is to increase
  # the temp space size for sparse multiplications
  # other optional arguments that typically describe the covariance function 
  # from mKrig are passed along in the list object$args
   np<- NNSize
   xObs<- mKrigObject$x
  
  # adjust grid if needed to include a margin of np+1 grid points beyond xObs
   
    marginInfo<- addMarginsGridList(xObs, gridList, np)

  # these are the slightly larger grids by adding margins.
    gridListNew<-  marginInfo$gridListNew
    offGridObject<- offGridWeights( xObs,
                                    gridListNew,
                                    mKrigObject = mKrigObject,
                                    np=np,
                                    giveWarnings = giveWarnings
    )
    cov.obj<- stationary.image.cov( setup=TRUE, 
                                    grid=gridListNew,
                                    cov.args= mKrigObject$args)
    return( 
       list(offGridObject = offGridObject,
                  cov.obj = cov.obj,
               marginInfo = marginInfo
            )
          )
  }
