%#
%# fields  is a package for analysis of spatial data written for
%# the R software environment.
%# Copyright (C) 2024 Colorado School of Mines
%# 1500 Illinois St., Golden, CO 80401
%# Contact: Douglas Nychka,  douglasnychka@gmail.edu,
%#
%# This program is free software; you can redistribute it and/or modify
%# it under the terms of the GNU General Public License as published by
%# the Free Software Foundation; either version 2 of the License, or
%# (at your option) any later version.
%# This program is distributed in the hope that it will be useful,
%# but WITHOUT ANY WARRANTY; without even the implied warranty of
%# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%# GNU General Public License for more details.
%#
%# You should have received a copy of the GNU General Public License
%# along with the R software environment if not, write to the Free Software
%# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%# or see http://www.r-project.org/Licenses/GPL-2
%##END HEADER
%##END HEADER
\name{LENSExample}
\alias{LENSExample}

\docType{data}
\title{
Surface temperature sensitivity from the NCAR Climate model.
}
\description{
Based on  the NCAR Earth System Model (CESM). 
the local change in surface temperature (C) expected by a 1 degree change in the global average. 
This statistic is referred as the climate sensitivity and was found by 
OLS regression on a 30 member ensemble of climate runs. 
}
%\usage{data(NorthAmericanRainfall)}
\format{
 The usual image format with components \code{x}, \code{y} and \code{z}.
 
\code{x} Equally spaced longitude values.

\code{y} Equally spaced latitude values.

\code{z} Climate sensitivity on a 288X192 grid. 

}

\source{
The CESM Large Ensemble Project, led by Clara Deser and Jennifer Kay,
is a publicly available set of climate model simulations intended for 
understanding  variability in the Earth's climate system and climate change.
See \url{https://www.cesm.ucar.edu/community-projects/lens} for background. 

All simulations are performed with the nominal 1-degree latitude/longitude version 
of the Community Earth System Model version 1 (CESM1) with CAM5.2 as its atmospheric 
component. Each member is subject to the same radiative forcing scenario
historical up to 2005 and RCP8.5 thereafter. Here the main component of these scenarios are 
the observed or projected greenhouse gas emisions leading to global warming. 
Each climate run begins from a slightly different 
initial atmospheric state, created by randomly perturbing temperatures at the level
of round-off error. Even this subtle difference due to nonlinear dynamics 
leads to climate simulations that can be assumed to be independent realizations. 

Here 30 members of the ensemble are summarized to give local sensitivity for the Northern Hemisphere
summer months, June, July, and August (aka JJA). 
This involves an linear regression on the summer mean time series for each model grid box with the
global summer average and for each ensemble member. The 30 slopes from the regressions for each 
grid box are averaged to get the sensitivity in this data set. 

See the fields package github repository for the R script that creates 
this data set from the model output. 

}
\examples{
data(LENSExample)
imagePlot( LENSExample)

# Nathan Lenssen's version
 
divMap(LENSExample$x, LENSExample$y, LENSExample$z,
      worldCol="white", lwd=.5, zCap=2, wrapVal=c(-180, 180))
}

\keyword{datasets}
