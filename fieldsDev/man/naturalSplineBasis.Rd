%#
%# fields  is a package for analysis of spatial data written for
%# the R software environment.
%# Copyright (C) 2026 Colorado School of Mines
%# 1500 Illinois St., Golden, CO 80401
%# Contact: Douglas Nychka,  douglasnychka@gmail.com,
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
\name{naturalSplineBasis}
\alias{naturalSplineBasis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Natural spline basis.
}
\description{
This simple function relies on the splines package and implements 
de Boor's algorithm to create a B-spline basis. The "natural boundary"
conditions
constrain higher order derivatives to be zero at the endpoints.
}
\usage{
naturalSplineBasis(sGrid, sKnots, degree = 3, derivative = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sGrid}{
  Values to evaluate the basis functions.
%%     ~~Describe \code{sGrid} here~~
}
  \item{sKnots}{
  Knot points. There will be as many basis functions as knot points, and they
  must be unique.
%%     ~~Describe \code{sKnots} here~~
}
  \item{degree}{
  Degree of spline; \code{degree = 3} gives a basis that is piecewise
  cubic polynomial at the knots. 
%%     ~~Describe \code{degree} here~~
}
  \item{derivative}{
  Order of derivative to evaluate at \code{sGrid}. Default is to 
  evaluate the basis functions.
%%     ~~Describe \code{derivative} here~~
}
}
\details{
This implements the famous and elegant De Boor algorithm to find the basis
functions. Following De Boor's analysis, natural boundary conditions are 
enforced by repeating knots
at the boundaries. With the natural boundary conditions, there are as many basis 
functions as knots. Although piecewise polynomials with the segments defined by
the knots, the functions are constrained to be smooth. For the cubic case, all 
basis functions will be continuous and have continuous first and second 
derivatives. At the end points (first and last knot), the functions will have 
third and second derivatives equal to zero. Thus, extrapolation beyond the 
end points is linear. For cubic natural splines, the support is over 5
knot values
away from the boundaries. 

The natural spline basis is also the basis for a smoothing spline and is 
distinct from more arbitrary choices such as regression splines. 

For the curious reader, one can back out the cubic coefficients for this basis 
by examining the derivatives in each interval -- and using your calculus!

%%  ~~ If necessary, more details than the description above ~~
}
\value{
If m is the length of \code{sGrid} and n is the length of \code{sKnots}, the 
returned value is an m X n matrix with i,j element being the j-th basis function
evaluated at \code{sGrid[i]}.

}

\author{
Douglas Nychka
%%  ~~who you are~~
}

\seealso{
\code{\link{Tps}}, \code{\link{Tps.cov}}, \code{\link{sreg}},
\code{\link{splint}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library( splines)
sGrid<- seq( 0,1, length.out=200)
sKnots<- seq( 0,1, length.out=8) # 8 knot points
BasisMatrix<- naturalSplineBasis(sGrid, sKnots)
fields.style()
matplot(sGrid, BasisMatrix, lty=1, type="l" )
rug( sKnots, lwd=3)
# Note symmetry but also different shapes near endpoints. 
}

\keyword{smooth} 

