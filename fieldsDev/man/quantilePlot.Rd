
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

\name{quantilePlot}
\alias{quantilePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile plots of general distributions and with uncertainty bounds
}
\description{
This function extends the usual normal QQ plot to other distributions and also
uses Monte Carlo sampling to get an idea of the finite sample uncertainty if the
sample was indeed from the target distribution. Uncertainty is characterized
by a simultaneous coverage envelope around the expected QQ plot. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
quantilePlot(y, quantileFunction = "qnorm", 
quantileArgs = NULL, alpha = 0.05, M = 200,
smoothQuantiles = FALSE, dfSmooth = 80, 
col = "thistle1",
lineCol = "thistle3", cex = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  The data 
%%     ~~Describe \code{y} here~~
}
  \item{quantileFunction}{
  The target distribution in the form of a quantile function. 
%%     ~~Describe \code{quantileFunction} here~~
}
  \item{quantileArgs}{
  Other arguments to be used in the quantile distribution function.
%%     ~~Describe \code{quantileArgs} here~~
}
  \item{alpha}{
  Coverage probability level. 
%%     ~~Describe \code{alpha} here~~
}
  \item{M}{
  Number of Monte Carlo draws to approximate the coverage envelope. 
%%     ~~Describe \code{M} here~~
}
  \item{smoothQuantiles}{
  If TRUE, smooth out the Monte Carlo results to improve the
  plotting. 
%%     ~~Describe \code{smoothQuantiles} here~~
}
  \item{dfSmooth}{
  Effective degrees of freedom to use for the quantile envelope smoothing. 
  Default is the minimum of 80 degrees of freedom or .95* sample size.
%%     ~~Describe \code{dfSmooth} here~~
}
  \item{col}{Fill color of the envelope.
%%     ~~Describe \code{col} here~~
}
  \item{lineCol}{
  Boundary color of the envelope.
%%     ~~Describe \code{lineCol} here~~
}
  \item{cex}{
  Character size for main plot.
%%     ~~Describe \code{cex} here~~
}
  \item{\dots}{
  Any additional arguments to pass to the scatterplot function (\code{plot})
%%     ~~Describe \code{\dots} here~~
}
}
\details{
If p= F(y) is the distribution for a random variable Y, then the quantile 
function is y = Q(p). A QQ plot plots a sorted sample,
\code{YS}, against \code{YTarget} with  \code{YTarget[i]} being 
 Q( (i- (1/2))/n) and n the sample size. 
 The coverage envelope is found by generating  M Monte Carlo samples of size n 
 from F. These are sorted and the mean, m[i], and standard deviation, s[i], are
 found for each 
 i. The coverage envelope is the value of B so that m[i] +/- B*s[i] contains 
 (1-alpha) percent of the samples. The smoothing option just smooths this
 envelope boundary slightly using a cubic smoothing spline
 for a cleaner figure.

}

\author{
Douglas Nychka %%  ~~who you are~~
}
\note{
We are not sure why this kind of function, with the uncertainty bounds, was not
part of base R. QQ plots have all kinds of hidden dependencies, and it is useful to 
take this into account. 
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{qqnorm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed( 123)
n<- 100
y<- rnorm(n)
quantilePlot( y)
# note the width of envelope at the ends!

set.seed( 222)
yE<- rexp(n)
# standardize
yEStandardized<- (yE- mean( yE))/ sd( yE)
quantilePlot( yEStandardized)
# suggests yE is not normally distributed N(0,1)

quantilePlot( yE,quantileFunction = "qexp" ,M=200)
# now within bounds.

quantilePlot( y, smoothQuantiles=TRUE, dfSmooth=20)
# removes some kinks

}

\keyword{hplot} 

