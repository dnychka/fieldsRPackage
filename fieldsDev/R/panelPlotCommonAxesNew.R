#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2026 Colorado School of Mines
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
# this is a script to create an integrated panel of plots using just the
# graphics primitives in base R graphics. ggplot and other packages may do 
# this as well but often it is difficult to modify the choices these higher 
# level panel functions make. Inevitably for a finished graphic I will need
# to tweak something and it is just as easy at that point to 
# build the complex plots from simple elements using base R. 
# DWN nychka@mines.edu

panelPlotCommonAxes<- function(m, n, xr, yr,
                               plotCode,
                               xlabels= "",
                               ylabels="",
                               outerXLabel=NULL,
                               outerYLabel=NULL,
                               PANESPACE=0,
                               LEGENDSPACE=c(0,0,0,0),
                               axes=TRUE,
                               addAxis=c(TRUE, TRUE, TRUE, TRUE),
                               alternateAxes=TRUE,
                               axisMargin=c( 3.1,3.1,3.1,3.1),
                               axisList=list( line=0),
                               axisLabelList=list(line=2),
                               outerLabelMargin= c(3.1,3.1,0,0),
                               outerLabelList=list(line=2),
                               colBg="grey90",
                               colGrid="grey98",
                               boxLwd= 1,
                               titleBoxHeight= 0,
                               titleBoxList= list( 
                               boxColor = "lightgreen",
                               boxBorder= "black",
                               lwd=boxLwd),
                               verbose=FALSE
                               ){
  
########################################################
### case specific
#########################################################

# following are panel specific here  this is a panel of plots with
# these are the axes ranges -- set these for
# each application 
# Here the ranges are just repeated if given singly
  
  if( !is.matrix( xr)){
    xr<- matrix( xr, n,2, byrow=TRUE)
  }
  if( !is.matrix( yr)){
    yr<- matrix( yr, m,2, byrow=TRUE)
  }
  if( length( xlabels)==1){
    xlabels<- rep( xlabels, n)
  }
  if( length( ylabels)==1){
    ylabels<- rep( ylabels, m)
  }
  if( verbose){
    cat("x ranges", fill=TRUE)
    print(xr)
    cat("x labels", fill=TRUE)
    print( xlabels)
    cat("y ranges", fill=TRUE)
    print( yr)
    cat("y labels", fill=TRUE)
    print(ylabels)
  }

# If you plan to add legends or colorbars 
# after drawing the plot use the LEGENDSPACE
# argument to reserve extra space around the 
# plot to do this.
# and see the function addLegend to 
# reset to make the plotting easier for the legends
# in [0,1]X [0,1] figure/user coordinates.
# You might  have to tinker with the legend spacing in most cases 
# we just do this by trial and error!
# set LEGENDSPACE  to zero -- the default if you are not including a legend outside of the 
# panel 

# Each expression is the graphics code for one of the plots in the panel. 
# if you are plotting data objects either pass them to this function or expect
#  them to be
# in your workspace.
# The code for each plot can be any graphics function or R code. 
# Below are just some simple examples 
# you can also hard code  plot function  into the loop below but this 
# makes the coding more complex.
# Note that some of these example plots use  multiple plot commands. 

########################################################
###  The following  stuff should not  need to be changed 
#########################################################


# nonspecific settings
# laying out canvas, panel, and each plot area (pane)
# Note that some extra room is added on top if a 
# title box is added.
# units are in character size (aka line argument in the functions)
#
par( 
mfrow=c(m,n),
mar=c(PANESPACE, PANESPACE, titleBoxHeight + PANESPACE, 0)
)
 
  
# outer label real estate
addOuterLabels<- ! (is.null(outerXLabel) &
                       is.null(outerYLabel))
outerLabelMargin<- c(
   ifelse(!is.null(outerXLabel), outerLabelMargin[1], 0),
   ifelse(!is.null(outerYLabel), outerLabelMargin[2], 0),
   ifelse(!is.null(outerXLabel), outerLabelMargin[3], 0),
   ifelse(!is.null(outerYLabel), outerLabelMargin[4], 0)
)

# Note these are vectors of 4 spacings
# c(bottom, left, top,right)

outerMargin<- outerLabelMargin + LEGENDSPACE + axisMargin +
  c(0,0, PANESPACE+ titleBoxHeight ,PANESPACE )
# set outer margins to make room for everything outside of the panel.
# this is pretty heavy-handed and the ggplot folks will cringe at this.
# However, it is obvious what is happening -- and  just base R graphics.
# panel spacing only added to the left and bottom sides so
# add at the end along top and right 

if( verbose){
  cat("legend margins", fill=TRUE)
  print( LEGENDSPACE)
  cat("outer label  margins", fill=TRUE)
  print( outerLabelMargin)
  cat("axis margins", fill=TRUE)
  print( axisMargin)
  cat("oma used in par", fill=TRUE)
  print(outerMargin )
}

par(oma=outerMargin)
#
########################################################
###  loop over plots this is done row by row
#########################################################

counter<- 0
for( I in 1:m ){
  for( J in 1:n){
    counter<- counter+1
    #blank canvas for plot I,J
    plot( xr[J,],yr[I,], xaxt = "n",
          yaxt = "n",type="n",
          bty="n")
  
    if( length(plotCode ) >1){
    thePlot<- plotCode[[counter]]
    }
    else{
      thePlot<- plotCode[[1]]
    }
    
    if(verbose){
      cat(" ",fill=TRUE )
      cat("counter: ",counter,  fill=TRUE )
      cat("I,J: ", I,J,      fill=TRUE )
      cat("Blank plot with  x range: ",xr[J,], fill=TRUE )
      cat("Blank plot with y range: ", yr[I,], fill=TRUE )
      cat(" the plot code for this pane: ",fill=TRUE)
      cat(as.character(thePlot), fill=TRUE)
    }

# skip this pane if NULL
    if(!is.null(thePlot)){
    # some default background setup for the I,J plot ...
    # draw background first
    xy<- par()$usr
    # here is a handy function for the grid we missed this function for about 25 years!
   
    rect(xy[1],xy[3],xy[2],xy[4],col=colBg)
    grid(col=colGrid, lty=1, lwd=.5)
    box(lwd=boxLwd )
    #
    
    
    # now add the title box with  a  background color 
    if( titleBoxHeight>0){
      do.call( "addTitleBox", c(list( boxHeight=titleBoxHeight),
                                titleBoxList)
                )
    }
    # box around the plot rectangle
    box(lwd=boxLwd)
    # here the plot code for plot I,J is run
    eval(thePlot)
    #
    }
########################################################
### Add axes 
########################################################
  # copy passed list for axis so it can be modified for a 
  # specific pane ( e.g. top row)
 
    
    if(axes){
      # add appropriate axis and axis label if axes is TRUE 
      # but only add axes  to  the plots in the  outer edges  
      # add more options to these axis functions  as needed through the list arguments
      # e.g. cex to increase text size 
    # figure out what to do 
     ### Top 
      addYAxis<- FALSE
      addXAxis<- FALSE
      axisListTmp <- axisList
      axisLabelListTmp<- axisLabelList
    if( I == 1 & addAxis[3]){ 
     
      addYAxis<- TRUE
      tickmarkLabel<- ifelse(alternateAxes, !(J%%2==0), TRUE)
      axisLabelTmp<- xlabels[J]
      sideTmp <- 3
      # add extra space due to the title box
       axisListTmp$line <- axisListTmp$line  + titleBoxHeight
       axisLabelListTmp$line<- axisLabelListTmp$line + titleBoxHeight
    }
    ### Bottom
    if( I == m & addAxis[1]){ 
      addYAxis<- TRUE
      tickmarkLabel<- ifelse(alternateAxes, !(J%%2==0), TRUE)
      axisLabelTmp<- xlabels[J]
      sideTmp <- 1
    }
      if( addYAxis){
        # axis with/without the labels
        do.call( "axis",c( 
           list( side=sideTmp), 
                  axisListTmp,
                 labels=tickmarkLabel)
        )
        # text along axis 
        do.call("mtext",c( list(text=axisLabelTmp,
                                side=sideTmp),
                           axisLabelListTmp)
        )
      }
    # reset the temporary list 
    axisListTmp <- axisList
    axisLabelTmp<- axisLabelList
    ################################################################
    ### Left 
    if( J == 1 & addAxis[2]){ 
      addXAxis<- TRUE
      tickmarkLabel<- ifelse(alternateAxes, !(I%%2==0), TRUE)
      sideTmp<- 2
      axisLabelTmp<- ylabels[I]
    }
    ### Right 
    if( J == n & addAxis[4]){ 
      addXAxis<- TRUE
      tickmarkLabel<- ifelse(alternateAxes, !(I%%2==0), TRUE)
      sideTmp<- 4
      axisLabelTmp<- ylabels[I]
    }
      if( addXAxis){
        
        # axis with/without the labels
        do.call( "axis",c( 
          list( side=sideTmp), 
                axisListTmp,
              labels=tickmarkLabel)
                 )
        # text along axis  
        do.call("mtext",c( list(text=axisLabelTmp,
                                side=sideTmp),
                           axisLabelList)
                )
      }
    ################################################################
      }
# 
}

# done looping through the plots now add stuff outside the panel in the 
# margins. 

########################################################
### Add overall axis titles
########################################################
# margin labels added on bottom and left edges  
# modify outerLabelList if you don't like the label placement/style
#
if( addOuterLabels){
  if(!is.null(outerXLabel)){
  do.call("mtext", 
      c(list( outer=TRUE,side=1, text= outerXLabel),outerLabelList)
      )
  }
  if(!is.null(outerYLabel)){
    do.call("mtext", 
      c(list( outer=TRUE, side=2,text= outerYLabel),outerLabelList)
      )
   
  }
 
}
# all done see comments above about LEGENDSPACE to add a color bar or legend.
}
}



