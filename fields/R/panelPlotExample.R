# this is a script to create an integrated panel of plots using just the
# graphics primitives in base R graphics. ggplot and other packages may do 
# this as well but often it is difficult to modify the choices these higher 
# level panel functions make. Inevitably for a finished graphics I will need
# to tweak something and it is just as easy at that point to 
# build the complex plots from simple elements using base R. 
# DWN nychka@mines.edu

panelPlotExample<- function(){
  
########################################################
### case specific
#########################################################

# following are panel specific here  this is a panel of plots with
# 3 rows  and 2 columns
  m<-3
  n<- 2
  theXLabel<- "Overall X label"
  theYLabel<- "Overall Y label"
# these are the axes ranges -- set these for
# each application 
# Here the ranges are just repeated
# change these if a the ranges for the rows or columns varies
#
  yr<- matrix(c(0,11)+100, n, 2, byrow=TRUE)
  xr<- matrix(c( 0,20 ), m, 2, byrow=TRUE)
  
# You will have to tinker with the legend spacing in most cases 
# we just do this by trail and error!
# set LEGENDSPACE  to zero if you are not including a legend outside of the 
# panel 
# LEGENDSPACE<- 0
  
  LEGENDSPACE<- 8
  

# Each expression are the graphics code for one of the plots in the panel. 
# if you plotting data objects either pass them to this function or expcet them to be
# in your workspace.
# The code for each plot can be any graphics function or R code. 
# Below are just some simple examples 
# you can also hard code  plot function  into the loop below but this 
# makes the coding more complex.
# Note that  some of these are example plots use  multiple plot commands. 
 
  
plotCode<- list(
  expression(points( (1:10)*2, 101:110, pch=4), 
             points( (1:5), 101:105, pch=19, col="blue")
             ),
  
  expression( points( (1:5), 101:105, pch=19, col="blue"),
              # cheat and put title within the plot
              title("title for 2nd plot", line=-1) 
             ),
  
  expression(image( 1:5*2, 104:108, 
                    outer(1:5, 106:110,"+"),
                    add=TRUE, 
                    col=tim.colors(), 
                    zlim=c( 105, 120))
             ,
             points(5, 107,pch=16, col="grey20"),
             title("multiple plot commands 
                    in this plot"
                   , line=-2)
             ),
  
  expression(points( (1:5)*2, 101:105, pch=14, cex=2, col="orange2"),
             title("(alpha = .5)", adj=.05, line=-1)
             ),
  NULL, # this means skip the 5th plot ...
  expression(lines( (3:10)*2, 103:110, lwd=4,col="green4") 
  )
)
 

########################################################
###  The following  stuff may not need to be changed 
#########################################################

# suggested sizes for readable figures
par(cex.axis = 1.5, cex.lab = 1.5, cex = 1.5, 
     lwd = 1.5, bg = "transparent", pch = 16)
# oranges2+ grey background = Lufthansa
palette(c("orange2", "green2", "blue2", "red1"))


# nonspecific settings
#laying out canvas, panel, and each plot area (pane)
# units are in chararcter size
PANESPACE<- .5 # set larger to make room for plot titles.
par( mfrow=c(m,n),
mar=c(PANESPACE, PANESPACE,0,0),
oma=c(5,5,3,3+LEGENDSPACE))
colBg<-"grey85"

verbose<- FALSE # handy switch for printing debug info

########################################################
###  loop over plots this is done row by row
#########################################################

counter<- 0
for( I in 1:m ){
  for( J in 1:n){
    counter<- counter+1
    
    #blank canvas for plot I,J
    plot( xr[I,],yr[J,], xaxt = "n", yaxt = "n",type="n", bty="n")
    thePlot<- plotCode[[counter]]
    
    if(verbose){
      cat(" ",fill=TRUE )
      cat("counter: ",counter,  fill=TRUE )
      cat("I,J: ", I,J,      fill=TRUE )
      cat("plot: ", as.character(thePlot), fill=TRUE)
    }

# skip this pane if NULL
    if(!is.null(thePlot)){
    # some default background setup for the I,J plot ...
    # draw background first
    xy<- par()$usr
    rect(xy[1],xy[3],xy[2],xy[4],col=colBg)
    #
    # here is a handy function for the grid we missed this function for about 25 years!
    grid(col="grey98", lty=1, lwd=.5)
    box(lwd=.5)
    
    # here the plot code for plot I,J is run
    eval(thePlot)
    #
    # add more plot stuff here if needed 
    # use the counter or I,J for special 
    # things for particular panes 
    # e.g.  
    # if( counter == 4){ text(5,10,"a funky label")}
    #
    # if( I=1 & J=2){
    #   title(" a special case title for just plot (1,2)")
    # }
    #
    }
    
# add axes but only to  the plots in the  outer edges  
# add more options to these axis functions  as needed
# e.g. cex to increase text size or label info 
    if( I==1){   axis( side=3)}
    if( I==m){    axis(side=1)}
    if( J ==1 ){  axis(side=2)}
    if( J == n ){ axis(side=4)}
#   
  }
}

# done looping through the plots now add stuff outside the panel in the 
# margins. 

########################################################
### Add overall axis titles
########################################################
# margin labels added on bottom and left edges  
# modify these if you don't like the label placement/style
 mtext(outer=TRUE, side=1, text = theXLabel, line=3)
 mtext(outer=TRUE, side=2, text = theYLabel, line=3)
 
#### legend placement
# We don't like this code -- but it works ...
# Creates a new canvas with plot coords 0,1 w/o erasing
# your beautiful graphic above 
# this means youcan put the legend and I guess
# any other annotation 
# anywhere -- keeping in mind the new plotting coordinates are 0,1
# so pretty easy to guess where things will go 
 
 par( mar = c(0, 0, 0, 0),
      oma=c(0,0,0,0),
      fig=c(0,1,0,1), 
      new = TRUE, # not sure what this does
      usr=c(0,1,0,1))
########################################################
### Begin figure specific placement of legend/color bar
########################################################
#
# usually you will have to adjust exactly where the legend/colorbar info goes
# easiest to just try different locations
# 
# first the more basic way 
# legend(
#        x=.85,y=.85, 
#        xpd=TRUE,
#        cex=2,
#        legend= c("IM", "IBD", "1R", "2R"),
#        pch = c(4, 2, 15, 19), col = 1:4
# )
#
#  more readable code for a legend  using inset argument
 legend("topright",
        horiz = FALSE,
        inset = c(.02,.05),
        cex=2,
        legend= c("IM", "IBD", "1R", "2R"),
        pch = c(4, 2, 15, 19), col = 1:4
        )

 # add a color bar in margin too
 # this needs the fields library
 # the small plot coordinates  determine position and size
 # relative to the full size of the figure region. 
 # see the imagePlot help file for more details
 # 
 
 imagePlot( legend.only=TRUE, col=tim.colors(),
            zlim=c( 105, 120),
            smallplot= c( .81,.81+.03,.2,.5),
            legend.lab="ColorBar Label",
            legend.line=4)
 
} 
 

 