###########################################################
#
# MSC 301: Introduction to Physical Oceanography
# Homework 01 Sample Code
#
# Author: mgrossi
# Date: 01 Sep 2018
#
###########################################################

# =========================================================
# QUESTION 2d:
# Carefully plot the temperature and salinity profiles from
# the accompanying data file, being sure to label all axes.

# First move into the directory (folder) in which the files
# are saved and then read in the csv file.

setwd('/path/to/wherever/you/saved/the/homework/files/')            # This is the FOLDER!
dat = read.csv(file='HW01_GoM_data.csv', header=TRUE, row.names=1)  # This is the FILE!

# Remember that by default read.csv() assumes that the first
# row in the file contains data. In our case, the first
# row contains column names, so we specify this by passing
# header=TRUE. Similarly, row.names=1 tells R that the
# first column are actually row names, not data.
#
# Now to make the plots!
#
# We need to tell R what to plot. Just like in Excel, where
# we would specify which column contains the x values and
# which column contains the y values, we need to do the
# same thing for the R function plot().
#
# We have saved our data to a variable we called 'dat'.
# Recall that there are two ways to index specific columns
# in R: we can either specify row and column using
# dat[row#,column#], or we can use a column name, such as
# dat$Depth_m. (Hint: Use colnames(dat) to see what the
# column names are, if you are unsure.)
#
# To make a temperature profile, we want temperature as a
# function of depth. For visualization purposes, we
# typically plot depth on the y axis with 0 at the top.
#
# In order to compare the profiles, we need the x axis
# ranges to be the same on both temperature plots. We can
# do this either by combining the two temperature profiles
# together and then using the range() function to extract
# the joint min and max, or by indexing the two
# temperature columns and using range() to find the max
# and min. See below to see what I mean.
#
# Tip: The function png() is used to save a plot as a png
# file. I'll demonstrate it here, along with some other
# handy options, for your future use.

# Temperature vs. Depth for profile 1

png(filename='P1_TvsD.png')                   # save plot to file (optional)
plot(x=dat$Prof1_Temp_degC, y=dat$Depth_m,    # x, y values
     xlim=range(dat[,c(2,4)]),                # set x axis range (see comments above)
     ylim=rev(range(dat$Depth_m)),            # set y axis range to have 0 at the top
     type='l', col='blue', lwd=3,             # plot as a blue line of width 3
     main='Profile 1: Temperature vs. Depth', # add a title
     xlab='Temp (deg C)', ylab='Depth (m)')   # set x, y axis labels
dev.off()                                     # write plot to file (used with png())

# png() opens a file and dev.off() closes the file
# (literally, 'turn the plotting device off.')
# Everything in between actually writes the file.
# Note that the file will be saved in your current
# working directory. If you're not sure what this
# is, use getwd().

# Now for profile 2:

png(filename='P2_TvsD.png')
plot(x=dat$Prof2_Temp_degC, y=dat$Depth_m,
     xlim=range(dat[,c(2,4)]),
     ylim=rev(range(dat$Depth_m)),
     type='l', col='blue', lwd=3,
     main='Profile 2: Temperature vs. Depth',
     xlab='Temp (deg C)', ylab='Depth (m)')
dev.off()

# Salinity vs. Depth

png(filename='P1_SvsD.png')
plot(x=dat$Prof1_Sal_psu, y=dat$Depth_m,
     xlim=range(dat[,c(3,5)]),
     ylim=rev(range(dat$Depth_m)),
     type='l', lwd=3, col='blue',
     main='Profile 1: Salinity vs. Depth',
     xlab='Salinity (psu)', ylab='Depth (m)')
dev.off()

png(filename='P2_SvsD.png')
plot(x=dat$Prof2_Sal_psu, y=dat$Depth_m,
     xlim=range(dat[,c(3,5)]),
     ylim=rev(range(dat$Depth_m)),
     type='l', lwd=3, col='blue',
     main='Profile 2: Salinity vs. Depth',
     xlab='Salinity (psu)', ylab='Depth (m)')

# There are always many ways to do the same thing! This is
# just one way. You may discover alternative options.
# You might, for example, make the depths negative instead
# of setting the ylim to get 0m at the top of the graph.
#
# Since we're comparing two profiles, it may be helpful to
# show both temperature profiles on one plot, and
# similarly, both salinity profiles on one plot. How might
# we do that in R?

# Temperature vs. Depth

plot( x=dat$Prof1_Temp_degC, y=dat$Depth_m,   # Plot one profile the same we
     xlim=range(dat[,c(2,4)]),                # did above. It doesn't matter
     ylim=rev(range(dat$Depth_m)),            # which one (profile 1 or 2)
     type='l', lwd=3, col='blue',             # you do first.
     main='Temperature vs. Depth',
     xlab='Temp (deg C)', ylab='Depth (m)')
points(y=dat$Depth_m, x=dat$Prof2_Temp_degC,  # Add the second profile using points()
       type='l', lwd=3, col='red')            # function with the same settings.
legend('topleft',                             # Add a legend. Notice we used different colors
       legend=c('Profile 1', 'Profile 2'),    # series labels
       col=c('blue', 'red'), lty=1, lwd=3)    # color and line type

# Salinity vs. Depth

plot(x=dat$Prof1_Sal_psu, y=dat$Depth_m,
     xlim=range(dat[,c(3,5)]),
     ylim=rev(range(dat$Depth_m)),
     type='l', lwd=3, col='blue',
     main='Salinity vs. Depth',
     xlab='Salinity (psu)', ylab='Depth (m)')
points(y=dat$Depth_m, x=dat$Prof2_Sal_psu,
       type='l', lwd=3, col='red')
legend('bottomleft', legend=c('Profile 1', 'Profile 2'),
       col=c('blue', 'red'), lty=1, lwd=3)

# =========================================================
# QUESTION 2e:
# Calculate the density profile at each station. Plot your
# results, being careful to label all axes.

# First load the equation of state function, and then use
# it to calculate the two density profiles.

source('EOS80.R')

# Some of you discovered another way of loading the function:
# Double-clicking the file opens it in RStudio as an R
# script, which you can then run. Both methods work. One way
# to be sure the function loaded is to look for water.density
# under "Functions" in your Global Environment in the upper
# right panel of RStudio.
#
# Remember that water.density() is a function, like plot().
# Just like a math function, for which one needs to provide
# numerical values for each variable in order to obtain the
# value of the function, we need to pass values to an R
# function in order to get anything out. water.density() has
# three variables: temp, sal, and pres. We need to supply a
# value (or a vector of values) for each of these variables
# in order to get an output.
#
# For convienence, let's define two new variables, one for
# each profile. What does the function water.density()
# produce?

dens1 = water.density(temp = dat$Prof1_Temp_degC, 
                      sal = dat$Prof1_Sal_psu, 
                      pres = dat$Depth_m)

dens2 = water.density(temp = dat$Prof2_Temp_degC,
                      sal = dat$Prof2_Sal_psu,
                      pres = dat$Depth_m)

# Make the plots, as above.

# Option 1: Individual plots

# Density vs. Depth
plot(x=dens1, y=dat$Depth_m,
     xlim=range(c(dens1, dens2)),
     ylim=rev(range(dat$Depth_m)),
     type='l', lwd=3, col='blue',
     main='Profile 1: Density vs. Depth',
     xlab='Density (kg/m^3)', ylab='Depth (m)')

plot(x=dens2, y=dat$Depth_m,
     xlim=range(c(dens1, dens2)),
     ylim=rev(range(dat$Depth_m)),
     type='l', lwd=3, col='blue',
     main='Profile 2: Density vs. Depth',
     xlab='Density (kg/m^3)', ylab='Depth (m)')

# Option 2: Both profiles on one plot

# Density vs. Depth
plot(x=dens1, y=dat$Depth_m,
     xlim=range(c(dens1, dens2)),
     ylim=rev(range(dat$Depth_m)),
     type='l', lwd=3, col='blue',
     main='Density vs. Depth',
     xlab='Density (kg/m^3)', ylab='Depth (m)')
points(x=dens2, y=dat$Depth_m,
       type='l', lwd=3, col='red')
legend('bottomleft', legend=c('Profile 1', 'Profile 2'),
       col=c('blue', 'red'), lty=1, lwd=3)

# That's all there is to it!
#
# Take some time to play around with these samples and
# become comfortable with setting plot parameters, indexing
# data, and extracting information from the data. Remember,
# practice makes perfect!
# =========================================================

