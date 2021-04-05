###############################################################################
#
# MSC 301: Introduction to Physical Oceanography
# Homework 05 Sample Code
#
# Author: mgrossi
# Date: 01 Nov 2018
#
###############################################################################

require(R.matlab)
require(fields)

setwd("~/Documents/TA/IntroPhysOce/HW03/")
dat <- readMat("hw_3_data.mat")

# =============================================================================
# QUESTION 2a: In homework 3, we calculated geostrophic velocities across a
# loop current eddy. Estimate the planetary vorticity at the center of this
# eddy.

# Save each variable in the data object 'dat' as unique R variables
# (note that this is neither necessary nor recommended practice for R code,
# but we'll do it here for simplicity and readability):

lat <- dat$lat
lon <- dat$lon
tm <- dat$t
ssh <- dat$ssh

# Visualize the data to locate the center of the eddy.

image.plot(t(ssh[,,1]), axes=FALSE,
           main="Sea Surface Height Anomalies",
           xlab="Longitude (deg)", ylab="Latitude (deg)",
           legend.lab="SSH anomalies (m)")
axis(1, at=seq(0,1,length=21), labels=seq(min(lon), max(lon), length=21))
axis(2, at=seq(0,1,length=14), labels=seq(min(lat), max(lat), length=14))
box()

# Planetary vorticity is given by:  f = 2*omega*sin(latitude)
# The center of the eddy is at roughly 25.5ºN.

omega <- 7.2921e-5 # rad/s
(f <- 2 * omega * sin((25.5*pi)/180)) # Converting to radians

# =============================================================================
# QUESTION 2a: Use the geostrophic velocities to estimate the relative
# vorticity of the loop current eddy (Hint: You will need two perpendicular
# transects across the entire eddy.)

# To calculate relative vorticity, we need to know dv/dx and du/dy. Note that
# the LCE is centered at roughly 25.5ºN, 87.5ºW, characterized by larger SSH
# values. Our transects should go through the center of the eddy.

# Calculate u and v in the same way as Homework 3. First, u across a
# north-south transect:

NStrans_lon_ind <- which(lon == (-87.5))
NStrans_lat_ind <- which(lat > 24 & lat < 27)

NStrans_lon <- lon[NStrans_lon_ind]
NStrans_lat <- lat[NStrans_lat_ind]

eta <- ssh[NStrans_lat_ind, NStrans_lon_ind, 1]

# Geostrophic balance in the y direction:  -g * d_eta/d_y = fu
g <- 9.81 # m/s
deta <- diff(eta) # m
dy <- diff(NStrans_lat) * 111e3 # m
u <- - (g/f) * (deta/dy)

# Now u across an east-west transect:

EWtrans_lat_ind <- which(lat == 25.5)
EWtrans_lon_ind <- which(lon > (-89) & lon < (-85))

EWtrans_lat <- lat[EWtrans_lat_ind]
EWtrans_lon <- lon[EWtrans_lon_ind]

eta <- ssh[EWtrans_lat_ind, EWtrans_lon_ind, 1]

# Geostrophic balance in the x direction is:  -g * d_eta/d_x = -fv
g <- 9.81 # m/s
deta <- diff(eta) # m
dx <- diff(EWtrans_lon) * 111e3 * cos((EWtrans_lat*pi)/180) # m
v <- (g/f) * (deta/dx)

# Relative vorticity is given by:  zeta = dv/dx - du/dy.

du <- diff(u)
dv <- diff(v)

# Since we have a uniform spatial grid, every dx is the same distance, so we
# can take the first element of dx and divide it into every dv. Same with du
# and dy.

dvdx <- dv / dx[1]
dudy <- du / dy[1]

# Remember that dvdx and dudy are still vectors. Let's use the mean of these
# to calculate zeta:

(zeta = mean(dvdx) - mean(dudy)) # 1/s

# =============================================================================
# QUESTION 2c: Assume that the eddy can be represented by a cylinder that
# extends to a depth of 1000 m. Calculate the potential vorticity (PV).

# Potential vorticity is defined as: (zeta + f)/h

h <- 1000 # m

(PV <- (zeta + f) / h)
