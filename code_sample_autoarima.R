# =========================================================================== #
# Author:  mgrossi
# Created: 13 Apr 2020
# Updated: 06 May 2023 - Fixed typos in comments
#
# The functions below are designed to automatically and efficiently fit ARIMA
# (AutoRegressive Integrated Moving Average) models to many (O(1000)) drifter
# time series. Forecasts from these traditional statistical regression methods 
# serve as performance benchmarks for any new machine learning based model we
# develop.
#
# These routines are called and used by a larger and more complicated script 
# that runs a simulation of an operational drifter deployment.
#
# auto.arima eliminates the need to conduct a manual parameter search for all
# time series (once the limits have been determined - see comments below) and
# implementing it with R's "apply" family of functions allows this to be done
# in a vectorized fashion. The original parent script compiles many such
# forecasts from run.autoarima (over the course of a simulated deployment) and
# writes them out to a netCDF file.
#
# =========================================================================== #

normalize <- function(x){
  # Return normalized dataset
  x_norm <- (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
  return(x_norm)
}

unnormalize <- function(x, mins, maxs){
  # Returns unnormalized dataset based on mins and maxs
  x_unnorm <- (x * (maxs - mins)) + mins
  return(x_unnorm)
}

get_var <- function(data, cols){
  # Get ARIMA output variable
  data.frame(matrix(unlist(data), ncol=length(cols), byrow=FALSE,
             dimnames=list(NULL, cols)))
}

plot_forecast <- function(predicts, targets, xmin=NA, xmax=NA, fname=NA){
  # Plot an ARIMA prediction object with target values. Writes to file if
  # 'fname' is provided.
  
  # Check that the class of the object passed to 'predicts' is 'forecast' and
  # convert if necessary.
  if(class(predicts) != 'forecast'){
    predicts <- structure(predicts, class='forecast')
  }
  
  # If fname is passed, create a file to be written
  if(!is.na(fname)){
    bitmap(fname, width=4, height=2, res=800, units='in', type='jpegcmyk')
  }
  
  # Get min and max from 'predicts' if not passed
  if(is.na(xmin)){
    xmin <- start(predicts$fitted)[1]
  }
  if(is.na(xmax)){
    xmax <- end(predicts$mean)[1]
  }
  
  # Plot the time series
  plot(predicts, main=predicts$method, xlim=c(xmin, xmax),
       ylim=range(c(predicts$upper, predicts$lower, predicts$fitted, targets)))
  points(ts(targets, start=end(predicts$fitted)[1]), type='l', lty='dashed')
  
  # Write file (if applicable)
  if(!is.na(fname)){
    dev.off()
  }
}

run.autoarima <- function(u, v, forecastDays=1){
  # =========================================================================== #
  # Fit ARIMA models to a collection of drifter velocity time series. This
  # function operates on the full data sets "u" and "v" (utilizing the "apply" 
  # function family for vectorization) to create unique ARIMA models for each 
  # time series (columns in "u" and "v"). These models are then run to create
  # forecasts of length "forecastDays".
  #
  # Inputs:
  #   u : data.frame containing zonal drifter velocities arranged with time
  #		    in rows and drifters in columns, with drifter ID as column names
  #   v : data.frame containing meridional drifter velocities arranged with
  #       time in rows and drifters in columns, with drifter ID as column
  #       names
  #   forecastDays : int, number of days over which to forecast drifters
  #
  # Output:
  #   Returns a list containing forecasts, model parameters, and other output
  #   information that may be of interest later on
  
  # Required packages
  require(forecast)
  require(lubridate)
  require(ncdf4)
  require(abind)
  require(plyr)
  
  # ======================================================================= #
  # Normalize data
  uNorm <- normalize(u)
  vNorm <- normalize(v)
  
  # Extract timestamps
  datetime <- ymd_hms(rownames(u), tz='GMT')
  
  # Number of observations per day contained in the data set
  obsPerDay <- max(table(as.Date(datetime)))
  
  # Prediction dimension: length of each forecast
  ptimeDim <- forecastDays * obsPerDay
  
  # Observation time step in hours
  timeInt <- as.numeric(difftime(datetime[2], datetime[1], units='hours'))

  # Exclude drifters containing less than 24 hrs worth of observations
  ind <- (colSums(!is.na(uNorm)) >= (24/timeInt)) & !is.na(uNorm[nrow(uNorm),])
  u_ss <- uNorm[,ind]
  v_ss <- vNorm[,ind]
  
  # Estimate an ARIMA model for each drifter
  # Note: These max.p and max.q values were determined by a full run with these
  # set to 96 and 48, respectively.
  # This is the most computationally expensive step, as it conducts parameter
  # searches for all available drifters
  uAR <- lapply(u_ss, function(x) auto.arima(x, max.p=16, max.q=12))
  vAR <- lapply(v_ss, function(x) auto.arima(x, max.p=16, max.q=12))

  # Apply respective models to the data
  uMapped <- mapply(function(uAR, u_ss) forecast(uAR, h=ptimeDim), uAR, u_ss)              
  vMapped <- mapply(function(vAR, v_ss) forecast(vAR, h=ptimeDim), vAR, v_ss)              

  # Issue forecasts for the drifters
  forecasts_u <- lapply(apply(uMapped, 2, list), function(x) x[[1]])
  forecasts_v <- lapply(apply(vMapped, 2, list), function(x) x[[1]])

  # Extract forecasts and standard errors
  # These are reformatted for ease of use later on.
  uPred <- get_var(lapply(forecasts_u, function(x) as.vector(x$mean)), 
                   cols=colnames(u_ss))
  vPred <- get_var(lapply(forecasts_v, function(x) as.vector(x$mean)), 
                   cols=colnames(v_ss))
  ulo80 <- get_var(lapply(forecasts_u, function(x)
                   as.vector(x$lower[,'80%'])), cols=colnames(u_ss))
  uhi80 <- get_var(lapply(forecasts_u, function(x)
                   as.vector(x$upper[,'80%'])), cols=colnames(u_ss))
  ulo95 <- get_var(lapply(forecasts_u, function(x)
                   as.vector(x$lower[,'95%'])), cols=colnames(u_ss))
  uhi95 <- get_var(lapply(forecasts_u, function(x)
                   as.vector(x$upper[,'95%'])), cols=colnames(u_ss))
  vlo80 <- get_var(lapply(forecasts_v, function(x)
                   as.vector(x$lower[,'80%'])), cols=colnames(v_ss))
  vhi80 <- get_var(lapply(forecasts_v, function(x)
                   as.vector(x$upper[,'80%'])), cols=colnames(v_ss))
  vlo95 <- get_var(lapply(forecasts_v, function(x)
                   as.vector(x$lower[,'95%'])), cols=colnames(v_ss))
  vhi95 <- get_var(lapply(forecasts_v, function(x)
                   as.vector(x$upper[,'95%'])), cols=colnames(v_ss))

  # Unnormalize output (mins and maxs are taken from original data)
  uPred <- unnormalize(uPred, mins=min(u[,ind], na.rm=TRUE),
                       maxs=max(u[,ind], na.rm=TRUE))
  vPred <- unnormalize(vPred, mins=min(v[,ind], na.rm=TRUE),
                       maxs=max(v[,ind], na.rm=TRUE))
  ulo80 <- unnormalize(ulo80, mins=min(u[,ind], na.rm=TRUE),
                       maxs=max(u[,ind], na.rm=TRUE))
  uhi80 <- unnormalize(uhi80, mins=min(u[,ind], na.rm=TRUE),
                       maxs=max(u[,ind], na.rm=TRUE))
  ulo95 <- unnormalize(ulo95, mins=min(u[,ind], na.rm=TRUE),
                       maxs=max(u[,ind], na.rm=TRUE))
  uhi95 <- unnormalize(uhi95, mins=min(u[,ind], na.rm=TRUE),
                       maxs=max(u[,ind], na.rm=TRUE))
  vlo80 <- unnormalize(vlo80, mins=min(v[,ind], na.rm=TRUE),
                       maxs=max(v[,ind], na.rm=TRUE))
  vhi80 <- unnormalize(vhi80, mins=min(v[,ind], na.rm=TRUE),
                       maxs=max(v[,ind], na.rm=TRUE))
  vlo95 <- unnormalize(vlo95, mins=min(v[,ind], na.rm=TRUE),
                       maxs=max(v[,ind], na.rm=TRUE))
  vhi95 <- unnormalize(vhi95, mins=min(v[,ind], na.rm=TRUE),
                       maxs=max(v[,ind], na.rm=TRUE))

  # Extract standard ARIMA model parameters
  upList <- lapply(uAR, function(x) as.numeric(arimaorder(x)['p']))
  udList <- lapply(uAR, function(x) as.numeric(arimaorder(x)['d']))
  uqList <- lapply(uAR, function(x) as.numeric(arimaorder(x)['q']))
  vpList <- lapply(vAR, function(x) as.numeric(arimaorder(x)['p']))
  vdList <- lapply(vAR, function(x) as.numeric(arimaorder(x)['d']))
  vqList <- lapply(vAR, function(x) as.numeric(arimaorder(x)['q']))
    
  # Extract seasonal ARIMA model parameters (fillVal=NA if not applicable)
  uPList <- lapply(uAR, function(x) as.numeric(arimaorder(x)['P']))
  uDList <- lapply(uAR, function(x) as.numeric(arimaorder(x)['D']))
  uQList <- lapply(uAR, function(x) as.numeric(arimaorder(x)['Q']))
  vPList <- lapply(vAR, function(x) as.numeric(arimaorder(x)['P']))
  vDList <- lapply(vAR, function(x) as.numeric(arimaorder(x)['D']))
  vQList <- lapply(vAR, function(x) as.numeric(arimaorder(x)['Q']))   
  uFreqList <- lapply(uAR, function(x) as.numeric(arimaorder(x)['freq']))
  vFreqList <- lapply(vAR, function(x) as.numeric(arimaorder(x)['freq']))

  # Compile results
  # (Note to self: consider replacing this with a hash (R dictionary) for
  # efficiency, if needed)
  output <- list(uPred, vPred, ulo80, uhi80, ulo95, uhi95, ulo80, vhi80,
                 vlo95, vhi95, upList, udList, uqList, vpList, vdList,
                 vqList, uPList, uDList, uQList, vPList, vDList, vQList,
                 uFreqList, vFreqList)
  names(output) <- c('uPred', 'vPred', 'ulo80', 'uhi80', 'ulo95', 'uhi95',
                     'ulo80', 'vhi80', 'vlo95', 'vhi95', 'upList', 'udList',
                     'uqList', 'vpList', 'vdList', 'vqList', 'uPList',
                     'uDList', 'uQList', 'vPList', 'vDList', 'vQList',
                     'uFreqList', 'vFreqList')
  return(output)
  }

# =========================================================================== #
