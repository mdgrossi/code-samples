
# R<~/round_3/scripts/qc_summary.R --vanilla --slave >> ~/round_3/out/qc_summary.out 2>&1 &

# =============================================================================
# Author:  Matt Grossi
# Created: 4 Jan 2015
# Updated: 26 Mar 2019 - uses ncdf4, since ncdf has been deprecated; added
#                        pressure inversion test and >=20 data points test
#          03 Apr 2019 - added "profile >10 m" QC test
#          09 Apr 2019 - revised to operate on netCDF files in their original
#                        directories, IDed using the directory name instead of
#                        "A", "I", and "P"
#          11 Apr 2019 - added max depth and number of bins columns
#          12 Apr 2019 - moved TooShallow test to argo_filtering.R
#          15 Apr 2019 - added simulated density profile QC checks, profile bin
#                        counts, and profile max/min depths
#          06 May 2023 - fixed spacing in comments, removed home directory
#                        references
#
# OVERVIEW
# The global Argo data set from the program’s conception to the present (2019)
# contains around 2 million CTD profiles. The goal is to use as many of these
# as possible, but some must be rejected based on a variety of data quality
# control (QC) tests. Instead of compiling all 2 million profiles and then
# removing those that fail the tests, one can first loop through all the files,
# extract the necessary metadata and relevant QC flags, and later use this
# information to filter the metadata and compile only the profiles that pass
# all QC tests.
#
# Although this requires looping through all files twice (once to extract the
# metadata, a second time to compile the profiles), this method is less
# memory-intensive than the alternative approach of extracting and compliling
# metadata, QC information, and data all at once before filtering.
#
# During this first loop through the files, the total number of profiles in
# each file is also added up. This sum represents the starting number of
# profiles before filtering out bad/rejected data.
# 
# =============================================================================
# Print header for ".out" file
message("DATA QC SUMMARY FOR ALL DOWNLOADED PROFILES CONTAINING REQ'D DATA")
message("Script: qc_summary.R")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
system("date")

# =============================================================================
# WORKFLOW & RATIONALE
#
# Packages
require(ncdf4)

# Functions
# Deal with entire profiles of NA (otherwise -Inf is returned)
na.min <- function(x) ifelse(!all(is.na(x)), min(x, na.rm=T), NA)
na.max <- function(x) ifelse(!all(is.na(x)), max(x, na.rm=T), NA)

# Calculate profile QC flags
flagCalc <- function(x){
	# In accordance with Argo QC convention, "good" profiles will be those
	# having QC flag of either 1, 2, 5, or 8 (see QC manual). Each profile
	# is graded based on the precentage of good observations.
    n = apply(x, 2, function(x) length(x[x %in% c(1, 2, 5, 8)]) /
                  length(na.omit(x[!x %in% 9]))) * 100
    
    flag = ifelse(n == 100, "A",
           ifelse(75 <= n & n <100, "B",
           ifelse(50 <= n & n < 75, "C",
           ifelse(25 <= n & n < 50, "D",
           ifelse(0 < n & n < 25, "E",
           ifelse(n == 0, "F", " "))))))
    
    return(flag)
}

# File list to loop through
setwd("~/argo/prof/all/")
basin <- c("atlantic_ocean", "indian_ocean", "pacific_ocean")
files <- NULL
for(b in basin){
    files <- c(files, paste(b, dir(paste(getwd(), b, sep="/"), pattern="*.nc"),
    							sep="/"))
    }
message(paste0(length(files), " files to process."))

# Object to store compiled metadata from all profiles (see explanations below)
qc <- NULL

# "Total number of profiles downloaded" counter: keep track of the (cumulative)
# sum of the number of profiles within each data file
numberOfProfiles <- 0

# Loop through downloaded files to extract information
for(fname in files){
	
# Full status update (for testing purposes)
    message(paste0("Starting file ", f, " of ", length(files)))
	nc <- nc_open(fname)

# Number of profiles contained in file
	profilesThisFile <- length(ncvar_get(nc, "JULD"))
	numberOfProfiles <- numberOfProfiles + profilesThisFile

# =============================================================================
# CHECK for required metadata, variables, and QC flags in netCDF file:
#		JULD					= Julian day of profile
#		JULD_QC					= Julian day QC flag
#		LATITUDE				= Latitude of profile
#		LONGITUDE				= Longitude of profile
#		POSITION_QC				= Position QC flag
#		DATA_MODE				= Real-time / Delayed-mode / Adjusted
#		DATA_STATE_INDICATOR	= Data processing stage
#		PROFILE_PRES_QC			= QC flag for entire pressure profile
#		PROFILE_TEMP_QC			= QC flag for entire temperature profile
#		PROFILE_PSAL_QC			= QC flag for entire salinity profile
#		PRES					= Pressure data
#		TEMP					= Temperature data
#		PSAL					= Salinity data
#		PRES_QC					= Pressure data QC flags
#		TEMP_QC					= Temperature data QC flags
#		PSAL_QC					= Salinity data QC flags
#
# Each file must contain all of the above variables and flags. If any are
# missing from the file, the entire file is rejected because all profiles
# within the file would be missing information. For example, early Argo floats
# did not measure salinity, so PSAL, PSAL_QC, and PROFILE_PSAL_QC are missing.
# Density cannot be calculated without salinity, so such profiles are useless
# here.
#
# NOTE: The number of files rejected can be determined afterwards by comparing
# the number of files in this list with the number of files downloaded.
	varNames <- c("JULD", "JULD_QC", "LATITUDE", "LONGITUDE", "POSITION_QC",
	             "DATA_MODE", "DATA_STATE_INDICATOR",
	             "PROFILE_PRES_QC", "PROFILE_TEMP_QC", "PROFILE_PSAL_QC",
	             "PRES", "TEMP", "PSAL", "PRES_QC", "TEMP_QC", "PSAL_QC")

# If the file contains all required information, continue with meatadata
# extraction.
if(sum(!varNames %in% names(nc$var))==0){

# =============================================================================
# EXTRACT data and metadata
    
# Data state and data mode
    dataState <- trimws(ncvar_get(nc, "DATA_STATE_INDICATOR"))
    dataMode <- as.character(unlist(strsplit(ncvar_get(nc,"DATA_MODE"),
                                             split="")))

# Data (real-time and adjusted) and replace fill value 99999 with NA
	pres <- as.matrix(ncvar_get(nc, "PRES"))
	pres <- replace(pres, pres==99999, NA)
	temp <- as.matrix(ncvar_get(nc, "TEMP"))
	temp <- replace(temp, temp==99999, NA)
	psal <- as.matrix(ncvar_get(nc, "PSAL"))
	psal <- replace(psal, psal==99999, NA)

	presAdj <- as.matrix(ncvar_get(nc, "PRES_ADJUSTED"))
	presAdj <- replace(presAdj, presAdj==99999, NA)
	tempAdj <- as.matrix(ncvar_get(nc, "TEMP_ADJUSTED"))
	tempAdj <- replace(tempAdj, tempAdj==99999, NA)
	psalAdj <- as.matrix(ncvar_get(nc, "PSAL_ADJUSTED"))
	psalAdj <- replace(psalAdj, psalAdj==99999, NA)
    
# =============================================================================
# EXTRACT QC flags for each data point (profile depth z)
# Convert from character string to numerical matrix
	presQC <- ncvar_get(nc, "PRES_QC")
	presQC <- matrix(as.numeric(unlist(strsplit(presQC, split=""))),
	                 ncol=length(presQC))
	tempQC <- ncvar_get(nc, "TEMP_QC")
	tempQC <- matrix(as.numeric(unlist(strsplit(tempQC, split=""))),
	                 ncol=length(tempQC))
	psalQC <- ncvar_get(nc, "PSAL_QC")
	psalQC <- matrix(as.numeric(unlist(strsplit(psalQC, split=""))),
	                 ncol=length(psalQC))
	
	presAdjQC <- ncvar_get(nc, "PRES_ADJUSTED_QC")
	presAdjQC <- matrix(as.numeric(unlist(strsplit(presAdjQC, split=""))),
	                    ncol=length(presAdjQC))
	tempAdjQC <- ncvar_get(nc, "TEMP_ADJUSTED_QC")
	tempAdjQC <- matrix(as.numeric(unlist(strsplit(tempAdjQC, split=""))),
	                    ncol=length(tempAdjQC))
	psalAdjQC <- ncvar_get(nc, "PSAL_ADJUSTED_QC")
	psalAdjQC <- matrix(as.numeric(unlist(strsplit(psalAdjQC, split=""))),
	                    ncol=length(psalAdjQC))
	
# Define parameter bounds and change QC flag to 4 ("bad data") for any
# measurements that fall outside these physical ranges:
#
#     Pressure:  >0 dbar*
#	  Temp:	     -2.5 to 40.0 ºC (from global climatology and Argo user manual)
#     Salinity:  20 to 41 psu (from global climatology - Argo uses [2,41])
#
# *Due to the introduction of deep profiling floats capable of diving to 6000 m,
# pressures greater than 2100 m can no longer be set as bad data if these
# deep profiles are to be retained. Instead, check for minimum depth only, and
# later mark as NA depths below 2100 m (data below this depth won't be used
# anyway.)
#
# While this should be done automatically by Argo teams before uploading data,
# a number of cases have been discovered for which this check had not been
# conducted. In these cases, the observations are outside of the physical
# limits yet the data are still flagged as good.
#
# Note: formerly "data out-of-bounds" QC filter
	presQC <- replace(presQC, pres < 0, 4)
	presAdjQC <- replace(presAdjQC, presAdj < 0, 4)
	tempQC <- replace(tempQC, temp < (-2.5) | temp > 40, 4)
	tempAdjQC <- replace(tempAdjQC, tempAdj < (-2.5) | tempAdj > 40, 4)
	psalQC <- replace(psalQC, psal < 20 | psal > 41, 4)
	psalAdjQC <- replace(psalAdjQC, psalAdj < 20 | psalAdj > 41, 4)

# Check for pressure inversion (i.e., pressure decreasing with depth) in both
# the real-time and adjusted pressure variables. If an inversion occurs between
# the first and second levels, make the first pressure reading (i.e., the
# surface reading, which is most prone to measurement error or noise due to
# exposure to the atmosphere) NA. Then check for monotonically increasing
# pressure.
	pres[1, which(apply(pres, 2, function(x) (x == cummin(x))[2]))] <- NA
	presQC[1, which(apply(pres, 2, function(x) (x == cummin(x))[2]))] <- 3
	presInv <- !apply(pres, 2, function(x) all(na.omit(x) == cummax(na.omit(x))))
	
	presAdj[1, which(apply(presAdj, 2, function(x) (x == cummin(x))[2]))] <- NA
	presAdjQC[1, which(apply(presAdj, 2, function(x) (x == cummin(x))[2]))] <- 3
	presAdjInv <- !apply(presAdj, 2, function(x) all(na.omit(x) == cummax(na.omit(x))))
	
# Remove any data points flagged as BAD according to "Reference table 2a:
# profile quality flag" in ARGO USER'S MANUAL  (i.e., NOT 1, 2, 5, or 8)
	pres <- replace(pres, presQC != 1 & presQC != 2 &
						  presQC != 5 & presQC != 8, NA)
	temp <- replace(temp, tempQC != 1 & tempQC != 2 &
						  tempQC != 5 & tempQC != 8, NA)
	psal <- replace(psal, psalQC != 1 & psalQC != 2 &
						  psalQC != 5 & psalQC != 8, NA)
	
	presAdj <- replace(presAdj, presAdjQC != 1 & presAdjQC != 2 &
								presAdjQC != 5 & presAdjQC != 8, NA)
	tempAdj <- replace(tempAdj, tempAdjQC != 1 & tempAdjQC != 2 &
								tempAdjQC != 5 & tempAdjQC != 8, NA)
	psalAdj <- replace(psalAdj, psalAdjQC != 1 & psalAdjQC != 2 & 
								psalAdjQC != 5 & psalAdjQC != 8, NA)

# Remove observations below 2100 m
	pres <- replace(pres, pres > 2100, NA)
	temp <- replace(temp, pres > 2100, NA)
	psal <- replace(psal, pres > 2100, NA)
	
	presAdj <- replace(presAdj, presAdj > 2100, NA)
	tempAdj <- replace(tempAdj, presAdj > 2100, NA)
	psalAdj <- replace(psalAdj, presAdj > 2100, NA)
	
# Extract provided profile QC flags
	profilePresQC <- as.character(unlist(strsplit(ncvar_get(nc,
	                              "PROFILE_PRES_QC"), split="")))
	profileTempQC <- as.character(unlist(strsplit(ncvar_get(nc,
	                              "PROFILE_TEMP_QC"), split="")))
	profilePsalQC <- as.character(unlist(strsplit(ncvar_get(nc,
	                              "PROFILE_PSAL_QC"), split="")))

# =============================================================================
# Observed quirk: In file P_20140109_prof.nc, variables "PROFILE_TEMP_QC" and
# "PROFILE_PSAL_QC" have length==169 for 216 profiles. These should match. No
# fill values or space holders are provided. The following steps address this
# fluky situation by replacing all profile QC flags with NA to be filtered out
# in a later routine.

if(length(profilePresQC) != profilesThisFile){
	profilePresQC <- rep(NA, profilesThisFile)}

if(length(profileTempQC) != profilesThisFile){
	profileTempQC <- rep(NA, profilesThisFile)}

if(length(profilePsalQC) != profilesThisFile){
	profilePsalQC <- rep(NA, profilesThisFile)}
	
# =============================================================================
# CALCULATE profile QC flags as percentage of GOOD data using data point flags
#
# Due to noted inconsistencies in the way the provided profile QC flags are
# calculated, new alpha flags that quantify the percentage of good data points
# within the profile are calculated here according to the methods described in 
# "Reference table 2a: profile quality flag" in ARGO USER'S MANUAL:
#
#		A: n = 100%	- All profile levels contain good data
#		B: 75% ≤ n < 100%
#		C: 50% ≤ n < 75%
#		D: 25% ≤ n < 50%
#		E: 0% ≤ n < 25%
#		F: n = 0% - No profile levels have good data
#
# Calculate the percentage of GOOD data in the profiles
    calculatedProfilePresQC <- flagCalc(presQC)
    calculatedProfileTempQC <- flagCalc(tempQC)
    calculatedProfilePsalQC <- flagCalc(psalQC)
    calculatedProfilePresAdjQC <- flagCalc(presAdjQC)
    calculatedProfileTempAdjQC <- flagCalc(tempAdjQC)
    calculatedProfilePsalAdjQC <- flagCalc(psalAdjQC)
    rm(list=c("presAdjQC", "tempAdjQC", "psalAdjQC",
              "presQC", "tempQC", "psalQC"))

# Simulated density: check for NA or NaN
# If pressure, temperature, and salinity readings do not overlap at every depth
# level in the profile, density will be NA at that depth. This step flags such
# a scenario. First, set a density QC flag to be the lowest grade (A through F)
# between the three state variables. Then change any flag to F if density
# cannot be calculated for the profile.
    densSim <- pres
    densSim[is.na(temp)] <- NA
    densSim[is.na(psal)] <- NA
    calculatedProfileDensQC <- apply(rbind(calculatedProfilePresQC,
                                           calculatedProfileTempQC,
                                           calculatedProfilePsalQC), 2, max)
    calculatedProfileDensQC[colSums(!is.na(densSim)) == 0] <- "F"
    
    densAdjSim <- presAdj
    densAdjSim[is.na(tempAdj)] <- NA
    densAdjSim[is.na(psalAdj)] <- NA
    calculatedProfileDensAdjQC <- apply(rbind(calculatedProfilePresAdjQC,
                                              calculatedProfileTempAdjQC,
                                              calculatedProfilePsalAdjQC),
                                              2, max)
    calculatedProfileDensAdjQC[colSums(!is.na(densSim)) == 0] <- "F"
    
# Minimum and maximum profile depths
# Note: this is done separately for real-time (rt) and delayed-mode (dm)
	rtMinDepth <- apply(as.matrix(pres), 2, na.min)
	dmMinDepth <- apply(as.matrix(presAdj), 2, na.min)
	profMinDepth <- rep(NA, profilesThisFile)
	profMinDepth[dataMode == "R"] <- rtMinDepth[dataMode == "R"]
	profMinDepth[dataMode=="D" | dataMode=="A"] <- dmMinDepth[dataMode=="D" |
	                                                          dataMode=="A"]
	rtMaxDepth <- apply(as.matrix(pres), 2, na.max)
    dmMaxDepth <- apply(as.matrix(presAdj), 2, na.max)
    profMaxDepth <- rep(NA, profilesThisFile)
    profMaxDepth[dataMode == "R"] <- rtMaxDepth[dataMode == "R"]
    profMaxDepth[dataMode == "D" | dataMode == "A"] <- dmMaxDepth[dataMode == "D" |
                                                              dataMode == "A"]

# Determine whether the profile has at least 10 data points
# (any fewer than this will be difficult to regress and not worth retaining)
    rtNumBins <- apply(rbind(colSums(!is.na(pres)),
                             colSums(!is.na(temp)),
                             colSums(!is.na(psal)),
                             colSums(!is.na(densSim))), 2, min)
    dmNumBins <- apply(rbind(colSums(!is.na(presAdj)),
                             colSums(!is.na(tempAdj)),
                             colSums(!is.na(psalAdj)),
                             colSums(!is.na(densAdjSim))), 2, min)
    profNumBins <- rep(NA, profilesThisFile)
    profNumBins[dataMode=="R"] <- rtNumBins[dataMode=="R"]
    profNumBins[dataMode=="D" | dataMode=="A"] <- dmNumBins[dataMode=="D" |
                                                            dataMode=="A"]
    lessThanTen <- profNumBins < 10
    
# Clean up vars no longer needed
    rm(list=c("pres", "presAdj", "temp", "tempAdj", "psal", "psalAdj",
              "densSim", "densAdjSim"))

# Time and time QC
    juld <- ncvar_get(nc, "JULD")
    juldQC <- as.numeric(unlist(strsplit(ncvar_get(nc, "JULD_QC"), split="")))

# Position and position QC
    lat <- ncvar_get(nc, "LATITUDE")
    lon <- ncvar_get(nc, "LONGITUDE")
    posQC <- as.numeric(unlist(strsplit(ncvar_get(nc, "POSITION_QC"), split="")))

# =============================================================================
# COMPILE profile info and QC flags from each netCDF file
	qc <- rbind(qc, data.frame(fileName = rep(fname, length(juld)),
	                           profInd = seq(1:length(juld)),
	                           rtMinDepth = rtMinDepth,
	                           dmMinDepth = dmMinDepth,
	                           profMinDepth = profMinDepth,
	                           rtMaxDepth = rtMaxDepth,
	                           dmMaxDepth = dmMaxDepth,
	                           profMaxDepth = profMaxDepth,
	                           rtNumBins = rtNumBins,
	                           dmNumBins = dmNumBins,
	                           profNumBins = profNumBins,
	                           juld = juld, juldQC = juldQC,
	                           lat = lat, lon = lon,
	                           posQC = posQC,
	                           dataMode = dataMode,
	                           dataState = dataState,
	                           origProfPresQC = profilePresQC,
	                           origProfTempQC = profileTempQC,
	                           origProfPsalQC = profilePsalQC,
	                           calcProfPresQC = calculatedProfilePresQC,
	                           calcProfTempQC = calculatedProfileTempQC,
	                           calcProfPsalQC = calculatedProfilePsalQC,
	                           calcProfDensQC = calculatedProfileDensQC,
	                           calcProfPresAdjQC = calculatedProfilePresAdjQC,
	                           calcProfTempAdjQC = calculatedProfileTempAdjQC,
	                           calcProfPsalAdjQC = calculatedProfilePsalAdjQC,
	                           calcProfDensAdjQC = calculatedProfileDensAdjQC,
	                           presInversion = presInv,
	                           presAdjInversion = presAdjInv,
	                           tooFewBins = lessThanTen,
	                           tooShallow = rep(FALSE, profilesThisFile),
	                           noSurfaceObs = rep(FALSE, profilesThisFile),
	                           stringsAsFactors = FALSE))
    
    rm(list=c("dataMode", "dataState", "juld", "juldQC", "lat", "lon", "posQC",
              "rtMinDepth", "dmMinDepth", "rtMaxDepth", "dmMaxDepth",
              "profMinDepth", "profMaxDepth", "rtNumBins", "dmNumBins",
              "profNumBins", "presInv", "presAdjInv", "lessThanTen",
              "profilePresQC", "profilePsalQC", "profileTempQC",
              "calculatedProfilePresQC", "calculatedProfileTempQC",
              "calculatedProfilePsalQC", "calculatedProfileDensQC",
              "calculatedProfilePresAdjQC", "calculatedProfileTempAdjQC",
              "calculatedProfilePsalAdjQC", "calculatedProfileDensAdjQC"))
	
	}	# Close if statement (check for all variables)

# Close the current netCDF file
	nc_close(nc)
	
# Print update
	if((f %% 1000) == 0){
	    message(paste0(f, " files complete. ", length(files)-f, " to go."))
	}    

	}  # End open file loop

# WRITE out files
message("Writing file 'meta_all_profiles.RData'")
save(qc, file="~/round_3/meta_all_profiles.RData")

# =============================================================================

# Print summary in ".out" file
message("SUMMARY:")
message(paste0("Total number of files downloaded: ", length(files)))
message(paste0("Total number of profiles downloaded: ", numberOfProfiles))
message(paste0("Total number of profiles containing all variables: ", nrow(qc)))
message("DONE")
system("date")

# =============================================================================
#			    					END
# =============================================================================
