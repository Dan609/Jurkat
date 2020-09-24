# Cell movement trajectory analysis by Dan Bobkov, 2020 # dan.bobkov@gmail.com
# Script desined to analyze tracks collected in ImageJ Manual tracking
# plugin from CQ1 cytometer derived images, recorded during 4 h experiment
# with Jurcat-HUVEC cells conditioned by co-cultivation and MbCD treatment.
# Images were collected from SPIONs and conrtol wells, to model the classical
# replicative senescense effects.
# If you use this script in your publications, please cite Bobkov D., Polyanskaya A., 
# Musorina A., Lomert E., Shabelnikov S., Poljanskaya G. Replicative senescence in MSCWJ-1 
# human umbilical cord mesenchymal stem cells is marked by characteristic changes in motility, 
# cytoskeletal organization, and RhoA localization. Molecular Biology Reports. 2020. 47(5):3867-3883. 
# doi: 10.1007/s11033-020-05476-6
# This script based on trajr package
# If you use trajr in your publications, please cite McLean DJ, Skowron Volponi MA. trajr:
# An R package for characterisation of animal trajectories. Ethology. 2018;12739.
# https://doi.org/10.1111/eth.12739.
#
#################################################
################## Set workspace  ###############
#################################################
# Load libraries-------------
library(PerformanceAnalytics)
library(GGally)
library(trajr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(plyr)
library(xlsx)
library(car)
library(rgl)
library(dunn.test)
library(DescTools)
library(FactoMineR)
library(factoextra)
library(Hmisc)
library(gplots)
library(ggsignif)
library(stringi)
library(PMCMRplus)
library(pca3d)
library(xtable)
library(DAAG)
library(cluster)
library(agricolae)
library(RVAideMemoire)
library(coin)
#################################################
################## Initialization ###############
#################################################

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE);
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)

file.copy(from=plots.png.paths, to="D:/YandexDisk/ExData/RAFTs/motility/R_script/png")

plots.png.detials <- file.info(plots.png.paths)
plots.png.detials <- plots.png.detials[order(plots.png.detials$mtime),]

sorted.png.names <- gsub(plots.dir.path, "D:/YandexDisk/ExData/RAFTs/motility/R_script/png", row.names(plots.png.detials), fixed=TRUE)
numbered.png.names <- paste0("D:/YandexDisk/ExData/RAFTs/motility/R_script/png", 1:length(sorted.png.names), ".png")

plots.dir.path


# Set names---------------

mynames <- c("track",
             "length",
             "distance",
             "straight",
             "square_displacement",
             "mean_speed",
             "sd_speed",
             "max_speed",
             "min_speed",
             "sinuosity",
             "emax",
             "DC",
             "SDDC",
             "mean_angle",
             "probe")

alltracks <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                      c("track",
                         "length",
                         "distance",
                         "straight",
                         "square_displacement",
                         "mean_speed",
                         "sd_speed",
                         "max_speed",
                         "min_speed",
                         "sinuosity",
                         "emax",
                         "DC",
                         "SDDC",
                         "mean_angle",
                         "probe"))
#tracks_W1J <- alltracks
tracks_W2H <- alltracks
tracks_W2J <- alltracks
tracks_W3H <- alltracks
#tracks_W4J <- alltracks
tracks_W5H <- alltracks
tracks_W5J <- alltracks
tracks_W6H <- alltracks
# tracks <- alltracks #?
# Remove outliers function------
###
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

#data$mean_speed <- remove_outliers(data$mean_speed)
#################################################
##### Trajectory Analisys Functions #############
#################################################
# Load Track analysis functions



traj_analysis_W2H <- function(input) {
  
  data <- read.csv(input)
  
  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))
  
  for (i in unique(data$Track)) {
    
    # Define x, y, and time coordinates
    coords <- data.frame(x = data$X[data$Track == i],
                         y = data$Y[data$Track == i],
                         # times = c(1:96))
                         timeCol = data$Slice[data$Track == i],
                         spatialUnits = "pixels", timeUnits = "hours")
    
    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 90/3600/3) #fps frames/3600/hours = 1/120
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 0.6579 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 0.4273504 / 1, "micrometer") # 234 px = 100 microns
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    # trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)
    
    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj), # D/L ratio
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           #mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           mean_angle = 0,
                           #probe = 'W1'
                           probe = 2
    )
    
    head(traj_params)
    
  }
  
  # print(traj_params)
  write.csv(traj_params, file = 'traj_W2H.csv')
  tracks <<- traj_params
  return(traj_params) } # 90 frames

traj_analysis_W2J <- function(input) {
  
  data <- read.csv(input)
  
  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))
  
  for (i in unique(data$Track)) {
    
    # Define x, y, and time coordinates
    coords <- data.frame(x = data$X[data$Track == i],
                         y = data$Y[data$Track == i],
                         # times = c(1:96))
                         timeCol = data$Slice[data$Track == i],
                         spatialUnits = "pixels", timeUnits = "hours")
    
    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 90/3600/3) #fps frames/3600/hours = 1/120
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 0.6579 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 0.4273504 / 1, "micrometer") # 234 px = 100 microns
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    # trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2) # xlim = grDevices::extendrange(x$x), ylim = grDevices::extendrange(x$y)
    
    # ,xlim = c(0,500), ylim = c(0,500)
    
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)
    
    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj), # D/L ratio
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           #mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           mean_angle = 0,
                           #probe = 'W1'
                           probe = 3
    )
    
    head(traj_params)
    
  }
  
  # print(traj_params)
  write.csv(traj_params, file = 'traj_W2J.csv')
  tracks <<- traj_params
  return(traj_params) } # 90 frames

traj_analysis_W3H <- function(input) {
  
  data <- read.csv(input)
  
  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))
  
  for (i in unique(data$Track)) {
    
    # Define x, y, and time coordinates
    coords <- data.frame(x = data$X[data$Track == i],
                         y = data$Y[data$Track == i],
                         # times = c(1:96))
                         timeCol = data$Slice[data$Track == i],
                         spatialUnits = "pixels", timeUnits = "hours")
    
    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 90/3600/3) #fps frames/3600/hours = 1/120
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 0.6579 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 0.4273504 / 1, "micrometer") # 234 px = 100 microns
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    # trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)
    
    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj), # D/L ratio
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           #mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           mean_angle = 0,
                           #probe = 'W1'
                           probe = 4
    )
    
    head(traj_params)
    
  }
  
  # print(traj_params)
  write.csv(traj_params, file = 'traj_W3H.csv')
  tracks <<- traj_params
  return(traj_params) } # 90 frames



traj_analysis_W5H <- function(input) {
  
  data <- read.csv(input)
  
  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))
  
  for (i in unique(data$Track)) {
    
    # Define x, y, and time coordinates
    coords <- data.frame(x = data$X[data$Track == i],
                         y = data$Y[data$Track == i],
                         # times = c(1:96))
                         timeCol = data$Slice[data$Track == i],
                         spatialUnits = "pixels", timeUnits = "hours")
    
    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 90/3600/3) #fps frames/3600/hours = 1/120
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 0.6579 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 0.4273504 / 1, "micrometer") # 234 px = 100 microns
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    # trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)
    
    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj), # D/L ratio
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           #mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           mean_angle = 0,
                           #probe = 'W1'
                           probe = 6
    )
    
    head(traj_params)
    
  }
  
  # print(traj_params)
  write.csv(traj_params, file = 'traj_W5H.csv')
  tracks <<- traj_params
  return(traj_params) } # 90 frames

traj_analysis_W5J <- function(input) {
  
  data <- read.csv(input)
  
  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))
  
  for (i in unique(data$Track)) {
    
    # Define x, y, and time coordinates
    coords <- data.frame(x = data$X[data$Track == i],
                         y = data$Y[data$Track == i],
                         # times = c(1:96))
                         timeCol = data$Slice[data$Track == i],
                         spatialUnits = "pixels", timeUnits = "hours")
    
    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 90/3600/3) #fps frames/3600/hours = 1/120
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 0.6579 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 0.4273504 / 1, "micrometer") # 234 px = 100 microns
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    # trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)
    
    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj), # D/L ratio
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           #mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           mean_angle = 0,
                           #probe = 'W1'
                           probe = 7
    )
    
    head(traj_params)
    
  }
  
  # print(traj_params)
  write.csv(traj_params, file = 'traj_W5J.csv')
  tracks <<- traj_params
  return(traj_params) } # 90 frames

traj_analysis_W6H <- function(input) {
  
  data <- read.csv(input)
  
  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))
  
  for (i in unique(data$Track)) {
    
    # Define x, y, and time coordinates
    coords <- data.frame(x = data$X[data$Track == i],
                         y = data$Y[data$Track == i],
                         # times = c(1:96))
                         timeCol = data$Slice[data$Track == i],
                         spatialUnits = "pixels", timeUnits = "hours")
    
    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 90/3600/3) #fps frames/3600/hours = 1/120
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 0.6579 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 0.4273504 / 1, "micrometer") # 234 px = 100 microns
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    # trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)
    
    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj), # D/L ratio
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           #mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           mean_angle = 0,
                           #probe = 'W1'
                           probe = 8
    )
    
    head(traj_params)
    
  }
  
  # print(traj_params)
  write.csv(traj_params, file = 'traj_W6H.csv')
  tracks <<- traj_params
  return(traj_params) } # 90 frames

# Choose dir function for mac -----------------
choose.dir <- function() {
  system("osascript -e 'tell app \"R\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
         intern = FALSE, ignore.stderr = TRUE)
  p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
  return(ifelse(length(p), p, NA))
}

file_list_W2H <- list.files(path = , choose.dir(),
                            pattern = "csv",
                            all.files = FALSE,
                            full.names = TRUE, recursive = TRUE,
                            ignore.case = FALSE, include.dirs = FALSE,
                            no.. = FALSE)

file_list_W2J <- list.files(path = , choose.dir(),
                            pattern = "csv",
                            all.files = FALSE,
                            full.names = TRUE, recursive = TRUE,
                            ignore.case = FALSE, include.dirs = FALSE,
                            no.. = FALSE)

file_list_W3H <- list.files(path = , choose.dir(),
                            pattern = "csv",
                            all.files = FALSE,
                            full.names = TRUE, recursive = TRUE,
                            ignore.case = FALSE, include.dirs = FALSE,
                            no.. = FALSE)


file_list_W5H <- list.files(path = , choose.dir(),
                            pattern = "csv",
                            all.files = FALSE,
                            full.names = TRUE, recursive = TRUE,
                            ignore.case = FALSE, include.dirs = FALSE,
                            no.. = FALSE)


file_list_W5J <- list.files(path = , choose.dir(),
                            pattern = "csv",
                            all.files = FALSE,
                            full.names = TRUE, recursive = TRUE,
                            ignore.case = FALSE, include.dirs = FALSE,
                            no.. = FALSE)

file_list_W6H <- list.files(path = , choose.dir(),
                            pattern = "csv",
                            all.files = FALSE,
                            full.names = TRUE, recursive = TRUE,
                            ignore.case = FALSE, include.dirs = FALSE,
                            no.. = FALSE)


test_list <- list.files(path = , choose.dir(),
                            pattern = "csv",
                            all.files = FALSE,
                            full.names = TRUE, recursive = TRUE,
                            ignore.case = FALSE, include.dirs = FALSE,
                            no.. = FALSE)

#################################################
##### Tracks importing from local dirs ##########
#################################################
# Choose dir ---------W1
#file_list_W1J <- list.files(path = , choose.dir(default = "",
#                                            caption = "Select folder"),
#                        pattern = "csv",
#                        all.files = FALSE,
#                        full.names = TRUE, recursive = TRUE,
#                        ignore.case = FALSE, include.dirs = FALSE,
#                        no.. = FALSE)

file_list_W2H <- list.files(path = , choose.dir(default = "",
                                            caption = "Select folder"),
                        pattern = "csv",
                        all.files = FALSE,
                        full.names = TRUE, recursive = TRUE,
                        ignore.case = FALSE, include.dirs = FALSE,
                        no.. = FALSE)

file_list_W2J <- list.files(path = , choose.dir(default = "",
                                               caption = "Select folder"),
                           pattern = "csv",
                           all.files = FALSE,
                           full.names = TRUE, recursive = TRUE,
                           ignore.case = FALSE, include.dirs = FALSE,
                           no.. = FALSE)

file_list_W3H <- list.files(path = , choose.dir(default = "",
                                               caption = "Select folder"),
                           pattern = "csv",
                           all.files = FALSE,
                           full.names = TRUE, recursive = TRUE,
                           ignore.case = FALSE, include.dirs = FALSE,
                           no.. = FALSE)

#file_list_W4J <- list.files(path = , choose.dir(default = "",
#                                               caption = "Select folder"),
#                           pattern = "csv",
#                           all.files = FALSE,
#                           full.names = TRUE, recursive = TRUE,
#                           ignore.case = FALSE, include.dirs = FALSE,
#                           no.. = FALSE)

file_list_W5H <- list.files(path = , choose.dir(default = "",
                                               caption = "Select folder"),
                           pattern = "csv",
                           all.files = FALSE,
                           full.names = TRUE, recursive = TRUE,
                           ignore.case = FALSE, include.dirs = FALSE,
                           no.. = FALSE)


file_list_W5J <- list.files(path = , choose.dir(default = "",
                                                caption = "Select folder"),
                            pattern = "csv",
                            all.files = FALSE,
                            full.names = TRUE, recursive = TRUE,
                            ignore.case = FALSE, include.dirs = FALSE,
                            no.. = FALSE)

file_list_W6H <- list.files(path = , choose.dir(default = "",
                                                caption = "Select folder"),
                            pattern = "csv",
                            all.files = FALSE,
                            full.names = TRUE, recursive = TRUE,
                            ignore.case = FALSE, include.dirs = FALSE,
                            no.. = FALSE)

#################################################
############# Call to function ##################
#################################################
for (file_name in test_list) {
  traj_analysis_W2H(file_name)
  tracks_W2H <- rbind(tracks_W2H, tracks)
}


# Start scan for W2H
for (file_name in file_list_W2H) {
  traj_analysis_W2H(file_name)
  tracks_W2H <- rbind(tracks_W2H, tracks)
}


# Start scan for W2J
for (file_name in file_list_W2J) {
  traj_analysis_W2J(file_name)
  tracks_W2J <- rbind(tracks_W2J, tracks)
}

# Start scan for W3H
for (file_name in file_list_W3H) {
  traj_analysis_W3H(file_name)
  tracks_W3H <- rbind(tracks_W3H, tracks)
}



# Start scan for W5H
for (file_name in file_list_W5H) {
  traj_analysis_W5H(file_name)
  tracks_W5H <- rbind(tracks_W5H, tracks)
}

# Start scan for W5J
for (file_name in file_list_W5J) {
  traj_analysis_W5J(file_name)
  tracks_W5J <- rbind(tracks_W5J, tracks)
}

# Start scan for W6H
for (file_name in file_list_W6H) {
  traj_analysis_W6H(file_name)
  tracks_W6H <- rbind(tracks_W6H, tracks)
}




#################################################
############# Prepare data ######################
#################################################
# Merge all tracks-----------------------
# alltracks <- tracks_W1

alltracks <- rbind(tracks_W6H, tracks_W3H, tracks_W2H, tracks_W5H, tracks_W2J, tracks_W5J ) # collect all tracks

summary(alltracks)
write.csv(alltracks, file = 'alltracks.csv') #save results
# Order probe levels--------------------
alltracks$probe <- as.factor(alltracks$probe)

alltracks$probe <- ordered(alltracks$probe,
                      levels = c("8", "4", "2", "6", "3", "7"))

# Set time to hours, remove tracks and mean_angle columns----------------
all.h <- alltracks
head(all.h)
all.h <- cbind(all.h[,c(6,7,8,9)]*3600, all.h[,c(2,3,4,5,10,11,12,13,15)])
data <- all.h
head(data)
tail(data)
summary(data)

# Density------------------------------ ------


data$cells <- 'NA'

data$cells[data$probe == 8] <- 'HUVEC EC medium'

data$cells[data$probe == 4] <- 'Huvec Jurkat medium'

data$cells[data$probe == 2] <- 'HUVEC control'

data$cells[data$probe == 6] <- 'HUVEC MbCD'

data$cells[data$probe == 3] <- 'Jurcat control'

data$cells[data$probe == 7] <- 'Jurcat MbCD'


data$line <- 'NA'

data$line[data$probe == 8] <- 'HUVEC'
data$line[data$probe == 4] <- 'HUVEC'
data$line[data$probe == 2] <- 'HUVEC'
data$line[data$probe == 6] <- 'HUVEC'

data$line[data$probe == 3] <- 'Jurcat'
data$line[data$probe == 7] <- 'Jurcat'


#


ggdensity(data, x = "mean_speed",
          add = "mean", rug = TRUE,
          color = "line", fill = "line",
          palette = "rainbow")


ggdensity(data[data$line=='Jurcat',], x = "mean_speed",
          add = "mean", rug = TRUE,
          color = "cells", fill = "cells",
          palette = "rainbow")

ggdensity(data[data$line=='HUVEC',], x = "mean_speed",
          add = "mean", rug = TRUE,
          color = "cells", fill = "cells",
          palette = "rainbow")
#


boxplot(mean_speed~probe, all.h)

boxplot(max_speed~probe, all.h)

boxplot(sinuosity~probe, all.h)

boxplot(straight~probe, all.h)

boxplot(DC~probe, all.h)

boxplot(emax~probe, all.h)



data$mean_speed <- remove_outliers(data$mean_speed)
boxplot(mean_speed~probe, data)
#################################################
######## Plot all correlations ##################
#################################################

# Plot all-in-one-------------------------
# plot(data)
# png()
ggpairs(data)

chart.Correlation(data[,1:12], histogram=TRUE, pch=19)

ggcorr(data, palette = "RdBu", label = TRUE)
# dev.off()

data$cells[data$probe == 8] <- 'HUVEC EC medium'

data$cells[data$probe == 4] <- 'Huvec Jurkat medium'

data$cells[data$probe == 2] <- 'HUVEC control'

data$cells[data$probe == 6] <- 'HUVEC MbCD'

data$cells[data$probe == 3] <- 'Jurcat control'

data$cells[data$probe == 7] <- 'Jurcat MbCD'

# Set X - axis names
CellSciGuylabs <- c("HUVEC\n(EC medium)", "HUVEC\n(Jurkat medium)", 
                    'HUVEC\nControl', 'HUVEC\nMbCD', 
                    'Jurcat\nControl', 'Jurcat\nMbCD')



# (01) Mean speed ----------------------
kruskal.test(data$mean_speed ~ data$probe)
compare_means(mean_speed ~ cells,  data = data, method = "wilcox.test") # pairwise comparisons

#compare_means(mean_speed ~ probe,  data = data, method = "t.test", p.adjust.method = "bonferroni")

compare_means(mean_speed ~ cells,  data = data, method = "t.test", p.adjust.method = "bonferroni")


kruskal.test(data$mean_speed ~ data$probe)


compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test")

# 1 Mean speed
ggplot(data, aes(x = probe, y = mean_speed)) +
  
  ggtitle("Mean Speed") +
  
  ylim(c(0, 500)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.12),
              cex = .99,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'Micrometers per hour',
       x = "", caption = "Kruskal-Wallis p-value = 1.675e-08") +
  
  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 90, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain")) +
  
  geom_boxplot(alpha = I(0.4), outlier.shape = NA, coef = 1.5, size = .8) +
  
  geom_signif(y_position = c(450),
              xmin = c(5),
              xmax = c(6),
              annotation = "****",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(250),
              xmin = c(1),
              xmax = c(2),
              annotation = "NS",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(490),
              xmin = c(4),
              xmax = c(6),
              annotation = "NS",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(420),
              xmin = c(3),
              xmax = c(5),
              annotation = "****",
              tip_length = 0.04) 



kruskal.test(data$DC ~ data$cells)


compare_means(DC ~ cells,  data = data, method = "wilcox.test")


compare_means(DC ~ cells,  data = data, method = "t.test", p.adjust.method = "bonferroni")

# 1 Mean speed
ggplot(data, aes(x = probe, y = DC)) +
  
  ggtitle("Directional Change") +
  
  ylim(c(0.1, .5)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.06),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'DC',
       x = "", caption = "Kruskal-Wallis p-value = 4.677e-14") +
  
  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 90, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain")) +
  
  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 1.5, size = .8) +
  
  geom_signif(y_position = c(.4),
              xmin = c(5),
              xmax = c(6),
              annotation = "NS",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(.490),
              xmin = c(4),
              xmax = c(6),
              annotation = "****",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(.440),
              xmin = c(3),
              xmax = c(5),
              annotation = "****",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(.4),
              xmin = c(1),
              xmax = c(2),
              annotation = "NS",
              tip_length = 0.04)



kruskal.test(data$length ~ data$cells)


compare_means(length ~ cells,  data = data, method = "wilcox.test")

compare_means(length ~ cells,  data = data, method = "t.test")


# 2 Length
ggplot(data, aes(x = probe, y = length)) +
  
  ggtitle("Length") +
  
  ylim(c(0, 1500)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.06),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'Micrometers',
       x = "", caption = "Kruskal-Wallis p-value = 1.993e-08") +
  
  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 90, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain")) +
  
  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 1.5, size = .8) +
  
  geom_signif(y_position = c(1000),
              xmin = c(1),
              xmax = c(2),
              annotation = "NS",
              tip_length = 0.04)  +
  
  geom_signif(y_position = c(1280),
              xmin = c(5),
              xmax = c(6),
              annotation = "****",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(1500),
              xmin = c(4),
              xmax = c(6),
              annotation = "NS",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(1350),
              xmin = c(3),
              xmax = c(5),
              annotation = "****",
              tip_length = 0.04)



kruskal.test(data$distance ~ data$cells)


compare_means(distance ~ cells,  data = data, method = "wilcox.test")

compare_means(distance ~ cells,  data = data, method = "t.test")


#
ggplot(data, aes(x = probe, y = distance)) +
  
  ggtitle("Distance") +
  
  ylim(c(0, 500)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.06),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'Micrometers',
       x = "", caption = "Kruskal-Wallis p-value = 0.01623") +
  
  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 90, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain")) +
  
  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 1.5, size = .8) +
  
  geom_signif(y_position = c(400),
              xmin = c(1),
              xmax = c(2),
              annotation = "NS",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(440),
              xmin = c(4),
              xmax = c(6),
              annotation = "**",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(360),
              xmin = c(3),
              xmax = c(5),
              annotation = "NS",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(500),
              xmin = c(6),
              xmax = c(5),
              annotation = "NS",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(500),
              xmin = c(3),
              xmax = c(4),
              annotation = "NS",
              tip_length = 0.04)


kruskal.test(data$max_speed ~ data$cells)


compare_means(max_speed ~ cells,  data = data, method = "wilcox.test")



# 2 Max speed
ggplot(data, aes(x = probe, y = max_speed)) +
  
  ggtitle("Max speed") +
  
  ylim(c(0, 2500)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.06),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'Micrometers per hour',
       x = "", caption = "Kruskal-Wallis p-value = 1.134e-06") +
  
  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 90, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain"))+
  
  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 1.5, size = .8) +
  
  geom_signif(y_position = c(1500),
              xmin = c(1),
              xmax = c(2),
              annotation = "NS",
              tip_length = 0.04)  +
  
  geom_signif(y_position = c(2220),
              xmin = c(5),
              xmax = c(6),
              annotation = "*",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(2450),
              xmin = c(4),
              xmax = c(6),
              annotation = "*",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(2000),
              xmin = c(3),
              xmax = c(5),
              annotation = "****",
              tip_length = 0.04)



  
compare_means(sinuosity ~ probe,  data = data, method = "t.test")

# 3 sinuosity

ggplot(all.h, aes(x = probe, y = sinuosity)) +
  
  ggtitle("Sinuosity") +
  
  ylim(c(0, 1)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.2),
              cex = .5,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'Sinuosity index',
       x = "Cell passage") +
  
  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 90, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain")) + 
  
  #  stat_summary(fun.data=data_summary, 
  #               geom="crossbar", 
  #               width=0.8) +
  
  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 1.5, size = .8)
  
  geom_signif(y_position = c(0.9),
              xmin = c(1),
              xmax = c(2),
              annotation = "ns",
              tip_length = 0.04)

  
  
  
  kruskal.test(data$emax ~ data$cells)
  
  
  compare_means(emax ~ cells,  data = data, method = "wilcox.test")

# 4 Emax
compare_means(emax ~ probe,  data = data, method = "t.test")

ggplot(data, aes(x = probe, y = emax)) +
  
  ggtitle("Straightness") +
  
  ylim(c(0, 7.5)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.06),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'E-max',
       x = "", caption = "Kruskal-Wallis p-value = 5.408e-13") +
  
  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 90, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain")) +
  
  #  stat_summary(fun.data=data_summary, 
  #               geom="crossbar", 
  #               width=0.8) +
  
  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 1.5, size = .8)+
  
  geom_signif(y_position = c(6),
              xmin = c(1),
              xmax = c(2),
              annotation = "NS",
              tip_length = 0.04)  +
  
  geom_signif(y_position = c(5),
              xmin = c(5),
              xmax = c(6),
              annotation = "NS",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(6),
              xmin = c(4),
              xmax = c(6),
              annotation = "****",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(6.9),
              xmin = c(3),
              xmax = c(5),
              annotation = "****",
              tip_length = 0.04)

# 5 DC

ggplot(all.h, aes(x = probe, y = DC)) +
  
  ggtitle("Directional change") +
  
  ylim(c(0.007, .048)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.2),
              cex = .5,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'DC',
       x = "Cell passage") +
  
  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain")) + 
  
  #  stat_summary(fun.data=data_summary, 
  #               geom="crossbar", 
  #               width=0.8) +
  
  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 0) +
  
  geom_signif(y_position = c(.045),
              xmin = c(1),
              xmax = c(2),
              annotation = "*",
              tip_length = 0.04) + 
  
  #  geom_signif(y_position = c(.049),
  #              xmin = c(1),
  #              xmax = c(4),
  #              annotation = "****",
  #              tip_length = 0.04) + 
  
  geom_signif(y_position = c(.042),
              xmin = c(2),
              xmax = c(3),
              annotation = "ns",
              tip_length = 0.04) 


#install.packages('agricolae')
Median.test(data$mean_speed, # ????????? ???????? ???????????
            data$probe,
            alpha = 0.05, 
            correct = TRUE,
            simulate.p.value = FALSE,
            group = TRUE, 
            main = NULL,
            console = TRUE)



# install.packages('RVAideMemoire')
mood.medtest(data$mean_speed, 
             data$probe,
             exact = FALSE)

mood.medtest(data$mean_speed, 
             data$probe,
             exact = TRUE)




#################################################
######## Investigate for normality ##############
#################################################

# QQ plot--------------------------
# ggqqplot(data$mean_speed, main = 'mean_speed')
ggqqplot(data$mean_speed, main = 'mean_speed')

ggqqplot(data[data$probe == '7',]$mean_speed, main = 'mean_speed')

#ggqqplot(data$sd_speed, main = 'sd_speed')
ggqqplot(data$max_speed, main = 'max_speed')
#ggqqplot(data$min_speed, main = 'min_speed')
ggqqplot(data$length, main = 'length')
ggqqplot(data$distance, main = 'distance')
ggqqplot(data$straight, main = 'straight')
ggqqplot(data$square_displacement, main = 'square_displacement')
ggqqplot(data$sinuosity, main = 'sinuosity')
ggqqplot(data$emax, main = 'emax')
ggqqplot(data$DC, main = 'DC')
ggqqplot(data$SDDC, main = 'SDDC')
ggqqplot(data[data$probe=='1',]$mean_speed, main = 'mean_speed')
ggqqplot(data[data$probe=='2',]$mean_speed, main = 'mean_speed')

# qq plot
#qqnorm(data$mean_speed, pch = 1, frame = FALSE)
#qqline(data$mean_speed, col = "steelblue", lwd = 2)
#qqplot(data$mean_speed)

# Shapiro test --------------------------
# Test for normality----------------------------

shapiro.test(data[data$probe=='2',]$mean_speed)
shapiro.test(data[data$probe=='3',]$mean_speed) # data is normally distributed
shapiro.test(data[data$probe=='7',]$mean_speed)




# Plot means
plotmeans(mean_speed ~ probe, data = data, frame = FALSE, 
          mean.labels=TRUE, connect=FALSE,
          n.label=TRUE, text.n.label="n = ",
          xlab = "Cells", ylab = "Micrometers per hour",
          main="Jurcat-HUVEC co-cultivation\n Mean Speed (Mean and 95% CI)", 
          ylim = c(0, 180))

plotmeans(mean_speed ~ probe, data = data, frame = FALSE, 
          mean.labels=FALSE, connect=FALSE,
          n.label=TRUE, text.n.label="n = ",
          xlab = "Cells", ylab = "Micrometers per hour",
          main="Jurcat-HUVEC co-cultivation\n Mean Speed (Mean and 95% CI)", 
          ylim = c(0, 180))

plotmeans(max_speed ~ probe, data = data, frame = FALSE, 
          mean.labels=FALSE, connect=TRUE,
          n.label=TRUE, text.n.label="n = ",
          xlab = "Cells", ylab = "Micrometers per hour",
          main="FetMSC 24 h Max Speed, 
          \nMean Plot with 95% CI", 
          ylim = c(0, 1200))

plotmeans(sinuosity ~ probe, data = data, frame = FALSE, 
          mean.labels=FALSE, connect=TRUE,
          n.label=TRUE, text.n.label="n = ",
          xlab = "Cells", ylab = "Sinuosity index",
          main="FetMSC 24 h Sinuosity, 
          \nMean Plot with 95% CI", 
          ylim = c(0, 1))

plotmeans(straight ~ probe, data = data, frame = FALSE, 
          mean.labels=FALSE, connect=TRUE,
          n.label=TRUE, text.n.label="n=",
          xlab = "Cells", ylab = "Straightness index",
          main="FetMSC 24 h Straightness, 
          \nMean Plot with 95% CI", 
          ylim = c(0, 1))

#################################################
### Compute mean and sd for track parameters ####
#################################################
# Mean and sd --------------------------

mean(data[data$probe=='3',]$mean_speed)
sd(data[data$probe=='3',]$mean_speed)

mean(data[data$probe=='7',]$mean_speed)
sd(data[data$probe=='7',]$mean_speed)

mean(data[data$probe=='3',]$max_speed)
sd(data[data$probe=='3',]$max_speed)


mean(data[data$probe=='7',]$max_speed)
sd(data[data$probe=='7',]$max_speed)


mean(data[data$probe=='3',]$distance)
sd(data[data$probe=='3',]$distance)


mean(data[data$probe=='7',]$distance)
sd(data[data$probe=='7',]$distance)


mean(data[data$probe=='3',]$length)
sd(data[data$probe=='3',]$length)


mean(data[data$probe=='7',]$length)
sd(data[data$probe=='7',]$length)


#

mean(data[data$probe=='3',]$DC)
sd(data[data$probe=='3',]$DC)


mean(data[data$probe=='7',]$DC)
sd(data[data$probe=='7',]$DC)

#
mean(data[data$probe=='3',]$emax)
sd(data[data$probe=='3',]$emax)

mean(data[data$probe=='7',]$emax)
sd(data[data$probe=='7',]$emax)



#  -----------
# compare means
compare_means(mean_speed ~ probe,  data = data, method = "t.test")
compare_means(straight ~ probe,  data = data, method = "t.test")


compare_means(mean_speed ~ probe,  data = data, method = "anova")
compare_means(mean_speed ~ probe,  data = data, method = "kruskal.test")
#write.xlsx(compare_means(mean_speed ~ probe,  data = data, method = "kruskal.test"),
#           file = 'kruskal.test.mean_speed.xlsx')

compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
compare_means(straight ~ probe,  data = data, method = "wilcox.test")

# (01) Mean speed ----------------------
kruskal.test(data$mean_speed ~ data$probe)
compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

qplot(probe, mean_speed, data = data,
      geom = c("jitter", "boxplot"),
      alpha = I(0.3), fill = probe,
      main = "Mean speed") +
  labs(y = 'Micrometers per hour',
       x = "Cell passage") +
  
  labs(title = "Mean speed", caption = "*Wilcoxon p-value") +
# xmin / xmax positions should match the x-axis labels' positions
geom_signif(y_position = c(170),
            xmin = c(1),
            xmax = c(2),
            annotation = "p = 0.036",
            tip_length = 0.04)


qplot(probe, mean_speed, data = data,
      geom = c("boxplot"),
      alpha = I(0.3), fill = probe,
      main = "Mean speed") +
  labs(y = 'Micrometers per hour',
       x = "Cell passage") +
  
  labs(title = "Mean speed", caption = "*Wilcoxon p-value") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(170),
              xmin = c(1),
              xmax = c(2),
              annotation = "p = 0.036",
              tip_length = 0.04)



# stat_compare_means(label.y = 3)

kruskal.test(data$max_speed ~ data$probe)
compare_means(max_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (02) Max speed ----------------------
qplot(probe, max_speed, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Max speed" ) +
  labs(y = 'Micrometers',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(550),
              xmin = c(1),
              xmax = c(2),
              annotation = "****",
              tip_length = 0.04) +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(500),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04) + theme_classic(base_size=14) + 
  
  # annotation_custom(grob) +
  
  labs(title = "Max speed",
       # subtitle = "",
       caption = "Kruskal-Wallis p-value < 2.2e-16")


kruskal.test(data$length ~ data$probe)
compare_means(length ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (03) Length ----------------------
qplot(probe, length, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Length" ) +
  labs(y = 'Micrometers',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(3500),
              xmin = c(1),
              xmax = c(2),
              annotation = "ns",
              tip_length = 0.04) +
  
  labs(title = "Length",
       # subtitle = "",
       caption = "Kruskal-Wallis p-value = 0.0557")

kruskal.test(data$distance ~ data$probe)
compare_means(distance ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (04) Distance ----------------------
qplot(probe, distance, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Distance" ) +
  labs(y = 'Micrometers',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1400),
              xmin = c(1),
              xmax = c(2),
              annotation = "NS",
              tip_length = 0.04) +

  
  # annotation_custom(grob) +
  
  labs(title = "Distance",
       # subtitle = "",
       caption = "Kruskal-Wallis p-value = 0.05576")


kruskal.test(data$straight ~ data$probe)
compare_means(straight ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (05) Straightness ----------------------
qplot(probe, straight, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Straightness", ) +
  labs(y = 'Straightness index',
       x = "Cell passage")  +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(1.1),
              xmin = c(1),
              xmax = c(2),
              annotation = "*",
              tip_length = 0.04) +
  
  theme_classic(base_size=14) + 
  
  # annotation_custom(grob) +
  
  labs(title = "Straightness",
       # subtitle = "",
       caption = "Kruskal-Wallis p-value = 0.01061")

kruskal.test(data$sinuosity ~ data$probe)
compare_means(sinuosity ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (6) Sinuosity------------------------
# Add jitter and change fill color by probe
qplot(probe, sinuosity, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Random search path tortuosity", ) +
  labs(y = 'Sinuosity index',
       x = "Cell passage") +
  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(0.7),
              xmin = c(1),
              xmax = c(2),
              annotation = "****",
              tip_length = 0.04) +
  
  theme_classic(base_size=14) + 
  
  # annotation_custom(grob) +
  
  labs(title = "Sinuosity",
       # subtitle = "",
       caption = "Kruskal-Wallis p-value < 2.2e-16")

kruskal.test(data$emax ~ data$probe)
compare_means(emax ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (7) Emax ------------------------
# Emax-a is a dimensionless, scale-independent measure of the maximum possible expected displacement
qplot(probe, emax, data = data,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Trajectory straightness index, E-max", ) +
  labs(y = 'emax',
       x = "Cell passage") +
  theme_classic(base_size=14) + 
  
  theme_classic(base_size=14) + 
  
  labs(title = "Straightness index, E-max",
       caption = "Kruskal-Wallis p-value = 0.04776")


kruskal.test(data$square_displacement ~ data$probe)
compare_means(square_displacement ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (9) square_displacement ------------------------
qplot(probe, square_displacement, data = data[data$probe != "p36",],
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "square_displacement", ) +
  labs(y = 'square_displacement',
       x = "Cell passage") +
  theme_classic(base_size=14) + 
  labs(title = "square_displacement",
       caption = "Kruskal-Wallis p-value < 2.2e-16")

kruskal.test(data$DC ~ data$probe)
compare_means(DC ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (10) DC ------------------------
qplot(probe, DC, data = data[data$probe != "p36",],
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "DC", ) +
  labs(y = 'DC',
       x = "Cell passage") +
  theme_classic(base_size=14) + 
  labs(title = "DC",
       caption = "Kruskal-Wallis p-value < 2.2e-16")

kruskal.test(data[data$probe != "p36",]$SDDC ~ data[data$probe != "p36",]$probe)
compare_means(SDDC ~ probe,  data = data[data$probe != "p36",], method = "wilcox.test") # pairwise comparisons
# (11) SDDC ------------------------
qplot(probe, SDDC, data = data[data$probe != "p36",],
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "SDDC", ) +
  labs(y = 'SDDC',
       x = "Cell passage") +
  theme_classic(base_size=14) + 
  labs(title = "SDDC",
       caption = "Kruskal-Wallis p-value < 2.2e-16")

kruskal.test(data$min_speed ~ data$probe)
compare_means(min_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons
# (12) min_speed ------------------------
qplot(probe, min_speed, data = data[data$probe != "p36",],
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "min_speed", ) +
  labs(y = 'min_speed',
       x = "Cell passage") +
  theme_classic(base_size=14) + 
  labs(title = "min_speed",
       caption = "Kruskal-Wallis p-value < 2.2e-16")


#################################################
################ Plot 3D images #################
#################################################

## Plot data in three factor space
scatter3d(x = data$mean_speed, y = data$sinuosity, z = data$straight, 
          groups = data$probe,
          xlab=deparse(substitute(mean_speed)), 
          ylab=deparse(substitute(Sinuosity)),
          zlab=deparse(substitute(Straightness)), 
          axis.scales=FALSE, axis.ticks=TRUE,
          surface=FALSE, 
          ellipsoid = TRUE,
          level=0.8, ellipsoid.alpha=0.5, id=FALSE,
          model.summary=FALSE)


#library(lattice)
#cloud(data$mean_speed ~ data$straight*data$sinuosity)

#################################################
########### Fit statistical models ##############
#################################################

# fit generalized linear models
fit <- glm(probe ~ mean_speed + sinuosity + straight, data,
           family = "binomial")
anova(fit, test='Chisq')
newobject1 <- anova(fit, test='Chisq')
print(xtable(newobject1, type = "latex"), file = "filename1.tex")

fit <- glm(probe ~ mean_speed + sinuosity + straight, data,
           family = "quasibinomial")
anova(fit, test='Chisq')
newobject2 <- anova(fit, test='Chisq')
print(xtable(newobject2, type = "latex"), file = "filename2.tex")

fit <- glm(probe ~ mean_speed + length + straight + sinuosity + 
             mean_speed * sinuosity + 
             mean_speed * length + 
             straight * sinuosity +
             mean_speed * straight * sinuosity , data,
           family = "binomial")

library(MuMIn)
options(na.action = "na.fail") # change the default "na.omit" to prevent models
# from being fitted to different datasets in
# case of missing values.
globalmodel <- lm(probe ~ length + straight + sinuosity, data = data)
combinations <- dredge(globalmodel)
print(combinations)
coefTable(combinations)

#write.xlsx(compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test"),
#           file = 'wilcox.test.mean_speed.xlsx')


####### unused ggplot elements
# Compute the analysis of variance------

res.aov <- aov(mean_speed ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)
# Straightness index
# Compute the analysis of variance
res.aov <- aov(straight ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)

# Compute the analysis of variance------
res.aov <- aov(sinuosity ~ probe, data = data) # One-way ANOVA
summary(res.aov)
TukeyHSD(res.aov)
plot(TukeyHSD(res.aov), las = 1)
kruskal.test(data$straight ~ data$probe)


compare_means(mean_speed ~ probe,  data = data, method = "kruskal.test")
write.xlsx(compare_means(mean_speed ~ probe,  data = data, method = "kruskal.test"),
           file = 'kruskal.test.mean_speed.xlsx')

compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

write.xlsx(compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test"),
           file = 'wilcox.test.mean_speed.xlsx')


# Compute mean and SD-----------------------------
df.summary.sinuosity <- group_by(data, probe) %>%
  summarise(
    sd = sd(sinuosity, na.rm = TRUE),
    sinuosity = mean(sinuosity)
  )

df.summary.sinuosity


# Perform pairwise comparisons -------------------------
compare_means(SDDC ~ probe,  data = data, method = "t.test")
compare_means(SDDC ~ probe,  data = data, method = "wilcox.test")

compare_means(sinuosity ~ probe,  data = data, method = "t.test")
compare_means(sinuosity ~ probe,  data = data, method = "wilcox.test")

write.csv(compare_means(SDDC ~ probe,  data = data, method = "t.test"), file="t.test.SDDC.csv")
write.csv(compare_means(SDDC ~ probe,  data = data, method = "wilcox.test"), file="wilcox.test.SDDC.csv")

write.csv(compare_means(sinuosity ~ probe,  data = data, method = "t.test"), file="t.test.sinuosity.csv")
write.csv(compare_means(sinuosity ~ probe,  data = data, method = "wilcox.test"), file="wilcox.test.sinuosity.csv")

compare_means(mean_speed ~ probe,  data = data, method = "wilcox.test") # pairwise comparisons

# Visualize: Specify the comparisons you want color---------------------
my_comparisons <- list( c("W1","W2"))

ggboxplot(data, x = "probe", y = "sinuosity")+ 
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 0)




#################################################
## Compute statistics and Plot results ##########
#################################################
# Summary -----------------------
summary(data[data$probe=='W1',])
summary(data[data$probe=='W2',])

# Kruskal test ---------
kruskal.test(data$mean_speed ~ data$probe)

shapiro.test(data[data$probe=='3',]$mean_speed) # ??????? ?????? ?????
shapiro.test(data[data$probe=='7',]$mean_speed) # 

# ????? ??????:
var.test(data$mean_speed ~ data$probe) # ????????? ?????? ????????? ?????????

# ???? ?????? ????????? ? ????????? ?? ??????????:
# t test
t.test(data$mean_speed ~ data$probe, var.equal = TRUE)

# ???? ??????? ?? ???????? - ???????? ???? ?????
wilcox.test(data[data$probe=='W1',]$mean_speed,
            data[data$probe=='W2',]$mean_speed,
            paired = FALSE) # ???? ?????




# function to produce summary statistics (mean and +/- sd), as required for ggplot2
data_summary <- function(x) {
  mu <- mean(x)
  sigma1 <- mu-sd(x)
  sigma2 <- mu+sd(x)
  return(c(y=mu,ymin=sigma1,ymax=sigma2))
}

# function for computing mean, DS, max and min values
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


wilcox.test(data[data$probe=='W1',]$mean_speed, 
            data[data$probe=='W2',]$mean_speed, 
            paired = FALSE)


###
glimpse(data)

shapiro.test(data[data$probe=='W1',]$mean_speed) # ??????? ?????? ?????
shapiro.test(data[data$probe=='W2',]$mean_speed) # 
# ????? ??????:
var.test(data[data$probe=='W1',]$mean_speed ~ 
         data[data$probe=='W2',]$mean_speed) # ????????? ?????? ????????? ?????????
# ???? ????????? ?? ??????????:
# t test
t.test(data[data$probe=='W1',]$mean_speed ~ 
       data[data$probe=='W2',]$mean_speed, var.equal = TRUE)
# ???? ? ?????? ??????????? ??????????,
# ?? ???????? ???????? ?????:
t.test(data[data$probe=='W1',]$mean_speed ~ 
       data[data$probe=='W2',]$mean_speed, var.equal = FALSE)


pairwise.t.test(data$mean_speed, data$probe, p.adj = "bonf")

pairwise.t.test(data$mean_speed, data$probe, p.adj = "holm")


compare_means(mean_speed ~ probe,
              data = data,
              method = "wilcox.test")

# plot 1 Mean +- SD
#ggsave("plot1.jpg", width = 12, height = 12, units = "cm")


ggplot(data, aes(x = probe, y = mean_speed)) +
  
  
  ylim(c(0, 160)) +
  
  theme_classic(base_size=14) +
  
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", 
               aes(width=0.7),
               position=position_dodge(width=1.5), size=0.8) +
  
  geom_jitter(position = position_jitter(.15),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'Micrometers per hour',
       x = "FetMSC 24h motility") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0,
                                   size = 12, face="bold",
                                   colour="black" ), 
        axis.text.y = element_text(angle = 0, hjust = 1,
                                   size = 11, face = 'bold',
                                   colour="black" )) +
  
  geom_signif(y_position = c(150),
              xmin = c(1),
              xmax = c(2),
              annotation = "p = 0.007", textsize = 5, family="serif",
              tip_length = 0.04, size = 1, face="bold")



#########

library(tidyverse)
library(gapminder)


data %>% 
  ggplot(aes(x=probe,y=mean_speed, fill=probe)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2)  +
  
  geom_signif(y_position = c(150),
              xmin = c(1),
              xmax = c(2),
              annotation = "p = 0.007", textsize = 4, family="serif",
              tip_length = 0.04, size = 1) +
  
  labs(y = 'Micrometers per hour',
       x = "FetMSC 24h motility")


ggplot(data, aes(x = probe, y = mean_speed))  + 
# Notched box plot with mean points
geom_boxplot(notch = TRUE, fill = "lightgray") +
  
stat_summary(fun.y = mean, geom = "point",
               shape = 18, size = 2.5, color = "#FC4E07") +
  
  geom_signif(y_position = c(150),
              xmin = c(1),
              xmax = c(2),
              annotation = "p = 0.007", textsize = 4, family="serif",
              tip_length = 0.04, size = 1) +
  
  labs(y = 'Micrometers per hour',
       x = "FetMSC 24h motility")





ggplot(data, aes(x = probe, y = mean_speed))  + 
  # Notched box plot with mean points
  geom_violin(trim = FALSE) + 
  stat_summary(
    fun.data = "mean_sdl",  fun.args = list(mult = 1), 
    geom = "pointrange", color = "black"
  ) +
  
  geom_signif(y_position = c(200),
              xmin = c(1),
              xmax = c(2),
              annotation = "p = 0.007", textsize = 4, family="serif",
              tip_length = 0.04, size = 1) +
  
  labs(y = 'Micrometers per hour',
       x = "FetMSC 24h motility")




############ C1-C6 column names
traj_analysis_W1 <- function(input) {
  
  data <- read.csv(input)
  
  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))
  
  for (i in unique(data$C1)) {
    
    # Define x, y, and time coordinates
    coords <- data.frame(x = data$C3[data$C1 == i],
                         y = data$C4[data$C1 == i],
                         # times = c(1:96))
                         timeCol = data$C2[data$C1 == i],
                         spatialUnits = "pixels", timeUnits = "hours")
    
    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 90/3600/3) #fps frames/3600/hours = 1/120
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 0.6579 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 0.4273504 / 1, "micrometer") # 234 px = 100 microns
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    # trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)
    
    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj), # D/L ratio
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           probe = 'W1')
    
    head(traj_params)
    
  }
  
  # print(traj_params)
  write.csv(traj_params, file = 'traj_W1.csv')
  tracks <<- traj_params
  return(traj_params) } # 90 frames

traj_analysis_W1J <- function(input) {
  
  data <- read.csv(input)
  
  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))
  
  for (i in unique(data$Track)) {
    
    # Define x, y, and time coordinates
    coords <- data.frame(x = data$X[data$Track == i],
                         y = data$Y[data$Track == i],
                         # times = c(1:96))
                         timeCol = data$Slice[data$Track == i],
                         spatialUnits = "pixels", timeUnits = "hours")
    
    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 90/3600/3) #fps frames/3600/hours = 1/120
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 0.6579 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 0.4273504 / 1, "micrometer") # 234 px = 100 microns
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    # trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)
    
    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj), # D/L ratio
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           #mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           mean_angle = 0,
                           #probe = 'W1'
                           probe = 1
    )
    
    head(traj_params)
    
  }
  
  # print(traj_params)
  write.csv(traj_params, file = 'traj_W1J.csv')
  tracks <<- traj_params
  return(traj_params) } # 90 frames


traj_analysis_W4J <- function(input) {

data <- read.csv(input)

traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                        c("track",
                          "length",
                          "distance",
                          "straight",
                          "square_displacement",
                          "mean_speed",
                          "sd_speed",
                          "max_speed",
                          "min_speed",
                          "sinuosity",
                          "emax",
                          "DC",
                          "SDDC",
                          "mean_angle",
                          "probe"))

for (i in unique(data$Track)) {
  
  # Define x, y, and time coordinates
  coords <- data.frame(x = data$X[data$Track == i],
                       y = data$Y[data$Track == i],
                       # times = c(1:96))
                       timeCol = data$Slice[data$Track == i],
                       spatialUnits = "pixels", timeUnits = "hours")
  
  trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 90/3600/3) #fps frames/3600/hours = 1/120
  TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
  # A 0.6579 object had length 1 pixels in the video, scale to micrometres
  trj <- TrajScale(trj, 0.4273504 / 1, "micrometer") # 234 px = 100 microns
  TrajGetUnits(trj) # Returns the spatial units of a trajectory
  TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
  TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
  # Rediscretization
  # The function TrajResampleTime linearly interpolates points along a trajectory
  # to create a new trajectory with fixed step time intervals.
  # trj <- TrajResampleTime(trj, 901)
  TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
  TrajGetFPS(trj)
  par(mar=c(5,5,5,5))
  # Plot it
  plot(trj, lwd = 2)
  points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
  # Trajectory analysis
  # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
  derivs <- TrajDerivatives(trj)
  
  traj_params <- add_row(traj_params,
                         track = i,
                         # total length of the trajectory
                         length = TrajLength(trj),
                         # straight-line distance from the start to the end of the trajectory
                         distance = TrajDistance(trj),
                         # Straightness index
                         straight = TrajStraightness(trj), # D/L ratio
                         # expected square displacement of a correlated random walk
                         square_displacement = TrajExpectedSquareDisplacement(trj),
                         # Measures of speed
                         mean_speed = mean(derivs$speed),
                         sd_speed = sd(derivs$speed),
                         max_speed = max(derivs$speed),
                         min_speed = min(derivs$speed),
                         # Measures of straightness
                         sinuosity = TrajSinuosity2(trj),
                         emax = TrajEmax(trj),
                         SDDC  =  sd(TrajDirectionalChange(trj)),
                         DC = mean(TrajDirectionalChange(trj)),
                         #mean_angle = TrajMeanVectorOfTurningAngles(trj),
                         mean_angle = 0,
                         #probe = 'W1'
                         probe = 5
  )
  
  head(traj_params)
  
}

# print(traj_params)
write.csv(traj_params, file = 'traj_W4J.csv')
tracks <<- traj_params
return(traj_params) } # 90 frames


# Start scan for W1J
#for (file_name in file_list_W1J) {
#  traj_analysis_W1J(file_name)
#  tracks_W1J <- rbind(tracks_W1J, tracks)
#}

tracks_W1J <- cbind(seq(20), matrix(0, 20, 13), rep(1, 20))
tracks_W1J <- as.data.frame(tracks_W1J)
colnames(tracks_W1J) <- mynames

# Start scan for W4J
#for (file_name in file_list_W4J) {
#  traj_analysis_W4J(file_name)
#  tracks_W4J <- rbind(tracks_W4J, tracks)
#}

tracks_W4J <- cbind(seq(20), matrix(0, 20, 13), rep(1, 20))
tracks_W4J <- as.data.frame(tracks_W4J)
colnames(tracks_W4J) <- mynames


rm(data)
rm(all.h)
rm(alltracks)
rm(tracks)
rm(tracks_W1J)
rm(tracks_W2H)
rm(tracks_W2J)
rm(tracks_W3H)
rm(tracks_W5H)
rm(tracks_W5J)
rm(tracks_W6H)

