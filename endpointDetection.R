# Author: Alex Searle-Barnes
# email: c.j.a.searle-barnes@soton.ac.uk
# GitHub: https://github.com/alexsb1
# date: 19 May 2021

# Description: This function imports a single dataframe containing one time resolved analysis
#              laser ablation mass spectrometry analysis of a foraminifera (or other carbonate shell),
#              then detects when the laser has burnt through the foraminifera test as a function of
#              change in 44Ca over time.

# Prerequisites: You must have already cleaned your data into a single data before using this function.
#                Your data should already be background corrected
#                This function will only return one endpoint per dataframe, so use with sapply or within
#                a loop if you have multiple endpoints to detect.

# Libraries
library("tidyverse")
library("smooth")

# Function to detect endpoint

endPoint <- function(df, dt = 10, smoothing = 5, scanRate = 6.4 , timeCol = "Time", Ca44 = "Ca44"){
  # Understanding inputs
  
  # df is your cleaned dataframe containing a single analysis
  
  # dt controls the change in time elapsed between data points (delta time).
  # An infinitesimal value does not produce a sharp dx/dt graph peak (max value)
  # for the diff() function value is used in the lag argument.
  # Use a lower value for a fast blast through of chamber wall.
  # Default = 10 
  
  # smoothing controls how much smoothing to apply over the dt period.
  # The bigger the number the more smoothing.
  # Default = 5
  
  # scanRate is the frequency of signal detection. The number of rows of data per second in your dataframe.
  # Default = 6.4
  
  # timeCol is the column title in your dataframe containing the time stamp in your TRA.
  # Default = "Time"
  
  # Ca44 is the column title in your dataframe containing the data 44Ca.
  # This could be any column of numerical data that you want to detect the endpoint,
  # not necessarily 44Ca or any particular elemental isotope.
  # Default = "Ca44"
  
  # Calculate first derivative (rate of change) to identify when the laser blasts through the chamber
  # Calculates the change of signal (cps) per change of time
  # Takes the absolute (modulus) value of this rate of change
  # Assigns a zero as first value and then the rate of change for each corresponding time interval.
  # A zero is used as a NA placeholder because there are n-1 difference values than n observations.
  
  # 44Ca is used to identify the time taken to blast through the chamber wall
  # smooth the signal before detecting the rate of change.
  # A large value of order causes more smoothing
  df$Ca44sma <- sma(df$Ca44, order = smoothing) %>% .$fitted %>% as.numeric()
  
  # This gives change as absolute
  Ca44dydt <- diff(df$Ca44sma, lag = dt) / diff(df$Time, lag = dt)
  Ca44dydt <- c(rep(0, times = dt), Ca44dydt)
  Ca44dydt <- (Ca44dydt - min(Ca44dydt)) / (max(Ca44dydt) - min(Ca44dydt)) #scaled to 0-1
  df$Ca44dydt <- Ca44dydt
  
  #The scaled values (between 0 - 1) of 44Ca are calculated here for plotting on a graph
  #for rate of change against time ($Time) comparison
  df$Ca44Scaled <- (df$Ca44sma - min(df$Ca44sma))/(max(df$Ca44sma) - min(df$Ca44sma))
 
  # A subset dataframe made of events when the dydt threshold is exceeded.
  # This could be interesting but not currently used.
#  exceededThresholdSubset <- subset(df, df$Ca44dydt >= dydtThreshold)
  
  # the number of seconds to remove after end point detection (endTime)
  tailSeconds <- dt/scanRate
  
  #Identify the time ($Time) in ms when the maximum rate of change in 44Ca occurs.
  # The subtraction of dt/x is in place to backtrack the zero dydt values at the start. This results in the endTime being before the sharp signal cutoff.
  # Notice that dt is in observations, but endTime in in seconds.
  endTime <- subset(df$Time, df$Ca44dydt == min(df$Ca44dydt)) - tailSeconds

  # A consequence of legacy code
  startTime <- min(df$Time)
  
  # Make a dataframe to append the outputs from this function to.
  # This dataframe will be returned outside this function
  dfReturn <- NULL
  
  dfReturn$df <- df # Your dataframe with only the rows that are between your startTime and endTime
  dfReturn$exceededThresholdSubset <- NULL #not currently used
  dfReturn$startTime <- startTime # The first time step in your analysis after this function
  dfReturn$endTime <- endTime # The last time step in your analysis after this function
  
  return(dfReturn)
  
}
