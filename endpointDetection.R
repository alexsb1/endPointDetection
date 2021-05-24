# Author: Alex Searle-Barnes
# email: c.j.a.searle-barnes@soton.ac.uk
# GitHub: https://github.com/alexsb1
# date: 19 May 2021

# Description: This function imports a single dataframe containing one time resolved analysis
#              laser ablation mass spectrometry analysis of a foraminifera (or other carbonate shell),
#              then detects when the laser has burnt through the foraminifera test as a function of
#              change in 44Ca over time.

# Prerequisites: You must have already cleaned your data into a single data before using this function.
#                This function will only return one endpoint per dataframe, so use with sapply or within
#                a loop if you have multiple endpoints to detect.

# Libraries
library("tidyverse")
library("smooth")

# Load example data
# data1 <- read.csv("exampleData/foram-72-shot-3.csv", header = TRUE)
# data2 <- read.csv("exampleData/foram-166-shot-7.csv", header = TRUE)
# data3 <- read.csv("exampleData/foram-174-shot-4.csv", header = TRUE)
# data4 <- read.csv("exampleData/coral1.csv", header = TRUE)
# data5 <- read.csv("exampleData/coral2.csv", header = TRUE)
# data6 <- read.csv("exampleData/coral3.csv", header = TRUE)
# data7 <- read.csv("exampleData/coral4.csv", header = TRUE)
# data8 <- read.csv("exampleData/coral5.csv", header = TRUE)
# data9 <- read.csv("exampleData/coral6.csv", header = TRUE)

# Function to detect endpoint

endPoint <- function(df, dt = 10, smoothing = 5, timeCol = "Time", signalCol = "Ca44", profile = "TRUE",  timeUnits = "seconds"){
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
  
  scanRate <- 1 / (df$Time[2] - df$Time[1])
  
  # timeCol is the column title in your dataframe containing the time stamp in your TRA.
  # Default = "Time"
  
  # signalCol is the column title in your dataframe containing the data 44Ca.
  # This could be any column of numerical data that you want to detect the endpoint,
  # not necessarily 44Ca or any particular elemental isotope.
  # Default = "Ca44"
  
  # Calculate first derivative (rate of change) to identify when the laser blasts through the chamber
  # Calculates the change of signal (cps) per change of time
  # Assigns a zero as first value and then the rate of change for each corresponding time interval.
  # A zero is used as a NA placeholder because there are n-1 difference values than n observations.
  
  # 44Ca is used to identify the time taken to blast through the chamber wall
  # smooth the signal before detecting the rate of change.
  # A large value of order causes more smoothing
  df$Ca44sma <- sma(df[, grep(signalCol, names(df))[1]], order = smoothing) %>% .$fitted %>% as.numeric()
  
  
  # This gives change as absolute
  Ca44dydt <- diff(df$Ca44sma, lag = dt) / diff(df[, grep(timeCol, names(df))[1]], lag = dt)
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
  tailSeconds <- round(tailSeconds, digits = 3) #not working with pipes
  
  #Identify the time ($Time) in ms when the maximum rate of change in 44Ca occurs.
  # The subtraction of dt/x is in place to backtrack the zero dydt values at the start. This results in the endTime being before the sharp signal cutoff.
  # Notice that dt is in observations, but endTime in in seconds.
  endTime <- subset(df[, grep(timeCol, names(df))[1]], df$Ca44dydt == min(df$Ca44dydt)) - tailSeconds
  
  # A consequence of legacy code
  startTime <- min(df[, grep(timeCol, names(df))[1]])
  maxTime <- max(df[, grep(timeCol, names(df))[1]])
  
  # Make a dataframe to append the outputs from this function to.
  # This dataframe will be returned outside this function
  dfReturn <- NULL
  
  dfReturn$df <- df # Your dataframe with only the rows that are between your startTime and endTime
  dfReturn$exceededThresholdSubset <- NULL #not currently used
  dfReturn$startTime <- startTime # The first time step in your analysis after this function
  dfReturn$endTime <- endTime # The last time step in your analysis after this function
  
  if(profile == "TRUE"){
    dfReturn$profile <- ggplot(df, aes(x=df[, grep(timeCol, names(df))[1]])) +
      annotate("rect", xmin = startTime - 10, xmax = startTime, ymin = -Inf, ymax = Inf, fill ="red", alpha = 0.5)+ #change -10 to a percentage
      annotate("rect", xmin = endTime, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.5)+
      geom_point(aes(y=Ca44Scaled, colour = "Signal")) +
      geom_line(aes(y=Ca44Scaled, colour = "Signal")) +
      geom_line(aes(y=Ca44dydt, colour = "dydt")) +
      geom_vline(xintercept = endTime, colour = "purple")+
      geom_label(x = startTime - 10 , y = median(abs(scale(df[, grep(signalCol, names(df))[1]], center = TRUE))), label = paste("TRA started at", startTime, timeUnits), size = 2, hjust = "left")+
      geom_label(x = endTime, y = mean(abs(scale(df[, grep(signalCol, names(df))[1]], center = TRUE))), label = paste("endTime \n (largest signal change - (dt/scanRate)) \n at", endTime, timeUnits), size = 2)+
      labs(y = paste("Scaled", signalCol, "signal and rate of change"),
           x = paste("Time elapsed in", timeUnits),
           subtitle = paste("With a smoothing of", smoothing, "observations, a dt of", dt, "observations and (dt/scanRate) of", tailSeconds, timeUnits)
      )+
      scale_x_continuous(limits = c(startTime-10, maxTime),
                         breaks = scales::pretty_breaks(n = 10))+
      scale_colour_manual("", 
                          breaks = c("Signal", "dydt"),
                          values = c("black", "blue")) +
      theme_bw()
    
    # Data points located within the red shading are removed in the returned dataframe as these are beyond 
    
    # The user can add a title to this plot after running this function by adding the following lines of code
    # dfReturn$profile + labs(title = "Title text")
    
  }
  
  return(dfReturn)
  
}
