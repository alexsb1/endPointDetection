# Laser ablation mass spectrometry automatic end point blast through detection

## Description
This function imports a single dataframe containing one time resolved analysis laser ablation mass spectrometry analysis of a foraminifera (or other carbonate shell), then detects when the laser has burnt through the foraminifera test as a function of change in 44Ca over time.

Calculate first derivative (rate of change) to identify when the laser blasts through the chamber
Calculates the change of signal (cps) per change of time
Takes the absolute (modulus) value of this rate of change
Assigns a zero as first value and then the rate of change for each corresponding time interval.
A zero is used as a NA placeholder because there are n-1 difference values than n observations.

44Ca is used to identify the time taken to blast through the chamber wall
smooth the signal before detecting the rate of change.
A large value of order causes more smoothing.

## Usage
`endPoint(df, dt = 10, smoothing = 5, scanrate = 6.4, timeCol = "Time", ca44 = "Ca44")`

## Arguments
**df** is your cleaned dataframe containing a single analysis
  
**dt** controls the change in time elapsed between data points (delta time).
An infinitesimal value does not produce a sharp dx/dt graph peak (max value)
for the diff() function value is used in the lag argument.
Use a lower value for a fast blast through of chamber wall.
Default = 10 

**smoothing** controls how much smoothing to apply over the dt period.
The bigger the number the more smoothing.
Default = 5

**scanRate** is the frequency of signal detection. The number of rows of data per second in your dataframe.
Default = 6.4

**timeCol** is the column title in your dataframe containing the time stamp in your TRA.
Default = "Time"

**Ca44** is the column title in your dataframe containing the data 44Ca.
This could be any column of numerical data that you want to detect the endpoint,
not necessarily 44Ca or any particular elemental isotope.
Default = "Ca44"

## Details
This function was designed for detecting when a laser had ablated through a foraminifera test and keeping only the outputted rows of a time resolved analysis that were relevant by returning a start time and end time of desired data.
The default values are based on outputs from an Agilent 8900 ICP-MS Triple Quad and a New Wave Research NWR193 laser ablation unit at the University of Southampton, Waterfront campus, National Oceanography Centre.

## Notes

## Examples


## Reference






