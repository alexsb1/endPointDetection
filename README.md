# Laser ablation mass spectrometry automatic end point blast through detection

## Description
This function imports a single dataframe containing one time resolved analysis laser ablation mass spectrometry analysis of a foraminifera (or other carbonate shell), then detects when the laser has burnt through the foraminifera test as a function of change in 44Ca over time.

Calculate first derivative (rate of change) to identify when the laser blasts through the chamber.
Calculates the change of signal (counts per second) per change of time.
A zero is used as a NA placeholder because there are n-1 difference values than n observations.

<sup>44</sup>Ca is used by default to identify the time taken to blast through the chamber wall
smooth the signal before detecting the rate of change.
A large value of order causes more smoothing.

## Usage
`endPoint(df, dt = 10, smoothing = 5, timeCol = "Time", signalCol = "Ca44", profile = "TRUE",  timeUnits = "seconds")`

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

**timeCol** is the column title in your dataframe containing the time stamp in your TRA.
Default = "Time"

**signalCol** is the column title in your dataframe containing the numerical data you want to sue to detect the endpoint.
Mostly this is "44Ca" containing a calcium isotope counts per second.
This could be any column of numerical data that you want to detect the endpoint, not necessarily <sup>44</sup>Ca or any particular elemental isotope.
Default = "Ca44"

**profile** is a visualisation of the endpoint detection mechanism in ggplot. This argument is logical.
Default = "TRUE"

**timeUnits** is the units that your time resolved analysis is measured in. This is the units of the timeCol.
This argument is a string and is only necessary if the argument `profile = "TRUE"`.
Default = "seconds"

## Details
This function was designed for detecting when a laser had ablated through a foraminifera test and keeping only the outputted rows of a time resolved analysis that were relevant by returning a start time and end time of desired data.
The default values are based on outputs from an Agilent 8900 ICP-MS Triple Quad and a New Wave Research NWR193 laser ablation unit at the University of Southampton, Waterfront campus, National Oceanography Centre, UK.

## Notes
The output of this function are returned in a new dataframe called dfReturn.

* `dfReturn$df` contains your dataframe with only the rows that are between your startTime and endTime.

* `dfReturn$startTime` contains the earliest time in your TRA as a numerical value.

* `dfReturn$endTime` contains the last timestep before the laser ablated through the carbonate shell in your TRA as a numerical value.

* `dfReturn$profile` contains visualisation of your TRA identifying where the laser ablated through the carbonate shell as a ggplot object. This is only available if a profile was generated using `profile = "TRUE"`.


Take note that parsing `profile = "TRUE"` substantially increases the time taken for this function to run.

The user can add a title to this plot after running this function by adding the following lines of code
`dfReturn$profile + labs(title = "Title text")`

## Examples

Example data is located at https://github.com/alexsb1/endpointDetection/

`data1 <- read.csv("exampleData/foram-72-shot-3.csv", header = TRUE)`

* The simplest usage, using all the default values. \
`endPoint(data1)`

* Specifying custom dt and smoothing arguments. \
`endPoint(data1, dt = 15, smoothing = 10)`

* Specifying a different dataframe column to use for detection by passing the "Mg24" dataframe heading as substitute for Ca44 and using default dt and smoothing. \n
`endPoint(data1, signalCol = "Mg24")`

* Specifying to not make a profile plot. This can save a substantial amount of time. \
_Note that if you use a time unit other than seconds you must specify these here._ \
`endPoint(data, profile = "FALSE", timeUnits = "seconds")`



## Reference

You should cite this function as

`Searle-Barnes C A, 2021, https://github.com/alexsb1/endPointDetection`

This is released under GNU General Public License v3.0




