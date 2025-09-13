# CPFU

This repository contains functions `CPF()`, `CBPF()`, `CBPF_intervals()` and `CBPF_normalized()` which create plots of conditional probability function (CPF) for air pollution studies. 

Conditional probability function is commonly used for directional analysis of air pollution data, which means that it provides connection between pollutant data and wind direction and speed. CPF is defined as:
```math
CPF = \frac{m_{\Delta\Theta}}{n_{\Delta\Theta}}
```
where $m_{\Delta\Theta}$ is the number of events from wind sector $\Delta\Theta$ when the pollutant exceeds previously set threshold, and $n_{\Delta\Theta}$ is the total number of events from wind sector $\Delta\Theta$. The threshold is usualy some given percentile of concentrations ($\mathrm{75^{th}}$ percentile is most often the choice). 

In the functions below, wind can be segmented into a given number of sectors according to direction and speed (in case of the bivariate version of the plot). CPF value is then calculated for each wind sector. 

The basic input for the functions is a `data.frame` with a column for the pollutant concentrations (the name for this column nneeds to be given), column for the wind speed, and column for the wind direction data. Other input is not mandatory. 

The user can choose the number of wind direction segments and (in case of bivariate plots, with functions `CBPF()`, `CBPF_intervals()`, `CBPF_normalized()`) the size of wind speed intervals. 

During the calm wind periods pollutant concentrations can be raised due to the local sources, which are not subjects of directional analysis (the sources can not be determined using the transport by wind), so the user often wants to omit events with the wind speed below some value. This wind speed limit is often 1 m/s. 

In order to perform the directional analysis, i. e. to estimate which wind directions contribute in air pollution, it is needed to define the CPF threshold value, which is the value that CPF would get if the given wind sector did not have any effect on pollutant concentrations. This value is defined as:
```math
CPF_t = \frac{100-perc}{100}
```
The chosen percentile is most often $\mathrm{75^{th}}$ so the CPF threshold is most often 0.25. In that case, every CPF value above 0.25 means that the wind from the given sector increases the pollution (relatively, compared to other wind sectors), while CPF value below 0.25 means that the wind from this sector decreases the pollution. 

But, an important question is whether the CPF value is ***significantly*** above/below the CPF threshold. To determine this, it is needed to estimate the uncertainties of CPF values. The functions contained in this repository calculate uncertainties as the confidence intervals using two methods:
- binomial ratio: confidence intervals are calculated using the function `BinomRatioCI()` from `DescTools` package, with the Koopman asymptotic score
- bootstraping: confidence intervals are calculated using the `boot` package, using the bootstrap percentile method

If the uncertainties are big for many sectors, it is possible that the number of events in many segments is very small, so using less segments (i. e. using bigger segments) could help in obtaining better results. 

## Testing data

File `mydata.csv` contains a `data.frame` with an example data that can be used for testing the functions. The file contains columns `ws`, `wd` and `PM` that represent wind speed (in knots), wind direction (in degrees) and particulate matter concentrations (in $\upmu \mathrm g\ \mathrm m^{-3}$) collected at the Adams, Colorado station during 2018, downloaded from the United States Environmental Protection Agency database (https://aqs.epa.gov/aqsweb/airdata/download_files.html#Raw). 




## Usage

The following packages should be loaded first: `boot`, `DescTools`, `ggplot2`.

The functions can be used simply by downloading the files `CPF.R`, `CBPF.R`, `CBPF_intervals.R`, and `CBPF_normalized.R` in the working directory.  After that, it is only needed to call the functions like this:
```R
source("CPF.R")
source("CBPF.R")
source("CBPF_intervals.R")
source("CBPF_normalized.R")
```
The functions are now in your workspace and you can use them with the following arguments:

```R
 CPF(
   data,
   sectors = 16,
   ws = "ws",
   wd = "wd",
   pollutant = NULL,
   percentile = 75,
   conf = 0.68,
   intervals = "br",
   R = 1000,
   remove_speed = 1,
   plot = TRUE
 )

 CBPF(
   data,
   sectors = 16,
   ws = "ws",
   wd = "wd",
   pollutant = NULL,
   percentile = 75,
   remove_speed = 1,
   speed_interval = 1,
   speed_unit = "m/s",
   plot = TRUE
 )
 
 CBPF_intervals(
   data,
   sectors = 16,
   ws = "ws",
   wd = "wd",
   pollutant = NULL,
   percentile = 75,
   conf = 0.68,
   intervals = "br",
   R = 1000,
   remove_speed = 1,
   speed_interval = 1,
   speed_unit = "m/s",
   plot = TRUE
 )

 CBPF_normalized(
   data,
   sectors = 16,
   ws = "ws",
   wd = "wd",
   pollutant = NULL,
   percentile = 75,
   conf = 0.68,
   intervals = "br",
   R = 1000,
   remove_speed = 1,
   speed_interval = 1,
   speed_unit = "m/s",
   plot = TRUE
 )
```

## Arguments

<dl>
  <dt>data</dt>
  <dd>A `data.frame` including wind direction, wind speed and pollutant concentration columns; mandatory</dd>
  <dt>sectors</dt>
  <dd>Number of wind direction sectors; default is 16</dd>
  <dt>ws</dt>
  <dd>Name for wind speed column in data; default is "ws"</dd>
  <dt>wd</dt>
  <dd>Name for wind direction column in data; default is "wd"</dd>
  <dt>pollutant</dt>
  <dd>Name for pollutant concentrations in data; mandatory</dd>
  <dt>percentile</dt>
  <dd>Percentile for threshold; default is 75Name for pollutant concentrations in data; mandatory</dd>
  <dt>remove_speed</dt>
  <dd>Remove speed below this value; default is 1</dd>
  <dt>speed_interval</dt>
  <dd>Interval for creating wind speed segments; default is 1</dd>
  <dt>speed_unit</dt>
  <dd>Measuring unit for wind speed; default is "m/s"</dd>
  <dt>conf</dt>
  <dd>Confidence level; default is 0.68</dd>
  <dt>intervals</dt>
  <dd>Method for interval calculation, c("br", "bootstrap"), "br" meaning binomial ratio; default is "br"</dd>
  <dt>R</dt>
  <dd>Number of bootstrap replicates; default is 1000</dd>
  <dt>plot</dt>
  <dd>If TRUE, a plot is created; if FALSE, a table is produced containing wind directions (wd), number of total observations (n), number of observations above threshold (m), CPF values, lower and upper limits for each wind segment; default is TRUENumber of bootstrap replicates; default is 1000</dd>
</dl>



## Examples

File `mydata.csv` from this repository can be used to test the functions. In this file, default names for wind speed and direction are already used, so the only input the user has to give is the pollutant name. Here, events with wind speed below 2 knots are also removed:
```R
# Test data:
mydata <- read.csv("mydata.csv")
CPF(mydata, pollutant = "PM", remove_speed = 2)
```
![Example of a CPF graph plotted with `CPF` function(https://github.com/mcargonja/CPFU/blob/db6cb2737a368e6f65470173fb2c758439ddab1d/example_CPF.jpg)
```R
# Using the bootstrap method for uncertainty intervals:
CPF(mydata, pollutant = "PM", remove_speed = 2, intervals = "bootstrap")

# Wind segmentation with different number of sectors:
CPF(mydata, pollutant = "PM", sectors = 32, remove_speed = 2)

# Creating a data.frame instead of plotting a graph:
CPF(mydata, pollutant = "PM", remove_speed = 2, plot = FALSE)

# Setting the confidence level at 95%:
CPF(mydata, pollutant = "PM", remove_speed = 2, conf = 95)
```
Bivariate CPF plots can be created using 3 different functions. The absic plot, without uncertainties, is created with `CBPF()`:
```R
# Basic usage of CBPF function:
CBPF(mydata, pollutant = "PM", remove_speed = 2, speed_interval = 4, speed_unit = "knots")
```
![Example of a bivariate CPF graph without uncertainties, plotted with `CBPF` function](https://github.com/mcargonja/CPFU/blob/c49077a3ccc9a772645b6b522a8b51a418712a98/example_CBPF.jpg)
The uncertainties are included in the graph created with `CBPF_intervals()` function:
```R
CBPF_intervals(mydata, pollutant = "PM", remove_speed = 2, speed_interval = 4, speed_unit = "knots")
```
![Example of a bivariate CPF graph with uncertainties, plotted with `CBPF\_intervals` function](https://github.com/mcargonja/CPFU/blob/b00790a591839d6ceebe813eec73cd6e85548d70/example_CBPF_intervals.jpg)
```R
# Confidence intervals from bootstrapping:
CBPF_intervals(mydata, pollutant = "PM", remove_speed = 2, speed_interval = 4, speed_unit = "knots", intervals = "bootstrap")

# Changing the number of bootstrap replicates:
CBPF_intervals(mydata, pollutant = "PM", remove_speed = 2, speed_interval = 4, speed_unit = "knots", intervals = "bootstrap", R = 2000)
```
```R
CBPF_normalized(mydata, pollutant = "PM", remove_speed = 2, speed_interval = 4, speed_unit = "knots")
```
![Example of a bivariate CPF graph normalized to uncertainties, plotted with `CBPF\_normalized` function](https://github.com/mcargonja/CPFU/blob/d632c282719f9dc0c83b3ac4190b7d0ea747c5d1/example_CBPF_normalized.jpg)
