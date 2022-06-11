# Queries for AWAPDaily, M. Kearney Feb 2020
#
# Tables: 1900-2020
#
# [id]
# [day]
# [rr] rainfall (mm), 1900-
# [tmax] maximum air temperature, 1911-
# [tmin] minimum air temperature, 1911-
# [vpr] vapour pressure hPa, 1971-
# [sol] global daily solar radiation MJ m-2 d-1, 1990-
#
# Table: latlon
#
# [id]
# [latitude]
# [longitude]


######################################################
######################################################
######################################################

# begin query for a specific site for a given period,
# site specified with Google Maps lookup
#### Load libraries and functions

# in this file is the function for cloud cover, run this line
source("2_Functions_cloud_2022.R")
install.packages("devtools")
install.packages("ncdf4")
devtools::install_github("mrke/NicheMapR")

library(NicheMapR)
library(dplyr)
library(chron)
library(RMySQL)

# Run this so that the cloud cover function from NicheMapR can work
get.global.climate("../Data/globalclimate") #downloads a data base and calls it global climate

memory.limit(size = 4000)

####### BEGIN USER INPUT #######################################
# Choose a csv file containing the species/individual, ReferenceNumber, latitude, longitude and months of activity
datab <- read.csv("../Data/5_Locations.csv")
head(datab)
sitesfile <- datab[,c("Spp", "Reg", 
                      "Latitude", "Longitude", 
                      "MonthsActivityALA", "MonthCollectionLabel")]
colnames(sitesfile)<-c("species", "reg", "lat", "lon", "months","months_Coll")

# Add the 0 before the months with only one digit:
sitesfile$months_Coll <- padz(sitesfile$months_Coll) # needs package chron

# The output directory is the same as the input file. This program will
# additionaly create a folder called output where results will go
output_dir <- dirname("../Output")

# Make sure this date starts at the start of the year (1st of January)
datestart_string <- "1/1/2010" # day, month, year

# Make sure this date ends at the end of the year (31st December)
datefinish_string <- "31/12/2020" # day, month, year
####### END USER INPUT #########################################

# close all database connections before opening a new one
killDbConnections()

sitesfile_example <- as.data.frame(sitesfile)[5:10, ] # example of some rows

source("./3_Get_data_2022.R")

results <- get_data(sitesfile, datestart_string, datefinish_string)
results$total_dataset
write.csv(results$total_dataset,paste(output_dir, "/2022ResultsEcolVar2.csv", sep = ""))
write.csv(results$check_missing,paste(output_dir, "/2022ResultsEcolVar_with_missing2.csv", sep = ""))
write.csv(results$output,paste(output_dir, "/results_everything.csv", sep = ""))
# close all database connections
killDbConnections()

plot_data <- function(output) {
  output$sol <- as.numeric(as.character(output$sol))
  # plot data
  dates <- seq.Date(as.Date(datestart), as.Date(datefinish), "days")

  par(mfrow = c(2, 3))
  output <- cbind(dates, output)
  plot(output$dates,
    output$rr,
    type = "h", xlab = "month", ylab = "rainfall (mm/d)"
  )
  plot(output$dates,
    output$tmax,
    type = "l",
    xlab = "month",
    ylab = expression("max air temp (" * degree * "C)")
  )
  plot(output$dates,
    output$tmin,
    type = "l",
    xlab = "month",
    ylab = expression("min air temp (" * degree * "C)")
  )
  plot(output$dates,
    output$tmax,
    ylim = c(min(output$tmin), max(output$tmax)),
    type = "l",
    xlab = "month",
    col = "2",
    ylab = expression("max air temp (" * degree * "C)")
  )
  points(output$dates,
    output$tmin,
    type = "l",
    xlab = "month",
    col = "4",
    ylab = expression("min air temp (" * degree * "C)")
  )
  plot(output$dates,
    output$vpr,
    type = "l", xlab = "month", ylab = "vapour pressure (hPa)"
  )
  plot(output$dates,
    output$sol,
    type = "l",
    xlab = "month",
    ylab = "solar radiation (MJ/m^2/d)"
  )
}