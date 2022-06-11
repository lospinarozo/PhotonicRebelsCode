new_db_connection <- function() {
  # open a new database connection to AWAP
  channel <- dbConnect(MySQL(),
    user = "general",
    password = "predecol",
    host = "115.146.93.180", # Updated server 2021
    dbname = "AWAPDaily",
    port = 3306
  )
  return(channel)
}

leapfix <- function(indata, yearlist, mult = 1) {
  leapyears <- seq(1900, 2060, 4)
  for (j in 1:length(yearlist)) {
    if (yearlist[j] %in% leapyears) {
      if (mult == 1) {
        data <- c(indata[1:59], indata[59], indata[60:365])
      } else {
        data <- c(indata[1:(59 * mult)], indata[(58 *
          mult + 1):(59 * mult)], indata[(59 * mult +
          1):(365 * mult)])
      }
    } else {
      data <- indata
    }
    if (j == 1) {
      alldata <- data
    } else {
      alldata <- c(alldata, data)
    }
  }
  return(alldata)
}

padz <- function(x, width = max(nchar(x)), fill = "0") {
  gsub(" ", fill, formatC(x, width = width))
}



get_cloud <- function(loc, yearlist, sol) {
  micro_clearsky <- micro_global(loc = loc, clearsky = 1, timeinterval = 365, solonly = 1) # compute clear sky for 365 days for a year
  clearskyrad <- micro_clearsky$metout[, c(1, 13)] # get the solar radiation column
  clearskysum <- aggregate(clearskyrad[, 2], by = list(clearskyrad[, 1]), FUN = sum)[, 2] # sum over all hours per day

  allclearsky <- leapfix(clearskysum, yearlist) # expand to a sequence of daily values across the time period, inserting leap years where necessary
  allclearsky <- (allclearsky * 3600) / 1e+06 # convert from J/h/m2 to MJ/day

  cloud <- (1 - sol / allclearsky) * 100

  mean_cloud <- mean(cloud)

  return(mean_cloud)
}



killDbConnections <- function() {
  all_cons <- dbListConnections(MySQL())

  print(all_cons)

  for (con in all_cons) {
    +dbDisconnect(con)
  }

  print(paste(length(all_cons), " connections killed."))
}