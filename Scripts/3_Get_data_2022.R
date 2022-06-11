### LOOP THROUGH SITES, AND FOR EACH SITE, THROUGH YEARS
### THIS FUNCTION NOW INCLUDES CLOUD COVER AND TAKING A DIFFERENT MONTH PER LINE
# loop through sites, make queries and collate data
get_data <- function(source_dataset, datestart_string, datefinish_string) {
    # process period of interest
    datestart <- strptime(datestart_string, "%d/%m/%Y") # convert to date format
    yearstart <- as.numeric(format(datestart, "%Y")) # get year start
    datefinish <- strptime(datefinish_string, "%d/%m/%Y") # convert to date format
    yearfinish <- as.numeric(format(datefinish, "%Y")) # yet year finish
    number_of_years <- yearfinish - yearstart
    years <- seq(yearstart, yearfinish, 1) # get sequence of years to do
    dates <- seq(datestart, datefinish + 3600, "DSTday") # sequence of dates
    juldaystart <- datestart$yday + 1 # get Julian day of year at start
    juldayfinish <- datefinish$yday + 1 # get Julian day of year at finish

    source_dataset$months <- (as.character(source_dataset$months))
    source_dataset$months_Coll <- as.character(source_dataset$months_Coll)
    
    source_dataset$lat <- as.numeric(source_dataset$lat)
    source_dataset$lon <- as.numeric(source_dataset$lon)

    # MAKE SURE DATA SETS ARE CLEAN EACH RUN#
    # This creates a template for one of the output tables.
    # This is so that data can be immediately added to this variable.
    # NOW INCLUDES CLOUD COVER
    # I included two sets because in case I need to compare the values for collection month vs. activity period
    
    
    individual_data <- data.frame(
        species = character(),
        reg = character(),
        lat = numeric(),
        lon = numeric(),
        
        # For activity period
        avg_temp_over_35 = numeric(),
        avg_max_temp = numeric(),
        avg_min_temp = numeric(),
        avg_sol = numeric(),
        avg_year_sol = numeric(), # Only available yearly
        avg_year_vpr = numeric(), # yearly, but will not use this one
        cloud_cover = numeric(),
        avg_rr = numeric(),
        avg_vpr = numeric(), # vapor pressure month
        
   
        # For collection month
        avg_temp_over_35_Coll = numeric(),
        avg_max_temp_Coll = numeric(),
        avg_min_temp_Coll = numeric(),
        avg_sol_Coll = numeric(),
        avg_rr_Coll = numeric(),
        avg_vpr_Coll = numeric() # vapor pressure month
        
        
    )

    # Another template for an output file.
    broken_records <- data.frame(
        species = character(),
        reg = character(),
        lat = numeric(),
        lon = numeric()
    )

    for (j in seq_len(length(source_dataset[, 1]))) {
        print(paste("Processing row ", j, " of ", length(source_dataset[, 1]), "...", sep = ""))

        # create a new database connection
        channel <- new_db_connection()

        # start loop through sites
        fail_flag <- 0

        species <- source_dataset[j, 1]
        reg <- source_dataset[j, 2]

        lat <- source_dataset[j, 3]
        lon <- source_dataset[j, 4]
        lat1 <- lat - 0.05
        lat2 <- lat + 0.05
        lon1 <- lon - 0.05
        lon2 <- lon + 0.05
        loc <- c(lon, lat) # Get data ready for cloud cover calculation

        months <- unlist(strsplit(as.character(source_dataset[j, "months"]),",")) # replaces as.character with strsplit to work with more than one month (activity period)
        months_Coll <- source_dataset[j, "months_Coll"]
        # Remember that: 
        # Month can not be written as 1 digit in the data frame we input, it has to be jan = 01
        # When working with more than one month the vector for each species should be 01,02,03 no space after the comma
       
        for (i in seq_len(length(years))) { # start loop through years
            # syntax for query
            if (length(years) == 1) { # doing a period within a year
                query <- paste("SELECT a.latitude, a.longitude, b.*
        FROM AWAPDaily.latlon as a
      , AWAPDaily.", years[i], " as b
      where (a.id = b.id) and (a.latitude between ",
                    lat1, " and ", lat2,
                    ") and (a.longitude between ", lon1, " and ", lon2,
                    ") and (b.day between ", juldaystart, " and ",
                    juldayfinish, ")
      order by b.day",
                    sep = ""
                )
            } else if (i == 1) { # doing first year, start at day requested
                query <- paste("SELECT a.latitude, a.longitude, b.*
          FROM AWAPDaily.latlon as a
          , AWAPDaily.", years[i], " as b
          where (a.id = b.id) and (a.latitude between ",
                    lat1, " and ", lat2, ") and (a.longitude between ",
                    lon1, " and ", lon2, ") and (b.day >= ", juldaystart, ")
          order by b.day",
                    sep = ""
                )
            } else if (i == length(years)) {
                # doing last year, only go up to last day requested
                query <- paste("SELECT a.latitude, a.longitude, b.*
            FROM AWAPDaily.latlon as a
            , AWAPDaily.", years[i], " as b
            where (a.id = b.id) and (a.latitude between ",
                    lat1, " and ", lat2, ") and (a.longitude between ",
                    lon1, " and ", lon2, ") and (b.day <= ", juldayfinish, ")
            order by b.day",
                    sep = ""
                )
            } else { # doing in between years, so get all data for this year
                query <- paste("SELECT a.latitude, a.longitude, b.*
        FROM AWAPDaily.latlon as a
        , AWAPDaily.", years[i], " as b
        where (a.id = b.id) and (a.latitude between ",
                    lat1, " and ", lat2, ") and (a.longitude between ",
                    lon1, " and ", lon2, ")
        order by b.day",
                    sep = ""
                )
            }

            response <- dbGetQuery(channel, query)

            # Detection of broken records

            if (nrow(response) == 0) {
                fail_flag <- 1
                broken_records <- rbind(
                    broken_records,
                    data.frame(species, reg, lat, lon)
                )
                break
            }

            # If the query worked properly, add it to the output table
            # This output table contains all the variables for each ind, but without the necessary processing. 
            # We will use it as a base to populate the empty individual data afterwards
            if (i == 1) {
                output <- cbind(reg, response) 
            } else {
                output <- rbind(output, cbind(reg, response))
            }
        } # end loop through years

        killDbConnections() # end the query

        # If the query failed, it will continue to next location
        # If the query worked:
        if (fail_flag != 1) {
            #### end loop to detect broken records

            # Add dates to the output
            output <- cbind(dates, output)

            sol <- as.numeric(output$sol) # making sure it is numeric for the clouds calculation below
            yearlist <- years
            cloud_cover <- get_cloud(loc, yearlist, sol) # adds the clouds from another code

            # Extracts only the summer records and processes them
            # In Laura's experiment: this part will work for all the months in the activity period.
            month_records <- output[format(output$dates, "%m") %in% months, ]
            print(head(month_records))
            avg_max_temp <- mean(as.numeric(month_records$tmax), na.rm = TRUE)
            avg_min_temp <- mean(month_records$tmin, na.rm = TRUE)
            avg_rr <- mean(month_records$rr, na.rm = TRUE)
            avg_vpr <- mean(month_records$vpr, na.rm = TRUE)
            month_records$sol <- sapply(month_records$sol,
                gsub,
                pattern = "\r", replacement = ""
            )
            month_records$sol <- as.numeric(month_records$sol)
            avg_sol <- mean(month_records$sol, na.rm = TRUE)
            avg_temp_over_35 <- nrow(month_records[month_records$tmax > 35, ]) / number_of_years
            
            
            output$sol <- sapply(output$sol, gsub, pattern = "\r", replacement = "")
            output$sol <- as.numeric(as.character(output$sol))
            avg_year_sol <- mean(output$sol, na.rm = TRUE) # Yearly, no need for extra modification
            avg_year_vpr <- mean(output$vpr, na.rm = TRUE) # Yearly, no need for extra modification

            
            #In Laura's experiment: I added this modification to also get the values for the collection month 
            
            month_records_Coll <- output[format(output$dates, "%m") %in% months_Coll, ]
            avg_max_temp_Coll <- mean(as.numeric(month_records_Coll$tmax), na.rm = TRUE)
            avg_min_temp_Coll <- mean(month_records_Coll$tmin, na.rm = TRUE)
            avg_rr_Coll <- mean(month_records_Coll$rr, na.rm = TRUE)
            avg_vpr_Coll <- mean(month_records_Coll$vpr, na.rm = TRUE)
            
            
            month_records_Coll$sol <- sapply(month_records_Coll$sol,
                                        gsub,
                                        pattern = "\r", replacement = "") # Sometimes we get an enter character in this variable. Eliminate this character
            month_records_Coll$sol <- as.numeric(month_records_Coll$sol) # And convert to numeric again
            avg_sol_Coll <- mean(month_records_Coll$sol, na.rm = TRUE)
            avg_temp_over_35_Coll <- nrow(month_records_Coll[month_records_Coll$tmax > 35, ]) / number_of_years
            
            
            # Create a data frame with one row with the results of one individual
            # below it will be added to the original empty data frame called individual data
            # Since it is a loop, one row for each individual will be added consecutively
            
            this_data <- data.frame(
                species,
                reg,
                lat,
                lon,
                avg_temp_over_35,
                avg_max_temp,
                avg_min_temp,
                avg_sol,
                avg_year_sol,
                avg_year_vpr,
                cloud_cover,
                avg_rr,
                avg_vpr,
                
                # Also add the ones for the collection month
                avg_temp_over_35_Coll, 
                avg_max_temp_Coll, 
                avg_min_temp_Coll,
                avg_sol_Coll,
                avg_rr_Coll,
                avg_vpr_Coll
                
            )

            individual_data <- rbind(individual_data, this_data) # Populate the empty data frame called individual data

            if (j == 1) {
                results <- output
            } else {
                results <- rbind(results, output)
            }
        }
    } # end loop through sites

    ###############

    # Final results:
    total_dataset <- individual_data

    # Check which ones didn't work:

    check_missing <- merge(total_dataset, source_dataset, "reg", all.y = TRUE)

    # Here paste the two things:
    res <- list(
        "total_dataset" = total_dataset,
        "check_missing" = check_missing
    )

    return(res)
}