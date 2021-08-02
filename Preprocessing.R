library(dplyr)
library(tidyverse)
library(lubridate) 
library(ggplot2)
library(survival)

options(tibble.print_min = 50)

# set current path to this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# sets local system language to English
Sys.setlocale("LC_ALL","English")

# clear the workspace
rm(list = ls())

# read the complete data
# download it from 
# https://data.melbourne.vic.gov.au/Transport/On-street-Car-Parking-Sensor-Data-2019/7pgd-bdf2
# (several gigabytes of data!)
data <- read.csv("Data/On-street_Car_Parking_Sensor_Data_-_2019.csv", header = TRUE)

colnames(data)[1] <- "DeviceId"

markers <- sort(unique(data$StreetMarker))
areas <- unique(data$Area)
streets <- unique(data$StreetName)

data <- as_tibble(data) %>% 
  mutate(ArrivalTime = as.POSIXct(as.character(ArrivalTime), format = "%m/%d/%Y %I:%M:%S %p", tz = "Australia/Melbourne"), 
         DepartureTime = as.POSIXct(as.character(DepartureTime), format = "%m/%d/%Y %I:%M:%S %p", tz = "Australia/Melbourne"),
         State = 1*(VehiclePresent == "false") + 2*(VehiclePresent == "true"),
         interval = interval(ArrivalTime, DepartureTime),
         DurationMinutes = as.numeric(as.character(DurationMinutes)),
         DurationMinutes2 = round(time_length(interval, unit = "minutes"), 2)) %>%
  arrange(StreetMarker, ArrivalTime) %>%
  relocate(State, .after = StreetMarker)

data.new <- NULL

for (m in 1:length(markers)){
  print(m)
  data.marker <- filter(data, StreetMarker == markers[m])
  data.marker <- data.marker %>% mutate(delete = 0)
  ind <- which(data.marker$DepartureTime[1:(nrow(data.marker)-1)] == data.marker$ArrivalTime[2:nrow(data.marker)] & 
                 data.marker$State[1:(nrow(data.marker)-1)] == data.marker$State[2:nrow(data.marker)] & 
                 data.marker$StreetMarker[1:(nrow(data.marker)-1)] == data.marker$StreetMarker[2:nrow(data.marker)])
  data.marker$delete[ind] <- 1
  i <- 1
  while (i <= length(ind)) {
    j <- 1
    while (is.element(ind[i]+j, ind)){
      j <- j+1
    }
    k <- ind[i]
    data.marker$ArrivalTime[k+j] <- data.marker$ArrivalTime[k]
    data.marker$DurationMinutes[k+j] <- sum(data.marker$DurationMinutes[k:(k+j)])
    data.marker$interval[k+j] <- interval(data.marker$ArrivalTime[k+j], data.marker$DepartureTime[k+j])
    data.marker$DurationMinutes2[k+j] <- round(time_length(data.marker$interval[k+j], unit = "minutes"), 2)
    i <- i+j
  }
  data.new <- bind_rows(data.new, data.marker)
}



data <- data.new %>% 
  filter(delete == 0) %>%
  mutate(d.start = yday(ArrivalTime),
         d.stop = yday(DepartureTime),
         h.start = hour(ArrivalTime),
         h.stop = hour(DepartureTime),
         m.start = 60*h.start + minute(ArrivalTime),
         m.stop = 60*h.stop + minute(DepartureTime),
         weekday = factor(weekdays(ArrivalTime, abbreviate = TRUE), levels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")),
         same.day = 1*(d.start == d.stop),
         irregular = 1*(ArrivalTime >= DepartureTime)) %>%
  relocate(DeviceId, .after = irregular)

saveRDS(data, "Data/data_2019_clean.rds")

###

# Melbourne parking lots in the CBD
# download: https://data.melbourne.vic.gov.au/Transport/On-street-Parking-Bays/crvt-b4kt
markers <- unique(data$StreetMarker)
parking.lots.old <- as_tibble(read.csv("layer_0_2019.csv", header = TRUE))
parking.lots <- as_tibble(read.csv("layer_0_2020.csv", header = TRUE))
missing.new <- markers[which(!is.element(markers, parking.lots$marker_id))]
parking.lots <- bind_rows(parking.lots, parking.lots.old[which(is.element(parking.lots.old$marker_id, missing.new)), ])
parking.lots$lon <- NA
parking.lots$lat <- NA
parking.lots$sensor <- factor(0, levels = c(0, 1))


for (i in 1:nrow(parking.lots)) {
  if (parking.lots$marker_id[i] != "" & is.element(parking.lots$marker_id[i], markers)){
    parking.lots$sensor[i] <- 1
  }
  polyg <- substring(gsub('.{3}$', '', parking.lots$the_geom[i]), first = 17)
  coord <- strsplit(polyg, ",")[[1]]
  lon.lat <- c(0, 0)
  for (j in 1:length(coord)) {
    coord[j] <- trimws(gsub("\\(|\\)", "", coord[j]), which = "both")
    lon.lat <- lon.lat + as.numeric(strsplit(coord[j], " +")[[1]])
  }
  lon.lat <- lon.lat/length(coord)
  parking.lots$lon[i] <- lon.lat[1]
  parking.lots$lat[i] <- lon.lat[2]
}

colnames(parking.lots)[2] <- "StreetMarker"

saveRDS(parking.lots, "Data/Parking_Lots.rds")