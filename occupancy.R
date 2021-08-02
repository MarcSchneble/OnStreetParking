markers <- levels(data$StreetMarker)
data$frac0 <- data$frac1 <- NA
for (i in 1:length(markers)) {
  print(i)
  ind_distance <- which(distance[i, ] <= 50 & distance[i, ] > 0)
  ind_marker <- which(data$StreetMarker == markers[i])
  if (length(ind_distance) == 0) {
    print(paste("A", i))
    data$frac0[ind_marker] <- 0
    data$frac1[ind_marker] <- 0
  } else {
    data_marker <- filter(data, StreetMarker == markers[i])
    data_distance_0 <- filter(data, StreetMarker %in% markers[ind_distance], State == 0)
    data_distance_1 <- filter(data, StreetMarker %in% markers[ind_distance], State == 1)
    frac0 <- frac1 <- rep(NA, length(ind_marker))
    for (j in 1:nrow(data_marker)) {
      frac0[j] <- length(which(data_marker$ArrivalTime[j] %within% data_distance_0$interval))/length(ind_distance)
      frac1[j] <- length(which(data_marker$ArrivalTime[j] %within% data_distance_1$interval))/length(ind_distance)
    }
    data$frac0[ind_marker] <- frac0
    data$frac1[ind_marker] <- frac1
  }
}
saveRDS(data, "data_frac.rds")
