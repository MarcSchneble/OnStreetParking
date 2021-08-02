options(tibble.print_min = 50)

# set current path to this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# sets local system language to English
Sys.setlocale("LC_ALL","English")

# clear the workspace
rm(list = ls())

library(ggmap)
library(dplyr)
library(lubridate) 
library(ggplot2)
library(survival)
library(grDevices)
library(survminer)
library(readr)
library(geonet)

# load data
data <- read_rds("Data/data_2019_clean.rds")
parking.lots <- read_rds("Data/Parking_Lots.rds")



# sensor or not

# add google API key
register_google("AIzaSyDX0CVJsIDBJF8NVFPAH84oLWPfvPa335Y")

theme_set(theme_bw() + theme(panel.grid = element_blank(), 
                             plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
                             legend.text = element_text(size = 20),
                             legend.title = element_text(size = 20),
                             axis.text = element_text(size = 20),
                             axis.title = element_text(size = 20)
))


map_Melbourne <- ggmap(get_map(
  location = c(144.963, -37.814), 
  maptype = "roadmap", zoom = 14)) + 
  geom_point(data = parking.lots, aes(x = lon, y = lat, color = sensor)) + 
  labs(x = "Longitude", y = "Latitude",
       color = "In-ground sensor") + 
  scale_color_hue(labels = c("No", "Yes")) + 
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  theme(legend.position = "bottom")

pdf(file = "Plots/ParkingBaysMelbourne.pdf", height = 10, width = 10)
print(map_Melbourne)
dev.off()

# duration not occupied

G <- readRDS("Data/network_geonet.rds")

data.free <- filter(data, State == 1, same.day == 1, irregular == 0, DurationMinutes2 > 0,
                    m.start >= 480, m.stop <= 1200)

data.free.grouped <- data.free %>% group_by(StreetMarker) %>%
  summarize(mean.free = mean(DurationMinutes2)) %>%
  ungroup() %>% arrange(mean.free)

data.free.grouped <- left_join(data.free.grouped, parking.lots, by = "StreetMarker") %>%
  select(StreetMarker, mean.free, lon, lat, sensor) %>%
  mutate(mean.free.cat = 1*(mean.free <= 15) + 2*(mean.free > 15 & mean.free <= 30) + 
           3*(mean.free > 30 & mean.free <= 45) + 4*(mean.free > 45 & mean.free <= 60) +
           5*(mean.free > 60),
         mean.free.cat = factor(mean.free.cat))

map_Melbourne.free <- ggmap(get_map(
  location = c(144.963, -37.814), 
  maptype = "roadmap", zoom = 15)) + 
  geom_segment(data = G$network$lins, aes(x = v1_x, y = v1_y, xend = v2_x, yend = v2_y), size = 0.8) +
  geom_segment(data = G$network$lins %>% filter(seg == 124 | seg == 143),
               aes(x = v1_x, y = v1_y, xend = v2_x, yend = v2_y), size = 5) +
  geom_point(data = data.free.grouped %>% filter(sensor == 1), aes(x = lon, y = lat, color = mean.free.cat)) +
  labs(x = "Longitude", y = "Latitude",
       color = "Average duration in \nstate 0 (in minutes)") +
  scale_color_hue(labels = c(expression(""<=" 15"), "15-30", "30-45", "45-60", "> 60")) + 
  guides(color = guide_legend(override.aes = list(size = 5), title.vjust = 3)) + 
  theme(legend.position = "bottom")

pdf(file = "Plots/DurationFree.pdf", height = 10, width = 10)
print(map_Melbourne.free)
dev.off()



# duration occupied

data.occ <- filter(data, State == 2, same.day == 1, irregular == 0, DurationMinutes2 > 0,
                   m.start >= 480, m.stop <= 1200)

data.occ.grouped <- data.occ %>% group_by(StreetMarker) %>%
  summarize(mean.occ = mean(DurationMinutes2)) %>%
  ungroup() %>% arrange(mean.occ)

data.occ.grouped <- left_join(data.occ.grouped, parking.lots, by = "StreetMarker") %>%
  select(StreetMarker, mean.occ, lon, lat, sensor) %>%
  mutate(mean.occ.cat = 1*(mean.occ <= 15) + 2*(mean.occ > 15 & mean.occ <= 30) + 
           3*(mean.occ > 30 & mean.occ <= 45) + 4*(mean.occ > 45 & mean.occ <= 60) +
           5*(mean.occ > 60),
         mean.occ.cat = factor(mean.occ.cat))

map_Melbourne.occ <- ggmap(get_map(
  location = c(144.963, -37.814), 
  maptype = "roadmap", zoom = 15)) + 
  geom_segment(data = G$network$lins, aes(x = v1_x, y = v1_y, xend = v2_x, yend = v2_y), size = 0.8) +
  geom_segment(data = G$network$lins %>% filter(seg == 124 | seg == 143),
               aes(x = v1_x, y = v1_y, xend = v2_x, yend = v2_y), size = 5) +
  geom_point(data = data.occ.grouped %>% filter(sensor == 1), aes(x = lon, y = lat, color = mean.occ.cat)) +
  labs(x = "Longitude", y = "Latitude",
       color = "Average duration in \nstate 1 (in minutes)") +
  scale_color_hue(labels = c(expression(""<=" 15"), "15-30", "30-45", "45-60", "> 60")) + 
  guides(color = guide_legend(override.aes = list(size = 5), title.vjust = 3)) + 
  theme(legend.position = "bottom") 


pdf(file = "Plots/DurationOccupied.pdf", height = 10, width = 10)
print(map_Melbourne.occ)
dev.off()


# kaplan meier estimates in lonsdale street between queen str and elizabeth str

data.lonsdale.west <- filter(data, StreetName == "LONSDALE STREET")#, 
                        # BetweenStreet1 == "QUEEN STREET" | BetweenStreet2 == "QUEEN STREET",
                        # BetweenStreet1 == "ELIZABETH STREET" | BetweenStreet2 == "ELIZABETH STREET")

data.lonsdale.west.free <- filter(data.lonsdale.west, State == 1, same.day == 1, irregular == 0, DurationMinutes2 > 0,
                             m.start >= 480, m.stop <= 1200)

data.lonsdale.west.occ <- filter(data.lonsdale.west, State == 2, same.day == 1, irregular == 0, DurationMinutes2 > 0,
                            m.start >= 480, m.stop <= 1200)


km.lonsdale.west.free <- survfit(formula = Surv(time = DurationMinutes2, 
                                           event = rep(1, nrow(data.lonsdale.west.free))) ~ weekday,
                            type = "kaplan-meier", data = data.lonsdale.west.free)
survplot.lonsdale.west.free <- ggsurvplot(fit = km.lonsdale.west.free, data = data.lonsdale.west.free, surv.median.line = "hv",
                                     conf.int = TRUE,
                                     xlab = "Availability in Minutes", 
                                     ylab = "Kaplan-Meier estimator of the duration in state 0", legend = c(0.90, 0.85),
                                     legend.title = "Weekday",
                                     legend.labs = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"),
                                     break.x.by = 5,
                                     xlim = c(0, 60),
                                     ggtheme=theme_bw() + theme(panel.grid = element_blank(), 
                                                                plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
                                                                legend.text = element_text(size = 20),
                                                                legend.title = element_text(size = 20),
                                                                axis.text = element_text(size = 20),
                                                                axis.title = element_text(size = 20)))
pdf(file = "Plots/AvailabilityLonsdaleWest.pdf", width = 10, height = 10, onefile = FALSE)
print(survplot.lonsdale.west.free)
dev.off()

km.lonsdale.west.occ <- survfit(formula = Surv(time = DurationMinutes2, 
                                           event = rep(1, nrow(data.lonsdale.west.occ))) ~ weekday,
                            type = "kaplan-meier", data = data.lonsdale.west.occ)
survplot.lonsdale.west.occ <- ggsurvplot(fit = km.lonsdale.west.occ, data = data.lonsdale.west.occ, surv.median.line = "hv",
                                     conf.int = TRUE,
                                     xlab = "Duration in Minutes", 
                                     ylab = "Kaplan Meier Estimator", legend = c(0.90, 0.85),
                                     title = "Parking duration in Lonsdale Street \nbetween Queen Street and Elizabeth Street",
                                     legend.title = "Weekday",
                                     legend.labs = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"),
                                     break.x.by = 5,
                                     xlim = c(0, 60),
                                     ggtheme=theme_bw() + theme(panel.grid = element_blank(), 
                                                                plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
                                                                legend.text = element_text(size = 20),
                                                                legend.title = element_text(size = 20),
                                                                axis.text = element_text(size = 20),
                                                                axis.title = element_text(size = 20))) 
pdf(file = "Plots/DurationLonsdaleWest.pdf", width = 10, height = 10, onefile = FALSE)
print(survplot.lonsdale.west.occ)
dev.off()


# kaplan meier estimates in lonsdale street between russel str and exhibition str

data.lonsdale.east <- filter(data, StreetName == "LONSDALE STREET",
                             BetweenStreet1 == "RUSSELL STREET" | BetweenStreet1 == "EXHIBITION STREET",
                             BetweenStreet2 == "EXHIBITION STREET" | BetweenStreet2 == "SPRING STREET")

data.lonsdale.east.free <- filter(data.lonsdale.east, State == 1, same.day == 1, irregular == 0, DurationMinutes2 > 0,
                                  m.start >= 480, m.stop <= 1200)

data.lonsdale.east.occ <- filter(data.lonsdale.east, State == 2, same.day == 1, irregular == 0, DurationMinutes2 > 0,
                                 m.start >= 480, m.stop <= 1200)


km.lonsdale.east.free <- survfit(formula = Surv(time = DurationMinutes2, 
                                                event = rep(1, nrow(data.lonsdale.east.free))) ~ SideOfStreetCode,
                                 type = "kaplan-meier", data = data.lonsdale.east.free)
survplot.lonsdale.east.free <- ggsurvplot(fit = km.lonsdale.east.free, data = data.lonsdale.east.free, surv.median.line = "hv",
                                          conf.int = TRUE,
                                          xlab = "Duration in state 0 (in minutes)", 
                                          ylab = "Kaplan-Meier estimator", legend = c(0.75, 0.85),
                                          legend.title = "Relative location of parking lot\non the street (side of street)",
                                          legend.labs = c("central", "north", "south"), 
                                          break.x.by = 5,
                                          xlim = c(0, 60),
                                          ggtheme=theme_bw() + theme(panel.grid = element_blank(), 
                                                                     plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
                                                                     legend.text = element_text(size = 20),
                                                                     legend.title = element_text(size = 20),
                                                                     axis.text = element_text(size = 20),
                                                                     axis.title = element_text(size = 20)))
pdf(file = "Plots/AvailabilityLonsdaleEast.pdf", width = 10, height = 10, onefile = FALSE)
print(survplot.lonsdale.east.free)
dev.off()

km.lonsdale.east.occ <- survfit(formula = Surv(time = DurationMinutes2, 
                                               event = rep(1, nrow(data.lonsdale.east.occ))) ~ SideOfStreetCode,
                                type = "kaplan-meier", data = data.lonsdale.east.occ)
survplot.lonsdale.east.occ <- ggsurvplot(fit = km.lonsdale.east.occ, data = data.lonsdale.east.occ, surv.median.line = "hv",
                                         conf.int = TRUE,
                                         xlab = "Duration in state 1 (in minutes)", 
                                         ylab = "Kaplan-Meier estimator", legend = c(0.75, 0.85),
                                         legend.title = "Relative location of parking lot\non the street (side of street)",
                                         legend.labs = c("central", "north", "south"),
                                         break.x.by = 5,
                                         xlim = c(0, 60),
                                         ggtheme=theme_bw() + theme(panel.grid = element_blank(), 
                                                                    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
                                                                    legend.text = element_text(size = 20),
                                                                    legend.title = element_text(size = 20),
                                                                    axis.text = element_text(size = 20),
                                                                    axis.title = element_text(size = 20))) 
pdf(file = "Plots/DurationLonsdaleEast.pdf", width = 10, height = 10, onefile = FALSE)
print(survplot.lonsdale.east.occ)
dev.off()

