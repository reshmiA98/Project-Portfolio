######################### Propogate the TLE Data ###############################
################################################################################
#Documentation
# https://cran.r-project.org/web/packages/asteRisk/vignettes/asteRisk.html
# this documentation is pretty step by step

#######Install packages#######
##install.packages("asteRisk")
##install.packages('asteRiskData',
##                 repos='https://rafael-ayala.github.io/drat/')
##install.packages("sf")
##install.packages("tidyverse")
library(asteRisk)
library(asteRiskData)
library(tidyr)
library(plotly)
library(lazyeval)
library(dplyr)

# Read the TLE data for Starlink satellites and BRIK_II
# Replace 'path_to_tle_file' with the actual path to your TLE files
all_tle_20231206_2039 <- readTLE("TLE_20231206_2039.txt")
print(all_tle_20231206_2039)

# Filter all TLEs containing "Starlink" in their name
starlink_tles <- Filter(function(x) grepl("STARLINK", x$objectName), all_tle_20231206_2039)
print(starlink_tles)
non_starlink_tles <- Filter(function(x) !grepl("STARLINK", x$objectName), all_tle_20231206_2039)
non_starlink_tles
non_starlink_leo<- Filter(function(x) x$meanMotion >11, non_starlink_tles)
length(non_starlink_leo)
# Filter the TLE for BRIK-II
brik_ii_tle <- Filter(function(x) grepl("BRIK-II", x$objectName), all_tle_20231206_2039)
print(brik_ii_tle)

# Find the index of the BRIK-II TLE in the original list
brik_ii_index <- which(sapply(all_tle_20231206_2039, function(x) grepl("BRIK-II", x$objectName)))

# Print the index
print(brik_ii_index)


# TLE number 3631 contains a state vectors for BRIK-II
brik_ii_tle <- all_tle_20231206_2039[[3631]]
brik_ii_tle

# Test a calculation for how many min to completed single orbit
# 0.0653615 is the inverse of mean motion or 0.06 days/orbit
# 0.0653615×60×24 min/orbit = 94.12056 min per single orbit
# or 94.12056 / 24 = 3.92169 revolutions per day
1/brik_ii_tle$meanMotion


########### build and select best model "sgp4" is selected
targetTimes <- seq(0, 720, by = 1)

brik_ii_results_position_matrix <- matrix(nrow = length(targetTimes), ncol = 3)
brik_ii_results_velocity_matrix <- matrix(nrow = length(targetTimes), ncol = 3)

for (i in 1:length(targetTimes)) {
  new_result <- sgdp4(n0 = brik_ii_tle$meanMotion * ((2 * pi)/(1440)), 
                      e0 = brik_ii_tle$eccentricity,
                      i0 = brik_ii_tle$inclination * pi/180, 
                      M0 = brik_ii_tle$meanAnomaly * pi/180, 
                      omega0 = brik_ii_tle$perigeeArgument * pi/180, 
                      OMEGA0 = brik_ii_tle$ascension * pi/180, 
                      Bstar = brik_ii_tle$Bstar, 
                      initialDateTime = brik_ii_tle$dateTime,
                      targetTime = targetTimes[i])
  brik_ii_results_position_matrix[i, ] <- new_result[[1]]
  brik_ii_results_velocity_matrix[i, ] <- new_result[[2]]
}
last_brik_ii_propagation <- new_result
brik_ii_results_position_matrix = cbind(brik_ii_results_position_matrix, targetTimes)
colnames(brik_ii_results_position_matrix) <- c("x", "y", "z", "time")
brik_ii_results_position_matrix

#### Predict ALL starlink satellites

########################## STARTING FROM HERE TAKES A TON OF PROCESSING POWER     #################
########################## SKIP AND READ "satellite_distances_720_1_20231206.csv" #################
########################## STARTING FROM THE NEXT LINES OF ##                     #################
##########################            TO SAVE PROCESSING POWER                    #####################

##Create empty lists that the prediction will eventually fill in
all_starlink_results_position_matrix <- list()
all_starlink_results_velocity_matrix <- list()

for (tle in starlink_tles) {
  results_position_matrix <- matrix(nrow = length(targetTimes), ncol = 3)
  results_velocity_matrix <- matrix(nrow = length(targetTimes), ncol = 3)
  
## Nested for-loop to loop through predictions for every satellite in the starlink_tle list  
  for (i in 1:length(targetTimes)) {
    new_result <- sgdp4(n0 = tle$meanMotion * ((2 * pi)/(1440)), 
                        e0 = tle$eccentricity,
                        i0 = tle$inclination * pi/180, 
                        M0 = tle$meanAnomaly * pi/180, 
                        omega0 = tle$perigeeArgument * pi/180, 
                        OMEGA0 = tle$ascension * pi/180, 
                        Bstar = tle$Bstar, 
                        initialDateTime = tle$dateTime,
                        targetTime = targetTimes[i])
    results_position_matrix[i, ] <- new_result[[1]]
    results_velocity_matrix[i, ] <- new_result[[2]]
  }
  all_starlink_results_position_matrix[[tle$objectName]] <- results_position_matrix
  all_starlink_results_velocity_matrix[[tle$objectName]] <- results_velocity_matrix
}

for (i in 1:length(all_starlink_results_position_matrix)) {
  satellite_matrix <- all_starlink_results_position_matrix[[i]]
  # Add targetTimes as a new column with the same value repeated
  satellite_matrix_with_time <- cbind(satellite_matrix, time = rep(targetTimes, 
  nrow(satellite_matrix)))
  # Update the matrix in the list
  all_starlink_results_position_matrix[[i]] <- satellite_matrix_with_time}

##Push the starlink predictions into a dataframe
all_starlink_positions_df <- data.frame()
for (i in 1:length(all_starlink_results_position_matrix)) {
  satellite_matrix <- all_starlink_results_position_matrix[[i]]
  satellite_df <- as.data.frame(satellite_matrix)
  names(satellite_df) <- c("x", "y", "z", "time")
  satellite_df$satellite <- names(all_starlink_results_position_matrix)[i]
  all_starlink_positions_df <- rbind(all_starlink_positions_df, satellite_df)
}

## Pushing brik-ii's prediction into a dateframe. Can't do the same way 
## for starlink's prediction
brik_ii_df <- as.data.frame(brik_ii_results_position_matrix)
names(brik_ii_df) <- c("x", "y", "z", "time")
brik_ii_df$satellite <- "BRIK-II"

##combining dataframes for euclidean distance purposes
combined_df <- rbind(brik_ii_df, all_starlink_positions_df)
combined_df

brik_ii_data <- combined_df[combined_df$satellite == "BRIK-II",]

## Doing a 3d euclidean distance for every starlink satellite against brik-ii
calculate_distances <- function(time) {
  # Ensure brik_ii_pos is a single row
   brik_ii_pos <- brik_ii_data[brik_ii_data$time == time, c("x", "y", "z")]
   if (nrow(brik_ii_pos) != 1) {
    stop("Multiple or no entries for BRIK-II at time ", time)
  }
  # Filter other satellites at the same time
  other_sats <- combined_df[combined_df$time == time & combined_df$satellite != "BRIK-II", ]
  other_sats$distance <- sqrt((other_sats$x - brik_ii_pos$x)^2 + 
                                (other_sats$y - brik_ii_pos$y)^2 + 
                                (other_sats$z - brik_ii_pos$z)^2)
  return(other_sats)
}
## Apply this function for each time point
 distances <- lapply(unique(combined_df$time), calculate_distances)
###### Combine the results into a single dataframe
 distance_df <- do.call(rbind, distances)
 distance_df

write.csv(distance_df, file="satellite_distances_720_1_20231206.csv")
########################## CONTINUE TO PROCESS FROM HERE #####################
########################## CONTINUE TO PROCESS FROM HERE #####################
########################## CONTINUE TO PROCESS FROM HERE #####################
########################## CONTINUE TO PROCESS FROM HERE #####################

# Let´s verify that the SDP4 algorithm was automatically chosen
last_brik_ii_propagation$algorithm

###Finding closest satellites, closest within 100km
filtered_df <- distance_df %>% 
  filter(distance <= 100) %>% 
  mutate(color = case_when(
    distance < 1 ~ "red",
    distance < 2 ~ "yellow",
    TRUE ~ "blue"
  ))
filtered_df
closest_satellites <- filtered_df %>% 
  group_by(satellite) %>% 
  summarise(closest_distance = min(distance)) %>% 
  arrange(closest_distance) %>% 
  head(5) %>% 
  pull(satellite)

## write csv that shows all possible conjunction events, only the events
write.csv(filtered_df, file="poss_conjunction_events.csv")

## grabbing all rows of the satellites that have potential conjunction events
satellites_within_100 <- distance_df %>%
  filter(distance <= 100) %>%
  distinct(satellite)

filtered_df <- distance_df %>%
  filter(satellite %in% satellites_within_100$satellite)
filtered_df

## Grabbing the top 3 satellites for animation purposes
brik_df <- brik_ii_df
brik_df$distance <- 0
brik_df
sat1_df <- filtered_df %>% filter(satellite == "STARLINK-30505")
sat2_df <- filtered_df %>% filter(satellite == "STARLINK-30405")
sat3_df <- filtered_df %>% filter(satellite == "STARLINK-30412")
sat4_df <- filtered_df %>% filter(satellite == "STARLINK-30292")
sat5_df <- filtered_df %>% filter(satellite == "STARLINK-30431")
sat6_df <- filtered_df %>% filter(satellite == "STARLINK-30414")
sat7_df <- filtered_df %>% filter(satellite == "STARLINK-30290")
sat8_df <- filtered_df %>% filter(satellite == "STARLINK-30392")
sat9_df <- filtered_df %>% filter(satellite == "STARLINK-30425")
sat10_df <- filtered_df %>% filter(satellite == "STARLINK-30272")

#################### 3d Plot the orbit for Brik-ii

# Function to create the accumulated dataframe for the animation
accumulate_by <- function(dat, var) {
  var <- lazyeval::f_eval(var, dat)
  lvls <- plotly:::getLevels(var)
  dats <- lapply(seq_along(lvls), function(x) {
    cbind(dat[var %in% lvls[seq(1, x)], ], frame = lvls[[x]])
  })
  dplyr::bind_rows(dats)
}

accumulated_df <- accumulate_by(as.data.frame(brik_ii_results_position_matrix), ~time)
accumulated_df
# Creating the orbit animation
sphere <- function(size, N_lat, N_lon) {
  theta <- seq(0, 2 * pi, length.out = N_lat)
  phi <- seq(0, pi, length.out = N_lon)
  
  # Coordinates for points on the sphere
  x0 <- size * outer(cos(theta), sin(phi))
  y0 <- size * outer(sin(theta), sin(phi))
  z0 <- size * outer(rep(1, N_lat), cos(phi))
  
  list(x = x0, y = y0, z = z0)
}

# Sphere size and plot
radius <- 6371 # Earth's radius in kilometers
N_lat <- 50    # Number of latitude points
N_lon <- 50    # Number of longitude points
coords <- sphere(radius, N_lat, N_lon)

custom_color <- 'lightgrey'

# Plotly surface plot with custom light grey color
fig <- plot_ly() %>%
  add_surface(x = coords$x, y = coords$y, z = coords$z,
              cauto = FALSE,
              cmin = -1,
              cmax = 1,
              colorscale = list(c(0, "lightgrey"), c(1, "lightgrey"))) %>%
  layout(scene = list(aspectmode = "cube"))

# Show plot
fig

orbit_animation <- plot_ly(data = accumulated_df, x = ~x, y = ~y, z = ~z, type = "scatter3d",
                           mode = "lines", opacity = 0.8, line = list(width = 6, color = ~time),
                           frame = ~frame) %>%
  animation_opts(frame = 100) %>%
  layout(scene = list(xaxis = list(range = c(-7000, 7000)),
                      yaxis = list(range = c(-7000, 7000)),
                      zaxis = list(range = c(-7000, 7000))))

# Add the sphere as a separate trace
orbit_animation <- orbit_animation %>%
  add_trace(x = coords$x, y = coords$y, z = coords$z, type = 'surface',
            cauto = FALSE, cmin = -1, cmax = 1,
            colorscale = list(c(0, "lightgrey"), c(1, "lightgrey"))) %>%
  layout(scene = list(xaxis = list(title = 'X (km)'),
                      yaxis = list(title = 'Y (km)'),
                      zaxis = list(title = 'Z (km)')),
         title = 'Satellite Orbits and Earth')

# Show the combined plot
orbit_animation

###-------------- TOP 4 SATELLITES -----------####

combined_sat_df <- bind_rows(brik_df, sat1_df, sat2_df, sat3_df,sat4_df, sat5_df,
                             sat6_df, sat7_df, sat8_df, sat9_df, sat10_df, .id = "satellite_id")

combined_sat_df <- combined_sat_df %>%
  mutate(color = ifelse(satellite == "BRIK-II", "blue", "red"))
combined_sat_df

##Cut down to 10 minute increments due to processing intensiveness

# Accumulate data for animation
accumulated_df <- accumulate_by(combined_sat_df, ~time)
accumulated_df
# Creating the orbit animation
orbit_animation2 <- plot_ly(data = accumulated_df, x = ~x, y = ~y, z = ~z, type = "scatter3d",
                           mode = "lines", opacity = 0.8, line = list(width = 3, color = ~color),
                           frame = ~frame) %>%
  animation_opts(frame = 100) %>%
  layout(scene = list(xaxis = list(range = c(-7000, 7000)),
                      yaxis = list(range = c(-7000, 7000)),
                      zaxis = list(range = c(-7000, 7000))))

# Add the sphere as a separate trace
orbit_animation2 <- orbit_animation2 %>%
  add_trace(x = coords$x, y = coords$y, z = coords$z, type = 'surface',
            cauto = FALSE, cmin = -1, cmax = 1,
            colorscale = list(c(0, "lightgrey"), c(1, "lightgrey"))) %>%
  layout(scene = list(xaxis = list(title = 'X (km)'),
                      yaxis = list(title = 'Y (km)'),
                      zaxis = list(title = 'Z (km)')),
         title = 'Satellite Orbits and Earth')

# Show the combined plot
orbit_animation2
###################----- TOP 4 SATELLITES DONE ######