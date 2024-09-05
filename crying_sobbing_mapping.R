<<<<<<< HEAD
#trying to make a map! 
#install.packages("leaflet")   # Install the 'leaflet' package if you haven't already
#install.packages('sp')
library(leaflet)              # Load the 'leaflet' package
library(ggplot2)
library(dplyr)
AT <- read.csv("AT_data.csv", stringsAsFactors = TRUE)
ATmini <- read.csv("Mini_AT_data.csv", stringsAsFactors = TRUE)

############# Test run with the Mini dataset 
deduplicated_ATmini <- distinct(ATmini, Lat_provenance, .keep_all = TRUE)
print(deduplicated_ATmini$Lat_provenance)
#############Actually collapsing rows from the big AT dataset 
deduplicated_ATbig <- distinct(AT, Lat_provenance, .keep_all = TRUE) 
print(deduplicated_ATbig$Lat_provenance)
##############

# Example data frame with latitude values
data <- data.frame(
  location = as.character(deduplicated_ATbig$Site_provenance),
  latitude = as.numeric(deduplicated_ATbig$Lat_provenance),
  longitude = as.numeric(deduplicated_ATbig$Lon_provenance)
)

map <- leaflet(data = data) %>%
  addTiles()   # Add default OpenStreetMap tiles as the base layer

# Add markers for each location using their latitudes and longitudes
map <- map %>%
  addMarkers(
    lng = ~longitude,    # Use the longitude column from your data frame
    lat = ~latitude,     # Use the latitude column from your data frame
    popup = ~location    # Display location names in the popups
  )

# Display the map
map





=======
#trying to make a map! 
#install.packages("leaflet")   # Install the 'leaflet' package if you haven't already
#install.packages('sp')
library(leaflet)              # Load the 'leaflet' package
library(ggplot2)
library(dplyr)
AT <- read.csv("AT_data.csv", stringsAsFactors = TRUE)
ATmini <- read.csv("Mini_AT_data.csv", stringsAsFactors = TRUE)

############# Test run with the Mini dataset 
deduplicated_ATmini <- distinct(ATmini, Lat_provenance, .keep_all = TRUE)
print(deduplicated_ATmini$Lat_provenance)
#############Actually collapsing rows from the big AT dataset 
deduplicated_ATbig <- distinct(AT, Lat_provenance, .keep_all = TRUE) 
print(deduplicated_ATbig$Lat_provenance)
##############

# Example data frame with latitude values
data <- data.frame(
  location = as.character(deduplicated_ATbig$Site_provenance),
  latitude = as.numeric(deduplicated_ATbig$Lat_provenance),
  longitude = as.numeric(deduplicated_ATbig$Lon_provenance)
)

map <- leaflet(data = data) %>%
  addTiles()   # Add default OpenStreetMap tiles as the base layer

# Add markers for each location using their latitudes and longitudes
map <- map %>%
  addMarkers(
    lng = ~longitude,    # Use the longitude column from your data frame
    lat = ~latitude,     # Use the latitude column from your data frame
    popup = ~location    # Display location names in the popups
  )

# Display the map
map





>>>>>>> 9c95ca3e6eea13d6f9f205838021ff0f9f6845bd
