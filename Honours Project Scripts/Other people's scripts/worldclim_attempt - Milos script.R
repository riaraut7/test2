# Install the necessary packages
install.packages("raster")
install.packages("sp")


# Load the necessary libraries
library(raster)
library(sp)
library(dplyr)

# Download WorldClim data
r <- getData("worldclim", var="bio", res=10)

# Subset lat-longs from AT dataset
df<- read.csv("working_standard_units_SE_with_filter_2_info.csv")
latlong<- data.frame(df$Lat_provenance, df$Long_provenance)
latlong<- na.omit(latlong)

latlong<-latlong %>% 
  rename(
    Latitude = df.Lat_provenance ,
    Longitude = df.Long_provenance
  )

latlong<- latlong[, c("Longitude", "Latitude")]

# Create a SpatialPoints object
points <- SpatialPoints(latlong, proj4string = r@crs)

# Extract the values
values <- extract(r, points)

# Combine the coordinates and values into a data frame
df <- cbind.data.frame(coordinates(points), values)

# Correct names for the "bio" columns
correct_names <- c("Longitude", "Latitude", "Annual Mean Temperature", "Mean Diurnal Range", "Isothermality", 
                   "Temperature Seasonality", "Max Temperature of Warmest Month", 
                   "Min Temperature of Coldest Month", "Temperature Annual Range", 
                   "Mean Temperature of Wettest Quarter", "Mean Temperature of Driest Quarter", 
                   "Mean Temperature of Warmest Quarter", "Mean Temperature of Coldest Quarter", 
                   "Annual Precipitation", "Precipitation of Wettest Month", 
                   "Precipitation of Driest Month", "Precipitation Seasonality", 
                   "Precipitation of Wettest Quarter", "Precipitation of Driest Quarter", 
                   "Precipitation of Warmest Quarter", "Precipitation of Coldest Quarter")

# Rename the columns
colnames(df) <- correct_names


# Convert WorldClim data to actual values
df <- df / 10

# Create dataframe without latlongs

df_climate<- df[, c("Annual Mean Temperature", "Mean Diurnal Range", "Isothermality", 
                    "Temperature Seasonality", "Max Temperature of Warmest Month", 
                    "Min Temperature of Coldest Month", "Temperature Annual Range", 
                    "Mean Temperature of Wettest Quarter", "Mean Temperature of Driest Quarter", 
                    "Mean Temperature of Warmest Quarter", "Mean Temperature of Coldest Quarter", 
                    "Annual Precipitation", "Precipitation of Wettest Month", 
                    "Precipitation of Driest Month", "Precipitation Seasonality", 
                    "Precipitation of Wettest Quarter", "Precipitation of Driest Quarter", 
                    "Precipitation of Warmest Quarter", "Precipitation of Coldest Quarter")]

#Remove NAs

df_climate<- na.omit(df_climate)

# Make a correlation plot with all of the climate variables

library(corrplot)
library(psych)
library(ggplot2)

df_climate_cor <- cor(df_climate)
pval <- corr.test(df_climate_cor, adjust="none")$p

png(height=1000, width=1000, file="climate_corrplot.png", type = "cairo")
climate_corrplot<- corrplot(df_climate_cor, method = "number", type = "lower", number.cex = 1, pch.cex = 2,tl.cex=1, insig = "pch", p.mat = pval)
dev.off()


