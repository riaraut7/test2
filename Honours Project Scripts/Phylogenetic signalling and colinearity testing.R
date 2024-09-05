#Blank 

# Load the necessary libraries
library(raster)
library(sp)
library(dplyr)
library(ggplot2)
library(plotly)
library(ggfortify)
library(geiger)

#### setting up your data #### 

#load dataset 
traits_data <- read.csv('phyloinfo_containing_dataset_new.csv', stringsAsFactors = TRUE)
mytree <- read.tree('your_working_tree.tre')
in_right_order <- read.csv('in_right_order_after_s.antarctici.csv', header = T)


traits_data <- traits_data %>% 
  select (Taxon, genus, family, av_breadth, av_topt, av_amax, av_lat, av_lon) %>% 
  filter(Taxon != 'Schistidium antarctici') %>% 
  filter(!is.na(av_lat))
rownames(traits_data) <- gsub(" ","_",traits_data$Taxon)

#### matching with phylogeny #### 
obj<-name.check(mytree,traits_data)
obj #Sort of!

#take out problematic species from your tree 
species_to_remove <- c("Aegilops_cylindrica", "Brassica_campestris", "Prosopis_juliflora", 'Atriplex_confertifolia',
                       'Helianthus_annuus', 'Solanum_lycopersicum', 'Atriplex_glabriuscula', 'Larix_decidua', 
                       'Solanum_tuberosum', 'Atriplex_hymenelytra', 'Oryza_sativa', 'Triticum_aestivum', 'Atriplex_sabulosa', 
                       'Phoenix_dactylifera')
# Remove the specified species from the tree
pruned_tree <- drop.tip(mytree, species_to_remove)

obj<-name.check(pruned_tree,traits_data)
obj #seems good! 
mytree=multi2di(pruned_tree)

#### Is there still temp signalling? #### 
#topt 
topt_trait <- traits_data[,5]
names(topt_trait)<-rownames(traits_data)

set.seed(0)
phylosig(mytree, topt_trait, method="lambda", test=TRUE, nsim=999)
## BEFORE taking out species with no lat values: 
# Phylogenetic signal lambda : 0.507889 
# logL(lambda) : -175.099 
# LR(lambda=0) : 4.19003 
# P-value (based on LR test) : 0.0406623 

## AFTER taking out species with no lat values: 
# Phylogenetic signal lambda : 0.507884 
# logL(lambda) : -137.012 
# LR(lambda=0) : 1.99581 
# P-value (based on LR test) : 0.157735


#### getting climatic data #### 
r <- getData("worldclim", var="bio", res=10)

###### STOPPED HERE 



#### testing for colinearity between climatic data and Topt #### 





#### random, ignore for now 
# Download WorldClim data 
r <- getData("worldclim", var="bio", res=10)


# Subset lat-longs from AT dataset


AT_df <- read.csv("final_curve_parameter_filter_2_info_11.3.24.csv")

#first, we extract some variables of interest
latlongsp<- data.frame(AT_df$Lat_provenance, AT_df$Long_provenance, AT_df$Taxon)
latlongsp<- na.omit(latlongsp)
latlongsp<-latlongsp %>% 
  rename(
    Latitude = AT_df.Lat_provenance ,
    Longitude = AT_df.Long_provenance, 
    Species = AT_df.Taxon
  )

latlong<- latlongsp[, c("Longitude", "Latitude")]

# Create a SpatialPoints object
points <- SpatialPoints(latlong, proj4string = r@crs)

# Extract the values
values <- extract(r, points)

# Combine the coordinates and values into a data frame, with a species column too! 
df <- cbind.data.frame(coordinates(points), values, latlongsp$Species)