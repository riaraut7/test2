

#let's load our libraries 
library(raster)
library(leaflet)            
library(ggplot2)
library(dplyr)
library(plotbiomes)

#load and clean data 
AT <- read.csv("final_curve_parameter_filter_2_info_11.3.24.csv", stringsAsFactors = TRUE)

#taking out the rows that don't have lat/lon values 
cleaned_data <- AT %>% 
  filter(is.na(june_failure))


#### Load worldclim variables all #### 
#loading the world clim data 
#NOTE: only need 1 and 12 (mean annual temp and precip) for whittakers, so you 
#don't need to run ALL of these lines 

fn_bio1 = 'wc10/wc2.1_10m_bio_1.tif' # BIO1 = Mean annual temp
fn_bio2 = 'wc10/wc2.1_10m_bio_2.tif'
fn_bio3 = 'wc10/wc2.1_10m_bio_3.tif' #isothermality 
fn_bio4 = 'wc10/wc2.1_10m_bio_4.tif'
fn_bio5 = 'wc10/wc2.1_10m_bio_5.tif'
fn_bio6 = 'wc10/wc2.1_10m_bio_6.tif'
fn_bio7 = 'wc10/wc2.1_10m_bio_7.tif'
fn_bio8 = 'wc10/wc2.1_10m_bio_8.tif'
fn_bio9 = 'wc10/wc2.1_10m_bio_9.tif'
fn_bio10 = 'wc10/wc2.1_10m_bio_10.tif' #BIO10 = mean warmest Q temp 
fn_bio11 = 'wc10/wc2.1_10m_bio_11.tif'
fn_bio12 = 'wc10/wc2.1_10m_bio_12.tif' 
fn_bio13 = 'wc10/wc2.1_10m_bio_13.tif'
fn_bio14 = 'wc10/wc2.1_10m_bio_14.tif'
fn_bio15 = 'wc10/wc2.1_10m_bio_15.tif'
fn_bio16 = 'wc10/wc2.1_10m_bio_16.tif'
fn_bio17 = 'wc10/wc2.1_10m_bio_17.tif'
fn_bio18 = 'wc10/wc2.1_10m_bio_18.tif'
fn_bio19 = 'wc10/wc2.1_10m_bio_19.tif'

bio1_raster = raster(fn_bio1) #annual mean temp 
bio2_raster = raster(fn_bio2) #mean diurnal range 
bio3_raster = raster(fn_bio3) #isothermality 
bio4_raster = raster(fn_bio4) #temp seasonality 
bio5_raster = raster(fn_bio5) #max temp of warmest month 
bio6_raster = raster(fn_bio6) #min temp of coldest month
bio7_raster = raster(fn_bio7) #temp annual range
bio8_raster = raster(fn_bio8) #mean temp of wettest quarter
bio9_raster = raster(fn_bio9) #mean temp of driest quarter
bio10_raster = raster(fn_bio10) #mmean temp of warmest quarter
bio11_raster = raster(fn_bio11) #mean temp of coldest quarter
bio12_raster = raster(fn_bio12) #annual precip 
bio13_raster = raster(fn_bio13) #precip of wettest month 
bio14_raster = raster(fn_bio14) #precip of driest month 
bio15_raster = raster(fn_bio15) #precipitation seasonality 
bio16_raster = raster(fn_bio16) #precip of wettest quarter
bio17_raster = raster(fn_bio17) #precip of driest quarter
bio18_raster = raster(fn_bio18) #preciptation of warmest quarter
bio19_raster = raster(fn_bio19) #precipitation of coldest quarter 


#### Choosing Specific Climatic variables of interest #### 
#variables needed for whittaker's space, and isothermality for interest 
cleaned_data$annual_temp = extract(bio1_raster, data.frame(x = cleaned_data$Long_provenance, 
                                                           y = cleaned_data$Lat_provenance))
cleaned_data$annual_precip = extract(bio12_raster, data.frame(x = cleaned_data$Long_provenance, 
                                                             y = cleaned_data$Lat_provenance))
cleaned_data$isothermality = extract(bio3_raster, data.frame(x = cleaned_data$Long_provenance, 
                                                              y = cleaned_data$Lat_provenance))

#looking at your distribution of isothermality out of curiosity
cleaned_isothermality_data <- cleaned_data %>% 
  filter(!is.na(isothermality))
max(cleaned_isothermality_data$isothermality)
min(cleaned_isothermality_data$isothermality)

#### Making the whittaker plot ####
cleaned_data <- cleaned_data %>% 
  filter(!is.na(Lat_provenance)) %>% 
  filter(!is.na(Long_provenance)) 

whittaker_base_plot()+
  theme_classic()+ 
  geom_point(data = cleaned_data, size =  3, alpha = 0.5, aes(x = annual_temp, y = annual_precip/10)) + 
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12))  
#print this 9x5.5 PDF 
class(whittaker_base_plot())


#unique latitude points (each dot in the whittaker's plot is a lat/long combo 
unique_counts_coordinate <- table(ATdata$Lat_provenance) 
print(nrow(unique_counts_coordinate)) #but you only have 124 unique latitude points 

