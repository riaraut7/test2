
#libraries 
library(dplyr)
library(tidyverse)

#let's load all your curves that passed filter 1 
passed_filter_1 <- read.csv('passed_filter_1.csv', stringsAsFactors = TRUE)

##### Passed filter 1 counts ####
#getting the first row for each curve 
first_rows <- passed_filter_1 %>%
  group_by(curveID) %>%
  slice(1)
nrow(first_rows) #859 

#let's replace all NA values with 'No error' 
first_rows <- first_rows %>%
  mutate(failure_status = if_else(is.na(failure_status), 'no error', failure_status))

#curves fewer than 5 points 
passed_number_points <- first_rows %>% 
  filter (failure_status != 'Rejected: fewer than 5 data points')  %>% 
  filter(failure_status != 'Rejected: fewer than two points before peak') %>% 
  filter(failure_status != 'Rejected: fewer than two points after peak') %>% 
  filter (failure_status != 'Odd units. Rejected: fewer than 5 data points') %>% 
  filter(failure_status != 'Odd units. Rejected: fewer than two points before peak') %>% 
  filter(failure_status != 'Odd units. Rejected: fewer than two points after peak') 
  
nrow(passed_number_points)


#ascending region too small 
passed_asc_desc <- passed_number_points %>% 
  filter (failure_status != 'Rejected: ascending region too small') %>% 
  filter (failure_status != 'Rejected: descending region too small')
nrow(passed_asc_desc)


#non-terrestrial count 
passed_terrestrial <- passed_asc_desc %>%
  filter(Habitat == "Terrestrial")
nrow(passed_terrestrial)

#in situ count 
passed_diurnal <- passed_terrestrial %>%
  filter(LI_6400_method != "In situ") %>%
  filter(General_method != "in situ")
nrow(passed_diurnal)



##### Passed filter 2 counts ####
passed_filter_2 <- read.csv('final_curve_parameter_filter_2_info_11.3.24.csv', stringsAsFactors = TRUE)
nrow(passed_filter_2)
passed_filter_2 <- passed_filter_2 %>%
  mutate(june_failure = if_else(is.na(june_failure), 'no error', june_failure))


passed_r2_test <- passed_filter_2 %>% 
  filter(june_failure != 'Warning: bad june fit.') 
nrow(passed_r2_test)

passed_SE_test <- passed_r2_test %>% 
  filter(june_failure != 'Warning: Topt has large SE.') %>% 
  filter(june_failure != 'Warning: Aopt has large SE.') %>% 
  filter(june_failure != 'Warning: breadth has large SE.') 
nrow(passed_SE_test)


#final remaining species
final_species_count <- passed_SE_test %>%
  group_by(Taxon) %>% 
  slice(1)
nrow(final_species_count) #64 final species 

##### Final curves latitude distribution #####
final_curves_noNA <- final_species_count %>% 
  filter(!is.na(Lat_provenance))

max(final_curves_noNA$Lat_provenance)
min(final_curves_noNA$Lat_provenance)

