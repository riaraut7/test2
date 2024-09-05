#downloading the right libraries first 
library(broom) 
#git trial 
library(tidyverse)
library(ggplot2)
library (nls2)
library(phytools)
library(patchwork)
library(geiger)

#load your data 
starting_dataset <- read.csv('final_curve_parameter_filter_2_info_11.3.24.csv', stringsAsFactors = TRUE) 

#### Cleaning the data and averaging by species #### 
#clean up your data 
cleaned_data <- starting_dataset %>% 
  filter(is.na(june_failure))
#and also let's turn all latitude values to absolute values: 
cleaned_data$abs_lat <- abs(cleaned_data$Lat_provenance)

#now let's make a species_averaged_data that has amax, topt, and breadth averaged per species 
species_averaged_data <- cleaned_data %>% 
  group_by(Taxon) %>% 
  summarise(average_breadth = mean(breadth_june), 
            average_topt = mean(topt_june), 
            average_amax = mean(amax_june), 
            Lat_provenance = mean(abs_lat)) #this is unideal, averaging the latitude. But for now, we'll keep it here, 
                                        #having acknolwedged its shortcomings 

#let's check how many rows (curves) don't have latitude values 
number_missing_curves <- sum(is.na(cleaned_data$abs_lat))
number_missing_curves #34 
nrow(cleaned_data) #276 
number_missing_species <- sum(is.na(species_averaged_data$Lat_provenance))
number_missing_species #14
nrow(species_averaged_data) #64 



#### JoAT checking assumptions of a linear model #### 
#let's see if the assumptions of a linnear model are met 
model <- lm(breadth_june ~ amax_june, data = cleaned_data)
residuals <- residuals(model)

histogram <- ggplot(data.frame(residuals), aes(x = residuals)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black", aes(y = ..density..), alpha = 0.5) +
  geom_density(color = "red") +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Density")
histogram

qqplot <- ggplot(data.frame(residuals), aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plot of Residuals")
qqplot


#### JoAT per curve ID  #### 

JoAT_per_curve_id <- ggplot (cleaned_data, aes(x = log(amax_june), y = breadth_june)) + 
  geom_point(alpha = 0.4, size = 3, color = 'yellowgreen') + 
  geom_smooth(method = 'lm', se = TRUE, color = 'seagreen', alpha = 0.5) + theme_classic() + 
  labs(x = expression(paste('Log-Transformed Optimal Assimilation Rate (µmol m'^-2, 's'^-1, ')')), 
       y = expression(paste('AT curve breadth (°C)'))) +  
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14)) + 
  geom_text(aes(x=0.5, y= 36, label = 'A'), size = 7)
print (JoAT_per_curve_id)

#getting values from your fitted linear regression 
model <- lm(breadth_june ~ log(amax_june), data = cleaned_data) #should always be y~x
summary(model)

#doing an actual correlation test 
# correlation_test <- cor.test(log(cleaned_data$amax_june), cleaned_data$breadth_june)
# print(correlation_test)

#### JoAT per species  #### 

JoAT_per_species <- ggplot (species_averaged_data, aes(x = log(average_amax), y = average_breadth)) + 
  geom_point(alpha = 0.4, size = 3, color = 'yellowgreen') + 
  geom_smooth(method = 'lm', se = TRUE, color = 'seagreen', alpha = 0.5) + theme_classic() + 
  labs(x = expression(paste('Log-Transformed Optimal Assimilation Rate (µmol m'^-2, 's'^-1, ')')), 
       y = expression(paste('AT curve breadth (°C)'))) + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14)) + 
  geom_text(aes(x=0.5, y= 30, label = 'B'), size = 7)
print (JoAT_per_species)

#correlation = 0.203 (this is the slope/strength of relationship), p = 0.110
model <- lm(average_breadth ~ average_amax, data = species_averaged_data)
summary(model)


# #let's just do a simple pearson's correlation test 
# correlation_test <- cor.test(log(species_averaged_data$average_amax), species_averaged_data$average_breadth)
# print(correlation_test)


#### JoAT PIC #### 
mytree <- read.tree('your_working_tree.tre')
in_right_order <- read.csv('in_right_order_after_s.antarctici.csv', header = T)
#let's make a matrix with PIC-corrected values of each trait of interest 
set.seed(0)
PIC_topt <- pic(in_right_order$av_topt, mytree)
PIC_amax <- pic(in_right_order$av_amax, mytree)
PIC_breadth <- pic(in_right_order$av_breadth, mytree)
PIC_amax_log <- pic(log(in_right_order$av_amax), mytree)
PIC_breadth_log <- pic(log(in_right_order$av_breadth), mytree)

pic_df <- data.frame(PIC_amax, PIC_topt, PIC_breadth, PIC_amax_log, PIC_breadth_log)

#and now let's graph it 
PIC_JoAT <- ggplot (pic_df, aes(x = PIC_amax_log, y = PIC_breadth)) + 
  geom_point(alpha = 0.4, size = 3, color = 'yellowgreen') + 
  geom_smooth(method = 'lm', se = TRUE, color = 'seagreen', alpha = 0.5) + 
  theme_classic() + 
  labs(x = expression(paste('PICs for Optimal Assimilation Rate (µmol m'^-2, 's'^-1, ')')), 
       y = expression(paste('PICs for AT curve breadth (°C)'))) + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14)) + 
  geom_text(aes(x=-1.1, y= 5.9, label = 'C'), size = 7)

print (PIC_JoAT) 

JoAT_lm <- lm(PIC_breadth ~ PIC_amax_log)
summary(JoAT_lm)

# correlation_test <- cor.test(PIC_amax_log, PIC_breadth)
# print(correlation_test) 


#### JoAT main Figure #### 
combined_plots <- JoAT_per_curve_id / JoAT_per_species / PIC_JoAT + plot_layout(ncol = 1)
combined_plots
#print this 12x6, pdf


#### Janzen checking assumptions of a linear model #### 
cleaned_no_lat_na <- cleaned_data %>% 
  filter (!is.na(cleaned_data$abs_lat))

model <- lm(breadth_june ~ abs_lat, data = cleaned_no_lat_na)
residuals <- residuals(model)
histogram <- ggplot(data.frame(residuals), aes(x = residuals)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black", aes(y = ..density..), alpha = 0.5) +
  geom_density(color = "red") +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Density")
histogram
qqplot <- ggplot(data.frame(residuals), aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plot of Residuals")
qqplot


#### Janzen per curve ID ####

cleaned_no_lat_na <- cleaned_data %>% 
  filter (!is.na(cleaned_data$abs_lat))

Janzen_per_curve_id <- ggplot (cleaned_no_lat_na, aes(x = abs_lat, y = breadth_june)) + 
  geom_point(alpha = 0.4, size = 3, color = 'yellowgreen') + 
  geom_smooth(method = 'lm', se = TRUE, color = 'seagreen', alpha = 0.5) + theme_classic() + 
  labs(x = expression(paste('Absolute global latitude (°)')), 
       y = expression(paste('AT curve breadth (°C)'))) + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14)) + 
  geom_text(aes(x= 8.9, y= 39, label = 'A'), size = 7)
print (Janzen_per_curve_id)


Janzens_but_tempt <- ggplot(cleaned_no_lat_na, aes(x= abs_lat, y = topt_june)) + 
  geom_point(alpha = 0.4) + 
  geom_smooth(method = 'lm', se = TRUE)

Janzens_but_tempt

##
model <- lm(breadth_june ~ abs_lat, data = cleaned_no_lat_na)
summary(model)
correlation_test <- cor.test(cleaned_no_lat_na$abs_lat, cleaned_no_lat_na$breadth_june)
print(correlation_test)

#### Janzen per species #### 

species_averaged_janzen <- cleaned_no_lat_na %>% 
  group_by(Taxon) %>% 
  summarise(average_breadth = mean(breadth_june), 
            average_topt = mean(topt_june), 
            average_amax = mean(amax_june), 
            Lat_provenance = mean(abs_lat))


Janzen_per_species <- ggplot (species_averaged_janzen, aes(x = Lat_provenance, y = average_breadth)) + 
  geom_point(alpha = 0.4, size = 3, color = 'yellowgreen') + 
  geom_smooth(method = 'lm', se = TRUE, color = 'seagreen', alpha = 0.5) + theme_classic() +
  labs(x = expression(paste('Absolute global latitude (°)')), 
       y = expression(paste('AT curve breadth (°C)'))) + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14)) + 
  geom_text(aes(x= 8.9, y= 39, label = 'B'), size = 7)

print (Janzen_per_species)
nrow(species_averaged_janzen)
correlation_test <- cor.test(species_averaged_janzen$Lat_provenance, species_averaged_janzen$average_breadth)
print(correlation_test)
model <- lm(average_breadth ~ Lat_provenance, data = species_averaged_janzen)
summary(model)

model$terms


#### Janzen PIC #### 

#first let's filter out the latitude na's 
in_right_order <- in_right_order %>% 
  filter(!is.na(av_lat))
rownames(in_right_order) <- gsub(" ","_",in_right_order$Taxon)

obj<-name.check(mytree,in_right_order)
obj #Sort of!

#take out problematic species from your tree 
species_to_remove <- c("Aegilops_cylindrica", "Brassica_campestris", "Prosopis_juliflora", 'Atriplex_confertifolia',
                       'Helianthus_annuus', 'Solanum_lycopersicum', 'Atriplex_glabriuscula', 'Larix_decidua', 
                       'Solanum_tuberosum', 'Atriplex_hymenelytra', 'Oryza_sativa', 'Triticum_aestivum', 'Atriplex_sabulosa', 
                       'Phoenix_dactylifera')
# Remove the specified species from the tree
pruned_tree <- drop.tip(mytree, species_to_remove)

obj<-name.check(pruned_tree,in_right_order)
obj #seems good! 
mytree=multi2di(pruned_tree)

set.seed(0)
PIC_topt <- pic(in_right_order$av_topt, mytree)
PIC_amax <- pic(in_right_order$av_amax, mytree)
PIC_breadth <- pic(in_right_order$av_breadth, mytree)
PIC_amax_log <- pic(log(in_right_order$av_amax), mytree)
PIC_breadth_log <- pic(log(in_right_order$av_breadth), mytree)
PIC_lat <- pic(in_right_order$av_lat, mytree)

pic_df <- data.frame(PIC_amax, PIC_topt, PIC_breadth, PIC_amax_log, PIC_breadth_log, PIC_lat)

#take out this weird outlier datapoint 
row_to_remove <- 11
pic_df <- pic_df[-row_to_remove, ]
#-23.39269559 for PIC_lat?
nrow(pic_df)

PIC_Janzen <- ggplot (pic_df, aes(x = PIC_lat, y = PIC_breadth)) + 
  geom_point(alpha = 0.4, size = 3, color = 'yellowgreen') + 
  geom_smooth(method = 'lm', se = TRUE, color = 'seagreen', alpha = 0.5) +   theme_classic() + 
  labs(x = expression(paste('PICs for absolute global latitude (°)')), 
       y = expression(paste('PICs for AT curve breadth (°C)'))) + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14)) + 
  geom_text(aes(x= -3.8, y= 4.7, label = 'C'), size = 7)
  

print (PIC_Janzen) 

Janzen_lm <- lm(PIC_breadth ~ PIC_lat)
summary(Janzen_lm)
correlation_test <- cor.test(PIC_lat, PIC_breadth)
print(correlation_test) 

#### Janzens Main figure #### 
combined_janzen <- Janzen_per_curve_id / Janzen_per_species / PIC_Janzen + plot_layout(ncol = 1) 
combined_janzen
#print this at 12x6

#### Janzen vertically logged - just for fun #####
#You tried to fit a logarithmic curve to this, just for fun 
Janzen_per_curve_id_logged <- ggplot (cleaned_no_lat_na, aes(x = abs_lat, y = breadth_june)) + 
  geom_point(alpha = 0.4, size = 3, color = 'yellowgreen') + 
  geom_smooth(method = 'lm', formula = y ~ log(x), se = TRUE, color = 'seagreen', alpha = 0.5) + theme_classic() + 
  labs(x = 'Absolute global latitude (°)', y = expression(paste('AT curve breadth (', Omega, ' )'))) + 
  ggtitle(expression(paste ('Plant Provenance Absolute Latitude(°) vs Curve Breadth (', Omega, ' ) Per Curve'))) + 
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), 
        axis.text = element_text(size = 16))

print (Janzen_per_curve_id_logged)
log_model <- lm(breadth_june ~ log(abs_lat), data = cleaned_data)
summary(log_model)
correlation_test <- cor.test(cleaned_data$abs_lat, log(cleaned_data$breadth_june))
print(correlation_test)







