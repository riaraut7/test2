# Testing Thomas vs June models 

#let's first load libraries and data 
library(broom)
library(tidyverse)
library(rTPC) 
library(nls2)
library(ggplot2)
library(nls.multstart)
library(patchwork)

#### Filtering your data #### 
AT_data <- read.csv('passed_filter_1.csv', stringsAsFactors = TRUE)
AT_data <- AT_data %>% 
  filter(is.na(failure_status)) %>% 
  filter(Habitat == 'Terrestrial') %>% 
  filter(!is.na(A_umol_m2_s)) #for here, we can just filter out the odd units 


#SIDE NOTE: just counting how many curves there are left # 
first_rows <- AT_data %>%
  group_by(curveID) %>%
  slice(1)
nrow(first_rows)

#### Setting up the June equation + errorcatch functions #### 
#1:The June model:  
June_2004 <- function (temp, aopt, topt, omega) 
{
  exponent <- (-1)*(((temp-topt)/omega)^2)
  p_t <- aopt*exp(exponent)
  return(p_t)
}

#errorcatch functions 
get_breadth_errorcatch <- function(x) {
  result <- tryCatch(
    {
      get_breadth(x)
    },
    error = function(e) {
      message("An error occurred: ", conditionMessage(e))
      return(NA)  # Store NA when an error occurs
    }
  )
  return(result)
}

get_omega_errorcatch <- function(x) {
  result <- tryCatch(
    {
      parameters <- coef(x) 
      parameters['omega']
    },
    error = function(e) {
      message("An error occurred: ", conditionMessage(e))
      return(NA)  # Store NA when an error occurs
    }
  )
  return(result)
}

get_pure_omega <- function(x) {
  this_is_omega <- coef(x)['omega']
  return (this_is_omega)
}

get_pure_topt <- function(x) {
  this_is_topt <- coef(x)['topt']
  return (this_is_topt)
}

get_pure_aopt_june <- function(x) {
  this_is_aopt <- coef(x)['aopt']
  return (this_is_aopt)
}

get_c_from_thomas <- function(x) {
  this_is_c <- coef(x)['c']
  return (this_is_c)
}


#### Fitting both models #### 
cvID <- unique(AT_data$curveID) 
data_out = data.frame() 
for (curve in cvID) {
  per_curve <- AT_data[AT_data$curveID == curve,] 
  print(curve) #prints the curve number 
  
  #try fitting thomas 
  mod = 'thomas_2012'
  start_vals <- get_start_vals(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'thomas_2012')
  low_lims <- get_lower_lims(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'thomas_2012')
  upper_lims <- get_upper_lims(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'thomas_2012')
  fit2 <- nls_multstart(A_umol_m2_s~thomas_2012(temp = Tleaf_C, a, b, c, tref),
                        data = per_curve,
                        iter = c(4,4,4,4),
                        start_lower = start_vals - 1,
                        start_upper = start_vals + 2,
                        lower = low_lims,
                        upper = upper_lims,
                        supp_errors = 'Y', 
                        convergence_count = FALSE)
  values_thomas <- calc_params(fit2)
  AIC_thomas <- AIC(fit2)
  
  #try fitting june 
  fit1 <- nls_multstart(A_umol_m2_s ~ June_2004(temp = Tleaf_C, aopt, topt, omega),
                        data = per_curve,
                        iter = 50,
                        start_lower = c(aopt = 0, topt = 0, omega = 2),
                        start_upper = c(aopt = 20, topt =50, omega =10),
                        lower = c(aopt = 0, topt = 0, omega = 0), 
                        upper = c(aopt = 90, topt = 100, omega = 100),
                        supp_errors = 'Y')
  values_june <- calc_params(fit1) 
  AIC_june <- AIC(fit1)

  #getting a list of variables per curve to combine into a dataset 

  each_AIC_june <- AIC_june[[1]]
  each_AIC_thomas <- AIC_thomas[[1]]
  taxon_list <- per_curve$Taxon
  each_taxon <- taxon_list[[1]]
  units_list <- per_curve$Units_other
  each_unit <- units_list[[1]]
  curve_id_list <- per_curve$curveID
  each_id <- curve_id_list[[1]]
  lat_list <- per_curve$Lat_provenance
  each_lat <- lat_list[[1]]
  long_list <- per_curve$Lon_provenance
  each_long <- long_list[[1]]
  data_out = bind_rows (data_out, data.frame(Curve_ID  = each_id, Taxon = each_taxon, 
                                             breadth_thomas = get_c_from_thomas(fit2),
                                             amax_thomas = get_rmax(fit2), 
                                             topt_thomas = get_topt(fit2), 
                                             breadth_june = get_omega_errorcatch(fit1), 
                                             amax_june = get_pure_aopt_june(fit1), 
                                             topt_june = get_pure_topt(fit1),
                                             Lat_prov = each_lat, 
                                             Long_prov = each_long, 
                                             AIC_june = AIC_june, 
                                             AIC_thomas = AIC_thomas))
}

t_J_modelled <- data_out
write.csv(t_J_modelled, file = 'thomas_june_comparison_dataset.csv')


#### START HERE Now let's compare the two functions visually #### 

t_J_modelled <- read.csv('thomas_june_comparison_dataset.csv', stringsAsFactors = TRUE)

breadth_compare <- ggplot(t_J_modelled, aes(x = breadth_june, y = breadth_thomas)) + 
  geom_point(alpha = 0.4, size = 3, color = '#CCCC33') + 
#  geom_abline(slope = 1, color = "#009933", linewidth = 1) + 
  geom_smooth(method = 'lm', se = TRUE, color = '#CC9900', alpha = 0.5) + theme_classic() + 
  labs(x = expression(paste('Breadth from the June model (°C)')), 
       y = expression(paste('Breadth from then Thomas model (°C)'))) +  
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14)) + 
  geom_text(aes(x=5, y= 120, label = 'A'), size = 7)
breadth_compare


t.test(t_J_modelled$breadth_june, t_J_modelled$breadth_thomas)
breadth_model <- lm(breadth_june ~ breadth_thomas, data = t_J_modelled)
breadth_model
summary(breadth_model) #curve is ~0.36 

aopt_compare <- ggplot(t_J_modelled, aes(x = amax_june, y = amax_thomas)) + 
  geom_point(alpha = 0.4, size = 3, color = '#CCCC33') + 
  #  geom_abline(slope = 1, color = "#009933", linewidth = 1) + 
  geom_smooth(method = 'lm', se = TRUE, color = '#CC9900', alpha = 0.5) + theme_classic() + 
  labs(x = expression(paste('A'[opt],  ' from the June model (µmol m'^-2, 's'^-1, ')')), 
       y = expression(paste('A'[opt],  ' from the Thomas model (µmol m'^-2, 's'^-1, ')'))) +  
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14)) + 
  geom_text(aes(x=2, y= 59, label = 'B'), size = 7)

aopt_compare #seems pretty 1:1 
t.test(t_J_modelled$amax_june, t_J_modelled$amax_thomas)
aopt_model <- lm(amax_june ~ amax_thomas, data = t_J_modelled)
summary(aopt_model) #curve is ~0.36 




topt_compare <- ggplot(t_J_modelled, aes(x = topt_june, y = topt_thomas)) + 
  geom_point(alpha = 0.4, size = 3, color = '#CCCC33') + 
#  geom_abline(slope = 1, color = "#009933", linewidth = 1) + 
  geom_smooth(method = 'lm', se = TRUE, color = '#CC9900', alpha = 0.5) + theme_classic() + 
  labs(x = expression(paste('T'[opt], ' from the June model (°C)')), 
       y = expression(paste('T'[opt], ' from the Thomas model (°C)'))) +  
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14)) + 
  geom_text(aes(x=9, y= 50, label = 'C'), size = 7)

topt_compare #seems pretty 1:1 
t.test(t_J_modelled$topt_june, t_J_modelled$topt_thomas)
topt_model <- lm(topt_june ~ topt_thomas, data = t_J_modelled)
summary(topt_model) 

#### let's print it all togetherrrr #### 
combined_plots <- breadth_compare / aopt_compare / topt_compare + plot_layout(ncol = 1)
combined_plots
#print this 12x6, pdf

#### Let's compare AIC values #### 
AIC_avjun <- mean(t_J_modelled$AIC_june)
AIC_avthom <- mean(t_J_modelled$AIC_thomas)
print(c(AIC_avjun, AIC_avthom))
# 59.28141 60.41899 


#### Do they show different results for hypothesis testing? ####

#Joat
Joat_june <- ggplot(t_J_modelled, aes(x = breadth_june, y = amax_june)) + 
  geom_point() + 
  geom_smooth(method = 'lm')
Joat_june_model <- lm(amax_june ~ breadth_june, t_J_modelled)
summary(Joat_june_model)

Joat_thomas <- ggplot(t_J_modelled, aes(x = breadth_thomas, y = amax_thomas)) + 
  geom_point() + 
  geom_smooth(method = 'lm') 
Joat_thomas
Joat_thomas_model <- lm(amax_thomas ~ breadth_june, t_J_modelled) 
summary(Joat_thomas_model)

#both give the same trend, both significant. Slightly different slopes 



