<<<<<<< HEAD
#New trying things out in rTPC

library(broom)
library(tidyverse)
library(rTPC) 
library(nls.multstart)
#RIA: Note that things works out/give you fewer errors if you have nls.multstart after broom and tidyverse 
#that's becuse you're dumb and have't figured out how to resolve the conflicts yet, but you can do that later 

ATdata <- read.csv("AT_data.csv", stringsAsFactors = TRUE)

#trying to fit and graph a gaussian on each curve 
Taxa <- unique(ATdata$Taxon) 
pdf("S-S-high fitted graphs_part2.pdf")
for (tx in Taxa) { 
  per_taxa <- ATdata[ATdata$Taxon == tx,] 
  ggplot(per_taxa, aes(Tleaf_C, A_umol_m2_s)) +
    geom_point() +
    theme_bw(base_size = 12) + facet_wrap(~Taxon, ncol = 2) + 
    labs(x = 'Temperature (ºC)',
         y = 'Metabolic rate',
         title = 'Respiration across temperatures in this species')
  
  mod = 'sharpeschoolhigh_1981'
  start_vals <- get_start_vals(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'sharpeschoolhigh_1981')
  
  
  fit <- nls_multstart(A_umol_m2_s~sharpeschoolhigh_1981(temp = Tleaf_C, r_tref,e,eh,th, tref = 15),
                       data = per_taxa,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  fit
  calc_params(fit)
  #if you want to round, just have 
  #calc_params(fit) %>% 
  #mutate_all(round, 2) #'2' being to 2 decimals 
  
  
  #And finally, if you wanna predict data: 
  new_data <- data.frame(temp = seq(min(per_taxa$Tleaf_C), max(per_taxa$Tleaf_C), 0.5))
  preds <- augment(fit, newdata2 = new_data)
  
  p <- ggplot(per_taxa, aes(Tleaf_C, A_umol_m2_s)) +
    geom_point() +
    geom_line(aes(x = Tleaf_C, y = .fitted), data = preds, col = 'green') +
    theme_bw(base_size = 12) + facet_wrap(~Taxon, ncol = 2)
    labs(x = 'Temperature (ºC)',
         y = 'Metabolic rate',
         title = 'aaaaah')
  print(p)
  #Ok so - turns out everything works fine until you get to the Q. macrocarpa curve? Idk why tho, wdym Tleaf_C not found? 
  #it's there the same as everything else? 
}
dev.off()

#*******************TRYING TO GET CTMAX, BREADTH, ETC FOR EACH SPECIES ***************************************** 
# a list of columns to have in your new CSV 
Species_ssh <- list() 
curve_breadth_ssh <- list() 
ct_max_ssh <- list() 
ct_min_ssh <- list() 
curve_breadth_g <- list() 
ct_max_g <- list() 
ct_min_g <- list() 
curve_breadth_lac <- list() 
ct_max_lac <- list() 
ct_min_lac <- list() 


#**************sharpe schoolfield high 
for (tx in Taxa) { 
  per_taxa <- ATdata[ATdata$Taxon == tx,] 
  # ggplot(per_taxa, aes(Tleaf_C, A_umol_m2_s)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) + facet_wrap(~Taxon, ncol = 2) + 
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Metabolic rate',
  #        title = 'Respiration across temperatures in this species')
 
  # one_species <- str(tx)  
  # # Species_ssh <- append(Species_ssh, list(one_species))
  # ct_max_ssh <- append (ct_max_ssh, list(values$ctmax))
  #^ I cannot for the love of god understand why this won't work PLEASE let me LIVE 
  
  mod = 'sharpeschoolhigh_1981'
  start_vals <- get_start_vals(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'sharpeschoolhigh_1981')
  
  fit <- nls_multstart(A_umol_m2_s~sharpeschoolhigh_1981(temp = Tleaf_C, r_tref,e,eh,th, tref = 15),
                       data = per_taxa,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  fit
  values <- calc_params(fit)

} 


#**************gaussian 
for (tx in Taxa) { 
  per_taxa <- ATdata[ATdata$Taxon == tx,] 

  mod = 'gaussian_1987'
  start_vals <- get_start_vals(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'gaussian_1987')
  low_lims <- get_lower_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'gaussian_1987')
  upper_lims <- get_upper_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'gaussian_1987')
  
  fit <- nls_multstart(A_umol_m2_s~gaussian_1987(temp = Tleaf_C, rmax, topt, a),
                       data = per_taxa,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  fit
  values <- calc_params(fit)
  
  
} 


#**************lactin2
for (tx in Taxa) { 
  per_taxa <- ATdata[ATdata$Taxon == tx,] 

  mod = 'lactin2_1995'
  start_vals <- get_start_vals(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'lactin2_1995')
  low_lims <- get_lower_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'lactin2_1995')
  upper_lims <- get_upper_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'lactin2_1995')
  
  fit <- nls_multstart(A_umol_m2_s~lactin2_1995(temp = Tleaf_C, a, b, tmax, delta_t),
                       data = per_taxa,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  fit
  values <- calc_params(fit)
  
} 




# Create a sample data frame
sharpe_high_data <- data.frame(Species = Species_ssh)
file_path <- "sharpe_high_processed_valued_part_3.csv"
write.csv(sharpe_high_data, file = file_path, row.names = FALSE)
#cat("CSV file has been created:", file_path, "\n")

#*************************************************************************************
#trying to fit all curves 
#*********************************************************************************** 
#write function to label ggplot2 panels
label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

#your data (chlorella_tpc in this case) is already loaded and shown 

# fit every model formulation in rTPC
d_fits <- nest(d2, data = c(Tleaf_C, A_umol_m2_s)) %>%
  mutate(beta = map(data, ~nls_multstart(A_umol_m2_s~beta_2012(temp = Tleaf_C, a, b, c, d, e),
                                         data = .x,
                                         iter = c(6,6,6,6,6),
                                         start_lower = get_start_vals(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'beta_2012') - 10,
                                         start_upper = get_start_vals(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'beta_2012') + 10,
                                         lower = get_lower_lims(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'beta_2012'),
                                         upper = get_upper_lims(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'beta_2012'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)),
         weibull = map(data, ~nls_multstart(rate~weibull_1995(temp = Tleaf_C, a,topt,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'weibull_1995') - 10,
                                            start_upper = get_start_vals(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'weibull_1995') + 10,
                                            lower = get_lower_lims(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'weibull_1995'),
                                            upper = get_upper_lims(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'weibull_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)))


# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:weibull)

# get parameters using tidy
params <- d_stack %>%
  mutate(., est = map(fit, tidy)) %>%
  select(-fit) %>%
  unnest(est)

# get predictions using augment
newdata <- tibble(temp = seq(min(d2$Tleaf_C), max(d2$Tleaf_C), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# plot
ggplot(d_preds, aes(Tleaf_C, A_umol_m2_s)) +
  geom_point(aes(Tleaf_C, A_umol_m2_s), d2) +
  geom_line(aes(Tleaf_C, .fitted), col = 'blue') +
  facet_wrap(~model_name, labeller = labeller(model_name = label_facets_num), scales = 'free', ncol = 5) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Fits of every model available in rTPC (yes lrf, no deutsch') +
  geom_hline(aes(yintercept = 0), linetype = 2)


#
ggplot(d2, aes(Tleaf_C, A_umol_m2_s)) +
  geom_point() +
  geom_line(aes(Tleaf_C, .fitted), preds, col = 'green') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures in A. tonduzii')

##Troubleshooting some things 
str(per_taxa)
str(preds)
=======
#New trying things out in rTPC

library(broom)
library(tidyverse)
library(rTPC) 
library(nls.multstart)
#RIA: Note that things works out/give you fewer errors if you have nls.multstart after broom and tidyverse 
#that's becuse you're dumb and have't figured out how to resolve the conflicts yet, but you can do that later 

ATdata <- read.csv("AT_data.csv", stringsAsFactors = TRUE)

#trying to fit and graph a gaussian on each curve 
Taxa <- unique(ATdata$Taxon) 
pdf("S-S-high fitted graphs_part2.pdf")
for (tx in Taxa) { 
  per_taxa <- ATdata[ATdata$Taxon == tx,] 
  ggplot(per_taxa, aes(Tleaf_C, A_umol_m2_s)) +
    geom_point() +
    theme_bw(base_size = 12) + facet_wrap(~Taxon, ncol = 2) + 
    labs(x = 'Temperature (ºC)',
         y = 'Metabolic rate',
         title = 'Respiration across temperatures in this species')
  
  mod = 'sharpeschoolhigh_1981'
  start_vals <- get_start_vals(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'sharpeschoolhigh_1981')
  
  
  fit <- nls_multstart(A_umol_m2_s~sharpeschoolhigh_1981(temp = Tleaf_C, r_tref,e,eh,th, tref = 15),
                       data = per_taxa,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  fit
  calc_params(fit)
  #if you want to round, just have 
  #calc_params(fit) %>% 
  #mutate_all(round, 2) #'2' being to 2 decimals 
  
  
  #And finally, if you wanna predict data: 
  new_data <- data.frame(temp = seq(min(per_taxa$Tleaf_C), max(per_taxa$Tleaf_C), 0.5))
  preds <- augment(fit, newdata2 = new_data)
  
  p <- ggplot(per_taxa, aes(Tleaf_C, A_umol_m2_s)) +
    geom_point() +
    geom_line(aes(x = Tleaf_C, y = .fitted), data = preds, col = 'green') +
    theme_bw(base_size = 12) + facet_wrap(~Taxon, ncol = 2)
    labs(x = 'Temperature (ºC)',
         y = 'Metabolic rate',
         title = 'aaaaah')
  print(p)
  #Ok so - turns out everything works fine until you get to the Q. macrocarpa curve? Idk why tho, wdym Tleaf_C not found? 
  #it's there the same as everything else? 
}
dev.off()

#*******************TRYING TO GET CTMAX, BREADTH, ETC FOR EACH SPECIES ***************************************** 
# a list of columns to have in your new CSV 
Species_ssh <- list() 
curve_breadth_ssh <- list() 
ct_max_ssh <- list() 
ct_min_ssh <- list() 
curve_breadth_g <- list() 
ct_max_g <- list() 
ct_min_g <- list() 
curve_breadth_lac <- list() 
ct_max_lac <- list() 
ct_min_lac <- list() 


#**************sharpe schoolfield high 
for (tx in Taxa) { 
  per_taxa <- ATdata[ATdata$Taxon == tx,] 
  # ggplot(per_taxa, aes(Tleaf_C, A_umol_m2_s)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) + facet_wrap(~Taxon, ncol = 2) + 
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Metabolic rate',
  #        title = 'Respiration across temperatures in this species')
 
  # one_species <- str(tx)  
  # # Species_ssh <- append(Species_ssh, list(one_species))
  # ct_max_ssh <- append (ct_max_ssh, list(values$ctmax))
  #^ I cannot for the love of god understand why this won't work PLEASE let me LIVE 
  
  mod = 'sharpeschoolhigh_1981'
  start_vals <- get_start_vals(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'sharpeschoolhigh_1981')
  
  fit <- nls_multstart(A_umol_m2_s~sharpeschoolhigh_1981(temp = Tleaf_C, r_tref,e,eh,th, tref = 15),
                       data = per_taxa,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  fit
  values <- calc_params(fit)

} 


#**************gaussian 
for (tx in Taxa) { 
  per_taxa <- ATdata[ATdata$Taxon == tx,] 

  mod = 'gaussian_1987'
  start_vals <- get_start_vals(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'gaussian_1987')
  low_lims <- get_lower_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'gaussian_1987')
  upper_lims <- get_upper_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'gaussian_1987')
  
  fit <- nls_multstart(A_umol_m2_s~gaussian_1987(temp = Tleaf_C, rmax, topt, a),
                       data = per_taxa,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  fit
  values <- calc_params(fit)
  
  
} 


#**************lactin2
for (tx in Taxa) { 
  per_taxa <- ATdata[ATdata$Taxon == tx,] 

  mod = 'lactin2_1995'
  start_vals <- get_start_vals(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'lactin2_1995')
  low_lims <- get_lower_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'lactin2_1995')
  upper_lims <- get_upper_lims(per_taxa$Tleaf_C, per_taxa$A_umol_m2_s, model_name = 'lactin2_1995')
  
  fit <- nls_multstart(A_umol_m2_s~lactin2_1995(temp = Tleaf_C, a, b, tmax, delta_t),
                       data = per_taxa,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  fit
  values <- calc_params(fit)
  
} 




# Create a sample data frame
sharpe_high_data <- data.frame(Species = Species_ssh)
file_path <- "sharpe_high_processed_valued_part_3.csv"
write.csv(sharpe_high_data, file = file_path, row.names = FALSE)
#cat("CSV file has been created:", file_path, "\n")

#*************************************************************************************
#trying to fit all curves 
#*********************************************************************************** 
#write function to label ggplot2 panels
label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

#your data (chlorella_tpc in this case) is already loaded and shown 

# fit every model formulation in rTPC
d_fits <- nest(d2, data = c(Tleaf_C, A_umol_m2_s)) %>%
  mutate(beta = map(data, ~nls_multstart(A_umol_m2_s~beta_2012(temp = Tleaf_C, a, b, c, d, e),
                                         data = .x,
                                         iter = c(6,6,6,6,6),
                                         start_lower = get_start_vals(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'beta_2012') - 10,
                                         start_upper = get_start_vals(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'beta_2012') + 10,
                                         lower = get_lower_lims(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'beta_2012'),
                                         upper = get_upper_lims(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'beta_2012'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)),
         weibull = map(data, ~nls_multstart(rate~weibull_1995(temp = Tleaf_C, a,topt,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'weibull_1995') - 10,
                                            start_upper = get_start_vals(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'weibull_1995') + 10,
                                            lower = get_lower_lims(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'weibull_1995'),
                                            upper = get_upper_lims(.x$Tleaf_C, .x$A_umol_m2_s, model_name = 'weibull_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)))


# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:weibull)

# get parameters using tidy
params <- d_stack %>%
  mutate(., est = map(fit, tidy)) %>%
  select(-fit) %>%
  unnest(est)

# get predictions using augment
newdata <- tibble(temp = seq(min(d2$Tleaf_C), max(d2$Tleaf_C), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# plot
ggplot(d_preds, aes(Tleaf_C, A_umol_m2_s)) +
  geom_point(aes(Tleaf_C, A_umol_m2_s), d2) +
  geom_line(aes(Tleaf_C, .fitted), col = 'blue') +
  facet_wrap(~model_name, labeller = labeller(model_name = label_facets_num), scales = 'free', ncol = 5) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Fits of every model available in rTPC (yes lrf, no deutsch') +
  geom_hline(aes(yintercept = 0), linetype = 2)


#
ggplot(d2, aes(Tleaf_C, A_umol_m2_s)) +
  geom_point() +
  geom_line(aes(Tleaf_C, .fitted), preds, col = 'green') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures in A. tonduzii')

##Troubleshooting some things 
str(per_taxa)
str(preds)
>>>>>>> 9c95ca3e6eea13d6f9f205838021ff0f9f6845bd
