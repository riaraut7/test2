<<<<<<< HEAD

#first load the important libraries (in the right order)
library(broom)
library(tidyverse)
library(rTPC) 
library(nls.multstart)

#introducing: The June et al. 2004 function!!! Good for getting breadths and Atmax's >:)) 
June_2004 <- function (temp, topt, omega, aopt) 
{
  p_t <- aopt*exp(-1*(temp-topt)/omega)^2 
  return(p_t)
}


#load your data and separate it per curve 
ATdata <- read.csv("AT_data.csv", stringsAsFactors = TRUE)
cvID <- unique(ATdata$curveID) 

#loading mini dataset for practice runs 
ATdata_mini <- read.csv('Mini_AT_data.csv', stringsAsFactors = TRUE)
cvID_mini <- unique(ATdata_mini$curveID)


#trying to get NA values and not just an error smh smh when I try to export values into a dataframe 
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


#Introduce the pdf and dataframe that you will be inputting your data into 
#pdf("all_curves_gaussian.pdf")
data_out = data.frame()

for (curve in cvID_mini) { 
  per_curve <- ATdata_mini[ATdata_mini$curveID == curve,] 
  print(curve) 
  Aumol_per_curve <- per_curve$A_umol_m2_s 
  #first you gotta check if your Assimilation values are in the right columns 
  if (is.na(Aumol_per_curve[[1]])) {
    #just fitting the model: 
    mod = 'gaussian_1987'
    start_vals <- get_start_vals(per_curve$Tleaf_C, per_curve$A_other_units, model_name = 'gaussian_1987')
    low_lims <- get_lower_lims(per_curve$Tleaf_C, per_curve$A_other_units, model_name = 'gaussian_1987')
    upper_lims <- get_upper_lims(per_curve$Tleaf_C, per_curve$A_other_units, model_name = 'gaussian_1987')
    
    fit <- nls_multstart(A_other_units~gaussian_1987(temp = Tleaf_C,rmax, topt,a),
                               data = per_curve,
                               iter = c(4,4,4),
                               start_lower = start_vals - 10,
                               start_upper = start_vals + 10,
                               lower = low_lims,
                               upper = upper_lims,
                               supp_errors = 'Y')
    
    fit
    values_other <- calc_params(fit)
    
    #if the model doesn't fit the data, just generate a graph with your plotted data. If it does,
    #then generate a graph with the fitted model. 
    
    if (is.na(values_other$e)) {
      #unmodelled_curves = bind_rows(unmodelled_curves,data.frame(species = per_curve$Taxon))
      p1 <- ggplot(per_curve, aes(Tleaf_C, A_other_units)) +
        geom_point() +
        theme_bw(base_size = 12) + facet_wrap(~curveID, ncol = 2) + 
        labs(x = 'Temperature (ºC)',
             y = 'Metabolic rate',
             title = 'Respiration across temperatures, unique units')
      print(p1)  
      #   print ('model didnt fit!')
      
    } else {
      #And finally, predicted a model and placing that over your raw data: 
      new_data <- data.frame(temp = seq(min(per_curve$Tleaf_C), max(per_curve$Tleaf_C), 0.5))
      preds <- augment(fit, newdata2 = new_data)
    
      p <- ggplot(per_curve, aes(Tleaf_C, A_other_units)) +
        geom_point() +
        geom_line(aes(x = Tleaf_C, y = .fitted), data = preds, col = 'darkgreen') +
        theme_bw(base_size = 12) + facet_wrap(~curveID, ncol = 2)
      labs(x = 'Temperature (ºC)',
           y = 'Metabolic rate',
           title = 'Anet across temperatures, unique units')
      print(p) 
    }
    
    #cat(curve, ' is other units') 
    
    #   next
  } else {
    #just fitting the model: 
    mod = 'gaussian_1987'
    start_vals <- get_start_vals(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'gaussian_1987')
    low_lims <- get_lower_lims(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'gaussian_1987')
    upper_lims <- get_upper_lims(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'gaussian_1987')
    
    fit <- nls_multstart(A_umol_m2_s~gaussian_1987(temp = Tleaf_C,rmax, topt,a),
                         data = per_curve,
                         iter = c(4,4,4),
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y')
    
    fit
    values <- calc_params(fit)
    
    #if the model doesn't fit the data, just generate a graph with your plotted data. If it does,
    #then generate a graph with the fitted model. 
    
    if (is.na(values$e)) {
      #unmodelled_curves = bind_rows(unmodelled_curves,data.frame(species = per_curve$Taxon))
      p1 <- ggplot(per_curve, aes(Tleaf_C, A_umol_m2_s)) +
        geom_point() +
        theme_bw(base_size = 12) + facet_wrap(~curveID, ncol = 2) + 
        labs(x = 'Temperature (ºC)',
             y = 'Metabolic rate',
             title = 'Respiration across temperatures')
      print(p1)  
      #   print ('model didnt fit!')
      
    } else {
      #And finally, predicted a model and placing that over your raw data: 
      new_data <- data.frame(temp = seq(min(per_curve$Tleaf_C), max(per_curve$Tleaf_C), 0.5))
      preds <- augment(fit, newdata2 = new_data)
      
      p <- ggplot(per_curve, aes(Tleaf_C, A_umol_m2_s)) +
        geom_point() +
        geom_line(aes(x = Tleaf_C, y = .fitted), data = preds, col = 'darkgreen') +
        theme_bw(base_size = 12) + facet_wrap(~curveID, ncol = 2)
      labs(x = 'Temperature (ºC)',
           y = 'Metabolic rate',
           title = 'Anet across temperatures')
      print(p) 
      
  data_out = bind_rows(data_out,data.frame(breadth = get_breadth_errorcatch(fit), 
                                           amax_1 = get_rmax(fit)))
      
    }
  }
  
}

data_out
dev.off()

gaussian_data <- data_out
file_path <- "mini_gaus_ver_3.csv"
write.csv(gaussian_data, file = file_path, row.names = FALSE)






# Test the function with different inputs
input_values <- c(5, 0)

for (input in input_values) {
  result <- calculate_result(input)
  print(paste("Result for input", input, "is:", result))
}



#trying to design a function for June_2004 
view(pawar_2018)

#Example 1: 
fit_linear_curve <- function(x, y) {
  # Create a linear model
  model <- lm(y ~ x)
  
  # Return the linear model
  return(model)
}

#Example 2: 
one_point = function(Asat, Ci, Tleaf) {
  
  Tleaf = Tleaf+273.15
  
  Gstar = Gstar_fn(Tleaf)
  Km = Km_fn(Tleaf)
  rho = rho_fn(Tleaf)
  
  Vcmax = Asat / ( (Ci-Gstar)/(Ci+Km) - 0.015*rho )
  
  return(Vcmax)
}

June_2004 <- function(leaf_temp, assimilation) {
  
  A_at_temp = Amax* e**(-((tleaf-topt)/omega)**2)
  return(final)
  
}




#Trying to get ctmax, min, breadth, etc from fitted models - 
#okay so - my plan was to get those values for each species (hence the for loop a bit further down) 
#and then add those values to a list within the loop 

data_out_2 = data.frame()

#generating ctmax, min, etc from a couple models !
for (curve in cvID_mini) { 
  #preliminary, seeing how your data looks in a graph: 
  per_curve <- ATdata_mini[ATdata_mini$curveID == curve,] 
  print(curve)
  
  #quadratic_2008
  mod = 'quadratic_2008'
  start_vals <- get_start_vals(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'quadratic_2008')
  low_lims <- get_lower_lims(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'quadratic_2008')
  upper_lims <- get_upper_lims(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'quadratic_2008')
  
  fit1 <- nls_multstart(A_umol_m2_s~quadratic_2008(temp = Tleaf_C, a, b, c),
                        data = per_curve,
                        iter = c(4,4,4),
                        start_lower = start_vals - 10,
                        start_upper = start_vals + 10,
                        lower = low_lims,
                        upper = upper_lims,
                        supp_errors = 'Y')
  
  #briere2_1999
  mod = 'briere2_1999'
  start_vals <- get_start_vals(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'briere2_1999')
  low_lims <- get_lower_lims(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'briere2_1999')
  upper_lims <- get_upper_lims(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'briere2_1999')
  
  fit2 <- nls_multstart(A_umol_m2_s~briere2_1999(temp = Tleaf_C, tmin, tmax, a, b),
                        data = per_curve,
                        iter = c(4,4,4,4),
                        start_lower = start_vals - 10,
                        start_upper = start_vals + 10,
                        lower = low_lims,
                        upper = upper_lims,
                        supp_errors = 'Y') 
  
  #Ratkowski 
  mod = 'ratkowsky_1983'
  start_vals <- get_start_vals(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'ratkowsky_1983')
  low_lims <- get_lower_lims(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'ratkowsky_1983')
  upper_lims <- get_upper_lims(per_curve$Tleaf_C, per_curve$A_umol_m2_s, model_name = 'ratkowsky_1983')
  
  fit3 <- nls_multstart(A_umol_m2_s~ratkowsky_1983(temp = Tleaf_C, tmin, tmax, a, b),
                        data = per_curve,
                        iter = c(4,4,4,4),
                        start_lower = start_vals - 10,
                        start_upper = start_vals + 10,
                        lower = low_lims,
                        upper = upper_lims,
                        supp_errors = 'Y') 
  
  
  data_out_2 = bind_rows(data_out_2,
                         data.frame(ctmin_quad = get_ctmin(fit1), ctmax_quad = get_ctmax(fit1), 
                                    ctmin_brie = get_ctmin(fit2), ctmax_brie = get_ctmax(fit2), 
                                    ctmin_ratk = get_ctmin(fit3), ctmax_ratk = get_ctmax(fit3)))  
  
  
} 
dev.off()
data_out_2 

ctmax_and_mins <- data_out_2
file_path <- "ctmax_and_mins_file_2.csv"
write.csv(ctmax_and_mins, file = file_path, row.names = FALSE)

warnings()
get_model_names()
=======
#New trying things out in rTPC

library(broom)
library(tidyverse)
library(rTPC) 
library(nls.multstart)
#RIA: Note that things works out/give you fewer errors if you have nls.multstart after broom and tidyverse 
#that's becuse you're dumb and haven't figured out how to resolve the conflicts yet, but you can do that later 

ATdata <- read.csv("AT_data.csv", stringsAsFactors = TRUE)

#trying to fit and graph a Sharpeschoolfield on each curve, and get a phf from that
Taxa <- unique(ATdata$Taxon) 
pdf("S-S-high fitted graphs_part2.pdf")
for (tx in Taxa) { 
  #preliminary, seeing how your data looks in a graph: 
  per_taxa <- ATdata[ATdata$Taxon == tx,] 
  ggplot(per_taxa, aes(Tleaf_C, A_umol_m2_s)) +
    geom_point() +
    theme_bw(base_size = 12) + facet_wrap(~Taxon, ncol = 2) + 
    labs(x = 'Temperature (ºC)',
         y = 'Metabolic rate',
         title = 'Respiration across temperatures')
  
  #just fitting the model: 
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
  
  
  #And finally, predicted a model and placing that over your raw data: 
  new_data <- data.frame(temp = seq(min(per_taxa$Tleaf_C), max(per_taxa$Tleaf_C), 0.5))
  preds <- augment(fit, newdata2 = new_data)
  
  p <- ggplot(per_taxa, aes(Tleaf_C, A_umol_m2_s)) +
    geom_point() +
    geom_line(aes(x = Tleaf_C, y = .fitted), data = preds, col = 'green') +
    theme_bw(base_size = 12) + facet_wrap(~Taxon, ncol = 2)
    labs(x = 'Temperature (ºC)',
         y = 'Metabolic rate',
         title = 'Respiration across temperatures')
  print(p)
  #Ok so - turns out everything works fine until you get to the Q. macrocarpa curve? Idk why tho,
  #Error is Tleaf_C not found? 
  #but species has the same columns as everything else? 
}
dev.off()





#Trying to get ctmax, min, breadth, etc from fitted models - 
#okay so - my plan was to get those values for each species (hence the for loop a bit further down) 
#and then add thoes values to a list within the loop 

# a list of columns to have in your new CSV 
Species_list_for_csv <- list() 
#sharpeschoolfield high 
curve_breadth_ssh <- list() 
ct_max_ssh <- list() 
ct_min_ssh <- list() 
#gaussian (unmodified)
curve_breadth_g <- list() 
ct_max_g <- list() 
ct_min_g <- list() 
#lactin 
curve_breadth_lac <- list() 
ct_max_lac <- list() 
ct_min_lac <- list() 


#generating ctmax, min, etc from sharpe schoolfield high 
for (tx in Taxa) { 
  per_taxa <- ATdata[ATdata$Taxon == tx,] 
  # ggplot(per_taxa, aes(Tleaf_C, A_umol_m2_s)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) + facet_wrap(~Taxon, ncol = 2) + 
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Metabolic rate',
  #        title = 'Respiration across temperatures in this species')
 
  one_species <- str(tx)  
  Species_ssh <- append(Species_ssh, list(one_species))
  Ct_max_ssh <- append (ct_max_ssh, list(values$ctmax))
  
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


#generating ctmax, min, etc from gaussian
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


#generating ctmax, min, etc from gaussian
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


#actually outputting the ctmax, etc values from each model into a csv: 
#well, just tried ssh so far, and then ran into an error 
#error:so the values you want are stored in the list, but when you try to generate a csv, each value 
#gets generate as its own column rather than as rows in one column (wtih each column being a list like Species)


sharpe_high_data <- data.frame(Species = Species_ssh)
file_path <- "sharpe_high_processed_valued_part_3.csv"
write.csv(sharpe_high_data, file = file_path, row.names = FALSE)




##Ignore this, this was just me playing around with a couple things: 
#also the code doesn't work because I deleted some things that you'd need for it to work 

#trying to fit all curves 
#write function to label ggplot2 panels
label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

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
