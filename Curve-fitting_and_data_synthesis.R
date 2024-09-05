
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
