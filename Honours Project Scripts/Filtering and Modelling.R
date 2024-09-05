#Ria Raut 
#let's load your packages first 
#library(nls2)
library(nls.multstart)
library(broom) 
library(tidyverse)
library(rTPC) 
library(scales)
library(gridExtra)
#library(nlstools)
library(beepr)

#### the june model  #### 
#The June model:  
June_2004 <- function (temp, aopt, topt, omega) {
  exponent <- (-1)*(((temp-topt)/omega)^2)
  p_t <- aopt*exp(exponent)
  return(p_t)
}


#### A series of preliminary errorcatch functions  ####
#the functions you've used: get_omega_errorcatch(fit1), get_topt(fit1), get_pure_aopt_june(fit1),

get_topt_errorcatch <- function(x) {
  result <- tryCatch(
    {
      get_topt(x)
    },
    error = function(e) {
      message("An error occurred with topt extraction: ", conditionMessage(e))
      return(NA)  # Store NA when an error occurs
    }
  )
  return(result)
}

get_omega_errorcatch <- function(x) {
  result <- tryCatch(
    {
      this_is_omega <- coef(x)['omega']
      return (this_is_omega)
    },
    error = function(e) {
      message("An error occured with omega extraction: ", conditionMessage(e))
      return(NA)  # Store NA when an error occurs
    }
  )
  return(result)
}

get_aopt_errorcatch <- function(x) {
  result <- tryCatch(
    {
      this_is_aopt <- coef(x)['aopt']
      return (this_is_aopt)

    },
    error = function(e) {
      message("An error occured with Aopt extraction: ", conditionMessage(e))
      return(NA)  # Store NA when an error occurs
    }
  )
  return(result)
}


#### filter 1 - add warnings to unfit curves #### 

filter_1 <- function(data) {
  
  filter_criteria = c() # Holder
  
  #new_dat = c() ###
  # Find metrics for filtering for each curve
  for (i in unique(data$curveID)) {
    
    cur_dat = subset(data, curveID == i) # Grab the current data set
    #resetting the 'odd' values each time so that when you have a curve with standard units, the curve doesn't 
    #fail just because the previous curve had odd units and had failed the criteria 
    nbefore_odd = 3
    nafter_odd = 3 
    Tbefore_odd = 5 
    Tafter_odd = 5
    #but if your curve actually has odd units, then you reset the nbefore_odd, etc values. 
    #this is possible because of this if statement below:  
    if (is.na(cur_dat$A_umol_m2_s[[1]])) { #i.e. if you don't have standard units for this curve 
      Tpk_odd = cur_dat$Tleaf_C[which(cur_dat$A_other_units == max(cur_dat$A_other_units)) [1]] #find peak w odd units 
      nbefore_odd = sum( cur_dat$Tleaf_C < Tpk_odd ) # Find number of points before peak
      nafter_odd = sum( cur_dat$Tleaf_C > Tpk_odd )  # Find number of points before peak
      Tbefore_odd = Tpk_odd-min(cur_dat$Tleaf_C)     # Find temperature span before peak
      Tafter_odd = max(cur_dat$Tleaf_C)-Tpk      # Find temperature span after peak 
    } else { #else, if you DO have standard units for this curve: 
      Tpk = cur_dat$Tleaf_C[which( cur_dat$A_umol_m2_s == max(cur_dat$A_umol_m2_s) )[1]] # Find peak
      #cur_dat = subset(cur_dat, Tleaf <= Tpk) ###
      nbefore = sum( cur_dat$Tleaf_C < Tpk ) # Find number of points before peak
      nafter = sum( cur_dat$Tleaf_C > Tpk )  # Find number of points before peak
      Tbefore = Tpk-min(cur_dat$Tleaf_C)     # Find temperature span before peak
      Tafter = max(cur_dat$Tleaf_C)-Tpk      # Find temperature span after peak 
      
    } #this bracket ends your if statement 
    filter_criteria = rbind(filter_criteria, data.frame(curveID = i, nbefore, nafter, Tbefore, Tafter, 
                                                        nbefore_odd, nafter_odd, Tbefore_odd, Tafter_odd)) 
    print(i)
  } #this bracket ends your for loop 

  # Merge filtering metrics with dataset for ease of filtering
  data = merge(data, filter_criteria, by="curveID")
  #make a new column to add your filter statuses in: 
  data$failure_status = NA
  
  # Reject curves if the are decreasing only, have too little T span, or too few points
  data$failure_status = ifelse(data$nbefore < 2, "Rejected: fewer than two points before peak", data$failure_status)
  data$failure_status = ifelse(data$nafter < 2, "Rejected: fewer than two points after peak", data$failure_status)
  data$failure_status = ifelse(data$Tbefore < 2, "Rejected: ascending region too small", data$failure_status)
  data$failure_status = ifelse(data$Tafter < 2, "Rejected: descending region too small", data$failure_status)
  data$failure_status = ifelse(data$nbefore+data$nafter < 5, "Rejected: fewer than 5 data points", data$failure_status)
  #same thing but for odd units: 
  data$failure_status = ifelse(data$nbefore_odd < 2, "Odd units. Rejected: fewer than two points before peak", data$failure_status)
  data$failure_status = ifelse(data$nafter_odd < 2, "Odd units.Rejected: fewer than two points after peak", data$failure_status)
  data$failure_status = ifelse(data$Tbefore_odd < 2, "Odd units.Rejected: ascending region too small", data$failure_status)
  data$failure_status = ifelse(data$Tafter_odd < 2, "Odd units.Rejected: descending region too small", data$failure_status)
  data$failure_status = ifelse(data$nbefore_odd+data$nafter_odd < 5, "Odd units. Rejected: fewer than 5 data points", data$failure_status)
  
  
  data_after_filter1 = data %>% select(-nbefore, -nafter, -Tbefore, -Tafter,
                                       -nbefore_odd, -nafter_odd, -Tbefore_odd, -Tafter_odd)  # Remove junk
  beep(1) #beeps when the function is done! 
  return(data_after_filter1)
}

#now let's call this function with our dataset: 
all_data <- read.csv('All_AT_data.csv', stringsAsFactors = TRUE)
passed_filter_1 <- filter_1(all_data)
#view(passed_filter_1)





#### What criteria is taking out most of the curves? ####
#We are counting to see how many curves have specific failure statuses. Just to examine our data and filtering! 

counter <- function (data) {
  fewer_than_2_before_peak <- 0 
  fewer_than_2_after_peak <- 0 
  ascending_region_too_small <- 0 
  descending_region_too_big <- 0 
  fewer_than_5_points <- 0 
  no_issue <- 0 
  yes_issue <- 0 
  total_count <- 0 
  
  cvID <- unique(data$curveID)
  for (cv in cvID) {
    total_count <- total_count + 1
    per_curve <- data[data$curveID == cv,] #take the dataframe for the first curve 
    curve_status <- per_curve[1,29] #take the failure_status from the first row for this curve  
    #print(curve_status)
    if (is.na(curve_status)) {
      no_issue <- no_issue + 1
    } else if (curve_status == 'Rejected: fewer than two points before peak' | curve_status == 'Odd units. Rejected: fewer than two points before peak') {
      fewer_than_2_before_peak <- fewer_than_2_before_peak + 1
      
    } else if (curve_status == 'Rejected: fewer than two points after peak' | curve_status == 'Odd units.Rejected: fewer than two points after peak') {
      fewer_than_2_after_peak <- fewer_than_2_after_peak + 1
      
    } else if (curve_status == 'Rejected: ascending region too small' | curve_status == 'Odd units.Rejected: ascending region too small') {
      ascending_region_too_small <- ascending_region_too_small + 1
      
    } else if (curve_status == 'Rejected: descending region too small' | curve_status == 'Odd units.Rejected: descending region too small') {
      descending_region_too_big <- descending_region_too_big + 1
      
    } else if (curve_status == 'Rejected: fewer than 5 data points' | curve_status == 'Odd units. Rejected: fewer than 5 data points') {
      fewer_than_5_points <- fewer_than_5_points + 1
    } else {
      yes_issue <- yes_issue + 1
    }
  }
  return (list(no_issue,                         #1
               fewer_than_2_before_peak,         #2
               fewer_than_2_after_peak,          #3
               ascending_region_too_small,       #4
               descending_region_too_big,        #5
               fewer_than_5_points,              #6
               yes_issue,                        #to catch some error messages that might show up. Theoretically, this should be a value of zero
               total_count))  
}

#let's run the function to see how many curves get filtered out at what stage
issue_check <- counter(passed_filter_1)
print(issue_check)


#### split and model all functions #### 

#We'll have a couple nested functions to make sure all spitting (by units) and modelling happens in one 
#go, and is easy to understand. 

# one_dataset <- passed_filter_1

split_and_model_all <- function(one_dataset) {
  
  #defining function 1
  innerfunction1 <- function(some_dataset) {
    filter_1_done <- some_dataset %>% 
      filter(is.na(failure_status), #taking out all curves that failed the first filter step 
             LI_6400_method != 'In situ', General_method != 'in situ',    #taking out all the in situ curves 
             Habitat == 'Terrestrial') #making sure to take out anything aquatic 
    #then make a dataset with just curves with odd units. As in, the A_umol_m2_s column must have 'NA'
    pf1_odd <- filter_1_done %>% 
      filter(is.na(A_umol_m2_s)) 
    #then make a dataset containing all the curves with the standard units. As long as the A_umol_m2_s column isn't empty! 
    pf1_standard <- subset(filter_1_done, A_umol_m2_s != 'NA')
    
    #then return two datasets, both sufficiently filtered: one with curves with standard units, 
    #one with curves of odd/miscellaneous units 
    return (list(pf1_standard, pf1_odd))
  }
  
  #calling function 1 
  innerfunction1_results <- innerfunction1(one_dataset)
  standard_units_dataset <- innerfunction1_results[[1]]
  odd_units_dataset <- innerfunction1_results[[2]]
  
  #now let's get to modelling: note that the next two function are the exact same, just deal with 
  #the different units of all the curves (which we split into odd and standard)
  #setting up blank datasets for the next 2 functions 
  data_out_standard = data.frame() 
  data_out_odd = data.frame() 

  cvID <- unique(standard_units_dataset$curveID) #make cvID a list of curve IDs (unique)

  #defining function 2 
  innerfunction2_standard <- function(data_to_fit) {
    listofplots_standard <- list() #this is a list for all the plots. Useful for generating a pdf of all the modelled plots! 
    for(cv in cvID){
      per_curve <- data_to_fit[data_to_fit$curveID == cv,] #grabbing data for your current curveID 
      #fitting the June model: 
      fit1 <- nls_multstart(A_umol_m2_s ~ June_2004(temp = Tleaf_C, aopt, topt, omega), 
                            data = per_curve,
                            iter = 50,
                            start_lower = c(aopt = 0, topt = 0, omega = 2),
                            start_upper = c(aopt = 20, topt =50, omega =10), #NOTE: you can have the lower and upper limits based on what your curves look like in your dataset
                            lower = c(aopt = 0, topt = 0, omega = 0), 
                            upper = c(aopt = 90, topt = 100, omega = 100),
                            supp_errors = 'Y')
      values <- calc_params(fit1)
      new_data1 <- data.frame(temp = seq(min(per_curve$Tleaf_C), max(per_curve$Tleaf_C)), length = 200) #I think the length arguement can help smoothen your curves? 
      preds1 <- augment(fit1, new_data = new_data1)
      print (cv) #just printing the curve ID to make sure each curve is going through/you're making progress
      
      #now let's actually graph the model + data 
      graph_main <- ggplot(per_curve, aes(Tleaf_C, A_umol_m2_s)) +
        geom_point() +
        geom_line(aes(x = Tleaf_C, y = .fitted), data = preds1, color = 'red')+ #this is the line graphing your preds 
        theme_bw(base_size = 12) + facet_wrap(~curveID, ncol = 2)
      labs(x = 'Temperature (ºC)',
           y = 'Metabolic rate')
      #let's store this graph in a list so that you can print it all in a pdf at the end 
      listofplots_standard [[cv]] <- graph_main
      
      #values to generate out to your csv. Taxon, units, curve ID  
      #Ria you named these backwards, so the 'each_taxon' variable, for example, is actually a list of taxons 
      #why are you like this, please fix it 
      taxon_list <- per_curve$Taxon
      each_taxon <- taxon_list[[1]]
      curve_id_list <- per_curve$curveID
      each_id <- curve_id_list[[1]]
      lat_list <- per_curve$Lat_provenance
      each_lat <- lat_list[[1]]
      long_list <- per_curve$Lon_provenance
      each_long <- long_list[[1]]
      
      
      #SE for each model (June_2004) 
      #variables: temp = Tleaf_C, aopt, topt, omega 
      Aopt_SE <- summary(fit1)$coefficients[, 'Std. Error']['aopt']
      each_Aopt_SE <- Aopt_SE[[1]]
      Topt_SE <- summary(fit1)$coefficients[, 'Std. Error']['topt'] 
      each_Topt_SE <- Topt_SE [[1]] 
      Omega_SE <- summary(fit1)$coefficients[, 'Std. Error']['omega']
      each_Omega_SE <- Omega_SE[[1]]
      rsq_june = 1 - sum( resid(fit1)^2)/ sum((per_curve$A_umol_m2_s - mean(per_curve$A_umol_m2_s))^2)
      rsq_june_list <- rsq_june[[1]]
      
      #this (below) is basically taking each curve's taxon and making it an object in a list. Then you're taking the 
      #first object in the list (the current held taxon) and adding that to your dataset, I think. Doing the same for 
      #all other properties too (like breadth). 
      data_out_standard = bind_rows (data_out_standard, data.frame(Curve_ID  = each_id, Taxon = each_taxon,
                                                                   breadth_june = get_omega_errorcatch(fit1),
                                                                   breadth_june_SE = each_Omega_SE,
                                                                   topt_june = get_topt_errorcatch(fit1),
                                                                   topt_june_SE = each_Topt_SE, 
                                                                   amax_june = get_aopt_errorcatch(fit1),
                                                                   amax_june_SE = each_Aopt_SE, 
                                                                   r_squared_june = rsq_june_list, 
                                                                   
                                                                   Lat_provenance = each_lat, 
                                                                   Long_provenance = each_long))
      
    } #close your for loop 
    pdf(file = 'graphs_standard_units_updated_again.pdf') #open a pdf
    print(listofplots_standard) #print all your standard unit graphs in said pdf
    dev.off() #close the pdf
    print(graph_main)
    return(data_out_standard) #and return the dataset with all modelled parameters 
  
  } #this bracket is to close your inner function 2
  
  #now let's set up inner function 3 the same way you did for 2  
  cvID_odd <- unique(odd_units_dataset$curveID)
  innerfunction3_odd <- function (data_to_fit) {
    listofplots_odd <- list()
    for(cv in cvID_odd){
      per_curve <- data_to_fit[data_to_fit$curveID == cv,]
      #June 
      fit1 <- nls_multstart(A_other_units ~ June_2004(temp = Tleaf_C, aopt, topt, omega),
                            data = per_curve,
                            iter = 50,
                            start_lower = c(aopt = 0, topt = 0, omega = 2),
                            start_upper = c(aopt = 500, topt =50, omega =10),
                            lower = c(aopt = 0, topt = 0, omega = 0),
                            upper = c(aopt = 900, topt = 100, omega = 100),
                            supp_errors = 'Y')
      values <- calc_params(fit1)
      new_data1 <- data.frame(temp = seq(min(per_curve$Tleaf_C), max(per_curve$Tleaf_C)), length = 200)
      preds1 <- augment(fit1, new_data = new_data1)
      print (cv)
      
      graph_main <- ggplot(per_curve, aes(Tleaf_C, A_other_units)) +
        geom_point() +
        geom_line(aes(x = Tleaf_C, y = .fitted), data = preds1, color = 'red')+
        theme_bw(base_size = 12) + facet_wrap(~curveID, ncol = 2)
      labs(x = 'Temperature (ºC)',
           y = 'Metabolic rate')
      listofplots_odd [[cv]] <- graph_main
      
      #values to generate out to your csv. Taxon, units, curve ID  
      taxon_list <- per_curve$Taxon
      each_taxon <- taxon_list[[1]]
      curve_id_list <- per_curve$curveID
      each_id <- curve_id_list[[1]]
      lat_list <- per_curve$Lat_provenance
      each_lat <- lat_list[[1]]
      long_list <- per_curve$Lon_provenance
      each_long <- long_list[[1]]
      
      #SE for fit1, June_2004 
      #varables: temp = Tleaf_C, aopt, topt, omega 
      Aopt_SE <- summary(fit1)$coefficients[, 'Std. Error']['aopt']
      each_Aopt_SE <- Aopt_SE[[1]]
      Topt_SE <- summary(fit1)$coefficients[, 'Std. Error']['topt'] 
      each_Topt_SE <- Topt_SE [[1]] 
      Omega_SE <- summary(fit1)$coefficients[, 'Std. Error']['omega']
      each_Omega_SE <- Omega_SE[[1]]
      rsq_june = 1 - sum( resid(fit1)^2)/ sum((per_curve$A_other_units - mean(per_curve$A_other_units))^2)
      rsq_june_list <- rsq_june[[1]]
      
      data_out_odd = bind_rows (data_out_odd, data.frame(Curve_ID  = each_id, Taxon = each_taxon,
                                                         breadth_june = get_omega_errorcatch(fit1),
                                                         breadth_june_SE = each_Omega_SE,
                                                         topt_june = get_topt_errorcatch(fit1),
                                                         topt_june_SE = each_Topt_SE, 
                                                         amax_june = get_aopt_errorcatch(fit1),
                                                         amax_june_SE = each_Aopt_SE, 
                                                         r_squared_june = rsq_june_list, 
                                                         Lat_provenance = each_lat, 
                                                         Long_provenance = each_long))
    } #close your for loop 
    pdf(file = 'graphs_odd_units_updated_again.pdf')
    print(listofplots_odd)
    dev.off()
    return(data_out_odd)
  } #close your inner function 3 
  
  #calling your two inner functions 
  #note that the function input has to match the input used to define cvID and cvID_odd (which you produce right before 
  #defining inner functions 1 and 2 respectively)
  data_out_standard <- innerfunction2_standard (standard_units_dataset)
  data_out_odd <- innerfunction3_odd (odd_units_dataset)
  #and you have the following outputs: (and also the 2 pdf's you made)
  beep(1) #beep when the dataset is ready and the pdfs have been printed! 
  return(list(data_out_standard, data_out_odd))
} #close your function 


#calling these functions  - you can delete this later 
#ncol(passed_filter_1)

two_dataset <- split_and_model_all(passed_filter_1) 
standard_units_updated <- two_dataset[[1]]
odd_units_updated <- two_dataset[[2]]

#let's save these generated datasets into csvs 
write.csv(standard_units_updated, 'standard_units_SE_updated_11.3.24.csv', row.names = FALSE)
write.csv(odd_units_updated, 'odd_units_SE_updated_11.3.24.csv', row.names = FALSE)


#### filter 2 warning flagging #### 
#after you've modelling the appropriate curves, you still have to run through the models and quality check 
#let's define a function for that 
filter_2_warningsystem <- function(dataset) {
  filter_criteria = c() # Holder
  for (i in unique(dataset$Curve_ID)) {
    #now, let's define a couple terms we can use to determine if a curve can fail or not 
    #ex. if breadth_error is negative, that means the SE is greater than the actual breadth 
    #same goes for topt and amax errors 
    curve_data = subset(dataset, Curve_ID == i) #getting the current dataset 
    r_2_error = curve_data$r_squared_june #just getting the r_squared/goodness of fit value 
    breadth_error = curve_data$breadth_june - curve_data$breadth_june_SE
    topt_error = curve_data$topt_june - curve_data$topt_june_SE
    amax_error = curve_data$amax_june  - curve_data$amax_june_SE

    filter_criteria = rbind(filter_criteria, data.frame(Curve_ID = i, r_2_error, 
                                                        breadth_error,
                                                        topt_error, amax_error)) 
    print(i)
  } #closing the for loop 
  dataset = merge(dataset, filter_criteria, by = 'Curve_ID') 
  
  #let's set up a failure status column in this dataset: 
  dataset$june_failure = NA 
  
  dataset$june_failure = ifelse(dataset$r_2_error <0.7, 'Warning: bad june fit.', dataset$june_failure) 
  dataset$june_failure = ifelse(dataset$breadth_error <0 , 'Warning: breadth has large SE.', dataset$june_failure)
  dataset$june_failure = ifelse(dataset$topt_error <0 , 'Warning: Topt has large SE.', dataset$june_failure)
  dataset$june_failure = ifelse(dataset$amax_error <0 , 'Warning: Amax has large SE.', dataset$june_failure)
  
  return(dataset)
} #close your function 

#calling this function for standard and odd units separately  
standard_after_filter_2 <- filter_2_warningsystem(standard_units_updated)
odd_after_filter_2 <- filter_2_warningsystem(odd_units_updated)

#let's combine these two into one dataset (combined vertically) 
final_combined <- rbind(standard_after_filter_2, odd_after_filter_2)
#NOTE: You have to keep in mind that if you are later comparing breadth vs aopt for these curves, you 
#CANNOT consider the odd unit amax values long with the standard unit curves!! 

write.csv(standard_after_filter_2, 'final_curve_parametesr_filter_2_info_11.3.24.csv', row.names = FALSE)

#and we've filtered and modelled all curves! Yay! 

