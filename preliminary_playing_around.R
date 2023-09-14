#Basics - Just playing with the data + GitHub 


library(ggplot2)
library(dplyr)
AT <- read.csv("AT_data.csv", stringsAsFactors = TRUE)

#Generating a PDF for P. Lanceolata curve (1 graph per sample)
pdf ("AT_T2.pdf")
ggplot(data = filter(AT, Taxon %in% c("Plantago lanceloata")), 
       aes(x= Tleaf_C, y = A_umol_m2_s))+ geom_point() + geom_smooth() + facet_wrap(Taxon~curveID, ncol = 2)
dev.off()

#Generating a PDF for all P.Lanceolata curves but on different pages per graph 
pdf ("AT_P.lanceloata.pdf")
for (g in unique(AT$curveID)) {
  p <- ggplot(subset((data = filter(AT, Taxon %in% c("Plantago lanceloata"))), curveID == g), 
    aes(x= Tleaf_C, y = A_umol_m2_s))+ geom_point() + geom_smooth() + 
    facet_wrap(Taxon~curveID, ncol = 2)
  print(p)
}
dev.off()  
#It gives you a pdf #, but that doesn't really indicate # of pages or anything like that 

#All species, 1 graph per species 
pdf("AT_Curve_Per_Species.pdf")
for (g in unique(AT$Taxon)) {
  p <- ggplot(subset(AT, Taxon == g), 
    aes(x= Tleaf_C, y = A_umol_m2_s))+ geom_point() + geom_smooth() + 
    facet_wrap(~Taxon, ncol = 2)
  print(p) 
}
dev.off()

#Comparing Maxima:
plgraphs <- ggplot(data = filter(AT, Taxon %in% c("Plantago lanceloata")), 
       aes(x= Tleaf_C, y = A_umol_m2_s))+ geom_point() + geom_smooth() + facet_wrap(Taxon~curveID, ncol = 2)

print (plgraphs)
max_values <- aggregate(A_umol_m2_s ~ Taxon, AT, max)
print(max_values)

#so this gives you the max asat value per species, and it arranges the species alphabetically 
#not sure how to find the max of each *curve* within a species - I guess you can filter out the species 
#beforehand, and then use that filtered data to find the max of each curve 

#Generating a PDF for each curves of all species: 
#pdf ("AT_All_Curve_Per_Page.pdf")
#for (g in unique(AT$curveID)) {
#  p <- ggplot(subset(AT, curveID == g), 
#              aes(x= Tleaf_C, y = A_umol_m2_s))+ geom_point() + geom_smooth() + 
#    facet_wrap(Taxon~curveID, ncol = 2)
#  print(p)
#}
#dev.off()  


#trying to fit a gaussian on a curve 

# load in ggplot
library(ggplot2)
library(nls.multstart)
library(rTPC)
library(tidyverse)

# subset for the first TPC curve
#data('AT') - #you've already loaded your data, so don't need this line  
d <- subset(AT, curveID == 694) #Schistidium antarctici, 694
d2 <- subset(AT, Taxon == 'Aristolochia tonduzii')
print(d2$A_umol_m2_s)
# get start values and fit model
start_vals <- get_start_vals(d2$Tleaf_C, d2$A_umol_m2_s, model_name = 'gaussian_1987')

# fit model
#d$modelweights <- NULL
mod <- nls.multstart::nls_multstart(A_umol_m2_s~gaussian_1987(temp = Tleaf_C,rmax, topt,a),
#mod <- nls.multstart::nls_multstart(rate~gaussian_1987(temp = Tleaf_C,rmax, topt,a),
                                    data = d2,
                                    iter = c(4,4,4),
                                    start_lower = start_vals - 10,
                                    start_upper = start_vals + 10,
                                    lower = get_lower_lims(d2$Tleaf_C, d2$A_umol_m2_s, model_name = 'gaussian_1987'),
                                    upper = get_upper_lims(d2$Tleaf_C, d2$A_umol_m2_s, model_name = 'gaussian_1987'),
                                    supp_errors = 'Y',
                                    convergence_count = FALSE)


# look at model fit
summary(mod)

# get predictions
preds <- data.frame(temp = seq(min(d2$Tleaf_C), max(d2$Tleaf_C), length.out = 100))
preds <- broom::augment(mod, newdatas = preds)
#I changed it from newdata to newdatas, because I think I have something else haved as newdata and that's interfering 
#with the function of this code 
# plot
ggplot(preds) +
  geom_point(aes(Tleaf_C, A_umol_m2_s), d2) +
  geom_line(aes(Tleaf_C, .fitted), col = 'blue') +
  theme_bw()

