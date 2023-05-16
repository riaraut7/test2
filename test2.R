#Trying out pushing and pulling code/playing with AT Data 

#radom thing you're trying 
library(ggplot2)
library(reshape2)
library(dplyr)

#Trying to get Data from AT 
AT <- read.csv("AT_data.csv", stringsAsFactors = TRUE)

library(ggplot2)
library(dplyr)
ggplot(data = filter(AT, Taxon %in% c("Plantago lanceloata", "Pinus contorta")), 
       aes(x= Tleaf_C, y = A_umol_m2_s))+ geom_point() + geom_smooth() + facet_wrap(Taxon~curveID, ncol = 2)

