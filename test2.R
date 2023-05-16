#Trying out pushing and pulling code/playing with AT Data 

#radom thing you're trying 
library(ggplot2)
library(dplyr)

#Generating a PDF for P. Lanceolata curve (1 graph per sample)
AT <- read.csv("AT_data.csv", stringsAsFactors = TRUE)
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
PDF("AT_Curve_Per_Species.pdf")
for (g in unique(AT$Taxon)) {
  p <- ggplot(subset(AT, Taxon == g), 
    aes(x= Tleaf_C, y = A_umol_m2_s))+ geom_point() + geom_smooth() + 
    facet_wrap(~Taxon, ncol = 2)
  print(p) 
}
dev.off()
##hmm this doesn't work, gives null device, 1 


#Generating a PDF for all curves of all species: 
#pdf ("AT_All_Curve_Per_Page.pdf")
#for (g in unique(AT$curveID)) {
#  p <- ggplot(subset(AT, curveID == g), 
#              aes(x= Tleaf_C, y = A_umol_m2_s))+ geom_point() + geom_smooth() + 
#    facet_wrap(Taxon~curveID, ncol = 2)
#  print(p)
#}
#dev.off()  


