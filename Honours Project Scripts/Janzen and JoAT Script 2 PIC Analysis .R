#PIC Analysis for both hypotheses 

#Load your libraries 
library(phytools)
library(ggplot2)


#Load your data 
mytree <- read.tree('your_working_tree.tre')
in_right_order <- read.csv('in_right_order_after_s.antarctici.csv', header = T)

#### Actual PIC Analysis and graphing JoAT #### 

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
  labs(x =expression(paste('PIC-Corrected Difference in Optimal Assimilation Rate (' *italic(A[opt]) * ')')),#'Optimal Assimilation Temperature, Aopt', 
       y = expression(paste('AT curve breadth (', Omega, ' )'))) +  
  ggtitle(expression('Optimal Assimilation Rate (' * italic(A[opt]) * ') vs Curve Breadth (' * Omega * ' ) Per Species ')) + 
  theme(axis.title.x = element_text(size = 16), axis.title.y =element_text(size = 16), 
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18))
print (PIC_JoAT) 

JoAT_lm <- lm(PIC_breadth ~ PIC_amax_log)
summary(JoAT_lm)
correlation_test <- cor.test(PIC_amax_log, PIC_breadth)
print(correlation_test) 

#### Actual PIC analysis and graphing, Janzen #### 

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
                       'Phoenix_dactylifera' )
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


PIC_Janzen <- ggplot (pic_df, aes(x = PIC_lat, y = PIC_breadth)) + 
  geom_point(alpha = 0.4, size = 3, color = 'yellowgreen') + 
  geom_smooth(method = 'lm', se = TRUE, color = 'seagreen', alpha = 0.5) +   theme_classic() + 
  labs(x = 'PIC-corrected Absolute global latitude (°)', y = expression(paste('PIC-corrected AT curve breadth (', Omega, ' )'))) + 
  ggtitle(expression(paste ('Plant Provenance Absolute Latitude(°) vs Curve Breadth (', Omega, ' ) Per Species'))) + 
  theme(axis.title.x = element_text(size = 16), axis.title.y =element_text(size = 16), 
        axis.text = element_text(size = 16), plot.title = element_text(size = 18))
print (PIC_Janzen) 

Janzen_lm <- lm(PIC_breadth ~ PIC_lat)
summary(Janzen_lm)
correlation_test <- cor.test(PIC_lat, PIC_breadth)
print(correlation_test) 




