
#In this script, you will take your AT curves that have all passed the second round of filtering, and then 
#obtain the species, genus, and family names for each represented species. Then, you plot that up in a 
#phylogenetic tree 

#first, just load your libraries 
library(dplyr)
library(ggtree)
library(phytools)
library(viridis)
library(ggplot2)

################################ START FROM HERE ####################################
#let's first load our data and tree 
Main_tree <- read.tree('your_working_tree.tre')
class(Main_tree) #class phylo 
Main_data <- read.csv('in_right_order_after_s.antarctici.csv', stringsAsFactors = FALSE)
class(Main_data) #class data frame. Let's turn it to a matrix 

################## topt ##########################################################
##Let's try this with Topt, the only one that shows signalling 
names(Main_data)
topt_main_data <-dplyr:: select(Main_data, Taxon, av_topt)
topt_main_data$Taxon <- gsub(' ', '_', topt_main_data$Taxon)
rownames(topt_main_data) <- topt_main_data$Taxon
topt_main_data <- topt_main_data %>% dplyr::select(-Taxon)
topt_main_maxtrix <- as.matrix(topt_main_data)


Main_tree$tip.label

p <- ggtree(Main_tree, branch.length = 'none', layout = "circular") %<+% 
  as.data.frame(topt_main_maxtrix[,1]) 
#p <- p + geom_tiplab(aes(label = label), size = 3, offset = 4)  #the [,1] you dn't need right now because you only have just one column anyway 
p2 <- gheatmap(p, as.data.frame(topt_main_maxtrix[,1]), width = 0.2, font.size = 0) +
  scale_fill_viridis(option="magma", discrete = F, 
                     name=expression(T[opt]~(phantom()^o*C)), 
                     direction = 1) + 
  theme(legend.key.size = unit(1,"line"), legend.position = c(1.07, 0.5),
        legend.title = element_text(hjust = 2)) 
p2 
#you printed the pdf as 5x5

################### for aopt ########################################################
names(Main_data)
aopt_main_data <-dplyr:: select(Main_data, Taxon, av_amax)
aopt_main_data$Taxon <- gsub(' ', '_', aopt_main_data$Taxon)
rownames(aopt_main_data) <- aopt_main_data$Taxon
aopt_main_data <- aopt_main_data %>% dplyr::select(-Taxon)
aopt_main_maxtrix <- as.matrix(aopt_main_data)

p <- ggtree(Main_tree, branch.length = 'none', layout = "circular") %<+% 
  as.data.frame(aopt_main_maxtrix[,1]) #the [,1] you dn't need right now because you only have just one column anyway 

p2_aopt <- gheatmap(p, as.data.frame(aopt_main_maxtrix[,1]), width = 0.2, font.size = 0) +
  scale_fill_viridis(option="mako", discrete = F, 
                     name = expression(atop("A"[opt], paste(µmol~m^-2~s^-1))), 
#                     name=expression(paste('A'[opt], '(µmol m'^-2, 's'^-1, ')')), 
                     direction = -1) + 
  theme(legend.key.size = unit(1,"line"), legend.position = c(1.1, 0.5),
        legend.title = element_text(hjust = 2)) 
p2_aopt

################################## AT breadth ########################################

names(Main_data)
bopt_main_data <-dplyr:: select(Main_data, Taxon, av_breadth)
bopt_main_data$Taxon <- gsub(' ', '_', bopt_main_data$Taxon)
rownames(bopt_main_data) <- bopt_main_data$Taxon
bopt_main_data <- bopt_main_data %>% dplyr::select(-Taxon)
bopt_main_maxtrix <- as.matrix(bopt_main_data)

p <- ggtree(Main_tree, branch.length = 'none', layout = "circular") %<+% 
  as.data.frame(bopt_main_maxtrix[,1]) #the [,1] you dn't need right now because you only have just one column anyway 

p2_breadth <- gheatmap(p, as.data.frame(bopt_main_maxtrix[,1]), width = 0.2, font.size = 0) +
  scale_fill_viridis(option="cividis", discrete = F, 
                     name=expression(paste( Omega,' (°C)')), 
                     direction = 1)  + 
  theme(legend.key.size = unit(1,"line"), legend.position = c(1.07, 0.5),
        legend.title = element_text(hjust = 1.0)) 
  
p2_breadth

### combining all into 1 figure 
library(patchwork)
p2 
p2_aopt
p2_breadth
combined_phylo <- p2 / p2_aopt / p2_breadth + plot_layout(ncol = 1)
combined_phylo

