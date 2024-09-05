#Excluding Cultivars and testing signalling 

#let's load packages 
library(dplyr)
library(ggplot2)
library(plotly)
library(ggfortify)
library(geiger)
library(phytools)
library(viridis)
library(ggtree)

#load data 
mytree <- read.tree('your_working_tree.tre')
traits_df <- read.csv('in_right_order_after_s.antarctici.csv', header = T)

#### excluding cultivars and prepping dataset + tree #### 
rownames(traits_df) <- gsub(" ","_",traits_df$Taxon)
names(traits_df)

#take this out of the dataset 
traits_df <- traits_df %>% 
  filter (crops.y.n. != 'y') %>% 
  select (Taxon, genus, family, av_breadth, av_topt, av_amax, av_lat, random)  
  
#now let's see where the mistmatches are 
obj_1 <- name.check(mytree, traits_df)
obj_1

#species to remove from our tree: 
species_to_remove <- c('Brassica_campestris', 'Cannabis_sativa', 
                       'Glycine_max', 'Malus_domestica', 
                       'Oryza_sativa', 'Phoenix_dactylifera',  
                       'Saccharum', 'Solanum_lycopersicum', 
                       'Solanum_tuberosum', 'Triticum_aestivum', 
                       'Vitis_vinifera', 'Zea_mays')
#take them out of the tree 
pruned_tree <- drop.tip(mytree, species_to_remove)
#write.tree(pruned_tree, file = 'pruned_non_domestic_trees.tre')

#check that tree and data match 
obj<-name.check(pruned_tree,traits_df)
obj #slay 
pruned_tree=multi2di(pruned_tree)

#### Checking for phylogenetic signal #### 
#breath 
breadth_trait <- traits_df[,4]
names(breadth_trait)<-rownames(traits_df)
#topt 
topt_trait <- traits_df[,5]
names(topt_trait)<-rownames(traits_df)
#aopt 
amax_trait <- traits_df[,6]
names(amax_trait)<-rownames(traits_df)

##Let's get?? Pegel's lambda? Whait would you ge a different value ever time you run it? 
phylosig(pruned_tree, breadth_trait, method="lambda", test=TRUE, nsim=999)
# Phylogenetic signal lambda : 7.36576e-05 
# logL(lambda) : -137.493 
# LR(lambda=0) : -0.00208048 
# P-value (based on LR test) : 1 

phylosig(pruned_tree, topt_trait, method="lambda", test=TRUE, nsim=999)
# Phylogenetic signal lambda : 0.454978 
# logL(lambda) : -144.242 
# LR(lambda=0) : 1.46496 
# P-value (based on LR test) : 0.226142 

phylosig(pruned_tree, amax_trait, method="lambda", test=TRUE, nsim=999)
# Phylogenetic signal lambda : 0.446303 
# logL(lambda) : -181.682 
# LR(lambda=0) : 2.61629 
# P-value (based on LR test) : 0.105772 


#### Let's visualize this signalling ####
################### topt ##########################################################

topt_main_data <-dplyr:: select(traits_df, Taxon, av_topt)
topt_main_data$Taxon <- gsub(' ', '_', topt_main_data$Taxon)
rownames(topt_main_data) <- topt_main_data$Taxon
topt_main_data <- topt_main_data %>% dplyr::select(-Taxon)
topt_main_maxtrix <- as.matrix(topt_main_data)

p <- ggtree(pruned_tree, branch.length = 'none', layout = "circular") %<+% 
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

################### aopt ########################################################

aopt_main_data <-dplyr:: select(traits_df, Taxon, av_amax)
aopt_main_data$Taxon <- gsub(' ', '_', aopt_main_data$Taxon)
rownames(aopt_main_data) <- aopt_main_data$Taxon
aopt_main_data <- aopt_main_data %>% dplyr::select(-Taxon)
aopt_main_maxtrix <- as.matrix(aopt_main_data)

p <- ggtree(pruned_tree, branch.length = 'none', layout = "circular") %<+% 
  as.data.frame(aopt_main_maxtrix[,1]) 
p2_aopt <- gheatmap(p, as.data.frame(aopt_main_maxtrix[,1]), width = 0.2, font.size = 0) +
  scale_fill_viridis(option="mako", discrete = F, 
                     name = expression(atop("A"[opt], paste(µmol~m^-2~s^-1))), 
                     #                     name=expression(paste('A'[opt], '(µmol m'^-2, 's'^-1, ')')), 
                     direction = -1) + 
  theme(legend.key.size = unit(1,"line"), legend.position = c(1.1, 0.5),
        legend.title = element_text(hjust = 2)) 
p2_aopt

################### breadth ########################################

bopt_main_data <-dplyr:: select(traits_df, Taxon, av_breadth)
bopt_main_data$Taxon <- gsub(' ', '_', bopt_main_data$Taxon)
rownames(bopt_main_data) <- bopt_main_data$Taxon
bopt_main_data <- bopt_main_data %>% dplyr::select(-Taxon)
bopt_main_maxtrix <- as.matrix(bopt_main_data)

p <- ggtree(pruned_tree, branch.length = 'none', layout = "circular") %<+% 
  as.data.frame(bopt_main_maxtrix[,1]) 
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




#### Making a random fig #### 

random_main_data <-dplyr:: select(traits_df, Taxon, random)
random_main_data$Taxon <- gsub(' ', '_', random_main_data$Taxon)
rownames(random_main_data) <- random_main_data$Taxon
random_main_data <- random_main_data %>% dplyr::select(-Taxon)
random_main_maxtrix <- as.matrix(random_main_data)

p <- ggtree(pruned_tree, branch.length = 'none', layout = "circular") %<+% 
  as.data.frame(random_main_maxtrix[,1]) 
p2_random <- gheatmap(p, as.data.frame(random_main_maxtrix[,1]), width = 0.2, font.size = 0) +
  scale_fill_viridis(option="viridis", discrete = F, 
                     name=expression(paste('Trait value')), 
                     direction = 1)  + 
  theme(legend.key.size = unit(1,"line"), legend.position = c(1.07, 0.5),
        legend.title = element_text(hjust = 1.0)) 

p2_random

