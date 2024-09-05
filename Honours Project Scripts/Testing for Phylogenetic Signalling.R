#generating a phylogeny and a tree 

#lets load your packages first 
library(stringr)
library(V.PhyloMaker)
library(ggtree)
library(ggplot2)
library(tidyverse) 
library(TNRS)  
library(beepr) 
library(taxize)  
library(geiger) 
library(ape)
#library(phytools)


#Here is your starting dataset with curve Id, species, breadth, topt, and lat/long values  
starting_dataset <- read.csv('final_curve_parameter_filter_2_info_11.3.24.csv', stringsAsFactors = TRUE)
starting_dataset <- starting_dataset %>%
  filter(is.na(june_failure))


#### Updating taxon names #### 

update_taxon = function(data) {
  
  data$Name_submitted = tolower(data$Taxon) # Set taxon names to all lower case - prevents problems
  query_data = unique(data$Name_submitted)  # Extract unique names to prevent too many queries
  resolved_taxa = TNRS(query_data)    # Query TNRS with species list, get correct names
  
  resolved_taxa[resolved_taxa$Name_matched == "[No match found]",]$Name_matched = str_to_sentence(resolved_taxa[resolved_taxa$Name_matched == "[No match found]",]$Name_submitted)
  
  data = resolved_taxa %>%                      # Grab dataframe with corrected and submitted names
    select(Name_submitted, Name_matched) %>%    # Select only the columns we're interested in 
    merge(data, by = "Name_submitted") %>%      # Merge results dataframe with corrected names dataframe
    mutate(Taxon = Name_matched) %>%            # Add new column Taxon with corrected name
    select(-Name_submitted, -Name_matched) %>%  # Remove extra columns leaving only the corrected Taxon
    arrange(Taxon)                            # Sort by curveID to put back in original order
  
  return(data)
}

updated_taxon_dataset <- update_taxon(starting_dataset)


####To see which species got taken out due to phylogenetic issues ####
starting_species <- starting_dataset %>% 
  group_by(Taxon) %>% 
  slice(1)
updated_species <- updated_taxon_dataset %>% 
  group_by(Taxon) %>% 
  slice(1)

class(updated_taxon_dataset)
class(starting_dataset)
#
number_starting_species <- as.matrix(starting_species$Taxon)
number_updated_species <- as.matrix(updated_species$Taxon)
number_updated_species <- rbind(number_updated_species, "value1")

species_comparison_df <- data.frame(Updated_species = number_updated_species, 
                                    Starting_species = number_starting_species)

write.csv(species_comparison_df, file = 'starting_vs_updated_species_list.csv')
#from this, you can confirm the two subspecies that get combined into 1 thing 

#view(updated_taxon_dataset) #and now all your species names have been updated 


#### Averaging values by species #### 
average_by_taxon <- function (data) {
  data <- data %>% subset(is.na(june_failure))
  data <- select(data, c(Taxon, topt_june, amax_june, breadth_june, Lat_provenance, Long_provenance)) 
  data <- data %>% 
    group_by(Taxon) %>% 
    summarize (av_topt = mean(topt_june), 
               av_amax = mean(amax_june), 
               av_breadth = mean(breadth_june), 
               av_lat = abs(mean(Lat_provenance)), 
               av_lon = abs(mean(Long_provenance)))
  data <- as.data.frame(data)
}

averaged_dataset <- average_by_taxon(updated_taxon_dataset) 


#### Adding a family and genus column and saving as a dataset #### 

extract_phylogenetic_info = function (data) {
  genus_list <- list () 
  family_list <- list()
  for (i in data$Taxon) {
    #errorcatch while getting a family 
    tryCatch({
      family_name <- tax_name(query = as.character(i), get = "family", db = "ncbi") #get the family name 
      print(family_name)
      #beep('coin')
    }, error = function(e) {
      beep('ping')
      family_name <- 'NA'
    })
    
    #errorcatch while getting a genus name 
    tryCatch({
      genus_name <- tax_name(query = as.character(i), get = "genus", db = "ncbi") #get the genus name 
      print(genus_name)
    }, error = function(e) {
      beep('mario')
      family_name <- 'NA'
    })
    genus_list <- append(genus_list, genus_name$genus)
    family_list <- append(family_list, family_name$family)

  }
  data$family <- family_list 
  data$genus <- genus_list 
  beep('coin')
  return(data)
}

#generate a dataset and save it as a csv for future use 
phyloinfo_containing_dataset <- extract_phylogenetic_info(averaged_dataset)
phyloinfo_containing_dataset <- as.matrix.data.frame(phyloinfo_containing_dataset)
write.csv(phyloinfo_containing_dataset, file = 'phyloinfo_containing_dataset_new.csv')



#### STARTING POINT FOR PRESAVED PHYLOGENY DATASET - Making the phylogeny #### 

phyloinfo_containing_dataset <- read.csv ('phyloinfo_containing_dataset.csv', stringsAsFactors = TRUE)
#and get a working dataset with all the columns in the right order 
in_right_order <- phyloinfo_containing_dataset %>% 
  select (Taxon, genus, family, av_breadth, av_topt, av_amax, av_lat)


#from your big dataset with trait data, let's generate a dataset you can make a phylogeny from 
phylogeny_producing_dataset <- function (some_dataset) { 
  #let's extract the species and family names 
  some_dataset %>% 
    select(Taxon, Family = family) %>% #capitalizing the family column 
    unique () -> step_1 
  #now let's get the genus name by breaking up the binomial 
  step_1$Taxon %>% 
    str_split_fixed(' ', 2) %>% 
    as.data.frame() %>% 
    rename (Genus = V1, Species = V2) %>% 
    bind_cols(step_1) -> step_2
  #let's rename things (non capitalize) and remove the 'species' column that jus has the specific epithet 
  step_2 %>% rename (genus = Genus, species = Taxon, family = Family) %>% 
    select(-Species) -> step_3 
  #let's get the columns in the right order 
  step_3 %>% select (species, genus, family) -> step_3 
  #and now let's change the speices names so there's underscores, not spaces 
  step_3$species <- gsub (' ', '_', step_3$species)
  
  #and now let's return this new, polished dataset 
  return (step_3)
  }

tree_generating_dataset <- phylogeny_producing_dataset(in_right_order)
#View(tree_generating_dataset) 

#### Making tree and quality checking  #### 

#making the tree 
set.seed(0)
tree <- phylo.maker(tree_generating_dataset, tree = GBOTB.extended, output.tree = T, r = 1)
tree3 <- tree$scenario.3 

in_right_order <- phyloinfo_containing_dataset %>% 
  select (Taxon, genus, family, av_breadth, av_topt, av_amax, av_lat)
rownames(in_right_order) <- gsub(" ","_",in_right_order$Taxon)

#checking that your tree and the dataset you're feeding in has the same species and such 
obj<-name.check(tree3,in_right_order)
obj #Sort of!

#for now, take out the problematic species that aren't in the tree 
in_right_order <- subset(in_right_order, 
                         in_right_order$Taxon != "Schistidium antarctici")
#let's make the tree dichotomous and rooted (and overall non-problematic)
mytree=multi2di(tree3)


#### Phylogenetic Signalling ##### 


# First, you need to define which trait you want to test and give names to each value according to species
names(in_right_order)
#breath 
breadth_trait <- in_right_order[,4]
names(breadth_trait)<-rownames(in_right_order)
#topt 
topt_trait <- in_right_order[,5]
names(topt_trait)<-rownames(in_right_order)
#aopt 
amax_trait <- in_right_order[,6]
names(amax_trait)<-rownames(in_right_order)

##Let's get?? Pegel's lambda? Whait would you ge a different value ever time you run it? 
phylosig(mytree, breadth_trait, method="lambda", test=TRUE, nsim=999)
# Phylogenetic signal lambda : 7.34969e-05 
# logL(lambda) : -177.529 
# LR(lambda=0) : -0.00295999 
# P-value (based on LR test) : 1 
phylosig(mytree, topt_trait, method="lambda", test=TRUE, nsim=999)
# Phylogenetic signal lambda : 0.507889 
# logL(lambda) : -175.099 
# LR(lambda=0) : 4.19003 
# P-value (based on LR test) : 0.0406623 
phylosig(mytree, amax_trait, method="lambda", test=TRUE, nsim=999)
# Phylogenetic signal lambda : 0.413986 
# logL(lambda) : -228.409 
# LR(lambda=0) : 2.2425 
# P-value (based on LR test) : 0.134264 





