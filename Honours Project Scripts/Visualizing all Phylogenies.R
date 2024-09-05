#Visualizing all phylogenies 
#Ria Raut 

#load libraries 
library(ggtree) 
library(tidyverse)
library(TNRS)
library(beepr)
library(taxize)
library(V.PhyloMaker)
library(geiger)

#### Visualizing Past-Filter phylogeny ####
Main_tree <- read.tree('your_working_tree.tre')

p2 <- ggtree(Main_tree)  + theme_tree() + 
  geom_tiplab( 
    color = 'black',                      
    offset = 5,
    size = 3) + 
  hexpand(-10, direction = -1) 
  
p2


#### Producing Phylogeny Prior to Filtering #### 
All_data <- read.csv('All_AT_data.csv', stringsAsFactors = TRUE)
first_rows <- All_data %>%
  group_by(Taxon) %>% 
  filter(Habitat == 'Terrestrial') %>% 
  slice(1) 


#let's first update taxon names 
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

#call function 
updated_taxons <- update_taxon(first_rows)


#Now we get genus and family names 
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

#call function 
sp_ge_fa_alldata <- extract_phylogenetic_info(updated_taxons)
in_right_order <- sp_ge_fa_alldata %>% 
  select (Taxon, genus, family)


##now let's organize the dataset so it's workable for phylomaker 
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

#and call: 
tree_generating_dataset <- phylogeny_producing_dataset(in_right_order)


#now input the dataset into phylomaker to generate a tree 
set.seed(0)
tree <- phylo.maker(tree_generating_dataset, tree = GBOTB.extended, output.tree = T, r = 1)
tree3 <- tree$scenario.3 
rownames(tree_generating_dataset) <- tree_generating_dataset$species

#let's check that tree and data match - important for future analysis/good practice 
check_congruency <- name.check(tree3, tree_generating_dataset)
check_congruency

#some bryophytes are not found on the tree, so let's take those out 
updated_species_dataset <- tree_generating_dataset %>% 
  filter(!(species %in% c("Schistidium_antarctici", 
                      "Ceratodon_purpureus", 
                      "Bryum_pseudotriquetrum"))) 

check_congruency <- name.check(tree3, updated_species_dataset)
check_congruency 
#seems ok now (tree and data match now)
#let's root the tree and make it dichotomous 
mytree=multi2di(tree3)
##WRITING THE TREE 
write.tree(mytree, file = 'full_phylogeny_tree.tre')

#### Visualizing full phylogeny  ####
p <- ggtree(mytree)  + theme_tree() + 
  geom_tiplab( 
    color = 'black',                      
    offset = 5,
    size = 2) + 
  hexpand(-2, direction = -1)

p



#### The Non-cultivar tree ### 
pruned_tree <- read.tree('pruned_non_domestic_trees.tre')

p3 <- ggtree(pruned_tree)  + theme_tree() + 
  geom_tiplab( 
    color = 'black',                      
    offset = 5,
    size = 3) + 
  hexpand(-10, direction = -1) 

p3

pruned_tree$Nnode
Ntip(pruned_tree)

