# Resolves bad taxon names with TNRS query

# This version is for the empirical dataset.

# Last updated 4 March 2021, Josef Garen

library(tidyverse)
library(TNRS)

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
    arrange(curveID)                            # Sort by curveID to put back in original order
  
  return(data)
}


