library(tidyverse)
library(here)

# define function to access the pbdb, returns an url
get_pbdb_url <- function(taxon){
  params <- paste(
    # select group
    paste("base_name=",taxon, sep = ""), 
    # only return occurrences identified to at least genus
    # level and lump multiple occurrences from same collection into a single occurrence
    "idreso=lump_genus", 
    # excludes occurrences with modifiers: 
    # aff. / cf. / ? / "" / informal / sensu lato. on genus name
    "idqual=genus_certain", 
    # only return extinct genera
    "extant=no", 
    # classext=taxonomic information, with taxon numbers;
    # ident=individual components of the taxonomic identification
    "show=classext,ident", 
    sep="&")
  
  # get url
  uri <- paste("https://paleobiodb.org/data1.2/occs/list.tsv?", params, sep="")
  
  uri
}

# get total number of species 

# all fossil groups we used in the analysis
nr_spec <- 
  list("arthropoda", "bivalvia", "cnidaria", "echinodermata", "foraminifera", "gastropoda", "reptilia") %>% 
  # for each group, get the pbdb url
  map(get_pbdb_url) %>% 
  # now use the url to download the data from the pbdb for each group
  map(read.delim, quote="") %>% 
  # select only the species name column for each group#
  map("species_name") %>% 
  # get the total number of unique species per group
  map_dbl(n_distinct) %>% 
  # simply sum all numbers of species up to get total number of species
  sum()
