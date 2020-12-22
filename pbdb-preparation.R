library(tidyverse)
library(divDyn)
library(here)



# functions ---------------------------------------------------------------


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


# make function with dataframe as argument
make_trends <- function(dfr) {
  
  dat_occ <- dfr %>% 
    # remove the NA's
    filter(order != "NO_ORDER_SPECIFIED", 
           family != "NO_FAMILY_SPECIFIED")
  
  # save higher taxonomy
  tax <- dat_occ %>% 
    select(phylum, class, order, family, genus) %>% 
    unique()
  
  # transform to range data
  dat_range <- fadlad(dat_occ, tax="genus", age= c("max_ma", "min_ma"))
  
  
  dat_range <- dat_range %>% 
    # remove singletons as they add noise 
    filter(FAD != LAD) %>% 
    # remove genera that range to the recent (bin 95)
    filter(LAD != 0) %>% 
    # assign genus names to a column
    add_column(genus = rownames(.)) %>% 
    # merge them together to get higher taxonomy
    left_join(tax, by = "genus") %>% 
    # order colums properly and make consistent column names
    select(phylum, class, order, family, genus, fad = FAD, lad = LAD)  
  
  
  ### bin data 
  
  
  # import stage data from Gradstein 2012
  gradstein <- read_csv(here("data/gradstein.csv"))
  
  # bin FAD and LAD of each genus to stages  
  dat_range$fad_bin <- NA
  dat_range$lad_bin <- NA
  
  for (i in 1:length(dat_range$genus)) {
    for (j in 1:length(gradstein$stg)) {
      ifelse(gradstein$bottom[j] >= dat_range$fad[i] & gradstein$top[j] < dat_range$fad[i],
             dat_range$fad_bin[i] <- gradstein$stg[j], NA)
      ifelse(gradstein$bottom[j] > dat_range$lad[i] & gradstein$top[j] <= dat_range$lad[i],
             dat_range$lad_bin[i] <- gradstein$stg[j], NA)
    }}
  
  # clean data
  dat_range <- dat_range %>% 
    # remove singletons as they add noise 
    filter(fad_bin != lad_bin) %>% 
    # remove extant genera (genera that range to bin 95)
    filter(lad_bin != 95) %>% 
    # subset to the same bins as we currently have for isotope data (isotemp) (>= 14)
    filter(fad_bin >= 14)
  
  
  #### Calculate multiple long-term temperature trends
  
  # read in preprocessed isotope data
  isotemp2 <- read_csv(here("data/isotemp_trends.csv"))
  
  # data for calculation
  dat_trends <- dat_range %>% 
    # subset
    select(genus, fad, lad) 
  
  # add bins with 0 
  namevector <- as.character(c(1:94))
  dat_trends[ , namevector] <- NA
  dat_trends <- dat_trends[, c(4:length(dat_trends[0,]))]
  
  # fill in with 0 for survival and 1 extinction
  for (i in 1:length(dat_trends[,1])) {
    dat_trends[i,dat_range[i,"fad_bin"]] <- 0 
    dat_trends[i,dat_range[i,"lad_bin"]] <- 1
    ifelse(dat_range[i,"fad_bin"]!= dat_range[i,"lad_bin"]-1 & 
             dat_range[i,"fad_bin"]!= dat_range[i,"lad_bin"],
           dat_trends[i,(dat_range[i,"fad_bin"]+1):(dat_range[i,"lad_bin"]-1)] <- 0,
           NA)
  }
  
  
  
  
  # Reorganise ranges through time data so we can bind it to temperature data
  dat_trends <- dat_trends %>% 
    as_tibble() %>%
    mutate(genus = dat_range$genus) %>% 
    group_by(genus) %>% 
    gather(-genus, key="bins", value="extinction", na.rm = T, factor_key = F) %>%
    arrange(desc(bins), .by_group = T) %>% 
    rename(stage = bins) %>% 
    mutate(stage = as.numeric(stage)) %>% 
    ungroup()
  
  
  
  # Now bind the two
  dat_trends <- full_join(dat_trends, isotemp2) %>% 
    # add higher taxonomy
    full_join(dat_range %>% select(phylum, class, order, family, genus)) %>% 
    # remove redundant colums from merging
    drop_na(stage) %>% 
    # order it properly
    select(phylum, class, order, family, genus, stage, age, extinction, temp, 
           change_prev, trend.st1:trend.st10)
  
  #### calculate trends based on taxonomy
  
  # for genera
  ori_bin <- dat_trends %>%
    select(genus, stage, age) %>%
    group_by(genus) %>%
    dplyr::slice(which.min(stage)) %>% 
    rename(ori.bin = stage, ori.age = age) %>% 
    ungroup()
  
  # bind the two
  dat_trends <- full_join(dat_trends, ori_bin) %>% 
    # add empty column for genus trends
    add_column(trend_genus = as.numeric(NA)) 
  
  # calculate the trend starting at origin of the focal family
  for (i in 1:length(dat_trends$genus)) {
    bin <- c(dat_trends$ori.bin[i], ifelse(dat_trends$stage[i]>dat_trends$ori.bin[i], 
                                           dat_trends$stage[i]-1, dat_trends$stage[i]))
    age <- c(isotemp2$age[isotemp2$stage %in% bin[1]:bin[2]])
    
    temp <- c(isotemp2$temp[isotemp2$stage %in% bin[1]:bin[2]])
    lin <- lm(temp~age)
    dat_trends$trend_genus[i] <- -lin$coefficients[2]
  }
  
  #remov ori age and bin
  dat_trends <- dat_trends %>% 
    select(-c(ori.bin, ori.age))
  
  # same for family
  ori_bin <- dat_trends %>%
    select(family, stage, age) %>%
    group_by(family) %>%
    dplyr::slice(which.min(stage)) %>% 
    rename(ori.bin = stage, ori.age = age) %>% 
    ungroup()
  
  # bind the two
  dat_trends <- full_join(dat_trends, ori_bin) %>% 
    # add empty column for genus trends
    add_column(trend_family = as.numeric(NA)) 
  
  # calculate the trend starting at origin of the focal family
  for (i in 1:length(dat_trends$family)) {
    bin <- c(dat_trends$ori.bin[i], ifelse(dat_trends$stage[i]>dat_trends$ori.bin[i], 
                                           dat_trends$stage[i]-1, dat_trends$stage[i]))
    age <- c(isotemp2$age[isotemp2$stage %in% bin[1]:bin[2]])
    
    temp <- c(isotemp2$temp[isotemp2$stage %in% bin[1]:bin[2]])
    lin <- lm(temp~age)
    dat_trends$trend_family[i] <- -lin$coefficients[2]
  }
  
  # remov ori age and bin
  dat_trends <- dat_trends %>% 
    select(-c(ori.bin, ori.age))
  
  # same for order
  ori_bin <- dat_trends %>%
    select(order, stage, age) %>%
    group_by(order) %>%
    dplyr::slice(which.min(stage)) %>% 
    rename(ori.bin = stage, ori.age = age) %>% 
    ungroup()
  
  # bind the two
  dat_trends <- full_join(dat_trends, ori_bin) %>% 
    # add empty column for genus trends
    add_column(trend_order = as.numeric(NA)) 
  
  # calculate the trend starting at origin of the focal family
  for (i in 1:length(dat_trends$order)) {
    bin <- c(dat_trends$ori.bin[i], ifelse(dat_trends$stage[i]>dat_trends$ori.bin[i], 
                                           dat_trends$stage[i]-1, dat_trends$stage[i]))
    age <- c(isotemp2$age[isotemp2$stage %in% bin[1]:bin[2]])
    
    temp <- c(isotemp2$temp[isotemp2$stage %in% bin[1]:bin[2]])
    lin <- lm(temp~age)
    dat_trends$trend_order[i] <- -lin$coefficients[2]
  }
  
  # remov ori age and bin
  dat_trends <- dat_trends %>% 
    select(-c(ori.bin, ori.age))
  
  dat_trends
}


# now we can access data via the pbdb api using the get_pbdb_url function and 
# then get a data frame based on this raw data with cleaned range data and 
# calculated long-term trends

# as we iterate this process for each fossil group, it is better to make a 
# function out of these processing steps: 

process_data <-  function(taxon) {
  # get url 
  uri <- get_pbdb_url(taxon = taxon)
  
  # download data
  pbdb_data <- read.delim(file=uri, quote="") 
  
  # get long-term trends
  trend_data <- make_trends(dfr = pbdb_data)
  
  output_list <- list(pbdb = pbdb_data, trend = trend_data)
  output_list
}


# download and clean data -------------------------------------------------


### arthropoda
arthropoda <- process_data(taxon = "arthropoda")

# save raw pbdb data once 
# write_csv(arthropoda$pbdb, path = here("data/raw-fossil-data/arthropoda.csv"))
# save trend data
# write_csv(arthropoda$trend, path = here("data/arthropoda_trends.csv"))


### bivalvia
bivalvia <- process_data(taxon = "bivalvia")

# save raw pbdb data once 
# write_csv(bivalvia$pbdb, path = here("data/raw-fossil-data/bivalvia.csv"))
# save trend data
# write_csv(bivalvia$trend, path = here("data/bivalvia_trends.csv"))


### cnidaria
cnidaria <- process_data(taxon = "cnidaria")

# save raw pbdb data once 
# write_csv(cnidaria$pbdb, path = here("data/raw-fossil-data/cnidaria.csv"))
# save trend data
# write_csv(cnidaria$trend, path = here("data/cnidaria_trends.csv"))


### echinodermata
echinodermata <- process_data(taxon = "echinodermata")

# save raw pbdb data once 
# write_csv(echinodermata$pbdb, path = here("data/raw-fossil-data/echinodermata.csv"))
# save trend data
# write_csv(echinodermata$trend, path = here("data/echinodermata_trends.csv"))


### gastropoda
gastropoda <- process_data(taxon = "gastropoda")

# save raw pbdb data once 
# write_csv(gastropoda$pbdb, path = here("data/raw-fossil-data/gastropoda.csv"))
# save trend data
# write_csv(gastropoda$trend, path = here("data/gastropoda_trends.csv"))


### reptilia
reptilia <- process_data(taxon = "reptilia")

# save raw pbdb data once 
# write_csv(reptilia$pbdb, path = here("data/raw-fossil-data/reptilia.csv"))
# save trend data
# write_csv(reptilia$trend, path = here("data/reptilia_trends.csv"))
