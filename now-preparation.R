library(tidyverse)
library(divDyn)
library(here)

# read in data
mammals <- read_csv2(here("data/raw-fossil-data/mammalia.csv")) 

# clean data
mammals <- filter(mammals, TAXON_STATUS != "family attrib of genus uncertain" &
                        TAXON_STATUS != "genus attrib of species uncertain" &
                        TAXON_STATUS != "taxonomic validity uncertain" &
                        ORDER != "incertae sedis" & ORDER != "Incertaesedis" & ORDER != "Indet" &
                        ORDER != "indet." & ORDER != "PROBOSCIDEA" & FAMILY != "div.indet." &
                        FAMILY != "fam." & FAMILY != "incertae sedis" & FAMILY != "Incertae sedis" &
                        FAMILY != "Incertaesedis" & FAMILY != "indet" & FAMILY != "Indet" &
                        FAMILY != "indet." & FAMILY != "indet" & FAMILY != "Indet" &
                        SUBFAMILY != "\\N" & GENUS != "gen." & GENUS != "Gen." &
                        GENUS != "Indet" & GENUS != "indet." & GENUS != "Indet." & GENUS != "Sp.")

# clean column names
mammals <- select(mammals, order = ORDER, family = FAMILY, subfamily = SUBFAMILY,
                  genus = GENUS, max_ma = MAX_AGE, min_ma = MIN_AGE) %>% 
  # change colum types to double
  mutate(max_ma = as.double(max_ma), min_ma = as.double(min_ma))

# save higher taxonomy
tax <- mammals %>% 
  select(order, family, subfamily, genus) %>% 
  unique()

# transform to range data
dat_range <- fadlad(mammals, tax="genus", age= c("max_ma", "min_ma"))


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
  select(order, family, subfamily, genus, fad = FAD, lad = LAD)  


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
  filter(fad_bin >= 14) %>% 
  as.data.frame()


#### Calculate multiple long-term temperature trends

# read in preprocessed isotope data
isotemp2 <- read_csv(here("data/isotemp_trends.csv"))

# data for calculation
dat_trends <- dat_range %>%
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
  full_join(dat_range %>% select(order, family, subfamily, genus)) %>% 
  # remove redundant colums from merging
  drop_na(genus) %>% 
  # order it properly
  select(order, family, subfamily, genus, stage, age, extinction, temp, 
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

# save trend data
# write_csv(dat_trends, path = here("data/mammalia_trends.csv"))

