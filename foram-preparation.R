library(tidyverse)
library(divDyn)
library(here)

# read in data
forams <- read_csv(here("data/raw-fossil-data/foraminifera.csv")) 

# data cleaning
dat_range <- forams %>% 
  # remove singletons as they add noise 
  filter(fad != lad) %>% 
  # remove genera that range to the recent (bin 95)
  filter(lad != 0) 

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
  full_join(dat_range %>% select(order, superfamily, family, genus)) %>% 
  # remove redundant colums from merging
  drop_na(genus) %>% 
  # order it properly
  select(order, superfamily, family, genus, stage, age, extinction, temp, 
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
# write_csv(dat_trends, path = here("data/foraminifera_trends.csv"))

