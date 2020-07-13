# load libraries
library(magrittr)
library(readr)
library(tidyr)
library(dplyr)
library(divDyn)

# using the coral data from pbdb as divdyn does not have higher taxonomic rank
###----
dino_pbdb <- read.csv("pbdb_dino.csv")

table(dino_pbdb$accepted_rank)

# clean it a bit
dino_pbdb <- subset(dino_pbdb, dino_pbdb$accepted_rank != "subgenus")

toOmit <- omit(dino_pbdb, bin="stg", tax="genus", om="ref", ref="reference_no")
dino_pbdb <- dino_pbdb[!toOmit,]

# subset the genera with brackets and assign them to the overarching genus
test <- gsub("\\s*\\([^\\)]+\\)","",as.character(dino_pbdb$genus))
dino_pbdb$genus <- test

# using basic bin lengths to change it to range data
dinos <- fadlad(dino_pbdb, tax="genus", age =c("max_ma", "min_ma"))

# transform values back to positive
# dinos$FAD <- -dinos$FAD
# dinos$LAD <- -dinos$LAD

dinos$genus <- rownames(dinos)

# merge them together to get taxonomy
dinos <- merge(dinos,dino_pbdb, by  = "genus")

# order them properly, but we remove the order as ornithischia got no order
dinos <- dinos %>%
  select(phylum, class, family, genus, FAD, LAD)  

# remove the duplicates
dinos <- dinos[!duplicated(dinos),]

# remove the NA's, but
dinos <- subset(dinos, dinos$class != "NO_CLASS_SPECIFIED" &
                  dinos$family !=  "NO_FAMILY_SPECIFIED")

dinos <- na.omit(dinos) # no NA's anymore, so good to go

#dinos$dur <- dinos$FAD - dinos$LAD
#mean(dinos$dur) #11.6516
#median(dinos$dur) #5.4165

#mean(corals$dur) #49.66307
#median(corals$dur) #35.5
# import data 
data(stages)
gradstein <- stages
rm(stages)
# as I have assigned the stage age from Gradstein 2012, we need to assign the same ages to the stage
# file
gradstein$bottom[c(21,50:55, 57:58, 71:74, 86:87, 90, 92, 95)] <- c(443.8, 259.8, 254.2, 252.2, 250.0,247.1, 241.5, 228.4, 209.5, 
                                                                    139.4, 133.9, 130.8, 126.3, 41.2, 37.8, 23.03, 11.63, 0.0118)

gradstein$top[21:length(gradstein$bottom)-1] <- gradstein$bottom[21:length(gradstein$bottom)]

# change colnames
colnames(dinos) <- c("Phylum", "Class", "Family", "Genus", "FAD.age", "LAD.age")

# prepare weizer & prokoph temperature data
isotemp <- read.csv(file="TimeSeriesUsed.csv", header=TRUE,row.names = 1) 

#assign age to isotope data
isotemp$age <- gradstein$mid[14:94]

#mean.dot<-ggplot(data=isotemp, aes(x=age, y=Temp))+
#  geom_line(colour="#65a3a4", size=1.5) + theme_classic()

# bin FAD and LAD of each genus  
dinos$FAD_bin <- NA
dinos$LAD_bin <- NA

for (i in 1:length(dinos$Genus)) {
  for (j in 1:length(gradstein$stg)) {
    ifelse(gradstein$bottom[j] >= dinos$FAD.age[i] & gradstein$top[j] < dinos$FAD.age[i],
           dinos$FAD_bin[i] <- gradstein$stg[j], NA)
    ifelse(gradstein$bottom[j] > dinos$LAD.age[i] & gradstein$top[j] <= dinos$LAD.age[i],
           dinos$LAD_bin[i] <- gradstein$stg[j], NA)
  }}

#remove singletons as they add noise 
#foram_sing <- subset(forams_nms, forams_nms$FAD_bin == forams_nms$LAD_bin) # 199 singletons
dinos <- subset(dinos, dinos$FAD_bin != dinos$LAD_bin)

# when forams range to the recent (bin 95), shift this bin to 94 to remove pull of the recent and 
# to enable temp-calculation
dinos$LAD_bin[dinos$LAD_bin==95] <- 94
dinos <- subset(dinos, dinos$FAD_bin != dinos$LAD_bin)

# subset to the same bins as we currently have for isotemp (>= 14)
dinos <- subset(dinos, dinos$FAD_bin >=14)

#write.table(dinos, file = "dinosauria.csv")

#detach("package:plyr", unload=TRUE)
# Calculate average Temperature and standard deviation per bin
# From 360 My regular time bins for every 0.5 mln years
# Add two columns with previous and next short-term change in Temperature, but take the lm
# as we do for the long term trend as well. We need to add (-) to get the real trend
# this time we will set the beginning long time trend 1 to the origination of the genus
# and long trend 2 to origination of the family
isotemp2<- isotemp %>%
  subset(age<=481.55000) %>%
  group_by(Stage) %>%
  summarise(Temp.mean = mean(Temp), n = n()) %>%
  arrange(desc(Stage))
isotemp2$change.prev <- NA
isotemp2$age <- gradstein$mid[94:14]


# do this one for the occurrence data??
for (i in unique(isotemp$Stage)) {
  dum1 <- filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-1, i)])
  dum2 <-  lm(Temp ~ age, data = dum1)
  isotemp2[isotemp2$Stage==i, "change.prev"] <- -dum2$coefficients[2]
}

# Equal the number of temperature time bins to occurrence time bins
isotemp2 %<>% subset(Stage>13 & Stage<95) %>%
  transform(bins = as.factor(Stage))  # We need bins as factor later on when binding with foram data

rm(dum1, dum2, i, j, dino_pbdb)

###----
# now calculate four trends 
###----
# build dummy for calculation
dumbo <- dinos
dinos_safe <- dinos[,c(4,7,8)]

# add bins with 0 and fill in with 1 for extinction and origination
namevector <- as.character(c(1:94))
dinos_safe[ , namevector] <- NA
dinos_safe <- dinos_safe[, c(4:length(dinos_safe[0,]))]

for (i in 1:length(dinos_safe[,1])) {
  dinos_safe[i,dumbo[i,7]] <- 0 
  dinos_safe[i,dumbo[i,8]] <- 1
  ifelse(dumbo[i,7]!= dumbo[i,8]-1 & dumbo[i,7]!= dumbo[i,8],
         dinos_safe[i,(dumbo[i,7]+1):(dumbo[i,8]-1)] <- 0, NA)
  ifelse(dumbo[i,8] == 94, dinos_safe[i,dumbo[i,8]] <- 0, 
         dinos_safe[i,dumbo[i,8]] <- 1)
}


# bind it again 
rownames(dinos_safe) <-  make.names(dumbo[,4], unique=TRUE)


# Reorganise ranges through time data so we can bind it to temperature data
dinos_safe %<>%  as_tibble() %>%
  mutate(Genus = rownames(dinos_safe)) %>%
  group_by(Genus) %>%
  gather(-Genus, key="bins", value="extinct", na.rm = T, factor_key = T) %>%
  arrange(desc(bins), .by_group = T) 


# Now bind the two
temp_dinos <- full_join(dinos_safe, isotemp2)

#add taxonomical levels
dumbo <- dumbo[, 2:4]

# Now bind the two
temp_dinos <- full_join(temp_dinos, dumbo)

temp_dinos <- temp_dinos %>%
  select(Class, Family, Genus, bins, age, extinct, Temp.mean, change.prev) %>%
  na.omit()



# Calculate long term changes in temp (1 stage to 10 stages)
# Use lm() to determine slope of the long-term  temperature trends
# If we want to investigate the interaction between long term change and short term,
# we need to exclude the short term bin from the long term (see Manuel's figures),
# and thus shift the long term lm result up by +1. This will be messy, I'm sorry
temp_dinos$bins <- as.numeric(temp_dinos$bins)

ori_bin <- temp_dinos %>%
  select(Genus, bins, age) %>%
  group_by(Genus) %>%
  dplyr::slice(which.min(bins)) 
colnames(ori_bin) <- c("Genus", "ori.bin", "ori.age" )

# Now bind the two
test2 <- full_join(temp_dinos, ori_bin)

# Set bins back to numeric
test2 %<>% transform(bins = as.numeric(as.character(bins)))
test2$trend1 <- NA
test2$trend2 <- NA
test2$trend3 <- NA

for (i in 1:length(test2$Genus)) {
  bin1 <- c(test2$ori.bin[i], ifelse(test2$bins[i]>test2$ori.bin[i], 
                                     test2$bins[i]-1, test2$bins[i]))
  age1 <- c(isotemp$age[isotemp$Stage %in% bin1[1]:bin1[2]])
  
  temp1 <-c(isotemp$Temp[isotemp$Stage %in% bin1[1]:bin1[2]])
  lin1 <- lm(temp1~age1)
  test2$trend1[i] <- -lin1$coefficients[2]
}



# replace the NA's with zeros as we have no trend there
#test2$trend1[is.na(test2$trend1)] <- 0

#remov ori age and bin
test2 <- test2[,-c(9:10)]

#now do the same for trend2
# first build a data frame where we save the origin of each family

dino_fam <- unique(dinos$Family)
fam_dat <- data.frame(matrix(ncol =2 , nrow = length(dino_fam)))
colnames(fam_dat) <- c("age.ori", "bin.ori")
rownames(fam_dat) <-  c(as.character(dino_fam))


for (i in 1:length(fam_dat$age.ori)) {
  fam_dat[i, 2] <- min(dinos$FAD_bin[dinos$Family == dino_fam[i]])
  fam_dat[i, 1] <- isotemp$age[isotemp$Stage == fam_dat[i, 2]]
}



# add a colum with the family names for binding
fam_dat$Family <- rownames(fam_dat)

# Now bind the two
test3 <- full_join(test2, fam_dat)

# now calculate the trend starting at origin of the specific family
for (i in 1:length(test3$Genus)) {
  bin1 <- c(test3$bin.ori[i], ifelse(test3$bins[i]>test3$bin.ori[i], 
                                     test3$bins[i]-1, test3$bins[i]))
  age1 <- c(isotemp$age[isotemp$Stage %in% bin1[1]:bin1[2]])
  
  temp1 <- c(isotemp$Temp[isotemp$Stage %in% bin1[1]:bin1[2]])
  lin1 <- lm(temp1~age1)
  test3$trend2[i] <- -lin1$coefficients[2]
}


# replace the NA's with zeros as we have no trend there
#test3$trend2[is.na(test3$trend2)] <- 0

#remov ori age and bin
test3 <- test3[,-c(12:13)]


#now do the same for trend3 class
# first build a data frame where we save the origin of each family

dino_cla <- unique(dinos$Class)
cla_dat <- data.frame(matrix(ncol =2 , nrow = length(dino_cla)))
colnames(cla_dat) <- c("age.ori", "bin.ori")
rownames(cla_dat) <-  c(as.character(dino_cla))


for (i in 1:length(cla_dat$age.ori)) {
  cla_dat[i, 2] <- min(dinos$FAD_bin[dinos$Class == dino_cla[i]])
  cla_dat[i, 1] <- isotemp$age[isotemp$Stage == cla_dat[i, 2]]
}



# add a colum with the superfamily names for binding
cla_dat$Class <- rownames(cla_dat)

# Now bind the two
test4 <- full_join(test3, cla_dat)

# now calculate the trend starting at origin of the specific Order
for (i in 1:length(test4$Genus)) {
  bin1 <- c(test4$bin.ori[i], ifelse(test4$bins[i]>test4$bin.ori[i], 
                                     test4$bins[i]-1, test4$bins[i]))
  age1 <- c(isotemp$age[isotemp$Stage %in% bin1[1]:bin1[2]])
  
  temp1 <- c(isotemp$Temp[isotemp$Stage %in% bin1[1]:bin1[2]])
  lin1 <- lm(temp1~age1)
  test4$trend3[i] <- -lin1$coefficients[2]
}

# replace the NA's with zeros as we have no trend there
#test4$trend3[is.na(test4$trend3)] <- 0


#remov ori age and bin
three_trends_dinos <- test4[,-c(12:13)]


# Set bins back to numeric
three_trends_dinos %<>% transform(bins = as.numeric(as.character(bins)))

# Cleanup
rm(dumbo, fam_dat, dinos, lin1, 
   ori_bin,  test2, test3, test4, age1, bin1, 
   i, dino_fam, namevector, temp1, temp_dinos, dinos_safe, toOmit, test, dino_cla, cla_dat)

#### -----
# Calculate long term changes in temp (1 stage to 10 stages)
#### ----
# Use lm() to determine slope of the long-term  temperature trends
# If we want to investigate the interaction between long term change and short term,
# we need to exclude the short term bin from the long term (see Manuel's figures),
# and thus shift the long term lm result up by +1.
for (i in unique(isotemp$Stage)) {
  sub1<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-1, i)])
  lin1<-lm(Temp~age, data=sub1)
  isotemp2[isotemp2$Stage==i+1, "trend.st1"]<- -lin1$coefficients[2]
  sub2<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-2, i)])
  lin2<-lm(Temp~age, data=sub2)
  isotemp2[isotemp2$Stage==i+1, "trend.st2"]<- -lin2$coefficients[2]
  sub3<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-3, i)])
  lin3<-lm(Temp~age, data=sub3)
  isotemp2[isotemp2$Stage==i+1, "trend.st3"]<- -lin3$coefficients[2]
  sub4<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-4, i)])
  lin4<-lm(Temp~age, data=sub4)
  isotemp2[isotemp2$Stage==i+1, "trend.st4"]<- -lin4$coefficients[2]
  sub5<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-5, i)])
  lin5<-lm(Temp~age, data=sub5)
  isotemp2[isotemp2$Stage==i+1, "trend.st5"]<- -lin5$coefficients[2]
  sub6<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-6, i)])
  lin6<-lm(Temp~age, data=sub6)
  isotemp2[isotemp2$Stage==i+1, "trend.st6"]<- -lin6$coefficients[2]
  sub7<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-7, i)])
  lin7<-lm(Temp~age, data=sub7)
  isotemp2[isotemp2$Stage==i+1, "trend.st7"]<- -lin7$coefficients[2]
  sub8<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-8, i)])
  lin8<-lm(Temp~age, data=sub8)
  isotemp2[isotemp2$Stage==i+1, "trend.st8"]<- -lin8$coefficients[2]
  sub9<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-9, i)])
  lin9<-lm(Temp~age, data=sub9)
  isotemp2[isotemp2$Stage==i+1, "trend.st9"]<- -lin9$coefficients[2]
  sub10<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-10, i)])
  lin10<-lm(Temp~age, data=sub10)
  isotemp2[isotemp2$Stage==i+1, "trend.st10"]<- -lin10$coefficients[2]
}

# Cleanup
rm(i, isotemp, sub1, lin1, sub2, lin2, sub3, lin3, sub4, lin4, sub5, lin5, 
   sub6, lin6, sub7, lin7, sub8, lin8, sub9, lin9, sub10, lin10)

# Set bins back to numeric
three_trends_dinos %<>% transform(bins = as.factor(as.character(bins)))

# Now bind the two
thirteen_trends_dinos <- full_join(three_trends_dinos, isotemp2)


# Set bins back to numeric
thirteen_trends_dinos %<>% transform(bins = as.numeric(as.character(bins)))


# remove redundant colums
thirteen_trends_dinos <- thirteen_trends_dinos[1:1604,]

rm(gradstein, isotemp2)


