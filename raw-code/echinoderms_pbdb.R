
# keep in mind that you have only isotope data starting at the ordovician!
#remove workspace
rm(list = ls())

library(magrittr)
library(readr)
library(tidyr)
library(dplyr)
library(divDyn)

# using the coral data from pbdb as divdyn does not have higher taxonomic rank
###----
echinoderms_pbdb <- read.csv("pbdb_echinoderms.csv")

table(echinoderms_pbdb$accepted_rank)

# clean it a bit
echinoderms_pbdb <- subset(echinoderms_pbdb, echinoderms_pbdb$accepted_rank != "subgenus")

toOmit <- omit(echinoderms_pbdb, bin="stg", tax="genus", om="ref", ref="reference_no")
echinoderms_pbdb <- echinoderms_pbdb[!toOmit,]

# subset the genera with brackets and assign them to the overarching genus
test <- gsub("\\s*\\([^\\)]+\\)","",as.character(echinoderms_pbdb$genus))
echinoderms_pbdb$genus <- test

# using basic bin lengths to change it to range data
echinoderms <- fadlad(echinoderms_pbdb, tax="genus", age =c("max_ma", "min_ma"))

# transform values back to positive
# echinoderms$FAD <- -echinoderms$FAD
# echinoderms$LAD <- -echinoderms$LAD

echinoderms$genus <- rownames(echinoderms)


# merge them together to get taxonomy
echinoderms <- merge(echinoderms,echinoderms_pbdb, by  = "genus")

# order them properly
echinoderms <- echinoderms %>%
  select(phylum, class, order, family, genus, FAD, LAD)  

# remove the duplicates
echinoderms <- echinoderms[!duplicated(echinoderms),]

# remove the NA's
echinoderms <- subset(echinoderms, echinoderms$class != "NO_CLASS_SPECIFIED" &
                        echinoderms$order != "NO_ORDER_SPECIFIED" &
                        echinoderms$family !=  "NO_FAMILY_SPECIFIED")

echinoderms <- na.omit(echinoderms) # no NA's anymore, so good to go

#echinoderms$dur <- echinoderms$FAD - echinoderms$LAD
#mean(echinoderms$dur) #47.2641
#median(echinoderms$dur) #27.4


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
colnames(echinoderms) <- c("Phylum", "Class", "Order", "Family", "Genus", "FAD.age", "LAD.age")

# prepare weizer & prokoph temperature data
isotemp <- read.csv(file="TimeSeriesUsed.csv", header=TRUE,row.names = 1) 

#assign age to isotope data
isotemp$age <- gradstein$mid[14:94]

#mean.dot<-ggplot(data=isotemp, aes(x=age, y=Temp))+
#  geom_line(colour="#65a3a4", size=1.5) + theme_classic()

# bin FAD and LAD of each genus  
echinoderms$FAD_bin <- NA
echinoderms$LAD_bin <- NA

for (i in 1:length(echinoderms$Genus)) {
  for (j in 1:length(gradstein$stg)) {
    ifelse(gradstein$bottom[j] >= echinoderms$FAD.age[i] & gradstein$top[j] < echinoderms$FAD.age[i],
           echinoderms$FAD_bin[i] <- gradstein$stg[j], NA)
    ifelse(gradstein$bottom[j] > echinoderms$LAD.age[i] & gradstein$top[j] <= echinoderms$LAD.age[i],
           echinoderms$LAD_bin[i] <- gradstein$stg[j], NA)
  }}

#remove singletons as they add noise 
#foram_sing <- subset(forams_nms, forams_nms$FAD_bin == forams_nms$LAD_bin) # 199 singletons
echinoderms <- subset(echinoderms, echinoderms$FAD_bin != echinoderms$LAD_bin)

# when forams range to the recent (bin 95), shift this bin to 94 to remove pull of the recent and 
# to enable temp-calculation
echinoderms$LAD_bin[echinoderms$LAD_bin==95] <- 94
echinoderms <- subset(echinoderms, echinoderms$FAD_bin != echinoderms$LAD_bin)

# subset to the same bins as we currently have for isotemp (>= 14)
echinoderms <- subset(echinoderms, echinoderms$FAD_bin >=14)

#write.table(echinoderms, file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/clean_data/echinodermata.csv")

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

rm(dum1, dum2, i, j, echinoderms_pbdb, test, toOmit)

###----
# now calculate four trends 
###----
# build dummy for calculation
dumbo <- echinoderms
echinoderms_safe <- echinoderms[,c(5,8,9)]

# add bins with 0 and fill in with 1 for extinction and origination
namevector <- as.character(c(1:94))
echinoderms_safe[ , namevector] <- NA
echinoderms_safe <- echinoderms_safe[, c(4:length(echinoderms_safe[0,]))]

for (i in 1:length(echinoderms_safe[,1])) {
  echinoderms_safe[i,dumbo[i,8]] <- 0 
  echinoderms_safe[i,dumbo[i,9]] <- 1
  ifelse(dumbo[i,8]!= dumbo[i,9]-1 & dumbo[i,8]!= dumbo[i,9],
         echinoderms_safe[i,(dumbo[i,8]+1):(dumbo[i,9]-1)] <- 0, NA)
  ifelse(dumbo[i,9] == 94, echinoderms_safe[i,dumbo[i,9]] <- 0, 
         echinoderms_safe[i,dumbo[i,9]] <- 1)
}


# bind it again 
rownames(echinoderms_safe) <-  make.names(dumbo[,5], unique=TRUE)


# Reorganise ranges through time data so we can bind it to temperature data
echinoderms_safe %<>%  as_tibble() %>%
  mutate(Genus = rownames(echinoderms_safe)) %>%
  group_by(Genus) %>%
  gather(-Genus, key="bins", value="extinct", na.rm = T, factor_key = T) %>%
  arrange(desc(bins), .by_group = T) 


# Now bind the two
temp_echinoderms <- full_join(echinoderms_safe, isotemp2)

#add taxonomical levels
dumbo <- dumbo[, 2:5]

# Now bind the two
temp_echinoderms <- full_join(temp_echinoderms, dumbo)

temp_echinoderms <- temp_echinoderms %>%
  select(Class, Order, Family, Genus, bins, age, extinct, Temp.mean, change.prev) %>%
  na.omit()



# Calculate long term changes in temp (1 stage to 10 stages)
# Use lm() to determine slope of the long-term  temperature trends
# If we want to investigate the interaction between long term change and short term,
# we need to exclude the short term bin from the long term (see Manuel's figures),
# and thus shift the long term lm result up by +1. This will be messy, I'm sorry
temp_echinoderms$bins <- as.numeric(temp_echinoderms$bins)
ori_bin <- temp_echinoderms %>%
  select(Genus, bins, age) %>%
  group_by(Genus) %>%
  dplyr::slice(which.min(bins)) 
colnames(ori_bin) <- c("Genus", "ori.bin", "ori.age" )

# Now bind the two
test2 <- full_join(temp_echinoderms, ori_bin)

# Set bins back to numeric
test2 %<>% transform(bins = as.numeric(as.character(bins)))
test2$trend1 <- NA
test2$trend2 <- NA
test2$trend3 <- NA
test2$trend4 <- NA

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
test2 <- test2[,-c(10:11)]

#now do the same for trend2
# first build a data frame where we save the origin of each family

ech_fam <- unique(echinoderms$Family)
fam_dat <- data.frame(matrix(ncol =2 , nrow = length(ech_fam)))
colnames(fam_dat) <- c("age.ori", "bin.ori")
rownames(fam_dat) <-  c(as.character(ech_fam))


for (i in 1:length(fam_dat$age.ori)) {
  fam_dat[i, 2] <- min(echinoderms$FAD_bin[echinoderms$Family == ech_fam[i]])
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
test3 <- test3[,-c(14:15)]

#now do the same for trend3 Order
# first build a data frame where we save the origin of each family

ech_ord <- unique(echinoderms$Order)
ord_dat <- data.frame(matrix(ncol =2 , nrow = length(ech_ord)))
colnames(ord_dat) <- c("age.ori", "bin.ori")
rownames(ord_dat) <-  c(as.character(ech_ord))


for (i in 1:length(ord_dat$age.ori)) {
  ord_dat[i, 2] <- min(echinoderms$FAD_bin[echinoderms$Order == ech_ord[i]])
  ord_dat[i, 1] <- isotemp$age[isotemp$Stage == ord_dat[i, 2]]
}



# add a colum with the superfamily names for binding
ord_dat$Order <- rownames(ord_dat)

# Now bind the two
test4 <- full_join(test3, ord_dat)

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
test4 <- test4[,-c(14:15)]

#now do the same for trend4 order
# first build a data frame where we save the origin of each family

ech_cla <- unique(echinoderms$Class)
cla_dat <- data.frame(matrix(ncol =2 , nrow = length(ech_cla)))
colnames(cla_dat) <- c("age.ori", "bin.ori")
rownames(cla_dat) <-  c(as.character(ech_cla))


for (i in 1:length(cla_dat$age.ori)) {
  cla_dat[i, 2] <- min(echinoderms$FAD_bin[echinoderms$Class == ech_cla[i]])
  cla_dat[i, 1] <- isotemp$age[isotemp$Stage == cla_dat[i, 2]]
}



# add a colum with the superfamily names for binding
cla_dat$Class <- rownames(cla_dat)

# Now bind the two
test5 <- full_join(test4, cla_dat)

# now calculate the trend starting at origin of the specific family
for (i in 1:length(test5$Genus)) {
  bin1 <- c(test5$bin.ori[i], ifelse(test5$bins[i]>test5$bin.ori[i], 
                                     test5$bins[i]-1, test5$bins[i]))
  age1 <- c(isotemp$age[isotemp$Stage %in% bin1[1]:bin1[2]])
  
  temp1 <-c(isotemp$Temp[isotemp$Stage %in% bin1[1]:bin1[2]])
  lin1 <- lm(temp1~age1)
  test5$trend4[i] <- -lin1$coefficients[2]
}


#remov ori age and bin
four_trends_echinoderms <- test5[,-c(14:15)]

# Set bins back to numeric
four_trends_echinoderms %<>% transform(bins = as.numeric(as.character(bins)))

# Cleanup
rm(dumbo, ech_cla, echinoderms, lin1, 
   ord_dat, ori_bin,  test2, test3, test4, test5, age1, bin1, 
   i, ech_fam, ech_ord, namevector, temp1, temp_echinoderms, echinoderms_safe, cla_dat, fam_dat)

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
four_trends_echinoderms %<>% transform(bins = as.factor(as.character(bins)))

# Now bind the two
fourteen_trends_echinoderms <- full_join(four_trends_echinoderms, isotemp2)


# Set bins back to numeric
fourteen_trends_echinoderms %<>% transform(bins = as.numeric(as.character(bins)))

# remove redundant colums
fourteen_trends_echinoderms <- fourteen_trends_echinoderms[-4835,]

rm(gradstein, isotemp2)
