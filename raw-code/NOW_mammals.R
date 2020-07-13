# load libraries
library(magrittr)
library(readr)
library(tidyr)
library(dplyr)
library(divDyn)


# using the mammal data from the NOW data base
mammals_raw <-read.csv2("NOW_mammals.csv")

# clean it
mammals_raw <- filter(mammals_raw, TAXON_STATUS != "family attrib of genus uncertain" &
                    TAXON_STATUS != "genus attrib of species uncertain" &
                    TAXON_STATUS != "taxonomic validity uncertain" &
                    ORDER != "incertae sedis" & ORDER != "Incertaesedis" & ORDER != "Indet" &
                    ORDER != "indet." & ORDER != "PROBOSCIDEA" & FAMILY != "div.indet." &
                    FAMILY != "fam." & FAMILY != "incertae sedis" & FAMILY != "Incertae sedis" &
                    FAMILY != "Incertaesedis" & FAMILY != "indet" & FAMILY != "Indet" &
                    FAMILY != "indet." & FAMILY != "indet" & FAMILY != "Indet" &
                    SUBFAMILY != "\\N" & GENUS != "gen." & GENUS != "Gen." &
                    GENUS != "Indet" & GENUS != "indet." & GENUS != "Indet." & GENUS != "Sp.")

mammals_raw <- select(mammals_raw, ORDER, FAMILY, SUBFAMILY, GENUS, MAX_AGE, MIN_AGE )

# change colnames
colnames(mammals_raw) <- c("order","family","subfamily","genus", "max_ma", "min_ma")

#mammals_raw %<>% transform(max_ma = as.numeric(as.character(max_ma)))%>%
#  transform(min_ma = as.numeric(as.character(min_ma)))

#mammals_raw$dur <- mammals_raw$max_ma - mammals_raw$min_ma

#mean(mammals_raw$dur) #1.619755
#median(mammals_raw$dur) #1.3


# change bins to numeric
mammals_raw %<>% transform(max_ma = as.numeric(as.character(max_ma)), min_ma = as.numeric(as.character(max_ma)))

# using basic bin lengths to change it to range data from divdyn
mammals <- fadlad(mammals_raw, tax="genus", age =c("max_ma","min_ma"))

# transform values back to positive
# mammals$FAD <- -mammals$FAD
# mammals$LAD <- -mammals$LAD

mammals$genus <- rownames(mammals)

# merge them together to get taxonomy
mammals_raw <- select(mammals_raw, order, family, subfamily, genus)
mammals <- merge(mammals, mammals_raw, by  = "genus", all.x = T)

# remove the duplicates
mammals <- mammals[!duplicated(mammals),]
# remove the genus duplicates
mammals <- mammals[!duplicated(mammals$genus),]

# order them properly
mammals <- mammals %>%
  select(order, family, subfamily, genus, FAD, LAD)  

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
colnames(mammals) <- c("Order", "Family","Subfamily", "Genus", "FAD.age", "LAD.age")

# prepare weizer & prokoph temperature data
isotemp <- read.csv(file="TimeSeriesUsed.csv", header=TRUE,row.names = 1) 

#assign age to isotope data
isotemp$age <- gradstein$mid[14:94]

#mean.dot<-ggplot(data=isotemp, aes(x=age, y=Temp))+
#  geom_line(colour="#65a3a4", size=1.5) + theme_classic()

# bin FAD and LAD of each genus  
mammals$FAD_bin <- NA
mammals$LAD_bin <- NA

for (i in 1:length(mammals$Genus)) {
  for (j in 1:length(gradstein$stg)) {
    ifelse(gradstein$bottom[j] >= mammals$FAD.age[i] & gradstein$top[j] < mammals$FAD.age[i],
           mammals$FAD_bin[i] <- gradstein$stg[j], NA)
    ifelse(gradstein$bottom[j] > mammals$LAD.age[i] & gradstein$top[j] <= mammals$LAD.age[i],
           mammals$LAD_bin[i] <- gradstein$stg[j], NA)
  }}

#remove singletons as they add noise 
#foram_sing <- subset(forams_nms, forams_nms$FAD_bin == forams_nms$LAD_bin) # 199 singletons
mammals <- subset(mammals, mammals$FAD_bin != mammals$LAD_bin)

# when forams range to the recent (bin 95), shift this bin to 94 to remove pull of the recent and 
# to enable temp-calculation
mammals$LAD_bin[mammals$LAD_bin==95] <- 94
mammals <- subset(mammals, mammals$FAD_bin != mammals$LAD_bin)

# subset to the same bins as we currently have for isotemp (>= 14)
mammals <- subset(mammals, mammals$FAD_bin >=14)

#write.table(mammals, file = "mammalia.csv")

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

rm(dum1, dum2, i, j)

###----
# now calculate four trends 
###----
# build dummy for calculation
dumbo <- mammals
mammals_safe <- mammals[,c(4,7,8)]

# add bins with 0 and fill in with 1 for extinction and origination
namevector <- as.character(c(1:94))
mammals_safe[ , namevector] <- NA
mammals_safe <- mammals_safe[, c(4:length(mammals_safe[0,]))]

for (i in 1:length(mammals_safe[,1])) {
  mammals_safe[i,dumbo[i,7]] <- 0 
  mammals_safe[i,dumbo[i,8]] <- 1
  ifelse(dumbo[i,7]!= dumbo[i,8]-1 & dumbo[i,7]!= dumbo[i,8],
         mammals_safe[i,(dumbo[i,7]+1):(dumbo[i,8]-1)] <- 0, NA)
  ifelse(dumbo[i,8] == 94, mammals_safe[i,dumbo[i,8]] <- 0, 
         mammals_safe[i,dumbo[i,8]] <- 1)
}


# bind it again 
rownames(mammals_safe) <-  make.names(dumbo[,4], unique=TRUE)


# Reorganise ranges through time data so we can bind it to temperature data
mammals_safe %<>%  as_tibble() %>%
  mutate(Genus = rownames(mammals_safe)) %>%
  group_by(Genus) %>%
  gather(-Genus, key="bins", value="extinct", na.rm = T, factor_key = T) %>%
  arrange(desc(bins), .by_group = T) 


# Now bind the two
temp_mammals <- full_join(mammals_safe, isotemp2)
temp_mammals <- temp_mammals[1:1877,]

#add taxonomical levels
dumbo <- dumbo[, 1:4]

# Now bind the two
temp_mammals <- full_join(temp_mammals, dumbo)

temp_mammals <- temp_mammals %>%
  select(Order, Family, Subfamily,  Genus, bins, age, extinct, Temp.mean, change.prev) %>%
  na.omit()



# Calculate long term changes in temp (1 stage to 10 stages)
# Use lm() to determine slope of the long-term  temperature trends
# If we want to investigate the interaction between long term change and short term,
# we need to exclude the short term bin from the long term (see Manuel's figures),
# and thus shift the long term lm result up by +1. This will be messy, I'm sorry
temp_mammals$bins <- as.numeric(temp_mammals$bins)

ori_bin <- temp_mammals %>%
  select(Genus, bins, age) %>%
  group_by(Genus) %>%
  dplyr::slice(which.min(bins)) 
colnames(ori_bin) <- c("Genus", "ori.bin", "ori.age" )

# Now bind the two
test2 <- full_join(temp_mammals, ori_bin)

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
# first build a data frame where we save the origin of each subfamily

for_sub <- unique(mammals$Subfamily)
sub_dat <- data.frame(matrix(ncol =2 , nrow = length(for_sub)))
colnames(sub_dat) <- c("age.ori", "bin.ori")
rownames(sub_dat) <-  c(as.character(for_sub))


for (i in 1:length(sub_dat$age.ori)) {
  sub_dat[i, 2] <- min(mammals$FAD_bin[mammals$Subfamily == for_sub[i]])
  sub_dat[i, 1] <- isotemp$age[isotemp$Stage == sub_dat[i, 2]]
}



# add a colum with the family names for binding
sub_dat$Subfamily <- rownames(sub_dat)

# Now bind the two
test3 <- full_join(test2, sub_dat)

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

#now do the same for trend3 Family
# first build a data frame where we save the origin of each family

for_fam <- unique(mammals$Family)
fam_dat <- data.frame(matrix(ncol =2 , nrow = length(for_fam)))
colnames(fam_dat) <- c("age.ori", "bin.ori")
rownames(fam_dat) <-  c(as.character(for_fam))


for (i in 1:length(fam_dat$age.ori)) {
  fam_dat[i, 2] <- min(mammals$FAD_bin[mammals$Family == for_fam[i]])
  fam_dat[i, 1] <- isotemp$age[isotemp$Stage == fam_dat[i, 2]]
}



# add a colum with the superfamily names for binding
fam_dat$Family <- rownames(fam_dat)

# Now bind the two
test4 <- full_join(test3, fam_dat)

# now calculate the trend starting at origin of the specific Order
for (i in 1:length(test4$Genus)) {
  bin1 <- c(test4$bin.ori[i], ifelse(test4$bins[i]>test4$bin.ori[i], 
                                     test4$bins[i]-1, test4$bins[i]))
  age1 <- c(isotemp$age[isotemp$Stage %in% bin1[1]:bin1[2]])
  
  temp1 <- c(isotemp$Temp[isotemp$Stage %in% bin1[1]:bin1[2]])
  lin1 <- lm(temp1~age1)
  test4$trend3[i] <- -lin1$coefficients[2]
}



#remov ori age and bin
test4 <- test4[,-c(14:15)]


#now do the same for trend4 order
# first build a data frame where we save the origin of each family

for_ord <- unique(mammals$Order)
ord_dat <- data.frame(matrix(ncol =2 , nrow = length(for_ord)))
colnames(ord_dat) <- c("age.ori", "bin.ori")
rownames(ord_dat) <-  c(as.character(for_ord))


for (i in 1:length(ord_dat$age.ori)) {
  ord_dat[i, 2] <- min(mammals$FAD_bin[mammals$Order == for_ord[i]])
  ord_dat[i, 1] <- isotemp$age[isotemp$Stage == ord_dat[i, 2]]
}



# add a colum with the superfamily names for binding
ord_dat$Order <- rownames(ord_dat)

# Now bind the two
test5 <- full_join(test4, ord_dat)

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
four_trends_mammals <- test5[,-c(14:15)]

# Set bins back to numeric
four_trends_mammals %<>% transform(bins = as.numeric(as.character(bins)))

# Cleanup
rm(dumbo, fam_dat, mammals_safe, lin1, 
   ord_dat, ori_bin, sub_dat,  test2, test3, test4, test5, age1, bin1, 
    i, for_fam, for_ord, for_sub, namevector, temp1, temp_mammals)

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

# Set bins to factor
four_trends_mammals %<>% transform(bins = as.factor(as.character(bins)))

# Now bind the two
ten_trends_mammals <- full_join(four_trends_mammals, isotemp2)
ten_trends_mammals <- ten_trends_mammals[1:1877,]

# Set bins back to numeric
ten_trends_mammals %<>% transform(bins = as.numeric(as.character(bins)))
#isotemp2 %<>% subset(Stage>25 & Stage<95) %>%
#  transform(bins = as.factor(Stage))

isotemp2 <- isotemp2[, -c(1,3)]


# Now bind the two
fourteen_trends_mammals <- full_join(four_trends_mammals, isotemp2)

# Set bins back to numeric
fourteen_trends_mammals %<>% transform(bins = as.numeric(as.character(bins)))

# remove redundant colums
fourteen_trends_mammals <- fourteen_trends_mammals[1:1877,]

rm(gradstein, isotemp2, mammals, mammals_raw)

