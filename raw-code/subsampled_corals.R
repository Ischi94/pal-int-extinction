# set working directory
setwd("C:/Users/gmath/Documents/Master/4.Semester/Masterarbeit") 

# keep in mind that you have only isotope data starting at the ordovician!
#remove workspace
rm(list = ls())


library(tidyr)
library(dplyr)
library(divDyn)

#function for occurrences into stage bins
binning <-function(x){gradstein$stg[x<gradstein$bottom  &  x>=gradstein$top]}

# import data 
gradstein <- read.csv("gradstein.csv")[,-1]


# download all bivalves occurrences from the paleobiology database
# this might take some time, as it is quite a big file
corals_pbdb <- read.csv("pbdb_corals.csv")

# clean it a bit
corals_pbdb <- subset(corals_pbdb, corals_pbdb$accepted_rank != "subgenus")

toOmit <- omit(corals_pbdb, bin="stg", tax="genus", om="ref", ref="reference_no")
corals_pbdb <- corals_pbdb[!toOmit,]

# subset the genera with brackets and assign them to the overarching genus
test <- gsub("\\s*\\([^\\)]+\\)","",as.character(corals_pbdb$genus))
corals_pbdb$genus <- test

#Add mean age
corals_pbdb$mean_ma<-(corals_pbdb$max_ma+corals_pbdb$min_ma)/2

#Add  stage number 
corals_pbdb$stg <-lapply(corals_pbdb$mean_ma,binning)

# Exclude unassigned bins
corals_pbdb <- subset(corals_pbdb, is.na(stg)==F)
corals_pbdb$stg <- as.integer(corals_pbdb$stg)

# subset to the same bins as we currently have for isotemp (>= 14)
corals_pbdb <- subset(corals_pbdb, corals_pbdb$stg >=14)

# SQS <- subsample(corals_pbdb,iter=1, q=0.4,tax="genus", bin="stg", 
#                 output="dist", type="sqs", FUN = fadlad)

CR <- subsample(corals_pbdb,iter=1, q=50,tax="genus", bin="stg", 
                output="dist", type="cr", FUN = fadlad)

corals <- CR$results[[1]]
corals <- corals[,-3]

corals$genus <- rownames(corals)


# merge them together to get taxonomy
corals <- merge(corals,corals_pbdb, by  = "genus")

# order them properly
corals <- corals %>%
  dplyr::select(phylum, class, order, family, genus, FAD, LAD)  

# remove the duplicates
corals <- corals[!duplicated(corals),]

# remove the NA's
corals <- subset(corals, corals$class != "NO_CLASS_SPECIFIED" &
                   corals$order != "NO_ORDER_SPECIFIED" &
                   corals$family !=  "NO_FAMILY_SPECIFIED")

corals <- na.omit(corals) # no NA's anymore, so good to go

#remove singletons as they add noise 
corals <- subset(corals, corals$FAD != corals$LAD)

# when forams range to the recent (bin 95), shift this bin to 94 to remove pull of the recent and 
# to enable temp-calculation
corals$LAD[corals$LAD==95] <- 94
corals <- subset(corals, corals$FAD != corals$LAD)

# subset to the same bins as we currently have for isotemp (>= 14)
corals <- subset(corals, corals$FAD >=14)

# change colnames
colnames(corals) <- c("Phylum", "Class", "Order", "Family", "Genus", "FAD", "LAD")


# prepare weizer & prokoph temperature data
isotemp <- read.csv(file="TimeSeriesUsed.csv", header=TRUE,row.names = 1) 

#assign age to isotope data
isotemp$age <- gradstein$mid[14:94]

# get isotemp2
isotemp2 <- read.csv("isotemp2.csv")[,-1]

# build dummy for calculation
dumbo <- corals
corals_safe <- corals[,c(5:7)]

# add bins with 0 and fill in with 1 for extinction and origination
namevector <- as.character(c(1:94))
corals_safe[ , namevector] <- NA
corals_safe <- corals_safe[, c(4:length(corals_safe[0,]))]

for (i in 1:length(corals_safe[,1])) {
  corals_safe[i,dumbo[i,6]] <- 0 
  corals_safe[i,dumbo[i,7]] <- 1
  ifelse(dumbo[i,6]!= dumbo[i,7]-1 & dumbo[i,6]!= dumbo[i,7],
         corals_safe[i,(dumbo[i,6]+1):(dumbo[i,7]-1)] <- 0, NA)
  ifelse(dumbo[i,7] == 94, corals_safe[i,dumbo[i,7]] <- 0, 
         corals_safe[i,dumbo[i,7]] <- 1)
}


# bind it again 
rownames(corals_safe) <-  make.names(dumbo[,5], unique=TRUE)


# Reorganise ranges through time data so we can bind it to temperature data
corals_safe %<>%  as_tibble() %>%
  mutate(Genus = rownames(corals_safe)) %>%
  group_by(Genus) %>%
  gather(-Genus, key="bins", value="extinct", na.rm = T, factor_key = T) %>%
  arrange(desc(bins), .by_group = T) 

corals_safe$bins <- as.integer(corals_safe$bins)

# Now bind the two
temp_corals <- full_join(corals_safe, isotemp2)

#add taxonomical levels
dumbo <- dumbo[, 2:5]

# Now bind the two
temp_corals <- full_join(temp_corals, dumbo)

temp_corals <- temp_corals %>%
  dplyr::select(Class, Order, Family, Genus, bins, age, extinct, Temp.mean, change.prev) %>%
  na.omit()



# Calculate long term changes in temp (1 stage to 10 stages)
# Use lm() to determine slope of the long-term  temperature trends
# If we want to investigate the interaction between long term change and short term,
# we need to exclude the short term bin from the long term (see Manuel's figures),
# and thus shift the long term lm result up by +1. This will be messy, I'm sorry
temp_corals$bins <- as.numeric(temp_corals$bins)

ori_bin <- temp_corals %>%
  dplyr::select(Genus, bins, age) %>%
  group_by(Genus) %>%
  dplyr::slice(which.min(bins)) 
colnames(ori_bin) <- c("Genus", "ori.bin", "ori.age" )

# Now bind the two
test2 <- full_join(temp_corals, ori_bin)

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
test2 <- test2[,-c(10:11)]

#now do the same for trend2
# first build a data frame where we save the origin of each family

biv_fam <- unique(corals$Family)
fam_dat <- data.frame(matrix(ncol =2 , nrow = length(biv_fam)))
colnames(fam_dat) <- c("age.ori", "bin.ori")
rownames(fam_dat) <-  c(as.character(biv_fam))


for (i in 1:length(fam_dat$age.ori)) {
  fam_dat[i, 2] <- min(corals$FAD[corals$Family == biv_fam[i]])
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
test3 <- test3[,-c(13:14)]

#now do the same for trend3 Order
# first build a data frame where we save the origin of each family

biv_ord <- unique(corals$Order)
ord_dat <- data.frame(matrix(ncol =2 , nrow = length(biv_ord)))
colnames(ord_dat) <- c("age.ori", "bin.ori")
rownames(ord_dat) <-  c(as.character(biv_ord))


for (i in 1:length(ord_dat$age.ori)) {
  ord_dat[i, 2] <- min(corals$FAD[corals$Order == biv_ord[i]])
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
test4 <- test4[,-c(13:14)]


#remov ori age and bin
three_trends_corals <- test4

# Set bins back to numeric
three_trends_corals %<>% transform(bins = as.numeric(as.character(bins)))

# Cleanup
rm(dumbo, fam_dat, corals, lin1, 
   ord_dat, ori_bin,  test2, test3, test4, age1, bin1, 
   i, biv_fam, biv_ord, namevector, temp1, temp_corals, corals_safe,
   test, toOmit, binning, SQS)

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
three_trends_corals %<>% transform(bins = as.factor(as.character(bins)))
isotemp2 %<>% transform(bins = as.factor(as.character(bins)))

# Now bind the two
thirteen_trends_corals <- full_join(three_trends_corals, isotemp2)


# Set bins back to numeric
thirteen_trends_corals %<>% transform(bins = as.numeric(as.character(bins)))

# remove redundant colums
thirteen_trends_corals <- thirteen_trends_corals[-c(1500:1502),]

rm(gradstein, isotemp2)

### calculate glmm ###

library(lme4)
library(geiger)
library(visreg)
library(ggplot2)

# replace the NA's with zeros as we have no trend there
thirteen_trends_corals$trend1[is.na(thirteen_trends_corals$trend1)] <- 0
thirteen_trends_corals$trend2[is.na(thirteen_trends_corals$trend2)] <- 0
thirteen_trends_corals$trend3[is.na(thirteen_trends_corals$trend3)] <- 0

# Split short term temperature change into warming and cooling:
thirteen_trends_corals$cooling<-ifelse(thirteen_trends_corals$change.prev<0, thirteen_trends_corals$change.prev, NA)
thirteen_trends_corals$warming<-ifelse(thirteen_trends_corals$change.prev>0, thirteen_trends_corals$change.prev, NA)

## Warming
# Iterate through each long term temperature change
vars = names(dplyr::select(thirteen_trends_corals, trend1:trend3, trend.st1:trend.st10)) 
modelsw4 = lapply(setNames(vars, vars), function(var) {
  form = paste("extinct~warming:", var, "+(1|Genus)")
  glmer(form, data=thirteen_trends_corals, family="binomial")
})

# Make data frame for model output
warm_fourteen_trends <- data.frame(model=names(dplyr::select(thirteen_trends_corals,
                                                             trend1:trend3, trend.st1:trend.st10)), 
                                   intercept=NA, interaction=NA, AIC=NA)

# Run loop to fill for.warm (coefficients and p-values)
for (i in warm_fourteen_trends$model) {
  sum <- summary(modelsw4[[i]])
  warm_fourteen_trends[warm_fourteen_trends$model==i, "intercept"]<-
    paste(round(sum$coefficients[1, 1], 3), sep=" ", "±",
          round(sum$coefficients[1, 2], 3),
          ifelse(sum$coefficients[1, 4]<0.001, "***", 
                 ifelse(sum$coefficients[1, 4]<0.01, "**", 
                        ifelse(sum$coefficients[1, 4]<0.05, "*", ""))))
  warm_fourteen_trends[warm_fourteen_trends$model==i, "interaction"]<-
    paste(round(sum$coefficients[2, 1], 3), sep=" ", "±",
          round(sum$coefficients[2, 2], 3),
          ifelse(sum$coefficients[2, 4]<0.001, "***", 
                 ifelse(sum$coefficients[2, 4]<0.01, "**", 
                        ifelse(sum$coefficients[2, 4]<0.05, "*", ""))))
  warm_fourteen_trends[warm_fourteen_trends$model==i, "AIC"]<-as.numeric(round(sum$AICtab[[1]], 1))
}

# Add column weith AIC weights
warm_fourteen_trends$dAIC<- as.numeric(round(aicw(warm_fourteen_trends$AIC)$delta, 1))
warm_fourteen_trends$AIC.weights<- as.numeric(signif(aicw(warm_fourteen_trends$AIC)$w, 3))

sum_warm <- summary(modelsw4[[which(warm_fourteen_trends$dAIC==0)]])

## Cooling
# Iterate through each long term temperature change
vars = names(dplyr::select(thirteen_trends_corals, trend1:trend3, trend.st1:trend.st10)) 
modelsc4 = lapply(setNames(vars, vars), function(var) {
  form = paste("extinct~cooling:", var, "+(1|Genus)")
  glmer(form, data=thirteen_trends_corals, family="binomial")
})

# Make data frame for model output
cool_fourteen_trends <-data.frame(model=names(dplyr::select(thirteen_trends_corals,
                                                            trend1:trend3, trend.st1:trend.st10)), 
                                  intercept=NA, interaction=NA, AIC=NA)


# Run loop to fill for.cool (coefficients and p-values)
for (i in cool_fourteen_trends$model) {
  sum <- summary(modelsc4[[i]])
  cool_fourteen_trends[cool_fourteen_trends$model==i, "intercept"]<-
    paste(round(sum$coefficients[1, 1], 3), sep=" ", "±",
          round(sum$coefficients[1, 2], 3),
          ifelse(sum$coefficients[1, 4]<0.001, "***", 
                 ifelse(sum$coefficients[1, 4]<0.01, "**", 
                        ifelse(sum$coefficients[1, 4]<0.05, "*", ""))))
  cool_fourteen_trends[cool_fourteen_trends$model==i, "interaction"]<-
    paste(round(sum$coefficients[2, 1], 3), sep=" ", "±",
          round(sum$coefficients[2, 2], 3),
          ifelse(sum$coefficients[2, 4]<0.001, "***", 
                 ifelse(sum$coefficients[2, 4]<0.01, "**", 
                        ifelse(sum$coefficients[2, 4]<0.05, "*", ""))))
  cool_fourteen_trends[cool_fourteen_trends$model==i, "AIC"]<-as.numeric(round(sum$AICtab[[1]], 1))
}

# Add column weith AIC weights
cool_fourteen_trends$dAIC<- as.numeric(round(aicw(cool_fourteen_trends$AIC)$delta, 1))
cool_fourteen_trends$AIC.weights<- as.numeric(signif(aicw(cool_fourteen_trends$AIC)$w, 3))

sum_cool <- summary(modelsc4[[which(cool_fourteen_trends$dAIC==0)]])

### now without long-term temperature ###

modelsw4 <- glmer("extinct~warming+(1|Genus)", data=thirteen_trends_corals, 
                  family="binomial")
sum_warm_wo <- summary(modelsw4)


modelsc4 <- glmer("extinct~cooling+(1|Genus)", data=thirteen_trends_corals, 
                  family="binomial")
sum_cool_wo <- summary(modelsc4)
