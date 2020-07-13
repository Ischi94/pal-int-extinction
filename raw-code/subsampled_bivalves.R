# Here we show how to download data for Bivalvia from the PBDB database
# and how to clean (remove duplicates, omit Pull of the Recent, 
# remove singletons) and process it (bin data to stages, transform 
# occurrence data to range data)

# subsequently we add temperature data and calculate short-term and long-
# term trends of temperature change as an input for the final GLMMs

# load libraries
library(tidyr)
library(dplyr)
library(divDyn)

#function for occurrences into stage bins
binning <-function(x){gradstein$stg[x<gradstein$bottom  &  x>=gradstein$top]}

# import data 
gradstein <- read.csv("gradstein.csv")[,-1]



# download all bivalves occurrences from the paleobiology database
# this might take some time, as it is quite a big file
bivalves_pbdb <- read.csv("pbdb_bivalvia.csv")

# clean it a bit
bivalves_pbdb <- subset(bivalves_pbdb, bivalves_pbdb$accepted_rank != "subgenus")

toOmit <- omit(bivalves_pbdb, bin="stg", tax="genus", om="ref", ref="reference_no")
bivalves_pbdb <- bivalves_pbdb[!toOmit,]

# subset the genera with brackets and assign them to the overarching genus
test <- gsub("\\s*\\([^\\)]+\\)","",as.character(bivalves_pbdb$genus))
bivalves_pbdb$genus <- test

#Add mean age
bivalves_pbdb$mean_ma<-(bivalves_pbdb$max_ma+bivalves_pbdb$min_ma)/2

#Add  stage number 
bivalves_pbdb$stg <-lapply(bivalves_pbdb$mean_ma,binning)

# Exclude unassigned bins
bivalves_pbdb <- subset(bivalves_pbdb, is.na(stg)==F)
bivalves_pbdb$stg <- as.integer(bivalves_pbdb$stg)

#SQS <- subsample(bivalves_pbdb,iter=1, q=0.4,tax="genus", bin="stg", 
#               output="dist", type="sqs", FUN = fadlad)

CR <- subsample(bivalves_pbdb,iter=1, q=50,tax="genus", bin="stg", 
               output="dist", type="cr", FUN = fadlad)

bivalves <- CR$results[[1]]
bivalves <- bivalves[,-3]

bivalves$genus <- rownames(bivalves)


# merge them together to get taxonomy
bivalves <- merge(bivalves,bivalves_pbdb, by  = "genus")

# order them properly
bivalves <- bivalves %>%
  dplyr::select(phylum, class, order, family, genus, FAD, LAD)  

# remove the duplicates
bivalves <- bivalves[!duplicated(bivalves),]

# remove the NA's
bivalves <- subset(bivalves, bivalves$class != "NO_CLASS_SPECIFIED" &
                     bivalves$order != "NO_ORDER_SPECIFIED" &
                     bivalves$family !=  "NO_FAMILY_SPECIFIED")

bivalves <- na.omit(bivalves) # no NA's anymore, so good to go

#remove singletons as they add noise 
bivalves <- subset(bivalves, bivalves$FAD != bivalves$LAD)

# when forams range to the recent (bin 95), shift this bin to 94 to remove pull of the recent and 
# to enable temp-calculation
bivalves$LAD[bivalves$LAD==95] <- 94
bivalves <- subset(bivalves, bivalves$FAD != bivalves$LAD)

# subset to the same bins as we currently have for isotemp (>= 14)
bivalves <- subset(bivalves, bivalves$FAD >=14)

# change colnames
colnames(bivalves) <- c("Phylum", "Class", "Order", "Family", "Genus", "FAD", "LAD")


# prepare weizer & prokoph temperature data
isotemp <- read.csv(file="TimeSeriesUsed.csv", header=TRUE,row.names = 1) 

#assign age to isotope data
isotemp$age <- gradstein$mid[14:94]

# get isotemp2
isotemp2 <- read.csv("isotemp2.csv")[,-1]

# build dummy for calculation
dumbo <- bivalves
bivalves_safe <- bivalves[,c(5:7)]

# add bins with 0 and fill in with 1 for extinction and origination
namevector <- as.character(c(1:94))
bivalves_safe[ , namevector] <- NA
bivalves_safe <- bivalves_safe[, c(4:length(bivalves_safe[0,]))]

for (i in 1:length(bivalves_safe[,1])) {
  bivalves_safe[i,dumbo[i,6]] <- 0 
  bivalves_safe[i,dumbo[i,7]] <- 1
  ifelse(dumbo[i,6]!= dumbo[i,7]-1 & dumbo[i,6]!= dumbo[i,7],
         bivalves_safe[i,(dumbo[i,6]+1):(dumbo[i,7]-1)] <- 0, NA)
  ifelse(dumbo[i,7] == 94, bivalves_safe[i,dumbo[i,7]] <- 0, 
         bivalves_safe[i,dumbo[i,7]] <- 1)
}


# bind it again 
rownames(bivalves_safe) <-  make.names(dumbo[,5], unique=TRUE)


# Reorganise ranges through time data so we can bind it to temperature data
bivalves_safe %<>%  as_tibble() %>%
  mutate(Genus = rownames(bivalves_safe)) %>%
  group_by(Genus) %>%
  gather(-Genus, key="bins", value="extinct", na.rm = T, factor_key = T) %>%
  arrange(desc(bins), .by_group = T) 

bivalves_safe$bins <- as.integer(bivalves_safe$bins)

# Now bind the two
temp_bivalves <- full_join(bivalves_safe, isotemp2)

#add taxonomical levels
dumbo <- dumbo[, 2:5]

# Now bind the two
temp_bivalves <- full_join(temp_bivalves, dumbo)

temp_bivalves <- temp_bivalves %>%
  dplyr::select(Class, Order, Family, Genus, bins, age, extinct, Temp.mean, change.prev) %>%
  na.omit()



# Calculate long term changes in temp (1 stage to 10 stages)
# Use lm() to determine slope of the long-term  temperature trends
# If we want to investigate the interaction between long term change and short term,
# we need to exclude the short term bin from the long term (see Manuel's figures),
# and thus shift the long term lm result up by +1. This will be messy, I'm sorry
temp_bivalves$bins <- as.numeric(temp_bivalves$bins)

ori_bin <- temp_bivalves %>%
  dplyr::select(Genus, bins, age) %>%
  group_by(Genus) %>%
  dplyr::slice(which.min(bins)) 
colnames(ori_bin) <- c("Genus", "ori.bin", "ori.age" )

# Now bind the two
test2 <- full_join(temp_bivalves, ori_bin)

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

biv_fam <- unique(bivalves$Family)
fam_dat <- data.frame(matrix(ncol =2 , nrow = length(biv_fam)))
colnames(fam_dat) <- c("age.ori", "bin.ori")
rownames(fam_dat) <-  c(as.character(biv_fam))


for (i in 1:length(fam_dat$age.ori)) {
  fam_dat[i, 2] <- min(bivalves$FAD[bivalves$Family == biv_fam[i]])
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

biv_ord <- unique(bivalves$Order)
ord_dat <- data.frame(matrix(ncol =2 , nrow = length(biv_ord)))
colnames(ord_dat) <- c("age.ori", "bin.ori")
rownames(ord_dat) <-  c(as.character(biv_ord))


for (i in 1:length(ord_dat$age.ori)) {
  ord_dat[i, 2] <- min(bivalves$FAD[bivalves$Order == biv_ord[i]])
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
three_trends_bivalves <- test4

# Set bins back to numeric
three_trends_bivalves %<>% transform(bins = as.numeric(as.character(bins)))

# Cleanup
rm(dumbo, fam_dat, bivalves, lin1, 
   ord_dat, ori_bin,  test2, test3, test4, age1, bin1, 
   i, biv_fam, biv_ord, namevector, temp1, temp_bivalves, bivalves_safe,
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
three_trends_bivalves %<>% transform(bins = as.factor(as.character(bins)))
isotemp2 %<>% transform(bins = as.factor(as.character(bins)))

# Now bind the two
thirteen_trends_bivalves <- full_join(three_trends_bivalves, isotemp2)


# Set bins back to numeric
thirteen_trends_bivalves %<>% transform(bins = as.numeric(as.character(bins)))

# remove redundant colums
thirteen_trends_bivalves <- thirteen_trends_bivalves[-2944,]

rm(gradstein, isotemp2)

### calculate glmm ###
# load libraries
library(lme4)
library(geiger)
library(visreg)
library(ggplot2)

# replace the NA's with zeros as we have no trend there
thirteen_trends_bivalves$trend1[is.na(thirteen_trends_bivalves$trend1)] <- 0
thirteen_trends_bivalves$trend2[is.na(thirteen_trends_bivalves$trend2)] <- 0
thirteen_trends_bivalves$trend3[is.na(thirteen_trends_bivalves$trend3)] <- 0

# Split short term temperature change into warming and cooling:
thirteen_trends_bivalves$cooling<-ifelse(thirteen_trends_bivalves$change.prev<0, thirteen_trends_bivalves$change.prev, NA)
thirteen_trends_bivalves$warming<-ifelse(thirteen_trends_bivalves$change.prev>0, thirteen_trends_bivalves$change.prev, NA)

## Warming
# Iterate through each long term temperature change
vars = names(dplyr::select(thirteen_trends_bivalves, trend1:trend3, trend.st1:trend.st10)) 
modelsw4 = lapply(setNames(vars, vars), function(var) {
  form = paste("extinct~warming:", var, "+(1|Genus)")
  glmer(form, data=thirteen_trends_bivalves, family="binomial")
})

# Make data frame for model output
warm_fourteen_trends <- data.frame(model=names(dplyr::select(thirteen_trends_bivalves,
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
vars = names(dplyr::select(thirteen_trends_bivalves, trend1:trend3, trend.st1:trend.st10)) 
modelsc4 = lapply(setNames(vars, vars), function(var) {
  form = paste("extinct~cooling:", var, "+(1|Genus)")
  glmer(form, data=thirteen_trends_bivalves, family="binomial")
})

# Make data frame for model output
cool_fourteen_trends <-data.frame(model=names(dplyr::select(thirteen_trends_bivalves,
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

modelsw4 <- glmer("extinct~warming+(1|Genus)", data=thirteen_trends_bivalves, 
                 family="binomial")
sum_warm_wo <- summary(modelsw4)


modelsc4 <- glmer("extinct~cooling+(1|Genus)", data=thirteen_trends_bivalves, 
                  family="binomial")
sum_cool_wo <- summary(modelsc4)

