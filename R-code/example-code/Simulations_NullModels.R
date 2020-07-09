# Here we show how to 

# 1 Simulate 100 data sets with 400 observations resembling empirical data, but with no dependency of extinction on 
 # paleoclimatic interaction
  # 2 Build various data sets with varying numbers of observations based on the simulated data
    # 3 Calculate long-term term temperature trends and short-term temperature change for each data set
      # 4 Calculate GLMM's and estimate effect size based on these simulated data sets
        # 5 Plot GLMM results to check for type I error rate of the analytical process
          # 6 Iterate analysis for data sets with 400 observations to obtain a null model distribution to compare 
           # our empirical results with
            # 6.1 Create 1000 data sets with 400 observations
              # 6.2 Calculate long-term term temperature trends and short-term temperature change for each data set
                # 6.3 Calculate GLMM's and estimate effect size based on these simulated data sets
                  # 6.4 Plot GLMM results of Null model



# set working directory to where files are stored using the "rstudioapi" package

# Getting the path of this script
current_path = rstudioapi::getActiveDocumentContext()$path 

# setting it as working directory 
setwd(dirname(current_path ))

# please cross-check that this is the path where you store the other 
# files used for this script
getwd()


##### load libraries ####
library(divDyn)
library(magrittr)
library(readr)
library(tidyr)
library(dplyr)
library(divDyn)
library(lme4)
library(geiger)
library(visreg)
library(ggplot2)
library(svMisc)





# 1 Simulate 100 data sets with 400 observations resembling empiri --------


# import empirical data set of genus durations for clades used in our analysis
durations <- read.csv("durations.csv")[, -1]

# import stage data 
data(stages)
gradstein <- stages
rm(stages)
# as we have used stage age from Gradstein 2012 in our analysis, we need to 
# assign the same ages to the stage file
gradstein$bottom[c(21,50:55, 57:58, 71:74, 86:87, 90, 92, 95)] <- c(443.8, 259.8, 254.2, 252.2, 250.0,247.1, 241.5, 228.4, 209.5, 
                                                                    139.4, 133.9, 130.8, 126.3, 41.2, 37.8, 23.03, 11.63, 0.0118)
gradstein$top[21:length(gradstein$bottom)-1] <- gradstein$bottom[21:length(gradstein$bottom)]

# create empty list to save generated data sets
raw_data <- list()

# start creating data sets. Create 100 data sets with 400 observations each. Durations
for (h in 1:100) {
  
  # create data frame for simulated data
  rand_genus <- data.frame(matrix(ncol = 3, nrow = 400))
  colnames(rand_genus) <-  c("Genus", "FAD.age", "LAD.age")
  nr <- c(1:400)
  rand_genus$Genus <- paste0("Genus", nr)
  
  # assign random FAD to each Order with random flat function runif(). We start in the ordovician at 488.3 myr
  rand_genus$FAD.age <- round(runif(400, min=-0.01, max=481.55), 2)
  
  # now we draw a number from the durations of all genera from our observed data and subtract it from the FAD
  # to get the FAD. Distributions of the generated datasets therefore simulate observed conditions 
  # (i.e. a log-normal distribution). 
  rand_genus$LAD.age <- rand_genus$FAD.age - sample(durations$gen_dur, 
                                                    400, replace = FALSE)
  
  # as we can't have future ranges, we replace every negative LAD with 0 
  # (= range to recent)
  rand_genus$LAD.age[rand_genus$LAD.age < 0] <- 0
  
  # let's group these into 80 Families as we have approximately 
  # 5 genera per family in the observed data (10379/2304)
  
  # get the Family names
  fam_names <- paste0("Family", 1:80)
  
  # and assign them to the data frame of genera
  fam_list <- rep(fam_names, 5) 
  rand_genus$Family <- fam_list
  
  # let's group these into 10 Order as we have approximately 
  # 40 genera per order in the observed data (10379/276)
  
  # get the Order names
  ord_names <- paste0("Order", 1:10)
  
  # and assign them to the data frame of genera
  ord_list <- rep(ord_names, 40) 
  rand_genus$Order <- ord_list
  
  # let's group these into 2 Classes as we have approximately 
  # 300 genera per order in the observed data (10379/33)
  
  # get the Order names
  cla_names <- paste0("Class", 1:2)
  
  # and assign them to the data frame of genera
  cla_list <- rep(cla_names, 200) 
  rand_genus$Class <- cla_list
  
  
  # bin FAD and LAD of each genus  
  rand_genus$FAD_bin <- NA
  rand_genus$LAD_bin <- NA
  
  for (i in 1:length(rand_genus$Genus)) {
    for (j in 1:length(gradstein$stg)) {
      ifelse(gradstein$bottom[j] >= rand_genus$FAD.age[i] & gradstein$top[j] < rand_genus$FAD.age[i],
             rand_genus$FAD_bin[i] <- gradstein$stg[j], NA)
      ifelse(gradstein$bottom[j] > rand_genus$LAD.age[i] & gradstein$top[j] <= rand_genus$LAD.age[i],
             rand_genus$LAD_bin[i] <- gradstein$stg[j], NA)
    }}
  
  raw_data[[h]] <- rand_genus
  
  progress(h, max.value = 100)
  if (h == 100) message("Done!")
}

# save(raw_data, file="raw_sim_data.RData")
# load("raw_sim_data.RData")

# Cleanup
rm(list=ls()[! ls() %in% c("raw_data", "gradstein")])




# 2 Build various data sets with varying numbers of observations b --------


# now we have 100 data sets with 400 genera each. We want to see the rate of false positives (type I error) of our 
# analytical process and how this rate changes with size of data sets/ number of observations. 
# For this purpose, we  build 12 data sets from our simulated data with a varying number of observations, 
# starting at 25 and ending at 3000 observations. 

# define number of observations per data sets
df_size <- c(25, 50, 100, 200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000)

# build empty list
rep_data <- list()

# create data sets using "sample_n" function

for (h in 1:100) {
  rep_data[[h]]  <-  sample_n(raw_data[[h]], df_size[1], replace = T)
  rep_data[[100+h]]  <-  sample_n(raw_data[[h]], df_size[2], replace = T)
  rep_data[[200+h]]  <-  sample_n(raw_data[[h]], df_size[3], replace = T)
  rep_data[[300+h]]  <-  sample_n(raw_data[[h]], df_size[4], replace = T)
  rep_data[[400+h]]  <-  sample_n(raw_data[[h]], df_size[5], replace = T)
  rep_data[[500+h]]  <-  sample_n(raw_data[[h]], df_size[6], replace = T)
  rep_data[[600+h]]  <-  sample_n(raw_data[[h]], df_size[7], replace = T)
  rep_data[[700+h]]  <-  sample_n(raw_data[[h]], df_size[8], replace = T)
  rep_data[[800+h]]  <-  sample_n(raw_data[[h]], df_size[9], replace = T)
  rep_data[[900+h]]  <-  sample_n(raw_data[[h]], df_size[10], replace = T)
  rep_data[[1000+h]]  <-  sample_n(raw_data[[h]], df_size[11], replace = T)
  rep_data[[1100+h]]  <-  sample_n(raw_data[[h]], df_size[12], replace = T)

  progress(h, max.value = 100)
  if (h == 100) message("Done!")
}




# 3 Calculate long-term term temperature trends and short-term tem --------


# this loop follows the procedure of data processing described in the file "data_preparation.R"
# for detailed comments, please refer to this file. 
# create empty list
finished_data <- list()

# 
for (h in 1:1200) {
  
  # load the data
  rand_genus <- subset(rep_data[[h]], rep_data[[h]]$FAD_bin != rep_data[[h]]$LAD_bin)
  rand_genus <- drop_na(rand_genus, c("FAD_bin"))
  rand_genus$LAD_bin[rand_genus$LAD_bin==95] <- 94

  # build dummy for calculation
  dumbo <- rand_genus
  rand_genus <- rand_genus[,c("Genus","FAD_bin","LAD_bin")]
  
  # add bins with 0 and fill in with 1 for extinction and origination
  namevector <- as.character(c(1:94))
  rand_genus[ , namevector] <- NA
  rand_genus <- rand_genus[, c(4:length(rand_genus[0,]))]
  
  for (i in 1:length(rand_genus[,1])) {
    rand_genus[i,dumbo[i,"FAD_bin"]] <- 0 
    rand_genus[i,dumbo[i,"LAD_bin"]] <- 1
    ifelse(dumbo[i,7]!= dumbo[i,"LAD_bin"]-1 & dumbo[i,"FAD_bin"]!= dumbo[i,"LAD_bin"],
           rand_genus[i,(dumbo[i,"FAD_bin"]+1):(dumbo[i,"LAD_bin"]-1)] <- 0, NA)
    ifelse(dumbo[i,"LAD_bin"] == 94, rand_genus[i,dumbo[i,"LAD_bin"]] <- 0, 
           rand_genus[i,dumbo[i,"LAD_bin"]] <- 1)
  }
  
  
  # bind with dummy to get higher taxonomy 
  rownames(rand_genus) <-  make.names(dumbo[,"Genus"], unique=TRUE)
  
  
  # Reorganise ranges through time data so we can bind it to temperature 
  # data
  rand_genus %<>%  as_tibble() %>%
    mutate(Genus = rownames(rand_genus)) %>%
    group_by(Genus) %>%
    gather(-Genus, key="bins", value="extinct", na.rm = T, factor_key = T) %>%
    arrange(desc(bins), .by_group = T) 
  
  
  # read in veizer and prokoph temperature data, already prepared. See "data_preparation.R"
  isotemp2 <- read.csv(file="isotemp2.csv", header=TRUE,row.names = 1) 
  isotemp2$bins <- as.factor(isotemp2$bins)
  
  # Now bind the two
  rand_genus <- full_join(rand_genus, isotemp2)
  
  
  #add taxonomical levels
  taxonomy <- dumbo[, c("Class","Order","Family", "Genus")]
  
  # Now bind the two
  rand_genus <- full_join(rand_genus, taxonomy)
  
  # reorder properly
  rand_genus <- rand_genus %>%
    select(Class, Order, Family, Genus, bins, age, extinct, Temp.mean, change.prev)
  
  
  
  # Calculate long term changes in temp (1 stage to 10 stages)
  # Use lm() to determine slope of the long-term  temperature trends
  # If we want to investigate the interaction between long term change and short term,
  # we need to exclude the short term bin from the long term,
  # and thus shift the long term lm result up by +1.
  ori_bin <- rand_genus %>%
    select(Genus, bins, age) %>%
    group_by(Genus) %>%
    dplyr::slice(which.min(bins)) 
  colnames(ori_bin) <- c("Genus", "ori.bin", "ori.age" )
  
  # Now bind the two
  rand_genus <- full_join(rand_genus, ori_bin)
  
  # Set bins back to numeric
  rand_genus %<>% transform(bins = as.numeric(as.character(bins)))
  rand_genus$trend1 <- NA
  rand_genus$trend2 <- NA
  rand_genus$trend3 <- NA
  rand_genus$trend4 <- NA
  
  # remove NA's
  rand_genus <- drop_na(rand_genus, c("Class","Order", "Family", "Genus",
                                      "bins", "age"))
  
  # calculate short-term temperature change
  for (i in 1:length(rand_genus$Genus)) {
    bin1 <- c(rand_genus$ori.bin[i], 
              ifelse(rand_genus$bins[i]>rand_genus$ori.bin[i], 
                     rand_genus$bins[i]-1, 
                     rand_genus$bins[i]))
    age1 <- c(isotemp2$age[isotemp2$Stage %in% bin1[1]:bin1[2]])
    
    temp1 <-c(isotemp2$Temp[isotemp2$Stage %in% bin1[1]:bin1[2]])
    lin1 <- lm(temp1~age1)
    rand_genus$trend1[i] <- -lin1$coefficients[2]
  }
  
  #remov ori age and bin
  rand_genus <- rand_genus[, !(names(rand_genus) %in% c("ori.bin","ori.age"))]
  
  #now do the same for trend2
  # first build a data frame where we save the origin of each family
  
  fam <- unique(taxonomy$Family)
  fam_dat <- data.frame(matrix(ncol =2 , nrow = length(fam)))
  colnames(fam_dat) <- c("ori.age", "ori.bin")
  rownames(fam_dat) <-  c(as.character(fam))
  
  
  for (i in 1:length(fam_dat$ori.age)) {
    fam_dat[i, 2] <- min(dumbo$FAD_bin[dumbo$Family == fam[i]])
    fam_dat[i, 1] <- isotemp2$age[isotemp2$Stage == fam_dat[i, 2]]
  }
  
  
  # add a colum with the family names for binding
  fam_dat$Family <- rownames(fam_dat)
  
  # Now bind the two
  rand_genus <- full_join(rand_genus, fam_dat)
  
  # now calculate the trend starting at origin of the specific family
  for (i in 1:length(rand_genus$Genus)) {
    bin1 <- c(rand_genus$ori.bin[i], 
              ifelse(rand_genus$bins[i]>rand_genus$ori.bin[i], 
                     rand_genus$bins[i]-1, 
                     rand_genus$bins[i]))
    age1 <- c(isotemp2$age[isotemp2$Stage %in% bin1[1]:bin1[2]])
    
    temp1 <-c(isotemp2$Temp[isotemp2$Stage %in% bin1[1]:bin1[2]])
    lin1 <- lm(temp1~age1)
    rand_genus$trend2[i] <- -lin1$coefficients[2]
  }
  
  
  #remov ori age and bin
  rand_genus <- rand_genus[,!(names(rand_genus) %in% c("ori.bin","ori.age"))]
  
  #now do the same for trend3 Order
  # first build a data frame where we save the origin of each Order
  
  ord <- unique(dumbo$Order)
  ord_dat <- data.frame(matrix(ncol =2 , nrow = length(ord)))
  colnames(ord_dat) <- c("ori.age", "ori.bin")
  rownames(ord_dat) <-  c(as.character(ord))
  
  
  for (i in 1:length(ord_dat$ori.age)) {
    ord_dat[i, 2] <- min(dumbo$FAD_bin[dumbo$Order == ord[i]])
    ord_dat[i, 1] <- isotemp2$age[isotemp2$Stage == ord_dat[i, 2]]
  }
  
  
  # add a colum with the Order names for binding
  ord_dat$Order <- rownames(ord_dat)
  
  # Now bind the two
  rand_genus <- full_join(rand_genus, ord_dat)
  
  # now calculate the trend starting at origin of the specific Order
  for (i in 1:length(rand_genus$Genus)) {
    bin1 <- c(rand_genus$ori.bin[i], 
              ifelse(rand_genus$bins[i]>rand_genus$ori.bin[i], 
                     rand_genus$bins[i]-1, rand_genus$bins[i]))
    age1 <- c(isotemp2$age[isotemp2$Stage %in% bin1[1]:bin1[2]])
    
    temp1 <-c(isotemp2$Temp[isotemp2$Stage %in% bin1[1]:bin1[2]])
    lin1 <- lm(temp1~age1)
    rand_genus$trend3[i] <- -lin1$coefficients[2]
  }
  
  
  #remov ori age and bin
  rand_genus <- rand_genus[,!(names(rand_genus) %in% c("ori.bin","ori.age"))]
  
  #now do the same for trend4 Class
  # first build a data frame where we save the origin of each Class
  
  cla <- unique(dumbo$Class)
  cla_dat <- data.frame(matrix(ncol =2 , nrow = length(cla)))
  colnames(cla_dat) <- c("ori.age", "ori.bin")
  rownames(cla_dat) <-  c(as.character(cla))
  
  
  for (i in 1:length(cla_dat$ori.age)) {
    cla_dat[i, 2] <- min(dumbo$FAD_bin[dumbo$Class == cla[i]])
    cla_dat[i, 1] <- isotemp2$age[isotemp2$Stage == cla_dat[i, 2]]
  }
  
  
  # add a colum with the Class names for binding
  cla_dat$Class <- rownames(cla_dat)
  
  # Now bind the two
  rand_genus <- full_join(rand_genus, cla_dat)
  
  # now calculate the trend starting at origin of the specific Class
  for (i in 1:length(rand_genus$Genus)) {
    bin1 <- c(rand_genus$ori.bin[i], 
              ifelse(rand_genus$bins[i]>rand_genus$ori.bin[i], 
                     rand_genus$bins[i]-1, rand_genus$bins[i]))
    age1 <- c(isotemp2$age[isotemp2$Stage %in% bin1[1]:bin1[2]])
    
    temp1 <-c(isotemp2$Temp[isotemp2$Stage %in% bin1[1]:bin1[2]])
    lin1 <- lm(temp1~age1)
    rand_genus$trend4[i] <- -lin1$coefficients[2]
  }
  
  #remov ori age and bin
  rand_genus <- rand_genus[,!(names(rand_genus) %in% c("ori.bin","ori.age"))]
  
  
  # load the temperature trends, already calculated. See "data_preparation.R"
  temp_trends <- read.csv("temp_trends.csv", header = T, row.names = 1)
  
  # Set bins to factor
  rand_genus %<>% transform(bins = as.factor(as.character(bins)))
  temp_trends %<>% transform(bins = as.factor(as.character(bins)))
  
  temp_trends <- temp_trends[, -c(1,3)]
  
  
  # Now bind the two
  fourteen_trends_rand <- full_join(rand_genus, temp_trends)
  
  # Set bins back to numeric
  fourteen_trends_rand %<>% transform(bins = as.numeric(as.character(bins)))
  
  
  finished_data[[h]] <- fourteen_trends_rand 
  
  progress(h, max.value = 1200)
  if (h == 1200) message("Done!")
}

# save(finished_data, file="sim_data_new_1200.RData")
# load("sim_data_new_1200.RData")

# Cleanup
rm(list=ls()[! ls() %in% c("finished_data", "df_size")])




# 4 Calculate GLMM's and estimate effect size based on these simul --------


#Build empty data frame for saving the results
result_rand <- data.frame(lower.CI.w = numeric(), estimate.w = numeric(), upper.CI.w =numeric(),
                          lower.CI.c = numeric(), estimate.c = numeric(), upper.CI.c =numeric())

# data sets with less than 25 observations have convergence issues and the loop will stop at these, so we start with 
# a minimum of 30 observations
# this analytic procedure follows the script "glmm_analysis.R". For detailled comments, please refer to this script. 

for (h in 30:1200) {
  
  # replace the NA's with zeros as we have no trend there
  finished_data[[h]]$trend1[is.na(finished_data[[h]]$trend1)] <- 0
  finished_data[[h]]$trend2[is.na(finished_data[[h]]$trend2)] <- 0
  finished_data[[h]]$trend3[is.na(finished_data[[h]]$trend3)] <- 0
  finished_data[[h]]$trend4[is.na(finished_data[[h]]$trend4)] <- 0
  
  # Split short term temperature change into warming and cooling:
  finished_data[[h]]$cooling<-ifelse(finished_data[[h]]$change.prev<0, finished_data[[h]]$change.prev, NA)
  finished_data[[h]]$warming<-ifelse(finished_data[[h]]$change.prev>0, finished_data[[h]]$change.prev, NA)
  
  ## Warming
  # Iterate through each long term temperature change
  vars = names(dplyr::select(finished_data[[h]], trend1:trend4, trend.st1:trend.st10)) 
  modelswr = lapply(setNames(vars, vars), function(var) {
    form = paste("extinct~warming:", var, "+(1|Genus)")
    glmer(form, data=finished_data[[h]], family="binomial")
  })
  
  # Make data frame for model output
  warm_fourteen_trends_rand <- data.frame(model=names(dplyr::select(finished_data[[h]], trend1:trend4, trend.st1:trend.st10)), 
                                          intercept=NA, 
                                          interaction=NA, AIC=NA)
  
  # Run loop to fill for.warm (coefficients and p-values)
  for (i in warm_fourteen_trends_rand$model) {
    sum <- summary(modelswr[[i]])
    warm_fourteen_trends_rand[warm_fourteen_trends_rand$model==i, "intercept"]<-
      paste(round(sum$coefficients[1, 1], 3), sep=" ", "±",
            round(sum$coefficients[1, 2], 3),
            ifelse(sum$coefficients[1, 4]<0.001, "***", 
                   ifelse(sum$coefficients[1, 4]<0.01, "**", 
                          ifelse(sum$coefficients[1, 4]<0.05, "*", ""))))
    warm_fourteen_trends_rand[warm_fourteen_trends_rand$model==i, "interaction"]<-
      paste(round(sum$coefficients[2, 1], 3), sep=" ", "±",
            round(sum$coefficients[2, 2], 3),
            ifelse(sum$coefficients[2, 4]<0.001, "***", 
                   ifelse(sum$coefficients[2, 4]<0.01, "**", 
                          ifelse(sum$coefficients[2, 4]<0.05, "*", ""))))
    warm_fourteen_trends_rand[warm_fourteen_trends_rand$model==i, "AIC"]<-as.numeric(round(sum$AICtab[[1]], 1))
  }
  
  #  Add column weith AIC weights
  warm_fourteen_trends_rand$dAIC<- as.numeric(round(aicw(warm_fourteen_trends_rand$AIC)$delta, 1))
  warm_fourteen_trends_rand$AIC.weights<- as.numeric(signif(aicw(warm_fourteen_trends_rand$AIC)$w, 3))
  
  
  # let's take the values of the visreg and split them into cooling and warming and compare the 
  # exctinction response of each by means of boxplots
  vis_out_wr <- visreg(modelswr[[which(warm_fourteen_trends_rand$dAIC==0)[1]]], "warming", scale="response", 
                       by=paste(warm_fourteen_trends_rand[which(warm_fourteen_trends_rand$dAIC==0)[1],1]), 
                       breaks = c(-0.1, 0, 0.1),
                       rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5), plot = F) 
  
  
  
  vis_out_wr$fit$ctrend <- ifelse(vis_out_wr$fit[2]<0,"Trend = Cooling", "Trend = Warming")
  
  wtestw <- wilcox.test(vis_out_wr$fit$visregFit[vis_out_wr$fit$ctrend=="Trend = Warming"],
                        vis_out_wr$fit$visregFit[vis_out_wr$fit$ctrend=="Trend = Cooling"], 
                        paired = F, conf.int = T)
  
  
  
  result_rand[h,1] <- wtestw$conf.int[1]
  result_rand[h,2] <- wtestw$estimate
  result_rand[h,3] <- wtestw$conf.int[2]
  
  
  ## Cooling
  # Iterate through each long term temperature change
  vars = names(dplyr::select(finished_data[[h]], trend1:trend4, trend.st1:trend.st10)) 
  modelscr = lapply(setNames(vars, vars), function(var) {
    form = paste("extinct~cooling:", var, "+(1|Genus)")
    glmer(form, data=finished_data[[h]], family="binomial")
  })
  
  # Make data frame for model output
  cool_fourteen_trends_rand <-data.frame(model=names(dplyr::select(finished_data[[h]], c(trend1:trend4, trend.st1:trend.st10))), 
                                         intercept=NA, 
                                         interaction=NA, AIC=NA)
  
  # Run loop to fill for.cool (coefficients and p-values)
  for (i in cool_fourteen_trends_rand$model) {
    sum <- summary(modelscr[[i]])
    cool_fourteen_trends_rand[cool_fourteen_trends_rand$model==i, "intercept"]<-
      paste(round(sum$coefficients[1, 1], 3), sep=" ", "±",
            round(sum$coefficients[1, 2], 3),
            ifelse(sum$coefficients[1, 4]<0.001, "***", 
                   ifelse(sum$coefficients[1, 4]<0.01, "**", 
                          ifelse(sum$coefficients[1, 4]<0.05, "*", ""))))
    cool_fourteen_trends_rand[cool_fourteen_trends_rand$model==i, "interaction"]<-
      paste(round(sum$coefficients[2, 1], 3), sep=" ", "±",
            round(sum$coefficients[2, 2], 3),
            ifelse(sum$coefficients[2, 4]<0.001, "***", 
                   ifelse(sum$coefficients[2, 4]<0.01, "**", 
                          ifelse(sum$coefficients[2, 4]<0.05, "*", ""))))
    cool_fourteen_trends_rand[cool_fourteen_trends_rand$model==i, "AIC"]<-as.numeric(round(sum$AICtab[[1]], 1))
  }
  
  # Add column weith AIC weights
  cool_fourteen_trends_rand$dAIC<- as.numeric(round(aicw(cool_fourteen_trends_rand$AIC)$delta, 1))
  cool_fourteen_trends_rand$AIC.weights<- as.numeric(signif(aicw(cool_fourteen_trends_rand$AIC)$w, 3))
  
  
  # let's take the values of the visreg and split them into cooling and warming and compare the 
  # exctinction response of each by means of boxplots
  vis_out_cr <- visreg(modelscr[[which(cool_fourteen_trends_rand$dAIC==0)[1]]], "cooling", scale="response", 
                       by=paste(cool_fourteen_trends_rand[which(cool_fourteen_trends_rand$dAIC==0)[1],1]), 
                       breaks = c(-0.1, 0, 0.1),
                       rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5), plot = F) 
  
  
  
  vis_out_cr$fit$ctrend <- ifelse(vis_out_cr$fit[2]<0,"Trend = Cooling", "Trend = Warming")
  
  
  wtestc <- wilcox.test(vis_out_cr$fit$visregFit[vis_out_cr$fit$ctrend=="Trend = Warming"],
                        vis_out_cr$fit$visregFit[vis_out_cr$fit$ctrend=="Trend = Cooling"], 
                        paired = F, conf.int = T)
  
  result_rand[h,4] <- wtestc$conf.int[1]
  result_rand[h,5] <- wtestc$estimate
  result_rand[h,6] <- wtestc$conf.int[2]
  
  progress(h, max.value = 1200)
}

# save(result_rand, file="result_rand__new_1200.RData")

# fun::shutdown()

# Cleanup
 rm(list=ls()[! ls() %in% c("finished_data", "result_rand")])

# load("result_rand__new_1200.RData")

 

# 5 Plot GLMM results to check for type I error rate of the analyt --------


# get number of observations
df_size <- c(25, 50, 100, 200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000)

# add those to the data set 
result_rand$trial <- c(rep("genus_25", 100), rep("genus_50", 100), rep("genus_100", 100), 
                           rep("genus_200", 100), rep("genus_400", 100), rep("genus_600", 100), 
                           rep("genus_800", 100), rep("genus_1000", 100), rep("genus_1500", 100),
                           rep("genus_2000", 100), rep("genus_2500", 100),
                           rep("genus_3000", 100))



# build an empty data frame to save results
result_rand2 <- matrix(nrow = 12, ncol = 7)
colnames(result_rand2) <- c("lower.CI.w", "estimate.w", "upper.CI.w", 
                                "lower.CI.c", "estimate.c", "upper.CI.c", "trial")
result_rand2 <- as.data.frame(result_rand2)

# add the number of trials in the empty dataframe
result_rand2$trial <- unique(result_rand$trial)


# loop through the results to get the highest variance values and the mean of the trials
for (i in unique(result_rand$trial)) {
  result_rand2[result_rand2$trial==i,1] <- min(result_rand$lower.CI.w[result_rand$trial==i])
  result_rand2[result_rand2$trial==i,2] <- mean(result_rand$estimate.w[result_rand$trial==i])
  result_rand2[result_rand2$trial==i,3] <- max(result_rand$upper.CI.w[result_rand$trial==i])
  result_rand2[result_rand2$trial==i,4] <- min(result_rand$lower.CI.c[result_rand$trial==i])
  result_rand2[result_rand2$trial==i,5] <- mean(result_rand$estimate.c[result_rand$trial==i])
  result_rand2[result_rand2$trial==i,6] <- max(result_rand$upper.CI.c[result_rand$trial==i])
}



# add the number of observations/ genera
result_rand2$observations <- df_size

# define theme
my_theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  title = element_text(size = 12),
  panel.border = element_blank(),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

ggplot(result_rand2, aes(df_size, estimate.w))+
  geom_vline(xintercept=0, linetype="dashed", color = "grey50")+
  geom_vline(xintercept=30, linetype="dashed", color = "grey50")+
  geom_vline(xintercept=400, linetype="dashed", color = "grey50")+
  geom_ribbon(data=result_rand2,aes(ymin=lower.CI.w,ymax=upper.CI.w),alpha=0.4, fill='coral2')+
  geom_ribbon(data=result_rand2,aes(ymin=lower.CI.c,ymax=upper.CI.c),alpha=0.4, fill= "#56B4E9")+
  geom_point(color='coral2', size=2.5, stroke=2)+
  geom_line(color='coral2', size=1.25)+
  geom_point(data=result_rand2, aes(df_size, estimate.c), color= "#56B4E9", size=2.5)+
  geom_line(data=result_rand2, aes(df_size, estimate.c), color= "#56B4E9", size=1.25)+
  my_theme+ 
  labs(x = "N (Observations)", y = "Change in P (Extinction)")+
  ylim(-1,1)

# save the plot to working directory
# ggsave("type_I_error_rate.pdf")


# 6 Iterate analysis for data set with 400 observations to obtain  --------


# as simulations have shown, type I error rate drops around 50 observations and remains low with increasing
# number of observations. Above 400 observations, GLMM's show no sign of singularity or convergence. 
# as each of our empirical data sets contain more than 400 observations, we use this number as a cutoff for our null
# model. We repeat the analytical steps for 1000 data sets, each containing simulated data resembling empirical data
# but having extinctions independent from paleoclimate interaction. 

# empty environment
rm(list = ls())

# import data 
durations <- read.csv("durations.csv")[, -1]

# stage data 
data(stages)
gradstein <- stages
rm(stages)
# as we have assigned the stage age from Gradstein 2012, we need to 
# assign the same ages to the stage file
gradstein$bottom[c(21,50:55, 57:58, 71:74, 86:87, 90, 92, 95)] <- c(443.8, 259.8, 254.2, 252.2, 250.0,247.1, 241.5, 228.4, 209.5, 
                                                                    139.4, 133.9, 130.8, 126.3, 41.2, 37.8, 23.03, 11.63, 0.0118)
gradstein$top[21:length(gradstein$bottom)-1] <- gradstein$bottom[21:length(gradstein$bottom)]




# 6.1 Create 1000 data sets with 400 observations -------------------------


# create empty list to save generated data sets
raw_data <- list()

# start creating data sets
for (h in 1:1000) {
  
  # create randomized data
  rand_genus <- data.frame(matrix(ncol = 3, nrow = 400))
  colnames(rand_genus) <-  c("Genus", "FAD.age", "LAD.age")
  nr <- c(1:400)
  rand_genus$Genus <- paste0("Genus", nr)
  
  # assign random FAD to each Order with random flat function runif(),  We start in the ordovician at 481.55 myr
  rand_genus$FAD.age <- round(runif(400, min=-0.01, max=481.55), 2)
  
  # now we draw a number from the durations of all genera from our observed data and subtract it from the FAD
  # to get the FAD. Distributions of the generated datasets therefore simulate observed conditions 
  # (i.e. a log-normal distribution).
  rand_genus$LAD.age <- rand_genus$FAD.age - sample(durations$gen_dur, 
                                                    400, replace = FALSE)
  
  # as we can't have future ranges, we replace every negative LAD with 0 
  # (= range to recent)
  rand_genus$LAD.age[rand_genus$LAD.age < 0] <- 0
  
  # let's group these into 80 Families as we have approximately 
  # 5 genera per family in the real data (10379/2304)
  
  # get the Family names
  fam_names <- paste0("Family", 1:80)
  
  # and assign them to the data frame of genera
  fam_list <- rep(fam_names, 5) 
  rand_genus$Family <- fam_list
  
  # let's group these into 10 Order as we have approximately 
  # 40 genera per order in the real data (10379/276)
  
  # get the Order names
  ord_names <- paste0("Order", 1:10)
  
  # and assign them to the data frame of genera
  ord_list <- rep(ord_names, 40) 
  rand_genus$Order <- ord_list
  
  # let's group these into 2 Classes as we have approximately 
  # 300 genera per order in the real data (10379/33)
  
  # get the Class names
  cla_names <- paste0("Class", 1:2)
  
  # and assign them to the data frame of genera
  cla_list <- rep(cla_names, 200) 
  rand_genus$Class <- cla_list
  
  
  # bin FAD and LAD of each genus  
  rand_genus$FAD_bin <- NA
  rand_genus$LAD_bin <- NA
  
  for (i in 1:length(rand_genus$Genus)) {
    for (j in 1:length(gradstein$stg)) {
      ifelse(gradstein$bottom[j] >= rand_genus$FAD.age[i] & gradstein$top[j] < rand_genus$FAD.age[i],
             rand_genus$FAD_bin[i] <- gradstein$stg[j], NA)
      ifelse(gradstein$bottom[j] > rand_genus$LAD.age[i] & gradstein$top[j] <= rand_genus$LAD.age[i],
             rand_genus$LAD_bin[i] <- gradstein$stg[j], NA)
    }}
  
  raw_data[[h]] <- rand_genus
  
  progress(h, max.value = 1000)
  if (h == 1000) message("Done!")
}

# save(raw_data, file="raw_sim_400.RData")
# load("raw_sim_400.RData")

# Cleanup
rm(list=ls()[! ls() %in% c("raw_data", "gradstein")])




# 6.2 Calculate long-term term temperature trends and short-term t --------



# create empty list
finished_data <- list()

for (h in 1:1000) {
  
  # load the data
  rand_genus <- raw_data[[h]] 
  rand_genus <- drop_na(rand_genus, c("FAD_bin"))
  rand_genus$LAD_bin[rand_genus$LAD_bin==95] <- 94
  

  # build dummy for calculation
  dumbo <- rand_genus
  rand_genus <- rand_genus[,c("Genus","FAD_bin","LAD_bin")]
  
  # add bins with 0 and fill in with 1 for extinction and origination
  namevector <- as.character(c(1:94))
  rand_genus[ , namevector] <- NA
  rand_genus <- rand_genus[, c(4:length(rand_genus[0,]))]
  
  for (i in 1:length(rand_genus[,1])) {
    rand_genus[i,dumbo[i,"FAD_bin"]] <- 0 
    rand_genus[i,dumbo[i,"LAD_bin"]] <- 1
    ifelse(dumbo[i,7]!= dumbo[i,"LAD_bin"]-1 & dumbo[i,"FAD_bin"]!= dumbo[i,"LAD_bin"],
           rand_genus[i,(dumbo[i,"FAD_bin"]+1):(dumbo[i,"LAD_bin"]-1)] <- 0, NA)
    ifelse(dumbo[i,"LAD_bin"] == 94, rand_genus[i,dumbo[i,"LAD_bin"]] <- 0, 
           rand_genus[i,dumbo[i,"LAD_bin"]] <- 1)
  }
  
  
  # bind it again 
  rownames(rand_genus) <-  make.names(dumbo[,"Genus"], unique=TRUE)
  
  
  # Reorganise ranges through time data so we can bind it to temperature 
  # data
  rand_genus %<>%  as_tibble() %>%
    mutate(Genus = rownames(rand_genus)) %>%
    group_by(Genus) %>%
    gather(-Genus, key="bins", value="extinct", na.rm = T, factor_key = T) %>%
    arrange(desc(bins), .by_group = T) 
  
  
  # read in veizer and prokoph temperature data, already prepared
  isotemp2 <- read.csv(file="isotemp2.csv", header=TRUE,row.names = 1) 
  isotemp2$bins <- as.factor(isotemp2$bins)
  
  # Now bind the two
  rand_genus <- full_join(rand_genus, isotemp2)
  
  
  #add taxonomical levels
  taxonomy <- dumbo[, c("Class","Order","Family", "Genus")]
  
  # Now bind the two
  rand_genus <- full_join(rand_genus, taxonomy)
  
  rand_genus <- rand_genus %>%
    select(Class, Order, Family, Genus, bins, age, extinct, Temp.mean, change.prev)
  
  
  
  # Calculate long term changes in temp (1 stage to 10 stages)
  # Use lm() to determine slope of the long-term  temperature trends
  # If we want to investigate the interaction between long term change and short term,
  # we need to exclude the short term bin from the long term (see Manuel's figures),
  # and thus shift the long term lm result up by +1. This will be messy, I'm sorry
  ori_bin <- rand_genus %>%
    select(Genus, bins, age) %>%
    group_by(Genus) %>%
    slice(which.min(bins)) 
  colnames(ori_bin) <- c("Genus", "ori.bin", "ori.age" )
  
  # Now bind the two
  rand_genus <- full_join(rand_genus, ori_bin)
  
  # Set bins back to numeric
  rand_genus %<>% transform(bins = as.numeric(as.character(bins)))
  rand_genus$trend1 <- NA
  rand_genus$trend2 <- NA
  rand_genus$trend3 <- NA
  rand_genus$trend4 <- NA
  
  rand_genus <- drop_na(rand_genus, c("Class","Order", "Family", "Genus",
                                      "bins", "age"))
  
  for (i in 1:length(rand_genus$Genus)) {
    bin1 <- c(rand_genus$ori.bin[i], 
              ifelse(rand_genus$bins[i]>rand_genus$ori.bin[i], 
                     rand_genus$bins[i]-1, 
                     rand_genus$bins[i]))
    age1 <- c(isotemp2$age[isotemp2$Stage %in% bin1[1]:bin1[2]])
    
    temp1 <-c(isotemp2$Temp[isotemp2$Stage %in% bin1[1]:bin1[2]])
    lin1 <- lm(temp1~age1)
    rand_genus$trend1[i] <- -lin1$coefficients[2]
  }
  
  #remov ori age and bin
  rand_genus <- rand_genus[, !(names(rand_genus) %in% c("ori.bin","ori.age"))]
  
  #now do the same for trend2
  # first build a data frame where we save the origin of each family
  
  fam <- unique(taxonomy$Family)
  fam_dat <- data.frame(matrix(ncol =2 , nrow = length(fam)))
  colnames(fam_dat) <- c("ori.age", "ori.bin")
  rownames(fam_dat) <-  c(as.character(fam))
  
  
  for (i in 1:length(fam_dat$ori.age)) {
    fam_dat[i, 2] <- min(dumbo$FAD_bin[dumbo$Family == fam[i]])
    fam_dat[i, 1] <- isotemp2$age[isotemp2$Stage == fam_dat[i, 2]]
  }
  
  
  # add a colum with the family names for binding
  fam_dat$Family <- rownames(fam_dat)
  
  # Now bind the two
  rand_genus <- full_join(rand_genus, fam_dat)
  
  # now calculate the trend starting at origin of the specific family
  for (i in 1:length(rand_genus$Genus)) {
    bin1 <- c(rand_genus$ori.bin[i], 
              ifelse(rand_genus$bins[i]>rand_genus$ori.bin[i], 
                     rand_genus$bins[i]-1, 
                     rand_genus$bins[i]))
    age1 <- c(isotemp2$age[isotemp2$Stage %in% bin1[1]:bin1[2]])
    
    temp1 <-c(isotemp2$Temp[isotemp2$Stage %in% bin1[1]:bin1[2]])
    lin1 <- lm(temp1~age1)
    rand_genus$trend2[i] <- -lin1$coefficients[2]
  }
  
  
  #remov ori age and bin
  rand_genus <- rand_genus[,!(names(rand_genus) %in% c("ori.bin","ori.age"))]
  
  #now do the same for trend3 Order
  # first build a data frame where we save the origin of each Order
  
  ord <- unique(dumbo$Order)
  ord_dat <- data.frame(matrix(ncol =2 , nrow = length(ord)))
  colnames(ord_dat) <- c("ori.age", "ori.bin")
  rownames(ord_dat) <-  c(as.character(ord))
  
  
  for (i in 1:length(ord_dat$ori.age)) {
    ord_dat[i, 2] <- min(dumbo$FAD_bin[dumbo$Order == ord[i]])
    ord_dat[i, 1] <- isotemp2$age[isotemp2$Stage == ord_dat[i, 2]]
  }
  
  
  # add a colum with the Order names for binding
  ord_dat$Order <- rownames(ord_dat)
  
  # Now bind the two
  rand_genus <- full_join(rand_genus, ord_dat)
  
  # now calculate the trend starting at origin of the specific Order
  for (i in 1:length(rand_genus$Genus)) {
    bin1 <- c(rand_genus$ori.bin[i], 
              ifelse(rand_genus$bins[i]>rand_genus$ori.bin[i], 
                     rand_genus$bins[i]-1, rand_genus$bins[i]))
    age1 <- c(isotemp2$age[isotemp2$Stage %in% bin1[1]:bin1[2]])
    
    temp1 <-c(isotemp2$Temp[isotemp2$Stage %in% bin1[1]:bin1[2]])
    lin1 <- lm(temp1~age1)
    rand_genus$trend3[i] <- -lin1$coefficients[2]
  }
  
  
  #remov ori age and bin
  rand_genus <- rand_genus[,!(names(rand_genus) %in% c("ori.bin","ori.age"))]
  
  #now do the same for trend4 Class
  # first build a data frame where we save the origin of each Class
  
  cla <- unique(dumbo$Class)
  cla_dat <- data.frame(matrix(ncol =2 , nrow = length(cla)))
  colnames(cla_dat) <- c("ori.age", "ori.bin")
  rownames(cla_dat) <-  c(as.character(cla))
  
  
  for (i in 1:length(cla_dat$ori.age)) {
    cla_dat[i, 2] <- min(dumbo$FAD_bin[dumbo$Class == cla[i]])
    cla_dat[i, 1] <- isotemp2$age[isotemp2$Stage == cla_dat[i, 2]]
  }
  
  
  # add a colum with the Class names for binding
  cla_dat$Class <- rownames(cla_dat)
  
  # Now bind the two
  rand_genus <- full_join(rand_genus, cla_dat)
  
  # now calculate the trend starting at origin of the specific Class
  for (i in 1:length(rand_genus$Genus)) {
    bin1 <- c(rand_genus$ori.bin[i], 
              ifelse(rand_genus$bins[i]>rand_genus$ori.bin[i], 
                     rand_genus$bins[i]-1, rand_genus$bins[i]))
    age1 <- c(isotemp2$age[isotemp2$Stage %in% bin1[1]:bin1[2]])
    
    temp1 <-c(isotemp2$Temp[isotemp2$Stage %in% bin1[1]:bin1[2]])
    lin1 <- lm(temp1~age1)
    rand_genus$trend4[i] <- -lin1$coefficients[2]
  }
  
  #remov ori age and bin
  rand_genus <- rand_genus[,!(names(rand_genus) %in% c("ori.bin","ori.age"))]
  
  
  # load the temperature trends, already calculated
  temp_trends <- read.csv("temp_trends.csv", header = T, row.names = 1)
  
  # Set bins to factor
  rand_genus %<>% transform(bins = as.factor(as.character(bins)))
  temp_trends %<>% transform(bins = as.factor(as.character(bins)))
  
  temp_trends <- temp_trends[, -c(1,3)]
  
  
  # Now bind the two
  fourteen_trends_rand <- full_join(rand_genus, temp_trends)
  
  # Set bins back to numeric
  fourteen_trends_rand %<>% transform(bins = as.numeric(as.character(bins)))
  
  
  finished_data[[h]] <- fourteen_trends_rand 
  
  progress(h, max.value = 1000)
  if (h == 1000) message("Done!")
}

# save(finished_data, file="sim_data_1000.RData")
# load("sim_data_1000.RData")

# Cleanup
rm(list=ls()[! ls() %in% c("finished_data")])




# 6.3 Calculate GLMM's and estimate effect size based on these sim --------


# built empty data frame
result_rand <- data.frame(lower.CI.w = numeric(), estimate.w = numeric(), upper.CI.w =numeric(),
                          lower.CI.c = numeric(), estimate.c = numeric(), upper.CI.c =numeric())

for (h in 1:1000) {
  
  fourteen_trends_rand <- finished_data[[h]]
  
  # replace the NA's with zeros as we have no trend there
  fourteen_trends_rand$trend1[is.na(fourteen_trends_rand$trend1)] <- 0
  fourteen_trends_rand$trend2[is.na(fourteen_trends_rand$trend2)] <- 0
  fourteen_trends_rand$trend3[is.na(fourteen_trends_rand$trend3)] <- 0
  fourteen_trends_rand$trend4[is.na(fourteen_trends_rand$trend4)] <- 0
  
  # Split short term temperature change into warming and cooling:
  fourteen_trends_rand$cooling<-ifelse(fourteen_trends_rand$change.prev<0, fourteen_trends_rand$change.prev, NA)
  fourteen_trends_rand$warming<-ifelse(fourteen_trends_rand$change.prev>0, fourteen_trends_rand$change.prev, NA)
  
  ## Warming
  # Iterate through each long term temperature change
  vars = names(dplyr::select(fourteen_trends_rand, trend1:trend4, trend.st1:trend.st10)) 
  modelswr = lapply(setNames(vars, vars), function(var) {
    form = paste("extinct~warming:", var, "+(1|Genus)")
    glmer(form, data=fourteen_trends_rand, family="binomial")
  })
  
  # Make data frame for model output
  warm_fourteen_trends_rand <- data.frame(model=names(dplyr::select(fourteen_trends_rand, trend1:trend4, trend.st1:trend.st10)), 
                                          intercept=NA, 
                                          interaction=NA, AIC=NA)
  
  # Run loop to fill for.warm (coefficients and p-values)
  for (i in warm_fourteen_trends_rand$model) {
    sum <- summary(modelswr[[i]])
    warm_fourteen_trends_rand[warm_fourteen_trends_rand$model==i, "intercept"]<-
      paste(round(sum$coefficients[1, 1], 3), sep=" ", "±",
            round(sum$coefficients[1, 2], 3),
            ifelse(sum$coefficients[1, 4]<0.001, "***", 
                   ifelse(sum$coefficients[1, 4]<0.01, "**", 
                          ifelse(sum$coefficients[1, 4]<0.05, "*", ""))))
    warm_fourteen_trends_rand[warm_fourteen_trends_rand$model==i, "interaction"]<-
      paste(round(sum$coefficients[2, 1], 3), sep=" ", "±",
            round(sum$coefficients[2, 2], 3),
            ifelse(sum$coefficients[2, 4]<0.001, "***", 
                   ifelse(sum$coefficients[2, 4]<0.01, "**", 
                          ifelse(sum$coefficients[2, 4]<0.05, "*", ""))))
    warm_fourteen_trends_rand[warm_fourteen_trends_rand$model==i, "AIC"]<-as.numeric(round(sum$AICtab[[1]], 1))
  }
  
  #  Add column weith AIC weights
  warm_fourteen_trends_rand$dAIC<- as.numeric(round(aicw(warm_fourteen_trends_rand$AIC)$delta, 1))
  warm_fourteen_trends_rand$AIC.weights<- as.numeric(signif(aicw(warm_fourteen_trends_rand$AIC)$w, 3))
  
  
  # let's take the values of the visreg and split them into cooling and warming and compare the 
  # exctinction response of each by means of boxplots
  vis_out_wr <- visreg(modelswr[[which(warm_fourteen_trends_rand$dAIC==0)[1]]], "warming", scale="response", 
                       by=paste(warm_fourteen_trends_rand[which(warm_fourteen_trends_rand$dAIC==0)[1],1]), 
                       breaks = c(-0.1, 0, 0.1),
                       rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5), plot = F) 
  
  
  
  vis_out_wr$fit$ctrend <- ifelse(vis_out_wr$fit[2]<0,"Trend = Cooling", "Trend = Warming")
  
  wtestw <- wilcox.test(vis_out_wr$fit$visregFit[vis_out_wr$fit$ctrend=="Trend = Warming"],
                        vis_out_wr$fit$visregFit[vis_out_wr$fit$ctrend=="Trend = Cooling"], 
                        paired = F, conf.int = T)
  
  
  
  result_rand[h,1] <- wtestw$conf.int[1]
  result_rand[h,2] <- wtestw$estimate
  result_rand[h,3] <- wtestw$conf.int[2]
  
  
  ## Cooling
  # Iterate through each long term temperature change
  vars = names(dplyr::select(fourteen_trends_rand, trend1:trend4, trend.st1:trend.st10)) 
  modelscr = lapply(setNames(vars, vars), function(var) {
    form = paste("extinct~cooling:", var, "+(1|Genus)")
    glmer(form, data=fourteen_trends_rand, family="binomial")
  })
  
  # Make data frame for model output
  cool_fourteen_trends_rand <-data.frame(model=names(dplyr::select(fourteen_trends_rand, c(trend1:trend4, trend.st1:trend.st10))), 
                                         intercept=NA, 
                                         interaction=NA, AIC=NA)
  
  # Run loop to fill for.cool (coefficients and p-values)
  for (i in cool_fourteen_trends_rand$model) {
    sum <- summary(modelscr[[i]])
    cool_fourteen_trends_rand[cool_fourteen_trends_rand$model==i, "intercept"]<-
      paste(round(sum$coefficients[1, 1], 3), sep=" ", "±",
            round(sum$coefficients[1, 2], 3),
            ifelse(sum$coefficients[1, 4]<0.001, "***", 
                   ifelse(sum$coefficients[1, 4]<0.01, "**", 
                          ifelse(sum$coefficients[1, 4]<0.05, "*", ""))))
    cool_fourteen_trends_rand[cool_fourteen_trends_rand$model==i, "interaction"]<-
      paste(round(sum$coefficients[2, 1], 3), sep=" ", "±",
            round(sum$coefficients[2, 2], 3),
            ifelse(sum$coefficients[2, 4]<0.001, "***", 
                   ifelse(sum$coefficients[2, 4]<0.01, "**", 
                          ifelse(sum$coefficients[2, 4]<0.05, "*", ""))))
    cool_fourteen_trends_rand[cool_fourteen_trends_rand$model==i, "AIC"]<-as.numeric(round(sum$AICtab[[1]], 1))
  }
  
  # Add column weith AIC weights
  cool_fourteen_trends_rand$dAIC<- as.numeric(round(aicw(cool_fourteen_trends_rand$AIC)$delta, 1))
  cool_fourteen_trends_rand$AIC.weights<- as.numeric(signif(aicw(cool_fourteen_trends_rand$AIC)$w, 3))
  
  
  # let's take the values of the visreg and split them into cooling and warming and compare the 
  # exctinction response of each by means of boxplots
  vis_out_cr <- visreg(modelscr[[which(cool_fourteen_trends_rand$dAIC==0)[1]]], "cooling", scale="response", 
                       by=paste(cool_fourteen_trends_rand[which(cool_fourteen_trends_rand$dAIC==0)[1],1]), 
                       breaks = c(-0.1, 0, 0.1),
                       rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5), plot = F) 
  
  
  
  vis_out_cr$fit$ctrend <- ifelse(vis_out_cr$fit[2]<0,"Trend = Cooling", "Trend = Warming")
  
  
  wtestc <- wilcox.test(vis_out_cr$fit$visregFit[vis_out_cr$fit$ctrend=="Trend = Warming"],
                        vis_out_cr$fit$visregFit[vis_out_cr$fit$ctrend=="Trend = Cooling"], 
                        paired = F, conf.int = T)
  
  result_rand[h,4] <- wtestc$conf.int[1]
  result_rand[h,5] <- wtestc$estimate
  result_rand[h,6] <- wtestc$conf.int[2]
  
  progress(h, max.value = 1000)
  if (h == 1000) message("Done!")
}

# save(result_rand, file="result_rand_1000.RData")
# write.csv(result_rand, file = "result_rand_1000_wo_seed.csv")

# fun::shutdown()

# Cleanup
rm(list=ls()[! ls() %in% c("finished_data", "result_rand")])




# 6.4 Plot GLMM results of Null model -------------------------------------


load("result_rand_1000.RData")

result_rand$trial <- rownames(result_rand)

# define colors
library(wesanderson)
col2 <- wes_palette("Darjeeling1")

my_theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  title = element_text(size = 12),
  panel.border = element_blank(),
  #  panel.background = element_rect(fill= "white", colour= "black"),
  panel.grid.major = element_line(colour= "grey80"),
  panel.grid.minor = element_blank())

# first warming
ggplot(result_rand, aes(trial, estimate.w))+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin= min(result_rand$lower.CI.w),
                ymax= max(result_rand$upper.CI.w)), fill = "darkgrey", 
            alpha = 0.1)+
  geom_hline(yintercept= min(result_rand$lower.CI.w), linetype="solid", 
             color= col2[3])+
  geom_hline(yintercept= max(result_rand$upper.CI.w), linetype="solid", 
             color= col2[3])+
  geom_hline(yintercept=0)+
  geom_point(color= "grey40", fill ="coral2", alpha = 0.8, size=2, pch= 21)+ 
  theme_classic()+
  ylim(-0.12,0.12)+
  scale_x_discrete(breaks=seq(0, 1000, 100))+
  labs(y = "Change in extinction risk [%]", x= "number of trial")+
  my_theme

# save it to working directory
ggsave("warm_null_model.pdf")

# then cooling
ggplot(result_rand, aes(trial, estimate.c))+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin= min(result_rand$lower.CI.c),
                ymax= max(result_rand$upper.CI.c)), fill = "darkgrey", 
            alpha = 0.1)+
  geom_hline(yintercept= min(result_rand$lower.CI.c), linetype="solid", 
             color= col2[3])+
  geom_hline(yintercept= max(result_rand$upper.CI.c), linetype="solid", 
             color= col2[3])+
  geom_hline(yintercept=0)+
  geom_point(color= "grey40", fill ="#56B4E9", alpha = 0.8, size=2, pch= 21)+
  theme_classic()+
  ylim(-0.12,0.12)+
  scale_x_discrete(breaks=seq(0, 1000, 100))+
  labs(y = "Change in extinction risk [%]", x= "number of trial")+
  my_theme

# save it to working directory
ggsave("cool_null_model.pdf")

