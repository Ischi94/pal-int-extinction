library(colorednoise)
library(lme4)
library(tidyverse)
library(geiger)
library(here)


# # import stage data  ----------------------------------------------------

# import stage data from Gradstein 2012 for binning
gradstein <- read_csv(here("data/gradstein.csv")) 

# Functions ---------------------------------------------------------------

### model_df

# take a model and produce a Data frame (tibble for nicer printing) with model 
# output and AIC, BIC, deltaAIC, AICweights and deltaBIC
model_df <- function(my_model) {
  # make data frame for model output
  df <-
    tibble(
      model = dat_trends %>% select(starts_with("trend.")) %>% names(),
      intercept = character(length(model)),
      interaction = character(length(model)),
      AIC = numeric(length(model)),
      BIC = numeric(length(model))
    )
  # Run loop to fill interaction_warm (coefficients and p-values)
  for (i in df$model) {
    sum <- summary(my_model[[i]])
    df[df$model == i, "intercept"] <-
      paste(
        round(sum$coefficients[1, 1], 2),
        sep = " ",
        "+-",
        round(sum$coefficients[1, 2], 2),
        ifelse(
          sum$coefficients[1, 4] < 0.001,
          "***",
          ifelse(
            sum$coefficients[1, 4] < 0.01,
            "**",
            ifelse(sum$coefficients[1, 4] < 0.05, "*", "")
          )
        )
      )
    df[df$model == i, "interaction"] <-
      paste(
        round(sum$coefficients[2, 1], 2),
        sep = " ",
        "+-",
        round(sum$coefficients[2, 2], 2),
        ifelse(
          sum$coefficients[2, 4] < 0.001,
          "***",
          ifelse(
            sum$coefficients[2, 4] < 0.01,
            "**",
            ifelse(sum$coefficients[2, 4] < 0.05, "*", "")
          )
        )
      )
    df[df$model == i, "AIC"] <- as.numeric(round(sum$AICtab[[1]], 1))
    df[df$model == i, "BIC"] <- as.numeric(round(sum$AICtab[[2]], 1))
  }
  
  # Add column weith AIC weights and deltaBIC
  df$dAIC <- as.numeric(round(aicw(df$AIC)$delta, 1))
  df$AICweights <- as.numeric(signif(aicw(df$AIC)$w, 3))
  df$dBIC <- as.numeric(round(aicw(df$BIC)$delta, 1))
  df
}

### simulate noise trends

# simulate trends with the same parameters as empirical data but autocorrelated, 
# choose between red and blue noise
sim_noise_trends <- function(noise, nr_observations) {
  # set type of autorrelation/ noise
  if(noise == "red") phi_val <- runif(1, min = 0, max = 1)
  if(noise == "blue") phi_val <- runif(1, min = -1, max = 0)
  
  ### autocorrelated temperature ###
  
  # generate autocorrelated temperature time series for 14:94 stages, include stage 
  # 95 as detrending will remove one stage
  temperature <- colored_noise(timesteps = length(95:14), 
                               mean = 0, sd = 5, phi = phi_val) %>% 
    # detrend temperature data
    diff(differences = 1)
  
  
  # bin to stages
  isotemp <- tibble(stage = 94:14, temp = temperature) %>% 
    # add 20 to all values to get reasonable outcome space for temp, this does not
    # change the autocorrelation
    mutate(temp = temp + 20) %>% 
    # add age from gradstein and init new column
    add_column(age = gradstein$mid[94:14], change_prev = double(length(94:14)))
  
  
  # set up new dfr
  isotemp_trends <- isotemp
  
  # calculate change prev using linear regression
  for (i in unique(isotemp$stage)) {
    dum1 <- filter(isotemp, isotemp$stage %in% isotemp$stage[between(isotemp$stage, i-1, i)])
    dum2 <-  lm(temp ~ age, data = dum1)
    isotemp_trends[isotemp_trends$stage == i, "change_prev"] <- -dum2$coefficients[2]
  }
  
  # arrange dfr and calculate lags
  isotemp_trends <- isotemp_trends %>% 
    arrange(desc(stage)) %>% 
    mutate(lag1 = lag(temp, order_by = stage), 
           lag2 = lag(lag1, order_by = stage),
           lag3 = lag(lag2, order_by = stage),
           lag4 = lag(lag3, order_by = stage),
           lag5 = lag(lag4, order_by = stage),
           lag6 = lag(lag5, order_by = stage),
           lag7 = lag(lag6, order_by = stage),
           lag8 = lag(lag7, order_by = stage),
           lag9 = lag(lag8, order_by = stage),
           lag10 = lag(lag9, order_by = stage)) 
  
  
  # calculate long-term trends based on lags
  for (i in unique(isotemp_trends$stage)) {
    sub1 <- filter(isotemp_trends,
                   isotemp_trends$stage %in% isotemp_trends$stage[between(isotemp_trends$stage, i -
                                                                            1, i)])
    lin1 <- lm(temp ~ age, data = sub1)
    isotemp_trends[isotemp_trends$stage == i + 1, "trend.st1"] <-
      -lin1$coefficients[2]
    sub2 <- filter(isotemp_trends,
                   isotemp_trends$stage %in% isotemp_trends$stage[between(isotemp_trends$stage, i -
                                                                            2, i)])
    lin2 <- lm(temp ~ age, data = sub2)
    isotemp_trends[isotemp_trends$stage == i + 1, "trend.st2"] <-
      -lin2$coefficients[2]
    sub3 <-
      filter(isotemp_trends,
             isotemp_trends$stage %in% isotemp_trends$stage[between(isotemp_trends$stage, i -
                                                                      3, i)])
    lin3 <- lm(temp ~ age, data = sub3)
    isotemp_trends[isotemp_trends$stage == i + 1, "trend.st3"] <-
      -lin3$coefficients[2]
    sub4 <-
      filter(isotemp_trends,
             isotemp_trends$stage %in% isotemp_trends$stage[between(isotemp_trends$stage, i -
                                                                      4, i)])
    lin4 <- lm(temp ~ age, data = sub4)
    isotemp_trends[isotemp_trends$stage == i + 1, "trend.st4"] <-
      -lin4$coefficients[2]
    sub5 <-
      filter(isotemp_trends,
             isotemp_trends$stage %in% isotemp_trends$stage[between(isotemp_trends$stage, i -
                                                                      5, i)])
    lin5 <- lm(temp ~ age, data = sub5)
    isotemp_trends[isotemp_trends$stage == i + 1, "trend.st5"] <-
      -lin5$coefficients[2]
    sub6 <-
      filter(isotemp_trends,
             isotemp_trends$stage %in% isotemp_trends$stage[between(isotemp_trends$stage, i -
                                                                      6, i)])
    lin6 <- lm(temp ~ age, data = sub6)
    isotemp_trends[isotemp_trends$stage == i + 1, "trend.st6"] <-
      -lin6$coefficients[2]
    sub7 <-
      filter(isotemp_trends,
             isotemp_trends$stage %in% isotemp_trends$stage[between(isotemp_trends$stage, i -
                                                                      7, i)])
    lin7 <- lm(temp ~ age, data = sub7)
    isotemp_trends[isotemp_trends$stage == i + 1, "trend.st7"] <-
      -lin7$coefficients[2]
    sub8 <-
      filter(isotemp_trends,
             isotemp_trends$stage %in% isotemp_trends$stage[between(isotemp_trends$stage, i -
                                                                      8, i)])
    lin8 <- lm(temp ~ age, data = sub8)
    isotemp_trends[isotemp_trends$stage == i + 1, "trend.st8"] <-
      -lin8$coefficients[2]
    sub9 <-
      filter(isotemp_trends,
             isotemp_trends$stage %in% isotemp_trends$stage[between(isotemp_trends$stage, i -
                                                                      9, i)])
    lin9 <- lm(temp ~ age, data = sub9)
    isotemp_trends[isotemp_trends$stage == i + 1, "trend.st9"] <-
      -lin9$coefficients[2]
    sub10 <-
      filter(isotemp_trends,
             isotemp_trends$stage %in% isotemp_trends$stage[between(isotemp_trends$stage, i -
                                                                      10, i)])
    lin10 <- lm(temp ~ age, data = sub10)
    isotemp_trends[isotemp_trends$stage == i + 1, "trend.st10"] <-
      -lin10$coefficients[2]
  }
  
  
  
  ### autocorrelated extinction ###
  
  # set number of observations
  nr_obs <- nr_observations
  
  # get autocorrelated extinction signal
  extinction <- colored_noise(timesteps = nr_obs, mean = 0, 
                              sd = 0.5, phi = phi_val) %>% 
    # detrend extinction data
    diff(differences = 1) 
  
  # get autocorrelated stages  
  stage <- colored_noise(timesteps = nr_obs, mean = 30, sd = 40, phi = phi_val) %>% 
    round(0)
  
  
  dat_trends <- extinction %>% 
    # convert to extinction signal, every value above 0.5 signals extinction
    # data is still autocorrelated
    as_tibble_col(column_name = "extinction") %>% 
    mutate(extinction = if_else(extinction < 0.5, 0, 1)) %>% 
    # assign genera based on autocorrelated extinction signal
    # set first value to 1 to allow proper cumsum counting
    add_row(extinction = 1, .before = 0) %>% 
    mutate(name_var = cumsum(extinction), # get genus names
           genus_dum = "genus") %>% # set dummy column
    # unite to get genus names
    unite(col = "genus", genus_dum:name_var)  %>% 
    # assign stages
    add_column(stage = stage) %>% 
    # get stages in right order, count backwards from extinction until origination
    group_by(genus) %>% 
    mutate(first_val = first(stage)) %>% 
    add_column(number = 1) %>% 
    mutate(ticker = cumsum(number), # get dummy variable
           ticker = ticker -1, 
           # get correct stage order
           stage = first_val - ticker) %>% 
    # remove dummy variable
    select(genus, extinction, stage) %>% 
    # remove singletons
    add_count(genus) %>% 
    filter(n != 1) %>% 
    # subset to stages where we got temperature information 14:94
    filter(stage >= 14 & stage <= 94) %>% 
    # remove dummy variables
    select(-n) %>% 
    # remove grouping
    ungroup()
  
  
  
  # Now bind range data with temperature data
  dat_trends <- full_join(dat_trends, isotemp_trends) %>% 
    # remove redundant colums from merging
    drop_na(genus) %>% 
    # order it properly
    select(genus, stage, extinction, temp, 
           change_prev, trend.st1:trend.st10)
  
  dat_trends
}






# simulate trends ---------------------------------------------------------

# the output of the simulations for null testing was 900 models
# as we want to compare the distributions of both simulation approaches, 
# we again simulate 900 models 
nr_models <- 900


# red noise simulations --------------------------------------------------


# build empty list
rep_data <- list()

# set progress bar
pb <- progress::progress_bar$new(total = nr_models) 

for (i in 1:nr_models) {
  # simulate trends
  # but suppress messages to see progress bar
  rep_data[[i]] <- suppressMessages(
    sim_noise_trends(noise = "red", nr_observations = 5000))
  
  # display progress
  pb$tick()
}

# save data for reproduction
save(rep_data, 
     file = here("data/simulation-results/simulated_rednoise_rep.RData"))


# select only those trends with a sufficient number of observations (<= 1000)

# count lengths
data_length <- sapply(rep_data, nrow)

# get boolean vector for subsetting
data_length_low <- data_length >= 1000

# subset
rep_data <- rep_data[data_length_low]

# calculate glmms ---------------------------------------------------------


# set up dfr for results
noise_results <- tibble(model = 1:length(rep_data), 
                        warm = as.double(0), warm_low = as.double(0), 
                        warm_high = as.double(0), 
                        cool = as.double(0), cool_low = as.double(0), 
                        cool_high = as.double(0)) 


# set progress bar
pb <- progress::progress_bar$new(total = length(rep_data)) 



for (i in 1:length(rep_data)) {
  
  # assign data
  dat_trends <- rep_data[[i]]
  
  # Iterate through each trend
  vars <- dat_trends %>% select(starts_with("trend.")) %>% names()
  
  interaction_model <- lapply(setNames(vars, vars), function(var) {
    form = paste("extinction~change_prev:", var, "+(stage|genus)")
    suppressMessages(glmer(form, data=dat_trends, family="binomial"))
  })
  
  # Make data frame for model output
  interaction_df <- model_df(interaction_model)
  
  # choose final model
  interaction_final <- interaction_model[[which(interaction_df$dAIC==0)[[1]]]]
  
  
  # make informed predictions based on a subset (palaeoclimate interaction)
  # type = response gives us probability instead of log Odds
  
  # first get the trends
  trend <- interaction_df$model[[which(interaction_df$dAIC==0)[[1]]]] 
  
  # predict palaeoclimate interaction
  # if too few data, assign average
  
  #  warming warming
  ww_raw <- dat_trends %>% filter(get(trend) >=0 & change_prev >= 0)
  if(length(ww_raw$genus) == 0) { 
    wc_raw <- dat_trends } else{
  ww_pred <- predict(interaction_final, newdata = ww_raw,
                     type = "response", allow.new.levels = TRUE)
  }
  
  
  #  warming cooling 
  wc_raw <- dat_trends %>% filter(get(trend) >=0 & change_prev <= 0)
  if(length(wc_raw$genus) == 0) { 
    wc_raw <- dat_trends } else{
      wc_pred <- predict(interaction_final, newdata = wc_raw,
                         type = "response", allow.new.levels = TRUE)
    }
  
  
  #  cooling cooling 
  cc_raw <- dat_trends %>% filter(get(trend) <= 0 & change_prev <= 0)
  if(length(cc_raw$genus) == 0) { 
    cc_raw <- dat_trends } else{
      cc_pred <- predict(interaction_final, newdata = cc_raw,
                         type = "response", allow.new.levels = TRUE)
    }
  
  
  
  #  cooling warming
  cw_raw <- dat_trends %>% filter(get(trend) <=0 & change_prev >= 0)
  if(length(cw_raw$genus) == 0) { 
    cw_raw <- dat_trends } else{
      cw_pred <- predict(interaction_final, newdata = cw_raw,
                         type = "response", allow.new.levels = TRUE)
    }

  
  
  
  # take the predictions and test whether they significantly are above or below the
  # baseline. This is a one sided Wilcoxon rank sum test 
  # (equivalent to the Mann-Whitney test)
  wilcox_warm <- wilcox.test(ww_pred, wc_pred, paired = F,  conf.int = T)
  wilcox_cool <- wilcox.test(cc_pred, cw_pred, paired = F,  conf.int = T)
  
  
  # assign results
  noise_results[i, "warm"] <- wilcox_warm$estimate
  noise_results[i, "warm_low"] <- wilcox_warm$conf.int[[1]]
  noise_results[i, "warm_high"] <- wilcox_warm$conf.int[[2]]
  noise_results[i, "cool"] <- wilcox_cool$estimate
  noise_results[i, "cool_low"] <- wilcox_cool$conf.int[[1]]
  noise_results[i, "cool_high"] <- wilcox_cool$conf.int[[2]]
  
  
  # display progress
  pb$tick()
}


# save data
write_csv(noise_results, 
          path = here("data/results/autocorrelation_red_noise.csv"))



# blue noise simulations --------------------------------------------------


# build empty list
rep_data <- list()

# set progress bar
pb <- progress::progress_bar$new(total = nr_models) 

for (i in 1:nr_models) {
  # simulate trends
  # but suppress messages to see progress bar
  rep_data[[i]] <- suppressMessages(
    sim_noise_trends(noise = "blue", nr_observations = 5000))
  
  # display progress
  pb$tick()
}

# save data for reproduction
save(rep_data, 
     file = here("data/simulation-results/simulated_bluenoise_rep.RData"))


# select only those trends with a sufficient number of observations (<= 1000)

# count lengths
data_length <- sapply(rep_data, nrow)

# get boolean vector for subsetting
data_length_low <- data_length >= 1000

# subset
rep_data <- rep_data[data_length_low]

# calculate glmms ---------------------------------------------------------


# set up dfr for results
noise_results <- tibble(model = 1:length(rep_data), 
                        warm = as.double(0), warm_low = as.double(0), 
                        warm_high = as.double(0), 
                        cool = as.double(0), cool_low = as.double(0), 
                        cool_high = as.double(0)) 


# set progress bar
pb <- progress::progress_bar$new(total = length(rep_data)) 



for (i in 1:length(rep_data)) {
  
  # assign data
  dat_trends <- rep_data[[i]]
  
  # Iterate through each trend
  vars <- dat_trends %>% select(starts_with("trend.")) %>% names()
  
  interaction_model <- lapply(setNames(vars, vars), function(var) {
    form = paste("extinction~change_prev:", var, "+(stage|genus)")
    suppressMessages(glmer(form, data=dat_trends, family="binomial"))
  })
  
  # Make data frame for model output
  interaction_df <- model_df(interaction_model)
  
  # choose final model
  interaction_final <- interaction_model[[which(interaction_df$dAIC==0)[[1]]]]
  
  
  # make informed predictions based on a subset (palaeoclimate interaction)
  # type = response gives us probability instead of log Odds
  
  # first get the trends
  trend <- interaction_df$model[[which(interaction_df$dAIC==0)[[1]]]] 
  
  # predict palaeoclimate interaction
  # if too few data, assign average
  
  #  warming warming
  ww_raw <- dat_trends %>% filter(get(trend) >=0 & change_prev >= 0)
  if(length(ww_raw$genus) == 0) { 
    wc_raw <- dat_trends } else{
      ww_pred <- predict(interaction_final, newdata = ww_raw,
                         type = "response", allow.new.levels = TRUE)
    }
  
  
  #  warming cooling 
  wc_raw <- dat_trends %>% filter(get(trend) >=0 & change_prev <= 0)
  if(length(wc_raw$genus) == 0) { 
    wc_raw <- dat_trends } else{
      wc_pred <- predict(interaction_final, newdata = wc_raw,
                         type = "response", allow.new.levels = TRUE)
    }
  
  
  #  cooling cooling 
  cc_raw <- dat_trends %>% filter(get(trend) <= 0 & change_prev <= 0)
  if(length(cc_raw$genus) == 0) { 
    cc_raw <- dat_trends } else{
      cc_pred <- predict(interaction_final, newdata = cc_raw,
                         type = "response", allow.new.levels = TRUE)
    }
  
  
  
  #  cooling warming
  cw_raw <- dat_trends %>% filter(get(trend) <=0 & change_prev >= 0)
  if(length(cw_raw$genus) == 0) { 
    cw_raw <- dat_trends } else{
      cw_pred <- predict(interaction_final, newdata = cw_raw,
                         type = "response", allow.new.levels = TRUE)
    }
  
  
  
  
  # take the predictions and test whether they significantly are above or below the
  # baseline. This is a one sided Wilcoxon rank sum test 
  # (equivalent to the Mann-Whitney test)
  wilcox_warm <- wilcox.test(ww_pred, wc_pred, paired = F,  conf.int = T)
  wilcox_cool <- wilcox.test(cc_pred, cw_pred, paired = F,  conf.int = T)
  
  
  # assign results
  noise_results[i, "warm"] <- wilcox_warm$estimate
  noise_results[i, "warm_low"] <- wilcox_warm$conf.int[[1]]
  noise_results[i, "warm_high"] <- wilcox_warm$conf.int[[2]]
  noise_results[i, "cool"] <- wilcox_cool$estimate
  noise_results[i, "cool_low"] <- wilcox_cool$conf.int[[1]]
  noise_results[i, "cool_high"] <- wilcox_cool$conf.int[[2]]
  
  
  # display progress
  pb$tick()
}


# save data
write_csv(noise_results, 
          path = here("data/results/autocorrelation_blue_noise.csv"))

