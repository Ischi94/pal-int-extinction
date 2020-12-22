# loading packages
library(tidyverse)
library(lme4)
library(geiger)
library(here)


# functions ---------------------------------------------------------------
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

###


### glmm_analysis

# calculate models and make output list
glmm_analysis <- function(dat_subset){
  
  # Iterate through each warming
  vars <- dat_subset %>% select(starts_with("trend.")) %>% names()
  interaction_model <- lapply(setNames(vars, vars), function(var) {
    form = paste("extinction~change_prev:", var, "+(stage|genus)")
    glmer(form, data=dat_subset, family="binomial")
  })
  
  
  # Make data frame for model output
  interaction_df <- model_df(interaction_model)
  
  # choose final model
  interaction_final <- interaction_model[[which(interaction_df$dAIC==0)[[1]]]]
  
  
  # make informed predictions based on a subset (palaeoclimate interaction)
  # type = response gives us probability instead of log Odds
  
  # first get the trends
  trend <- interaction_df$model[[which(interaction_df$dAIC==0)[[1]]]]
  
  
  #  warming warming
  ww_raw <- dat_subset %>% filter(get(trend) >=0 & change_prev >= 0)
  ww_pred <- predict(interaction_final, newdata = ww_raw,
                     type = "response", allow.new.levels = TRUE)
  
  
  #  warming cooling 
  wc_raw <- dat_subset %>% filter(get(trend) >=0 & change_prev <= 0)
  wc_pred <- predict(interaction_final, newdata = wc_raw,
                     type = "response", allow.new.levels = TRUE)
  
  
  #  cooling cooling 
  cc_raw <- dat_subset %>% filter(get(trend) <= 0 & change_prev <= 0)
  cc_pred <- predict(interaction_final, newdata = cc_raw,
                     type = "response", allow.new.levels = TRUE)
  
  
  #  cooling warming
  cw_raw <- dat_subset %>% filter(get(trend) <=0 & change_prev >= 0)
  cw_pred <- predict(interaction_final, newdata = cw_raw,
                     type = "response", allow.new.levels = TRUE)
  
  
  
  
  
  # take the predictions and test whether they significantly are above or below the
  # baseline. This is a one sided Wilcoxon rank sum test 
  # (equivalent to the Mann-Whitney test)
  wilcox_warm <- wilcox.test(ww_pred, wc_pred, paired = F,  conf.int = T)
  wilcox_cool <- wilcox.test(cc_pred, cw_pred, paired = F,  conf.int = T)
  
  # make output vector and return it
  output_vct <- c(warm = wilcox_warm$estimate, 
                  cool = wilcox_cool$estimate)
  
  output_vct
}

#####



# read in data ------------------------------------------------------------


# import empirical data set of genus durations for clades used in our analysis
all_data_trends <- read_csv(here("data/all_data_trends.csv"), guess_max = 1e5)

# calculate durations
durations <- all_data_trends %>% 
  group_by(genus) %>% 
  # calculate durations
  mutate(min_age = min(age), max_age = max(age), 
         duration = max_age - min_age) %>% 
  ungroup() %>% 
  # remove duplicates
  distinct(genus, .keep_all = TRUE) %>% 
  # select only the durations as vector
  pull(duration)

# import stage data 
gradstein <- read_csv(here("data/gradstein.csv"))


# create temperature independant data -------------------------------------

# set up data frame for simulated data
# maximum number of observations
n_max <- 3000

sim_data <- tibble(
  # assign genus names 
  genus = paste0("Genus", 1:n_max), 
  # simulate random lad ages from uniform distribution
  fad_age = runif(n_max, min = -0.01, max = 481.55) %>% 
                     round(digits = 2)) %>% 
  # draw from the durations of all genera from our observed data 
  # and subtract it from the fad to get lad
  mutate(lad_age = fad_age - sample(durations, n_max, replace = TRUE)) %>% 
  # as we can't have future ranges, we replace every negative LAD with 0 
  # (= range to recent)
  mutate(lad_age = if_else(lad_age < 0, 0, lad_age))

  
# group into higher taxonomy

# first let's take a look at the empirical data

# number of genera (mean) per order
all_data_trends %>% 
  distinct(genus, .keep_all = TRUE) %>% 
  count(order) %>% 
  summarise(mean_order = mean(n)) # circa 40

# number of genera (mean) per family
all_data_trends %>% 
  distinct(genus, .keep_all = TRUE) %>% 
  count(family) %>% 
  summarise(mean_order = mean(n)) # circa 5

# group the simulated data based on these empirical numbers
sim_data <- sim_data %>% 
  add_column(order = paste0("Order", 1:(n_max/40)) %>% rep(each = 40), 
             family = paste0("Family", 1:(n_max/5)) %>% rep(each = 5), 
             .before = "genus")

# bin simulated data 
sim_data$fad_bin <- NA
sim_data$lad_bin <- NA

for (i in 1:length(sim_data$genus)) {
  for (j in 1:length(gradstein$stg)) {
    ifelse(gradstein$bottom[j] >= sim_data$fad_age[i] & gradstein$top[j] < sim_data$fad_age[i],
           sim_data$fad_bin[i] <- gradstein$stg[j], NA)
    ifelse(gradstein$bottom[j] > sim_data$lad_age[i] & gradstein$top[j] <= sim_data$lad_age[i],
           sim_data$lad_bin[i] <- gradstein$stg[j], NA)
  }}

# clean data
sim_data <- sim_data %>% 
  # remove singletons as they add noise 
  filter(fad_bin != lad_bin) %>% 
  # remove extant genera (genera that range to bin 95)
  filter(lad_bin != 95) %>% 
  # subset to the same bins as we currently have for isotope data (isotemp) (>= 14)
  filter(fad_bin >= 14) %>% 
  as.data.frame()





# calculate trends --------------------------------------------------------

#### Calculate multiple long-term temperature trends

# read in preprocessed isotope data
isotemp2 <- read_csv(here("data/isotemp_trends.csv"))

# data for calculation
dat_trends <- sim_data %>% 
  # subset
  select(genus, fad = fad_age, lad = lad_age) 

# add bins with 0 
namevector <- as.character(c(1:94))
dat_trends[ , namevector] <- NA
dat_trends <- dat_trends[, c(4:length(dat_trends[0,]))]

# fill in with 0 for survival and 1 extinction
for (i in 1:length(dat_trends[,1])) {
  dat_trends[i,sim_data[i,"fad_bin"]] <- 0 
  dat_trends[i,sim_data[i,"lad_bin"]] <- 1
  ifelse(sim_data[i,"fad_bin"]!= sim_data[i,"lad_bin"]-1 & 
           sim_data[i,"fad_bin"]!= sim_data[i,"lad_bin"],
         dat_trends[i,(sim_data[i,"fad_bin"]+1):(sim_data[i,"lad_bin"]-1)] <- 0,
         NA)
}


# Reorganise ranges through time data so we can bind it to temperature data
dat_trends <- dat_trends %>% 
  as_tibble() %>%
  mutate(genus = sim_data$genus) %>% 
  group_by(genus) %>% 
  gather(-genus, key="bins", value="extinction", na.rm = T, factor_key = F) %>%
  arrange(desc(bins), .by_group = T) %>% 
  rename(stage = bins) %>% 
  mutate(stage = as.numeric(stage)) %>% 
  ungroup()



# Now bind the two
dat_trends <- full_join(dat_trends, isotemp2) %>% 
  # add higher taxonomy
  full_join(sim_data %>% select(order, family, genus)) %>% 
  # remove redundant colums from merging
  drop_na(stage) %>% 
  # order it properly
  select(order, family, genus, stage, age, extinction, temp, 
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

# save data for reproduction of simulations
save(dat_trends, file = here("data/simulation-results/simulated_trends.RData"))


# we want to see the rate of false positives (type I error) of our 
# analytical process and how this rate changes with size of data sets/
# number of observations. for this purpose, we  build data sets from our 
# simulated data with a varying number of observations. 
all_data_trends %>% 
  count(phylum)

# empirical data has a minimum number of observation at 4,632 (echinodermata) 
# and a maximum of 28,819. simulations will cover the same potential space, 
# starting at 1,000 and ending at 30,000 observations. for each observation size, 
# we generate 100 data sets to see the absolute range. 

# define number of observations per data sets
df_size <- c(1000, 2000, 3000, 4000, 6000, 8000, 10000, 20000, 30000)

# build empty list
rep_data <- list()

# create data sets using "sample_n" function, but before set progress bar
pb <- progress::progress_bar$new(total = 100) 

for (h in 1:100) {
  rep_data[[h]]  <-  sample_n(dat_trends, df_size[1], replace = T)
  rep_data[[100+h]]  <-  sample_n(dat_trends, df_size[2], replace = T)
  rep_data[[200+h]]  <-  sample_n(dat_trends, df_size[3], replace = T)
  rep_data[[300+h]]  <-  sample_n(dat_trends, df_size[4], replace = T)
  rep_data[[400+h]]  <-  sample_n(dat_trends, df_size[5], replace = T)
  rep_data[[500+h]]  <-  sample_n(dat_trends, df_size[6], replace = T)
  rep_data[[600+h]]  <-  sample_n(dat_trends, df_size[6], replace = T)
  rep_data[[700+h]]  <-  sample_n(dat_trends, df_size[6], replace = T)
  rep_data[[800+h]]  <-  sample_n(dat_trends, df_size[6], replace = T)
  
  # display progress
  pb$tick()
}

# save data for reproduction of simulations
save(rep_data, file = here("data/simulation-results/simulated_trends_rep.RData"))

# calculate effect size based on glmms ------------------------------------

#Build empty data frame for saving the results
sim_result <- tibble(nr_observations = rep(df_size, each = 100), 
                     effect_warm = numeric(length(df_size)*100), 
                     effect_cool = numeric(length(df_size)*100)) %>% 
  as.data.frame()
  
  
# data sets with less than 1000 observations have convergence issues and the loop 
# will stop at these, so we start with a minimum of 1000 observations
# we can suppress the message to see the progress bar
# set progress bar
pb <- progress::progress_bar$new(total = length(rep_data)) 

# loop through subsets
for (i in 1:length(rep_data)) {
  
  # assign data subset
  dat_subset <- rep_data[[i]]
  
  # get effect size for both warm and cool 
  effect_size <- suppressMessages(glmm_analysis(dat_subset = dat_subset))
  
  # save warm
  sim_result[i, 2] <- effect_size[1]
  # save cool 
  sim_result[i, 3] <- effect_size[2]
  
  # display progress
  pb$tick()
}


# save data
write_csv(sim_result, path = here("data/results/simulation_results.csv"))


# plot the results --------------------------------------------------------



# calculate mean and median
sim_result_av <- sim_result %>% 
  group_by(nr_observations) %>% 
  summarise_all(list(mean, median)) %>% 
  rename(warm_mean = effect_warm_fn1, cool_mean = effect_cool_fn1, 
         warm_median = effect_warm_fn2, cool_median = effect_cool_fn2)


### plot mean with 95 confidence interval
sim_result_mean <- sim_result %>% 
  # get 95% CI
  pivot_longer(-nr_observations,
               names_to = "pal_int", values_to = "effect_size") %>% 
  group_by(nr_observations, pal_int) %>% 
  # get confidence intervals from standard deviation 
  summarise(mean_effect = mean(effect_size),
            sd_effect = sd(effect_size),
            n_effect = n()) %>%
  mutate(se_effect = sd_effect / sqrt(n_effect),
         ci_low = mean_effect - qt(1 - (0.05 / 2), n_effect - 1) * se_effect,
         ci_high = mean_effect + qt(1 - (0.05 / 2), n_effect - 1) * se_effect) %>% 
  ungroup() %>% 
  select(- c(sd_effect, n_effect, se_effect)) %>% 
  # transform to wider format
  pivot_wider(id_cols = nr_observations, names_from = pal_int, 
              values_from = c(mean_effect, ci_low, ci_high)) %>% 
  # rename properly
  rename(cool = mean_effect_effect_cool, warm = mean_effect_effect_warm, 
         ci_low_cool = ci_low_effect_cool, ci_low_warm = ci_low_effect_warm, 
         ci_high_cool = ci_high_effect_cool, 
         ci_high_warm = ci_high_effect_warm) %>% 
  # plot it
  ggplot(aes(x = nr_observations)) +
  # add reference line
  geom_hline(yintercept = 0, colour = "orange", linetype = "dashed") +
  # confidence intervals
  geom_ribbon(aes(ymin = ci_low_warm, ymax = ci_high_warm), 
              fill =  "coral2", alpha = 0.3) +
  geom_ribbon(aes(ymin = ci_low_cool, ymax = ci_high_cool), 
              fill =  "#56B4E9", alpha = 0.3) +
  # means
  geom_point(aes(y = warm), colour = "grey30", 
             shape = 21,  fill = "coral2", size = 1.75, stroke = 1) +
  geom_point(aes(y = cool), colour = "grey30", 
             shape = 21,  fill = "#56B4E9", size = 1.75, stroke = 1) +
  # zoom out
  coord_cartesian(ylim = c(-0.02, 0.02)) +
  # add labs
  labs(y = "Simulated change in extinction risk [%]", 
       x = "Number of observations") +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.ticks.x = element_line())

ggsave(sim_result_mean, 
       filename = here("figures/extended-data/simulation_mean.png"), 
       width = 183, height = 89, units = "mm", dpi = 300)


# all data points, for warm
sim_range_warm <- ggplot(sim_result) +
  # highlight range
  geom_ribbon(aes(x = seq(from = 0, to = 32000, length.out = 900),
                  ymin = min(sim_result$effect_warm), 
                  ymax = max(sim_result$effect_warm)), 
              fill = "grey70", alpha = 0.5) +
  # add reference lines
  geom_hline(yintercept = c(min(sim_result$effect_warm), 
                            max(sim_result$effect_warm)), 
             colour = "orange") +
  geom_point(aes(x = nr_observations, y = effect_warm), colour = "grey30", 
             shape = 21,  fill = "coral2", size = 1.75, stroke = 1, 
             alpha = 0.3) +
  scale_y_continuous(limits = c(-0.06, 0.06)) +
  labs(y = "Simulated change in extinction risk [%]", 
       x = "Number of observations") +
  coord_cartesian(expand = FALSE) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.ticks.x = element_line())


# for cool
sim_range_cool <- ggplot(sim_result) +
  # highlight range
  geom_ribbon(aes(x = seq(from = 0, to = 32000, length.out = 900),
                  ymin = min(sim_result$effect_cool), 
                  ymax = max(sim_result$effect_cool)), 
              fill = "grey70", alpha = 0.5) +
  # add reference lines
  geom_hline(yintercept = c(min(sim_result$effect_cool), 
                            max(sim_result$effect_cool)), 
             colour = "orange") +
  geom_point(aes(x = nr_observations, y = effect_cool), colour = "grey30", 
             shape = 21,  fill = "#56B4E9", size = 1.75, stroke = 1, 
             alpha = 0.3) +
  scale_y_continuous(limits = c(-0.06, 0.06)) +
  labs(y = "Simulated change in extinction risk [%]", 
       x = "Number of observations") +
  coord_cartesian(expand = FALSE) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.ticks.x = element_line())

# load package patchwork to combine plots
library(patchwork)

sim_range_comb <- (sim_range_warm + labs(y = NULL, x = NULL)) / sim_range_cool +
  plot_annotation(tag_levels = c('a', 'b'))

# save plot
ggsave(sim_range_comb, filename = here("figures/raw/simulation_range.pdf"), 
       width = 183, height = 89, units = "mm")


