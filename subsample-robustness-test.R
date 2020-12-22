library(tidyverse)
library(divDyn)
library(lme4)
library(geiger)
library(here)


# read in data ------------------------------------------------------------

# import data 
# stages
gradstein <- read_csv(here("data/gradstein.csv"))

# read in preprocessed isotope data
isotemp2 <- read_csv(here("data/isotemp_trends.csv"))



# functions ---------------------------------------------------------------


### binning 

# function for occurrences into stage bins
binning <- function(x){gradstein$stg[x<gradstein$bottom  &  x>=gradstein$top]}


### subsample_pbdb

# function that returns  shareholder quorum and classical rarefaction subsampled 
# data in list format
subsample_pbdb <- function(pbdb_raw){
  # bin the data to stages using the binning (function)
  pbdb_raw <- pbdb_raw %>% 
    # add mean age
    mutate(mean_ma = (max_ma + min_ma)/2) %>% 
    # add stage number
    mutate(stg = map_dbl(mean_ma, binning)) %>% 
    # Exclude unassigned bins
    drop_na(stg) 
  
  # shareholder quorum
  pbdb_sqs <- subsample(pbdb_raw,iter = 1, q = 0.4, tax = "genus", bin = "stg",
                        output = "dist", type = "sqs", FUN = fadlad)
  
  pbdb_sqs <- pbdb_sqs$results %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "genus")
  
  # classical rarefaction
  pbdb_cr <- subsample(pbdb_raw,iter = 1, q = 50, tax = "genus", bin = "stg",
                       output = "dist", type = "cr", FUN = fadlad)
  
  pbdb_cr <- pbdb_cr$results %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "genus")
  
  # merge together with raw data (but keep only those from the subsets) 
  # to get taxonomy
  
  # shareholder
  pbdb_sqs <- pbdb_sqs %>% 
    # merge
    left_join(pbdb_raw, by = "genus") %>% 
    # remove potential duplicates from merging
    distinct()
  
  # rarefaction
  pbdb_cr <- pbdb_cr %>% 
    # merge
    left_join(pbdb_raw, by = "genus") %>% 
    # remove potential duplicates from merging
    distinct()
  
  # make list
  output_ls <- list(sqs = pbdb_sqs, cr = pbdb_cr)
  
  output_ls
}

### make trends 

# function that takes a subsampled dfr as argument and returns a processes file, 
# ready for glmm analysis
make_trends <- function(dfr) {
  dat_occ <- dfr %>% 
    # remove the NA's
    filter(order != "NO_ORDER_SPECIFIED", 
           family != "NO_FAMILY_SPECIFIED")
  
  # save higher taxonomy
  tax <- dat_occ %>% 
    select(phylum, order, family, genus) %>% 
    unique()
  
  # transform to range data
  dat_range <- fadlad(dat_occ, tax="genus", age= c("max_ma", "min_ma"))
  
  
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
    select(phylum, order, family, genus, fad = FAD, lad = LAD)  
  
  
  ### bin data 
  
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
    filter(fad_bin >= 14)
  
  
  #### Calculate multiple long-term temperature trends
  
  # data for calculation
  dat_trends <- dat_range %>% 
    # subset
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
    full_join(dat_range %>% select(phylum, order, family, genus)) %>% 
    # remove redundant colums from merging
    drop_na(stage) %>% 
    # order it properly
    select(phylum, order, family, genus, stage, age, extinction, temp, 
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
  
  dat_trends
}


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


### glmm_analysis

# calculate models and make output list
glmm_analysis <- function(dat_trends){
  
  # Iterate through each warming
  vars <- dat_trends %>% select(starts_with("trend.")) %>% names()
  interaction_model <- lapply(setNames(vars, vars), function(var) {
    form = paste("extinction~change_prev:", var, "+(stage|genus)")
    glmer(form, data=dat_trends, family="binomial")
  })
  
  
  # Make data frame for model output
  interaction_df <- model_df(interaction_model)
  
  # choose final model
  interaction_final <- interaction_model[[which(interaction_df$dAIC==0)]]
  
  
  # make informed predictions based on a subset (palaeoclimate interaction)
  # type = response gives us probability instead of log Odds
  
  # first get the trends
  trend <- interaction_df$model[[which(interaction_df$dAIC==0)]] 
  
  
  #  warming warming
  ww_raw <- dat_trends %>% filter(get(trend) >=0 & change_prev >= 0)
  ww_pred <- predict(interaction_final, newdata = ww_raw,
                     type = "response", allow.new.levels = TRUE)
  
  
  #  warming cooling 
  wc_raw <- dat_trends %>% filter(get(trend) >=0 & change_prev <= 0)
  wc_pred <- predict(interaction_final, newdata = wc_raw,
                     type = "response", allow.new.levels = TRUE)
  
  
  #  cooling cooling 
  cc_raw <- dat_trends %>% filter(get(trend) <= 0 & change_prev <= 0)
  cc_pred <- predict(interaction_final, newdata = cc_raw,
                     type = "response", allow.new.levels = TRUE)
  
  
  #  cooling warming
  cw_raw <- dat_trends %>% filter(get(trend) <=0 & change_prev >= 0)
  cw_pred <- predict(interaction_final, newdata = cw_raw,
                     type = "response", allow.new.levels = TRUE)
  
  
  
  
  
  # take the predictions and test whether they significantly are above or below the
  # baseline. This is a one sided Wilcoxon rank sum test 
  # (equivalent to the Mann-Whitney test)
  wilcox_warm <- wilcox.test(ww_pred, wc_pred, paired = F,  conf.int = T)
  wilcox_cool <- wilcox.test(cc_pred, cw_pred, paired = F,  conf.int = T)
  
  # make output list
  output_df <- list(model = interaction_model, 
                    final_trend = interaction_final, 
                    wilcox_warm = wilcox_warm, 
                    wilcox_cool = wilcox_cool)
  output_df
}

### get_criterion

# function that takes a model, summarises it and then returns a named vector with
# AIC and BIC values
get_criterion <- function(my_model){
  # make summary
  model_smr <- my_model %>% summary()
  
  # return AIC and BIC
  c(aic = model_smr$AICtab[[1]], bic = model_smr$AICtab[[2]])
}


### criterion_baseline

# function taking a processed dfr named dat_trends and returning the AIC and BIC
# in a named vector for a glm model based on short-term change only
criterion_baseline <- function(dat_trends){
  # short term only
  my_model <- glmer("extinction~change_prev+(stage|genus)", 
                    data = dat_trends, family = "binomial")
  
  # get aic and bic values
  get_criterion(my_model = my_model)
}


# process data ------------------------------------------------------------

# read raw data for both bivalvia and reptilia
bivalvia <- read_csv(here("data/raw-fossil-data/bivalvia.csv"), guess_max = 1e5)

reptilia <- read_csv(here("data/raw-fossil-data/reptilia.csv"), guess_max = 1e5)

### subsampling

# bivalvia
bivalvia_sub <- subsample_pbdb(pbdb_raw = bivalvia)

# reptilia
reptilia_sub <- subsample_pbdb(pbdb_raw = reptilia)

### make trends from subsamples

## bivalvia

# sqs
bivalvia_sqs_trends <- make_trends(dfr = bivalvia_sub$sqs)

# cr
bivalvia_cr_trends <- make_trends(dfr = bivalvia_sub$cr)

## reptilia

# sqs
reptilia_sqs_trends <- make_trends(dfr = reptilia_sub$sqs)

# cr
reptilia_cr_trends <- make_trends(dfr = reptilia_sub$cr)


### glmms

# set up dfr to save results
model_comparison_sub <- tibble(taxon = character(0), 
                               type = character(0), 
                               aic_int = numeric(0), 
                               bic_int = numeric(0),
                               aic_sonly = numeric(0), 
                               bic_sonly = numeric(0))

## bivalvia

# sqs
dat_trends <- bivalvia_sqs_trends
bivalvia_sqs_res <- glmm_analysis(dat_trends = dat_trends)

# get AIC and BIC values for interaction model 
criterion_int <- get_criterion(my_model = bivalvia_sqs_res$final_trend)

# get AIC and BIC values for short-term change only (sonly) model
criterion_sonly <- criterion_baseline(dat_trends = dat_trends)

# assign results to dataframe
model_comparison_sub <- model_comparison_sub %>% 
  add_row(taxon = "Bivalvia", type = "sqs", 
          aic_int = criterion_int[1], bic_int = criterion_int[2], 
          aic_sonly = criterion_sonly[1], bic_sonly = criterion_sonly[2])

# cr
dat_trends <- bivalvia_cr_trends
bivalvia_cr_res <- glmm_analysis(dat_trends = dat_trends)

# get AIC and BIC values for interaction model 
criterion_int <- get_criterion(my_model = bivalvia_cr_res$final_trend)

# get AIC and BIC values for short-term change only (sonly) model
criterion_sonly <- criterion_baseline(dat_trends = dat_trends)

# assign results to dataframe
model_comparison_sub <- model_comparison_sub %>% 
  add_row(taxon = "Bivalvia", type = "cr", 
          aic_int = criterion_int[1], bic_int = criterion_int[2], 
          aic_sonly = criterion_sonly[1], bic_sonly = criterion_sonly[2])


## reptilia

# sqs
dat_trends <- reptilia_sqs_trends
reptilia_sqs_res <- glmm_analysis(dat_trends = dat_trends)

# get AIC and BIC values for interaction model 
criterion_int <- get_criterion(my_model = reptilia_sqs_res$final_trend)

# get AIC and BIC values for short-term change only (sonly) model
criterion_sonly <- criterion_baseline(dat_trends = dat_trends)

# assign results to dataframe
model_comparison_sub <- model_comparison_sub %>% 
  add_row(taxon = "Reptilia", type = "sqs", 
          aic_int = criterion_int[1], bic_int = criterion_int[2], 
          aic_sonly = criterion_sonly[1], bic_sonly = criterion_sonly[2])

# cr
dat_trends <- reptilia_cr_trends
reptilia_cr_res <- glmm_analysis(dat_trends = dat_trends)

# get AIC and BIC values for interaction model 
criterion_int <- get_criterion(my_model = reptilia_cr_res$final_trend)

# get AIC and BIC values for short-term change only (sonly) model
criterion_sonly <- criterion_baseline(dat_trends = dat_trends)

# assign results to dataframe
model_comparison_sub <- model_comparison_sub %>% 
  add_row(taxon = "Reptilia", type = "cr", 
          aic_int = criterion_int[1], bic_int = criterion_int[2], 
          aic_sonly = criterion_sonly[1], bic_sonly = criterion_sonly[2])


# save data
write_csv(model_comparison_sub, 
          path = here("data/results/model_comparison_subsampling.csv"))

# evaluation and plotting -------------------------------------------------

# pivot it to wider format
model_comparison_long <- model_comparison_sub %>% 
  pivot_longer(cols = -c(taxon, type), names_to = "crit_type", values_to = "crit_value")


# aic only
model_comparison_sub_aic <- model_comparison_long %>% 
  # select aic values
  filter(str_detect(crit_type, "aic")) %>% 
  # add x dummy variable 
  add_column(x_value = rep(c(100, 0), 4)) %>% 
  # add column for offsetting of labels
  mutate(crit_label = if_else(crit_type == "aic_int", 
                              crit_value - 500, crit_value + 500)) %>% 
  # reorder factor
  mutate(crit_type = fct_relevel(crit_type, "aic_sonly")) %>% 
  # start plotting
  ggplot(aes(x = x_value, y = crit_value)) +
  # facetting
  facet_wrap(vars(taxon, type)) +
  # add lines
  geom_line(size = 0.55, colour = "grey25") + 
  # # add points 
  geom_point(aes(fill = crit_type), shape = 21, 
             colour = "grey10", size = 1.75, stroke = 0.6) +
  # add aic numbers/ labelling
  geom_text(aes(label = round(crit_value), y = crit_label),
            size = 2.5, colour = "grey10") +
  # increase column width (white space)
  scale_x_continuous(limits = c(-20, 120)) +
  # relabel
  labs(y = "AIC value", x = NULL) +
  # relabel legend and change colour of dots
  scale_fill_manual(name = NULL, labels = c("change only", "change & trend"), 
                    values = c("grey40", "grey85")) +
  # set theme and modify
  theme_minimal() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2),
        panel.grid.minor.y = element_blank(),
        legend.position = c(0.5, 0.5), 
        panel.background = element_rect(fill = NA, color = "grey98"), 
        legend.background = element_rect(fill = "white", colour = "white"))

# save plot
ggsave(model_comparison_sub_aic, 
       filename = here("figures/raw/model_comparison_subsampling_aic.pdf"), 
       width = 183, height = 89, units = "mm")
