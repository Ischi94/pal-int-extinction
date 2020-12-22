# loading packages
library(tidyverse)
library(lme4)
library(geiger)
library(here)


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

###


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

###

# Calculate GLMM's ----------------------------------------------------------------

### arthropoda 

# load data
dat_trends <- read_csv(here("data/arthropoda_trends.csv"))

# get output list based on glmm
arthropoda_model <- glmm_analysis(dat_trends = dat_trends)

# save data
save(arthropoda_model, file = here("data/model-output/arthropoda_model.Rds"))



### bivalvia 

# load data
dat_trends <- read_csv(here("data/bivalvia_trends.csv"))

# get output list based on glmm
bivalvia_model <- glmm_analysis(dat_trends = dat_trends)

# save data
save(bivalvia_model, file = here("data/model-output/bivalvia_model.Rds"))



### cnidaria 

# load data
dat_trends <- read_csv(here("data/cnidaria_trends.csv"))

# get output list based on glmm
cnidaria_model <- glmm_analysis(dat_trends = dat_trends)

# save data
save(cnidaria_model, file = here("data/model-output/cnidaria_model.Rds"))



### echinodermata 

# load data
dat_trends <- read_csv(here("data/echinodermata_trends.csv"))

# get output list based on glmm
echinodermata_model <- glmm_analysis(dat_trends = dat_trends)

# save data
save(echinodermata_model, file = here("data/model-output/echinodermata_model.Rds"))



### foraminifera 

# load data
dat_trends <- read_csv(here("data/foraminifera_trends.csv"))

# get output list based on glmm
foraminifera_model <- glmm_analysis(dat_trends = dat_trends)

# save data
save(foraminifera_model, file = here("data/model-output/foraminifera_model.Rds"))



### gastropoda 

# load data
dat_trends <- read_csv(here("data/gastropoda_trends.csv"))

# get output list based on glmm
gastropoda_model <- glmm_analysis(dat_trends = dat_trends)

# save data
save(gastropoda_model, file = here("data/model-output/gastropoda_model.Rds"))



### mammalia 

# load data
dat_trends <- read_csv(here("data/mammalia_trends.csv"))

# get output list based on glmm
mammalia_model <- glmm_analysis(dat_trends = dat_trends)

# save data
save(mammalia_model, file = here("data/model-output/mammalia_model.Rds"))



### reptilia 

# load data
dat_trends <- read_csv(here("data/reptilia_trends.csv"))

# get output list based on glmm
reptilia_model <- glmm_analysis(dat_trends = dat_trends)

# save data
save(reptilia_model, file = here("data/model-output/reptilia_model.Rds"))

