library(tidyverse)
library(lme4)
library(here)
library(geiger)
library(MuMIn)

## read in preprocessed data for all taxa combined
all_trends <- read_csv(here("data/all_data_trends.csv"), guess_max = 1e5)




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
  interaction_final <- interaction_model[[which(interaction_df$dAIC==0)[[1]]]]
  
  # return final model
  interaction_final
}


## calculate model with stages covering mass extinctions
dat_trends <- all_trends
with_ext <- glmm_analysis(dat_trends = dat_trends)
  

# save the model output for mass extinctions
save(with_ext, file = here("data/model-output/all_trends_model.Rds"))

## now the same model without stages covering mass extinctions
# omit stages covering mass extinctions
dat_trends <- all_trends %>% 
  filter(stage != 20 & stage != 35 & stage != 51 & stage != 58 & stage != 81)

# calculate model
without_ext <- glmm_analysis(dat_trends = dat_trends)

# save the model output without mass extinctions
save(without_ext, 
     file = here("data/model-output/all_trends_without_massextinction_model.Rds"))

## calculate pseudo r-squared
with_ext_r2 <- r.squaredGLMM(with_ext) 
without_ext_r2 <- r.squaredGLMM(without_ext)


# build data frame with conditional (both fixed and random effects) theoretical R2 
r2_df <- tibble(type = c("mass", "no_mass"), 
                pseudo_rsquared = c(with_ext_r2[1,2], without_ext_r2[1,2]))

# save it
write_csv(r2_df, path = here("data/results/r2_mass_extinctions.csv"))

# calculate percentage increase
perc_incr <-  round((0.11674929 - 0.11070678) / 0.11070678 * 100, digits = 1)
  
  
# plot it
mass_extinctions_plot <- ggplot(r2_df) +
  geom_segment(aes(x = type, xend = type,
                   y = 0, yend = pseudo_rsquared), 
               color="grey30", size = 0.8) +
  geom_point(aes(x = type, y = pseudo_rsquared, fill = type), shape = 21,
             colour = "grey20", size = 4, stroke = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("grey40", "grey80")) +
  scale_y_continuous(limits = c(0, 0.151), expand = c(0,0)) +
  scale_x_discrete(labels = c("with mass extinctions", 
                              "without mass extinctions")) +
  annotate(geom = "text", label = paste("+", perc_incr, "%"),
           x = 2.25, y = 0.117, size = 2.5) +
  labs(x = NULL, y = bquote("Pseudo" ~R^2)) +
  theme_classic(base_size = 7) +
  theme(panel.grid.major.y = element_line(), 
        axis.line.y = element_blank(), 
        axis.ticks.length = unit(0, units = "mm"))


# save plot
ggsave(mass_extinctions_plot, 
       filename = here("figures/extended-data/mass_extinctions.jpg"), 
       width = 89*1.25, height = 89, units = "mm", dpi = 400)

