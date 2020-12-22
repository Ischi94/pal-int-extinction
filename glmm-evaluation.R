# loading packages
library(tidyverse)
library(lme4)
library(here)



# functions ---------------------------------------------------------------

# function that takes a model, summarises it and then returns a named vector with
# AIC and BIC values
get_criterion <- function(my_model){
  # make summary
  model_smr <- my_model %>% summary()
  
  # return AIC and BIC
  c(aic = model_smr$AICtab[[1]], bic = model_smr$AICtab[[2]])
}


# function taking a processed dfr named dat_trends and returning the AIC and BIC
# in a named vector for a glm model based on short-term change only
criterion_baseline <- function(dat_trends){
  # short term only
  my_model <- glmer("extinction~change_prev+(stage|genus)", 
                       data = dat_trends %>% select(starts_with("trend.")), 
                    family = "binomial")
  
  # get aic and bic values
  get_criterion(my_model = my_model)
}



# model comparison --------------------------------------------------------

# set up dfr to save results
model_comparison <- tibble(taxon = character(0), 
                           aic_int = numeric(0), 
                           bic_int = numeric(0),
                           aic_sonly = numeric(0), 
                           bic_sonly = numeric(0))

### arthropoda

# load interaction model list
load(file = here("data/model-output/arthropoda_model.Rds"))

# load trend data
dat_trends <- read_csv(here("data/arthropoda_trends.csv"))

# get AIC and BIC values for interaction model 
criterion_int <- get_criterion(my_model = arthropoda_model$final_trend)

# get AIC and BIC values for short-term change only (sonly) model
criterion_sonly <- criterion_baseline(dat_trends = dat_trends)

# assign results to dataframe
model_comparison <- model_comparison %>% 
  add_row(taxon = "Arthropoda", 
          aic_int = criterion_int[1], bic_int = criterion_int[2], 
          aic_sonly = criterion_sonly[1], bic_sonly = criterion_sonly[2])


### bivalvia

# load interaction model list
load(file = here("data/model-output/bivalvia_model.Rds"))

# load trend data
dat_trends <- read_csv(here("data/bivalvia_trends.csv"))

# get AIC and BIC values for interaction model 
criterion_int <- get_criterion(my_model = bivalvia_model$final_trend)

# get AIC and BIC values for short-term change only (sonly) model
criterion_sonly <- criterion_baseline(dat_trends = dat_trends)

# assign results to dataframe
model_comparison <- model_comparison %>% 
  add_row(taxon = "Bivalvia", 
          aic_int = criterion_int[1], bic_int = criterion_int[2], 
          aic_sonly = criterion_sonly[1], bic_sonly = criterion_sonly[2])



### cnidaria

# load interaction model list
load(file = here("data/model-output/cnidaria_model.Rds"))

# load trend data
dat_trends <- read_csv(here("data/cnidaria_trends.csv"))

# get AIC and BIC values for interaction model 
criterion_int <- get_criterion(my_model = cnidaria_model$final_trend)

# get AIC and BIC values for short-term change only (sonly) model
criterion_sonly <- criterion_baseline(dat_trends = dat_trends)

# assign results to dataframe
model_comparison <- model_comparison %>% 
  add_row(taxon = "Cnidaria", 
          aic_int = criterion_int[1], bic_int = criterion_int[2], 
          aic_sonly = criterion_sonly[1], bic_sonly = criterion_sonly[2])



### echinodermata

# load interaction model list
load(file = here("data/model-output/echinodermata_model.Rds"))

# load trend data
dat_trends <- read_csv(here("data/echinodermata_trends.csv"))

# get AIC and BIC values for interaction model 
criterion_int <- get_criterion(my_model = echinodermata_model$final_trend)

# get AIC and BIC values for short-term change only (sonly) model
criterion_sonly <- criterion_baseline(dat_trends = dat_trends)

# assign results to dataframe
model_comparison <- model_comparison %>% 
  add_row(taxon = "Echinodermata", 
          aic_int = criterion_int[1], bic_int = criterion_int[2], 
          aic_sonly = criterion_sonly[1], bic_sonly = criterion_sonly[2])



### foraminifera

# load interaction model list
load(file = here("data/model-output/foraminifera_model.Rds"))

# load trend data
dat_trends <- read_csv(here("data/foraminifera_trends.csv"))

# get AIC and BIC values for interaction model 
criterion_int <- get_criterion(my_model = foraminifera_model$final_trend)

# get AIC and BIC values for short-term change only (sonly) model
criterion_sonly <- criterion_baseline(dat_trends = dat_trends)

# assign results to dataframe
model_comparison <- model_comparison %>% 
  add_row(taxon = "Foraminifera", 
          aic_int = criterion_int[1], bic_int = criterion_int[2], 
          aic_sonly = criterion_sonly[1], bic_sonly = criterion_sonly[2])



### gastropoda

# load interaction model list
load(file = here("data/model-output/gastropoda_model.Rds"))

# load trend data
dat_trends <- read_csv(here("data/gastropoda_trends.csv"))

# get AIC and BIC values for interaction model 
criterion_int <- get_criterion(my_model = gastropoda_model$final_trend)

# get AIC and BIC values for short-term change only (sonly) model
criterion_sonly <- criterion_baseline(dat_trends = dat_trends)

# assign results to dataframe
model_comparison <- model_comparison %>% 
  add_row(taxon = "Gastropoda", 
          aic_int = criterion_int[1], bic_int = criterion_int[2], 
          aic_sonly = criterion_sonly[1], bic_sonly = criterion_sonly[2])



### mammalia

# load interaction model list
load(file = here("data/model-output/mammalia_model.Rds"))

# load trend data
dat_trends <- read_csv(here("data/mammalia_trends.csv"))

# get AIC and BIC values for interaction model 
criterion_int <- get_criterion(my_model = mammalia_model$final_trend)

# get AIC and BIC values for short-term change only (sonly) model
criterion_sonly <- criterion_baseline(dat_trends = dat_trends)

# assign results to dataframe
model_comparison <- model_comparison %>% 
  add_row(taxon = "Mammalia", 
          aic_int = criterion_int[1], bic_int = criterion_int[2], 
          aic_sonly = criterion_sonly[1], bic_sonly = criterion_sonly[2])



### reptilia

# load interaction model list
load(file = here("data/model-output/reptilia_model.Rds"))

# load trend data
dat_trends <- read_csv(here("data/reptilia_trends.csv"))

# get AIC and BIC values for interaction model 
criterion_int <- get_criterion(my_model = reptilia_model$final_trend)

# get AIC and BIC values for short-term change only (sonly) model
criterion_sonly <- criterion_baseline(dat_trends = dat_trends)

# assign results to dataframe
model_comparison <- model_comparison %>% 
  add_row(taxon = "Reptilia", 
          aic_int = criterion_int[1], bic_int = criterion_int[2], 
          aic_sonly = criterion_sonly[1], bic_sonly = criterion_sonly[2])

# save data
write_csv(model_comparison, path = here("data/results/model_comparison.csv"))


model_comparison_aic <- model_comparison %>% 
  # aic percentage change
  mutate(perc_change = ((aic_sonly/ aic_int)-1)*100) %>% 
  ggplot(aes(y = fct_reorder(taxon, desc(taxon)))) +
  # add lines
  geom_segment(aes(x = 0, xend = perc_change, yend = taxon), 
               colour = "grey", size = 1) +
  geom_vline(xintercept = 0, colour = "grey25") +
  # add points 
  geom_point(aes(x = perc_change, fill = perc_change >= 0), 
             colour = "grey25", shape = 21,  size = 2.5, stroke = 0.6,
             show.legend = FALSE, alpha = 0.85) +
  # set colours
  scale_fill_manual(values = c("coral2", "#56B4E9")) +
  # labels
  labs(x = "AIC change [%]                                  ", y = NULL) +
  # limits
  xlim(c(-2.75, 5)) +
  expand_limits(y = c(0, 9)) +
  # add arrow
  geom_segment(x = - 0.05, xend = - 0.5, y = 8.75, yend = 8.75,
               lineend = "round", linejoin = "round", size = 0.5, 
               arrow = arrow(length = unit(0.075, "inches")),
               colour = "coral2", alpha = 0.3) +
  # and second arrow
  geom_segment(x = 0.05, xend = 0.5, y = 8.75, yend = 8.75,
               lineend = "round", linejoin = "round", size = 0.5, 
               arrow = arrow(length = unit(0.075, "inches")),
               colour = "#56B4E9", alpha = 0.3) +
  # add text for decrease
  annotate(geom = "text", label = "deterioration", 
           x = - 1.3, y = 8.81, size = 3.5, colour = "grey40") +
  # and for increase
  annotate(geom = "text", label = "improvement", 
           x = 1.35, y = 8.81, size = 3.5, colour = "grey40") +
  # set theme and modify
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5, linetype = "dashed"), 
        axis.title.x = element_text(hjust = 1))

# save plot
ggsave(model_comparison_aic, 
       filename = here("figures/model_comparison_aic.pdf"), 
       width = 88, height = 89, units = "mm", dpi = 400)





# effect size -------------------------------------------------------------

# read in the dataframe with results from the Mann-Whitney tests. results where 
# extracted from the saved lists. 
effect_intensity <- read_csv(here("data/results/effect_intensity.csv"))

# read in the data from simulations to get range of null distribution
sim_results <- read_csv(here("data/results/simulation_results.csv"))

# calculate max and min of simulated change in extinction risk
sim_range <- sim_results %>% 
  pivot_longer(-nr_observations) %>% 
  summarise(min_val = min(value), 
            max_val = max(value)) %>% 
  as.double()



# plot it
effect_intensity_plot <- 
  ggplot(effect_intensity, aes(y = fct_reorder(taxon, desc(taxon)))) +
  # add range of simulations/ null distribution
  geom_rect(xmin = sim_range[1], xmax = sim_range[2], ymin = 0, ymax = 8.5, 
            alpha = 0.05, fill = "grey70") +
  # add reference line
  geom_vline(xintercept = 0, colour = "orange") +
  # add 95% CI for warm
  geom_errorbarh(aes(xmin = warm_low, xmax = warm_high), colour = "coral2", 
                 height = 0.2, size = 0.75, alpha = 0.8) +
  # add 95% CI for cool
  geom_errorbarh(aes(xmin = cool_low, xmax = cool_high), colour = "#56B4E9", 
                 height = 0.2, size = 0.75, alpha = 0.8) +
  # add dots for warm
  geom_point(aes(x = warm), shape = 21, 
             colour = "grey30", fill = "coral2", size = 1.75, stroke = 1) +
  # add dots for cool
  geom_point(aes(x = cool), shape = 21, 
             colour = "grey30", fill = "#56B4E9", size = 1.75, stroke = 1) +
  # add arrow
  geom_segment(x = - 0.005, xend = - 0.04, y = 8.75, yend = 8.75,
    lineend = "round", linejoin = "round", size = 0.5, 
    arrow = arrow(length = unit(0.075, "inches")),
    colour = "grey40") +
  # and second arrow
  geom_segment(x = 0.005, xend = 0.04, y = 8.75, yend = 8.75,
               lineend = "round", linejoin = "round", size = 0.5, 
               arrow = arrow(length = unit(0.075, "inches")),
               colour = "grey40") +
  # add text for decrease
  annotate(geom = "text", label = "decreased risk", 
           x = - 0.13, y = 8.85, size = 3.5, colour = "grey20") +
  # and for increase
  annotate(geom = "text", label = "increased risk", 
           x = 0.13, y = 8.85, size = 3.5, colour = "grey20") +
  scale_x_continuous(limits = c(-0.2, 0.5)) +
  scale_y_discrete(expand = expansion(mult = c(0.06, 0.14))) +
  labs(y = NULL, x = "Change in extinction risk [%]") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(linetype = "dashed"), 
        panel.grid.major.x = element_line(colour = "grey97"),
        axis.line.y = element_blank(), 
        axis.ticks.length.y = unit(0, "mm"))

  
ggsave(effect_intensity_plot, filename = here("figures/raw/effect_intensity.pdf"), 
       width = 183, height = 89, units = "mm")



# duration regression -----------------------------------------------------

# read in all data trends combined
all_data_trends <- read_csv(here("data/all_data_trends.csv"), guess_max = 1e5)

# calculate durations
durations <- all_data_trends %>% 
  group_by(genus) %>% 
  # calculate durations
  mutate(min_age = min(age), max_age = max(age), 
         duration = max_age - min_age) %>% 
  ungroup() %>% 
  select(phylum:age, duration) %>% 
  # remove duplicates
  distinct(genus, .keep_all = TRUE) %>% 
  # make dummy column for grouping
  mutate(my_group = if_else(!is.na(phylum), phylum, class), 
         my_group2 = if_else(my_group == "Mollusca", class, my_group),
         my_group3 = if_else(my_group2 == "Chordata", "Reptilia", my_group2)) %>% 
  select(phylum:duration, taxon = my_group3) %>% 
  group_by(taxon) %>% 
  summarise(median_duration = median(duration), 
            mean_duration = mean(duration))

# merge with effect estimates
regression_data <- durations %>% 
  full_join(effect_intensity) %>% 
  select(taxon:warm, cool)

# pivot to longer format
# mean
regression_data_mean <- regression_data %>% 
  select(-median_duration) %>% 
  pivot_longer(cols = c(warm, cool), names_to = "pal_int", values_to = "effect")

# median
regression_data_median <- regression_data %>% 
  select(-mean_duration) %>% 
  pivot_longer(cols = c(warm, cool), names_to = "pal_int", values_to = "effect")


# plot it
# mean
# get r squared
mean_r2 <- lm(effect ~ mean_duration, data = regression_data_mean) %>% 
  summary() %>% 
  .$adj.r.squared %>% 
  round(digits = 2)

regression_plot_mean <- 
  ggplot(regression_data_mean, aes(x = mean_duration, y = effect)) +
  # add regression line
  geom_smooth(method = "lm", colour = "orange", fill = "grey70", 
              size = 1.5) +
  # add points
  geom_point(aes(fill = pal_int), shape = 21, show.legend = F,
             colour = "grey30", size = 1.75, stroke = 1) +
  scale_fill_manual(values = c("#56B4E9", "coral2")) +
  annotate("text", label = paste("R^2 == ", mean_r2),
           x = 15, y = 0.43, parse = TRUE) +
  labs(x = "Mean duration of genera [myr]", y = "Change in extinction risk [%]") +
  scale_y_continuous(limits = c(-0.1, 0.43)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())


# median
# get r squared
median_r2 <- lm(effect ~ median_duration, data = regression_data_median) %>% 
  summary() %>% 
  .$adj.r.squared %>% 
  round(digits = 2)

regression_plot_median <- 
  ggplot(regression_data_median, aes(x = median_duration, y = effect)) +
  # add regression line
  geom_smooth(method = "lm", colour = "orange", fill = "grey70", 
              size = 1.5) +
  # add points
  geom_point(aes(fill = pal_int), shape = 21, show.legend = F, 
             colour = "grey30", size = 1.75, stroke = 1) +
  scale_fill_manual(values = c("#56B4E9", "coral2")) +
  annotate("text", label = paste("R^2 == ", median_r2),
           x = 12, y = 0.43, parse = TRUE) +
  labs(x = "Median duration of genera [myr]", y = "Change in extinction risk [%]") +
  scale_y_continuous(limits = c(-0.1, 0.43)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())


# save plots
ggsave(regression_plot_mean, 
       filename = here("figures/raw/regression_plot_mean.pdf"), 
       width = 183, height = 89, units = "mm")

ggsave(regression_plot_median, 
       filename = here("figures/raw/regression_plot_median.pdf"), 
       width = 183, height = 89, units = "mm")

