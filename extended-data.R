library(tidyverse)
library(here) 
library(flextable)
library(officer)
library(deeptime) 
library(divDyn) 
library(patchwork)


# tables ------------------------------------------------------------------



# change in extinction risk  ----------------------------------------------


# load data
ext_risk <- read_csv(here("data/results/effect_intensity.csv"))

# preprocess it
ext_risk_fxt <- ext_risk %>%
  select(Taxon = taxon, warm_low, warm, warm_high, 
         cool_low, cool, cool_high) %>% 
  # round properly
  mutate(across(where(is.numeric), round, 3)) %>% 
  # make a flextable
  flextable() %>% 
  # highlight warm vs cool
  bg(j = c("warm_low", "warm", "warm_high"), 
     bg = "#FFD2A2", part = "body") %>% 
  bg(j = c("cool_low", "cool", "cool_high"), 
     bg = "#B3E3FF", part = "body") %>% 
  # set column names
  compose(j = "warm_low", part = "header", 
          value = as_paragraph("Lower CI")) %>% 
  compose(j = "warm", part = "header", 
          value = as_paragraph("Estimate")) %>% 
  compose(j = "warm_high", part = "header", 
          value = as_paragraph("Upper CI")) %>% 
  compose(j = "cool_low", part = "header", 
          value = as_paragraph("Lower CI")) %>% 
  compose(j = "cool", part = "header", 
          value = as_paragraph("Estimate")) %>% 
  compose(j = "cool_high", part = "header", 
          value = as_paragraph("Upper CI")) %>%
  # add header
  add_header_row(values = c("", "", "Warming", "",  "", "Cooling", ""), 
                 top = TRUE) %>% 
  # bold headers
  bold(part = "header") %>% 
  # italic taxon column
  italic(j = "Taxon") %>%
  autofit()


# open docx-file and add flextable
my_doc <- read_docx() %>% 
  body_add_flextable(ext_risk_fxt)



# simulation range --------------------------------------------------------

# load data
sim_result <- read_csv(here("data/results/simulation_results.csv"),
                       guess_max = 1e5)

# make flextable

sim_result_fxt <- sim_result %>% 
  # calculate mean and median
  group_by(nr_observations) %>% 
  summarise_all(list(mean, median)) %>% 
  rename(warm_mean = effect_warm_fn1, cool_mean = effect_cool_fn1, 
         warm_median = effect_warm_fn2, cool_median = effect_cool_fn2) %>% 
  # round properly
  mutate(across(c(warm_mean:cool_median), round, 3)) %>% 
  # make a flextable
  flextable() %>% 
  # highlight warm vs cool
  bg(j = c("warm_mean", "warm_median"), 
     bg = "#FFD2A2", part = "body") %>% 
  bg(j = c("cool_mean", "cool_median"), 
     bg = "#B3E3FF", part = "body") %>% 
  # set column names
  compose(j = "nr_observations", part = "header", 
          value = as_paragraph("# Observations")) %>% 
  compose(j = "warm_mean", part = "header", 
          value = as_paragraph("Mean")) %>% 
  compose(j = "cool_mean", part = "header", 
          value = as_paragraph("Mean")) %>% 
  compose(j = "warm_median", part = "header", 
          value = as_paragraph("Median")) %>% 
  compose(j = "cool_median", part = "header", 
          value = as_paragraph("Median")) %>% 
  # merge header
  merge_h(part = "header") %>%  
  # bold headers
  bold(part = "header") %>% 
  # add subheader
  add_header_row(values = c("", "warming", "cooling", "warming",  "cooling"), 
                 top = FALSE) %>% 
  # add border
  hline(i = 2, j = 2:5, border = fp_border(color = "black"), part = "header") %>% 
  autofit()


# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(sim_result_fxt, pos = "after")


# model comparison --------------------------------------------------------


# load data 
model_comparison <- read_csv(here("data/results/model_comparison.csv"))


# make flextable
model_comparison_fxt <- model_comparison %>% 
  # round properly
  mutate(across(where(is.numeric), round, 0)) %>% 
  # reorder
  select(taxon, aic_sonly, aic_int, bic_sonly, bic_int) %>% 
  # make a flextable
  flextable() %>% 
  # highlight warm vs cool
  bg(j = c("aic_int", "bic_int"), 
     bg = "grey60", part = "body") %>% 
  bg(j = c("aic_sonly", "bic_sonly"), 
     bg = "grey85", part = "body") %>% 
  # set column names
  compose(j = "taxon", part = "header", 
          value = as_paragraph("Taxon")) %>% 
  compose(j = "aic_int", part = "header", 
          value = as_paragraph("AIC")) %>% 
  compose(j = "bic_int", part = "header", 
          value = as_paragraph("BIC")) %>% 
  compose(j = "aic_sonly", part = "header", 
          value = as_paragraph("AIC")) %>% 
  compose(j = "bic_sonly", part = "header", 
          value = as_paragraph("BIC")) %>% 
  # bold headers
  bold(part = "header") %>% 
  # set second header
  add_header_row(values = c("", "Trend only", "Change & Trend", 
                            "Trend only", "Change & Trend"), top = FALSE) %>%
  # italic taxon column
  italic(j = "taxon") %>%
  # merge headers
  merge_h(part = "header") %>% 
  # align headers
  align(align = "center", part = "header") %>% 
  # add border
  hline(i = 2, j = 2:5, border = fp_border(color = "black"), part = "header") %>% 
  autofit()


# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(model_comparison_fxt, pos = "after")





# autocorrelation Durbin-Watson -------------------------------------------


# load data
autocorrelation <- read_csv(here("data/results/autocorrelation_results.csv"))

# make flextable
autocorrelation_fxt <- autocorrelation %>% 
  # add column for linerange
  add_column(estimate = .$durbin_watson) %>% 
  # rename
  select(Taxon = taxon, "Durbin-Watson statistic" = durbin_watson, estimate) %>% 
  # round 
  mutate(across(where(is.numeric), round, 2)) %>% 
  flextable() %>% 
  # add linerange
  compose(j = 3, value = as_paragraph(
    linerange(value = estimate, max = max(estimate), 
              stickcol = "grey40")),
    part = "body") %>% 
  # set column names
  compose(j = "estimate", part = "header", 
          value = as_paragraph("Durbin-Watson statistic")) %>% 
  # bold headers
  bold(part = "header") %>% 
  # merge header cells
  merge_h(part = "header") %>% 
  # italic taxon column
  italic(j = "Taxon") %>%
  autofit()
  
# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(autocorrelation_fxt, pos = "after")


# number of genera --------------------------------------------------------


# load cleaned data
all_data_trends <- read_csv(here("data/all_data_trends.csv"), guess_max = 1e5)


# calculate number of genera within fossil groups
nr_genera <- all_data_trends %>% 
  # make dummy column for grouping
  mutate(my_group = if_else(!is.na(phylum), phylum, class), 
         my_group2 = if_else(my_group == "Mollusca", class, my_group),
         my_group3 = if_else(my_group2 == "Chordata", "Reptilia", my_group2)) %>% 
  select(phylum:genus, taxon = my_group3) %>% 
  # remove duplicates
  distinct(genus, .keep_all = TRUE) %>% 
  # count number of genera within higher groups
  count(taxon)


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

# merge durations and nr_genera
data_info_fxt <- nr_genera %>% 
  full_join(durations) %>% 
  # rename
  rename("Taxon" = taxon, "# Genera" = n, 
         "Median duration [myr]" = median_duration,
         "Mean duration [myr]" = mean_duration) %>% 
  mutate(across(where(is.numeric), round, 2)) %>% 
  # make a flextable
  flextable() %>% 
  # bold headers
  bold(part = "header") %>% 
  # italic taxon column
  italic(j = "Taxon") %>% 
  autofit()


# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(data_info_fxt, pos = "after")






# model summaries ---------------------------------------------------------


## functions

# make own function to tidy glmm output in broom style
tidy_glmm <- function(model_smr, Taxon){
  # get coefficients
  model_smr$coefficients %>% 
    # convert to tibble
    as.data.frame() %>% 
    rownames_to_column(var = "Interaction") %>% 
    as_tibble() %>% 
    # remove intercept term
    filter(Interaction == Interaction[c(FALSE, TRUE)]) %>% 
    # add information criterion
    add_column(AIC = model_smr$AICtab[1], 
               BIC = model_smr$AICtab[2]) %>% 
    # add taxon identifier
    add_column(Taxon = Taxon, .before = "Interaction") %>% 
    # remove z value
    select(-"z value")
}

### load final models and make nice summary

## arthropoda
load(file = here("data/model-output/arthropoda_model.Rds"))

# get traditional model summary
arthropoda_smr <- arthropoda_model$final_trend %>% 
  summary()

# tidy it
arthropoda <- tidy_glmm(arthropoda_smr, "Arthropoda")


## bivalvia
load(file = here("data/model-output/bivalvia_model.Rds"))

# get traditional model summary
bivalvia_smr <- bivalvia_model$final_trend %>% 
  summary()

# tidy it
bivalvia <- tidy_glmm(bivalvia_smr, "Bivalvia")


## cnidaria
load(file = here("data/model-output/cnidaria_model.Rds"))

# get traditional model summary
cnidaria_smr <- cnidaria_model$final_trend %>% 
  summary()

# tidy it
cnidaria <- tidy_glmm(cnidaria_smr, "Cnidaria")


## echinodermata
load(file = here("data/model-output/echinodermata_model.Rds"))

# get traditional model summary
echinodermata_smr <- echinodermata_model$final_trend %>% 
  summary()

# tidy it
echinodermata <- tidy_glmm(echinodermata_smr, "Echinodermata")


## foraminifera
load(file = here("data/model-output/foraminifera_model.Rds"))

# get traditional model summary
foraminifera_smr <- foraminifera_model$final_trend %>% 
  summary()

# tidy it
foraminifera <- tidy_glmm(foraminifera_smr, "Foraminifera")


## gastropoda
load(file = here("data/model-output/gastropoda_model.Rds"))

# get traditional model summary
gastropoda_smr <- gastropoda_model$final_trend %>% 
  summary()

# tidy it
gastropoda <- tidy_glmm(gastropoda_smr, "Gastropoda")


## mammalia
load(file = here("data/model-output/mammalia_model.Rds"))

# get traditional model summary
mammalia_smr <- mammalia_model$final_trend %>% 
  summary()

# tidy it
mammalia <- tidy_glmm(mammalia_smr, "Mammalia")


## reptilia
load(file = here("data/model-output/reptilia_model.Rds"))

# get traditional model summary
reptilia_smr <- reptilia_model$final_trend %>% 
  summary()

# tidy it
reptilia <- tidy_glmm(reptilia_smr, "Reptilia")


### merge models together
model_summaries <- arthropoda %>% 
  full_join(bivalvia) %>% 
  full_join(cnidaria) %>% 
  full_join(echinodermata) %>% 
  full_join(foraminifera) %>% 
  full_join(gastropoda) %>% 
  full_join(mammalia) %>% 
  full_join(reptilia)

# save model summaries
write_csv(model_summaries, path = here("data/results/model_summaries.csv"))

# make flextable
model_summaries_fxt <- model_summaries %>% 
  # get trend
  separate(col = Interaction, sep = ":", into = c(NA, "Ttrend")) %>% 
  # round digits
  mutate(across(c(Estimate, `Std. Error`), round, 2)) %>% 
  mutate(`Pr(>|z|)` = round(`Pr(>|z|)`, 4)) %>% 
  mutate(across(c(AIC, BIC), round, 0)) %>% 
  # make a flextable
  flextable() %>% 
  # bold headers
  bold(part = "header") %>% 
  # italic taxon column
  italic(j = "Taxon") %>%
  autofit()


# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(model_summaries_fxt, pos = "after")




# mean and sd for predictions ---------------------------------------------


## function to return the distribution for each palaeoclimate interaction type
palint_dist <- function(interaction_final, dat_trends, Taxon){
  
  #  warming warming
  ww_raw <- dat_trends %>% filter(trend_genus >=0 & change_prev >= 0)
  ww_pred <- predict(interaction_final, newdata = ww_raw,
                     type = "response", allow.new.levels = TRUE)
  
  #  warming cooling 
  wc_raw <- dat_trends %>% filter(trend_genus >=0 & change_prev <= 0)
  wc_pred <- predict(interaction_final, newdata = wc_raw,
                     type = "response", allow.new.levels = TRUE)
  
  
  #  cooling cooling 
  cc_raw <- dat_trends %>% filter(trend_genus <= 0 & change_prev <= 0)
  cc_pred <- predict(interaction_final, newdata = cc_raw,
                     type = "response", allow.new.levels = TRUE)
  
  
  #  cooling warming
  cw_raw <- dat_trends %>% filter(trend_genus <=0 & change_prev >= 0)
  cw_pred <- predict(interaction_final, newdata = cw_raw,
                     type = "response", allow.new.levels = TRUE)
  
  
  tibble(taxon = Taxon, 
         ww = list(ww_pred), cw = list(cw_pred), 
         wc = list(wc_pred), cc = list(cc_pred))
}


## arthropoda

# load final model 
interaction_final <- arthropoda_model$final_trend

# load raw trend data
dat_trends <- read_csv(here("data/arthropoda_trends.csv"), guess_max = 1e5)

# get distributions
arthropoda_dist <- palint_dist(interaction_final, dat_trends, Taxon = "Arthropoda")


## bivalvia

# load final model 
interaction_final <- bivalvia_model$final_trend

# load raw trend data
dat_trends <- read_csv(here("data/bivalvia_trends.csv"), guess_max = 1e5)

# get distributions
bivalvia_dist <- palint_dist(interaction_final, dat_trends, Taxon = "Bivalvia")



## cnidaria

# load final model 
interaction_final <- cnidaria_model$final_trend

# load raw trend data
dat_trends <- read_csv(here("data/cnidaria_trends.csv"), guess_max = 1e5)

# get distributions
cnidaria_dist <- palint_dist(interaction_final, dat_trends, Taxon = "Cnidaria")



## echinodermata

# load final model 
interaction_final <- echinodermata_model$final_trend

# load raw trend data
dat_trends <- read_csv(here("data/echinodermata_trends.csv"), guess_max = 1e5)

# get distributions
echinodermata_dist <- palint_dist(interaction_final, dat_trends, 
                                  Taxon = "Echinodermata")



## foraminifera

# load final model 
interaction_final <- foraminifera_model$final_trend

# load raw trend data
dat_trends <- read_csv(here("data/foraminifera_trends.csv"), guess_max = 1e5)

# get distributions
foraminifera_dist <- palint_dist(interaction_final, dat_trends, 
                                  Taxon = "Foraminifera")


## gastropoda

# load final model 
interaction_final <- gastropoda_model$final_trend

# load raw trend data
dat_trends <- read_csv(here("data/gastropoda_trends.csv"), guess_max = 1e5)

# get distributions
gastropoda_dist <- palint_dist(interaction_final, dat_trends, 
                                 Taxon = "Gastropoda")


## mammalia

# load final model 
interaction_final <- mammalia_model$final_trend

# load raw trend data
dat_trends <- read_csv(here("data/mammalia_trends.csv"), guess_max = 1e5)

# get distributions
mammalia_dist <- palint_dist(interaction_final, dat_trends, 
                               Taxon = "Mammalia")


## reptilia

# load final model 
interaction_final <- reptilia_model$final_trend

# load raw trend data
dat_trends <- read_csv(here("data/reptilia_trends.csv"), guess_max = 1e5)

# get distributions
reptilia_dist <- palint_dist(interaction_final, dat_trends, 
                             Taxon = "Reptilia")


## merge together
all_dist <- arthropoda_dist %>% 
  full_join(bivalvia_dist) %>% 
  full_join(cnidaria_dist) %>% 
  full_join(echinodermata_dist) %>% 
  full_join(foraminifera_dist) %>% 
  full_join(gastropoda_dist) %>% 
  full_join(mammalia_dist) %>% 
  full_join(reptilia_dist) 


## calculate mean and sd 
all_dist_smr <- all_dist %>% 
  # transform to longer format
  pivot_longer(-taxon, names_to = "pal_int", values_to = "ext_risk") %>% 
  # unlist values
  unnest(cols = c(ext_risk)) %>% 
  # get mean and sd per taxon and pal_int
  group_by(taxon, pal_int) %>% 
  summarise(mean_risk = mean(ext_risk), 
            sd_risk = sd(ext_risk)) %>% 
  ungroup() %>% 
  add_column("Tchange" = rep(c("cool", "cool", "warm", "warm"), 8), 
             "Ttrend" = rep(c("cool", "warm", "cool", "warm"), 8)) %>% 
  select(Taxon = taxon, Ttrend, Tchange, Mean = mean_risk, Sd = sd_risk) 


all_dist_fxt <- all_dist_smr %>% 
  mutate(across(where(is.numeric), round, 2)) %>% 
  # make a flextable
  flextable() %>% 
  # bold headers
  bold(part = "header") %>% 
  # italic taxon column
  italic(j = "Taxon") %>%
  # vertical merging
  merge_v(j = "Taxon") %>% 
  # add horizontal reference lines
  hline(i = c(4, 8, 12, 16, 20, 24, 28, 31), part = "body",
        border = fp_border(color = "grey50", width = 1)) %>% 
  # add background colour
  bg(j = "Ttrend", i =  rep(c(TRUE, FALSE), 16), 
     bg = "#B3E3FF", part = "body") %>% 
  bg(j = "Ttrend", i =  rep(c(FALSE, TRUE), 16), 
     bg = "#FFD2A2", part = "body") %>% 
  bg(j = "Tchange", i =  rep(c(TRUE, TRUE, FALSE, FALSE), 8), 
     bg = "#B3E3FF", part = "body") %>% 
  bg(j = "Tchange", i =  rep(c(FALSE, FALSE, TRUE, TRUE), 8), 
     bg = "#FFD2A2", part = "body") %>%
  fix_border_issues() %>% 
  autofit() 
  

# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(all_dist_fxt, pos = "after")




# autocorrelation noise ---------------------------------------------------


## read in data
red_noise <- read_csv(here("data/results/autocorrelation_red_noise.csv"))

blue_noise <- read_csv(here("data/results/autocorrelation_blue_noise.csv"))


## extract estimates only and render to longer format
# red noise
red_noise <- red_noise %>% 
  select(model, warm, cool) %>% 
  pivot_longer(-model, names_to = NULL, values_to = "estimate") %>% 
  add_column(type = "red_noise") %>% 
  select(-model)

# blue noise
blue_noise <- blue_noise %>% 
  select(model, warm, cool) %>% 
  pivot_longer(-model, names_to = NULL, values_to = "estimate") %>% 
  add_column(type = "blue_noise") %>% 
  select(-model)

# null model 
null_noise <- sim_result %>% 
  pivot_longer(-nr_observations, names_to = NULL, values_to = "estimate") %>% 
  add_column(type = "null_model") %>% 
  select(-nr_observations)

# combine in one dfr and add simulation null model results
auto_noise <- blue_noise %>% 
  full_join(red_noise) %>% 
  full_join(null_noise)



# calculate summaries
auto_noise_fxt <- auto_noise %>% 
  # calculate mean, sd and 95 % Confidence Intervals
  group_by(type) %>% 
  summarise(ci = list(mean_cl_normal(estimate)),
            Sd = sd(estimate)) %>% 
  unnest(cols = c(ci)) %>% 
  # set better names
  rename(Simulation = type, Mean = y, 
         "Lower 95% CI" = ymin, "Upper 95% CI" = ymax) %>% 
  # get it in the right order
  arrange(desc(Sd)) %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  # set better value names
  mutate(Simulation = c("Null Model", "Blue Noise", "Red Noise")) %>% 
  # transform to flextable
  flextable() %>% 
  # bold headers
  bold(part = "header") %>% 
  # italic taxon column
  italic(j = "Simulation") %>%
  autofit() 


# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(auto_noise_fxt, pos = "after")





# render word file ----------------------------------------------------------

# convert to word file/ add input to empty docx
print(my_doc, target = here("tables/extended_tables.docx"))



# figures -----------------------------------------------------------------



# distributions of effect  ------------------------------------------------

# use same data as all_dist_fxt table
all_dist_plot <- all_dist %>% 
  # transform to longer format
  pivot_longer(-taxon, names_to = "pal_int", values_to = "ext_risk") %>% 
  # unlist values
  unnest(cols = c(ext_risk)) %>% 
  mutate(Ttrend = str_sub(pal_int, 1, 1), 
         Ttrend = if_else(Ttrend == "c", "Ttrend = cooling", "Ttrend = warming"), 
         Tchange = str_sub(pal_int, 2, 2), 
         Tchange = if_else(Tchange == "c", "cooling", "warming")) %>% 
  
  ggplot() +
  geom_density(aes(ext_risk, fill = Tchange), colour = "grey40", size = 0.1, 
               alpha = 0.5) + 
  facet_grid(rows = vars(Ttrend), cols = vars(taxon)) + 
  scale_fill_manual(values = c("#56B4E9", "coral2")) +
  coord_cartesian(xlim = c(0,1), ylim = c(0, 30)) + 
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  labs(x = "Extinction risk [%]") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_line(colour ="grey20"))

  
# save plot 
ggsave(all_dist_plot, filename = here("figures/raw/palint_distribution.pdf"), 
       width = 183, height = 89, units = "mm")



# autocorrelation ---------------------------------------------------------

# use data from autocorrelation table
auto_noise_plot <- ggplot(auto_noise) +
  # add reference line
  geom_vline(xintercept = 0, colour = "orange", size = 0.5) +
  geom_histogram(aes(estimate, fill = type),
                 binwidth = 0.001) +
  scale_fill_manual(values = c("#56B4E9", "grey", "coral2"), 
                    name = "Simulation", 
                    labels = c("Blue Noise", "Null Model", "Red Noise")) +
  coord_cartesian(xlim = c(-0.21, 0.1)) +
  labs(x = "Simulated change in extinction risk [%]", y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(linetype = "dashed"), 
        panel.grid.minor.y = element_blank(), 
        legend.position = c(0.2, 0.8), 
        legend.background = element_rect(fill = "white", colour = "grey40"), 
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"))
  
# save plot 
ggsave(auto_noise_plot, 
       filename = here("figures/extended-data/autocorrelation_distributions.png"), 
       width = 183, height = 89, units = "mm", dpi = 300)




# number of genera --------------------------------------------------------


## function for occurrences into stage bins
# stages
gradstein <- read_csv(here("data/gradstein.csv"))

binning <- function(x){gradstein$stg[x<gradstein$bottom  &  x>=gradstein$top]}


## read in raw fossil data/ occurrence data
all_raw_data <- read_csv(here("data/raw-fossil-data/all_raw_data.csv"))


## assign grouping variable
all_raw_data <- all_raw_data %>% 
  # make dummy column for grouping
  mutate(my_group = if_else(!is.na(phylum), phylum, class), 
         my_group2 = if_else(my_group == "Mollusca", class, my_group),
         my_group3 = if_else(my_group2 == "Chordata", "Reptilia", my_group2)) %>% 
  select(taxon = my_group3, genus,max_ma:min_ma) 


## bin data 
all_raw_data <- all_raw_data %>% 
  # add mean age
  mutate(mean_ma = (max_ma + min_ma)/2) %>% 
  filter(mean_ma <= 540) %>% 
  # add stage number
  mutate(stg = map_dbl(mean_ma, binning)) %>% 
  # Exclude unassigned bins
  drop_na(stg) 


## calculate number of genera per stage
nr_genus_plot <- all_raw_data %>% 
  count(taxon, stg, name = "nr_genus") %>% 
  left_join(gradstein %>% select(stg, mid)) %>% 
  select(taxon, age = mid, stg, nr_genus) %>% 
  ggplot() +
  geom_area(aes(x = age, y = nr_genus), fill = "grey89", colour = "grey15") +
  scale_x_reverse(limits = c(540, 0)) +
  coord_geo(size = 2.5, height = unit(0.75, "lines"), alpha = 2/3) +
  labs(x = "Age [myr]", y = "Number of genera per stage") +
  facet_wrap( ~ taxon, scales = "free", ncol = 2) +
  theme_minimal(base_size = 7) +
  theme(axis.ticks.x = element_line(),
        axis.ticks.length.x = unit(1, units = "mm"))
  
# save plot 
ggsave(nr_genus_plot, 
       filename = here("figures/extended-data/number_of_genera.jpg"), 
       width = 180, height = 185, units = "mm", dpi = 400)



### calculate diversity metrics

## remove foraminifera as we have no occurrence data for them, only ranges
taxon <- all_raw_data %>% 
  filter(taxon != "Foraminifera") %>% 
  distinct(taxon) %>% pull()

## preallocate empty list
metrics <- vector("list", length = length(taxon))

## give names to the list for subsetting
names(metrics) <- taxon

## calculate metric per phyla and save it in the list as a tibble
## note the additional phylum column, which is necessary for collapsing of the list later
for (i in taxon){
  metrics[[i]] <- all_raw_data %>% 
    # go through phyla
    filter(taxon == i) %>% 
    # calculate metrics
    divDyn(., tax = "genus", bin = "stg") %>% 
    as_tibble() %>% 
    # add phylum names and geologic ages calculated with stages file
    mutate(taxon = i, age = gradstein$mid[1:nrow(.)]) %>% 
    # order it properly
    select(taxon, age, stg, everything())  
}

## collapse or flatten the list to a dataframe for facetting  
metrics <- bind_rows(metrics)



## plot diversity
csib_div <- ggplot(metrics) +
  geom_area(aes(x = age, y = divCSIB), fill = "grey89", colour = "grey15") +
  scale_x_reverse(limits = c(540, 0)) +
  coord_geo(size = 2.5, height = unit(0.75, "lines"), alpha = 2/3) +
  labs(x = "Age [myr]", y = "Corrected sampled-in-bin diversity") +
  facet_wrap( ~ taxon, scales = "free", ncol = 2) +
  theme_minimal() +
  theme(axis.ticks.x = element_line(),
    axis.ticks.length.x = unit(1, units = "mm"))

# save plot 
ggsave(csib_div, 
       filename = here("figures/extended-data/sampled_in_bin.png"), 
       width = 183, height = 183, units = "mm", dpi = 300)


## plot number of taxa going extinct
nr_extinct <- ggplot(metrics) +
  geom_area(aes(x = age, y = tExt), fill = "grey89", colour = "grey15") +
  scale_x_reverse(limits = c(540, 0)) +
  coord_geo(size = 2.5, height = unit(0.75, "lines"), alpha = 2/3) +
  labs(x = "Age [myr]", y = "Number of taxa getting extinct") +
  facet_wrap( ~ taxon, scales = "free", ncol = 2) +
  theme_minimal() +
  theme(axis.ticks.x = element_line(),
        axis.ticks.length.x = unit(1, units = "mm"))


# save plot 
ggsave(nr_extinct, 
       filename = here("figures/extended-data/nr_extinct.png"), 
       width = 183, height = 183, units = "mm", dpi = 300)


## plot extinction rate
extinct_rate <- ggplot(metrics) +
  geom_step(aes(x = age, y = extPC), colour = "grey15") +
  scale_x_reverse(limits = c(540, 0)) +
  coord_geo(size = 2.5, height = unit(0.75, "lines"), alpha = 2/3) +
  labs(x = "Age [myr]", y = "Per capita extinction rate") +
  facet_wrap( ~ taxon, scales = "free", ncol = 2) +
  theme_minimal() +
  theme(axis.ticks.x = element_line(),
        axis.ticks.length.x = unit(1, units = "mm"))


# save plot 
ggsave(extinct_rate, 
       filename = here("figures/extended-data/ectinction_rate.png"), 
       width = 183, height = 183, units = "mm", dpi = 300)




# effect vs memory --------------------------------------------------------


## calculate mean stage duration
mean_dur <- gradstein %>% 
  summarise(mean_dur = mean(dur)) %>% 
  pull()


## get temporal memory
temp_mem <- tibble(taxon = c("Arthropoda", "Bivalvia", "Cnidaria", "Echinodermata",
                 "Foraminifera", "Gastropoda", "Mammalia", "Reptilia"), 
                 trend = c(2, 4, 4, 9, 2, 5, 10, 3)) %>% 
  # calculate temporal memory 
  mutate(temp_mem = trend * mean_dur)


## transform ext_risk to longer format
ext_risk_long <- ext_risk %>% 
  pivot_longer(cols = c(warm, cool), names_to = "type", values_to = "pal_int") %>% 
  select(taxon, type, pal_int) %>% 
  # merge with durations 
  left_join(durations) %>% 
  # merge with temporal memory
  left_join(temp_mem) %>% 
  select(-trend)

## median duration

# get r squared
median_r2 <- lm(temp_mem ~ median_duration, data = ext_risk_long) %>% 
  summary() %>% 
  .$adj.r.squared %>% 
  round(digits = 2)

# plot it
regression_plot_median <- 
  ggplot(ext_risk_long, aes(x = median_duration, y = temp_mem)) +
  # add regression line
  geom_smooth(method = "lm", colour = "orange", fill = "grey70", 
              size = 1) +
  # add points
  geom_point(fill = "grey60", shape = 21, show.legend = F,
             colour = "grey30", size = 2.2, stroke = 1) +
  annotate("text", label = paste("R^2 == ", median_r2),
           x = 22, y = 50, parse = TRUE, size = 2.5) +
  coord_cartesian(ylim = c(0, 63)) +
  labs(x = "Median duration of genera [myr]", y = "Temporal Memory [myr]") +
  theme_minimal(base_size = 7) +
  theme(panel.grid.minor = element_blank())



## mean duration

# get r squared
mean_r2 <- lm(temp_mem ~ mean_duration, data = ext_risk_long) %>% 
  summary() %>% 
  .$adj.r.squared %>% 
  round(digits = 2)

# plot it
regression_plot_mean <- 
  ggplot(ext_risk_long, aes(x = mean_duration, y = temp_mem)) +
  # add regression line
  geom_smooth(method = "lm", colour = "orange", fill = "grey70", 
              size = 1) +
  # add points
  geom_point(fill = "grey60", shape = 21, show.legend = F,
             colour = "grey30", size = 2.2, stroke = 1) +
  annotate("text", label = paste("R^2 == ", mean_r2),
           x = 36.5, y = 50, parse = TRUE, size = 2.5) +
  coord_cartesian(ylim = c(0, 63)) +
  labs(x = "Mean duration of genera [myr]", y = "Temporal Memory [myr]") +
  theme_minimal(base_size = 7) +
  theme(panel.grid.minor = element_blank())


## combine plots
temp_mem_comb <- regression_plot_median / regression_plot_mean +
  plot_annotation(tag_levels = 'a')

# save plot 
ggsave(temp_mem_comb, 
       filename = here("figures/extended-data/temporal_memory_regression.jpg"), 
       width = 180, height = 185, units = "mm", dpi = 300)
