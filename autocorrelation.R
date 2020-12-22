# load packages
library(tidyverse)
library(here)
library(DHARMa) 


# functions ---------------------------------------------------------------

### temp_ac

# wrapper function to calculate temporal autocorrelation per model
temp_ac <- function(mymodel){
  # simulate residuals
  mymodel_res <- simulateResiduals(mymodel)
  # calculate autocorrelation based on stages
  mymodel_ac <- testTemporalAutocorrelation(mymodel_res, time = 14:94, plot = FALSE)
  mymodel_ac$statistic
}
# load data ---------------------------------------------------------------

# arthropoda
load(here("data/model-output/arthropoda_model.Rds"))
arthropoda <- arthropoda_model$final_trend

# bivalvia
load(here("data/model-output/bivalvia_model.Rds"))
bivalvia <- bivalvia_model$final_trend

# cnidaria
load(here("data/model-output/cnidaria_model.Rds"))
cnidaria <- cnidaria_model$final_trend

# echinodermata
load(here("data/model-output/echinodermata_model.Rds"))
echinodermata <- echinodermata_model$final_trend

# foraminifera
load(here("data/model-output/foraminifera_model.Rds"))
foraminifera <- foraminifera_model$final_trend

# gastropoda
load(here("data/model-output/gastropoda_model.Rds"))
gastropoda <- gastropoda_model$final_trend

# mammalia
load(here("data/model-output/mammalia_model.Rds"))
mammalia <- mammalia_model$final_trend

# reptilia
load(here("data/model-output/reptilia_model.Rds"))
reptilia <- reptilia_model$final_trend



# calculate autocorrelation -----------------------------------------------

# build dataframe for model output
df_ac <- tibble(taxon = c("Arthropoda", "Bivalvia", "Cnidaria", "Echinodermata", 
                          "Foraminifera", "Gastropoda", "Mammalia", "Reptilia"),
                durbin_watson = double(8))

# calculate Durbin-Watson statistic and assign to dataframe
df_ac[1 , 2] <- temp_ac(arthropoda)
df_ac[2 , 2] <- temp_ac(bivalvia)
df_ac[3 , 2] <- temp_ac(cnidaria)
df_ac[4 , 2] <- temp_ac(echinodermata)
df_ac[5 , 2] <- temp_ac(foraminifera)
df_ac[6 , 2] <- temp_ac(gastropoda)
df_ac[7 , 2] <- temp_ac(mammalia)
df_ac[8 , 2] <- temp_ac(reptilia)


# save data
write_csv(df_ac, path = here("data/results/autocorrelation_results.csv"))


# plotting
df_ac_plot <- df_ac %>%
  # plot it
  ggplot(aes(x = durbin_watson, y = fct_reorder(taxon, desc(taxon)))) +
  geom_vline(xintercept = 2, colour = "orange") +
  geom_vline(xintercept = c(1, 3), linetype = "dotted") +
  geom_point(pch = 21, size = 2.6, colour = "grey10", fill = "grey70") +
  coord_cartesian(xlim = c(0, 4)) + 
  labs(y = NULL, x = "Durbin-Watson statistic") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank())



ggsave(df_ac_plot, filename = here("figures/extended-data/autocorrelation.png"), 
       width = 183, height = 89, units = "mm", dpi = 300)


