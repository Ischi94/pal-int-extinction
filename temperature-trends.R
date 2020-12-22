# load packages
library(tidyverse)
library(here)


# import stage data from Gradstein 2012
gradstein <- read_csv(here("data/gradstein.csv"))

# prepare weizer & prokoph temperature data.
# We are loading the file which was already processed as described
# in the methods paragraph and as described by Reddin et al. 2017
isotemp <- 
  # load file
  read_csv(file=here("data/TimeSeriesUsed.csv")) %>% 
  # remove redundant columns and make clean names
  select(stage = Stage, temp = Temp) %>% 
  # add numeric age and empty column for subsequent calculation  
  add_column(age = gradstein$mid[14:94], change_prev = NA) %>% 
  mutate(change_prev = as.double(change_prev))

# take a glimpse
ggplot(data=isotemp, aes(x = age, y = temp))+
  geom_line(colour = "grey30", size = 1) + 
  scale_x_reverse() +
  labs(x = "Age [myr]", y = "Temperature [°C]") +
  theme_classic()

# set up new dfr
isotemp_trends <- isotemp

# fill in the values for short-term change using the lm as we do for the long term trend as well
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

# save data
write_csv(isotemp_trends, path = here("data/isotemp_trends.csv"))
