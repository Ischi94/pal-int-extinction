# set working directory
setwd("C:/Users/gmath/Documents/4.Semester/Masterarbeit") 

#remove workspace
rm(list = ls())

# Load range through and temperature data
source("17dinos_and_birds_pbdb.R")



#find outliers at the short term temp change
#outliers = boxplot(fourteen_trends$change.prev, plot=FALSE)$out
#Extract the outliers from the original data frame
#fourteen_trends <- fourteen_trends[!fourteen_trends$change.prev %in% outliers,]
#######################################

library(ggplot2)
library(lme4)
library(visreg)
library(geiger)
library(stargazer)
#remove.packages("plyr")
#######################################

### we start with four trends ###
###----
# replace the NA's with zeros as we have no trend there
thirteen_trends_dinos$trend1[is.na(thirteen_trends_dinos$trend1)] <- 0
thirteen_trends_dinos$trend2[is.na(thirteen_trends_dinos$trend2)] <- 0
thirteen_trends_dinos$trend3[is.na(thirteen_trends_dinos$trend3)] <- 0

# Split short term temperature change into warming and cooling:
thirteen_trends_dinos$cooling<-ifelse(thirteen_trends_dinos$change.prev<0, thirteen_trends_dinos$change.prev, NA)
thirteen_trends_dinos$warming<-ifelse(thirteen_trends_dinos$change.prev>0, thirteen_trends_dinos$change.prev, NA)

## Warming
# Iterate through each long term temperature change
vars = names(dplyr::select(thirteen_trends_dinos, trend1:trend3, trend.st1:trend.st10)) 
modelsw4 = lapply(setNames(vars, vars), function(var) {
  form = paste("extinct~warming:", var, "+(1|Genus)")
  glmer(form, data=thirteen_trends_dinos, family="binomial")
})

# Make data frame for model output
warm_fourteen_trends <- data.frame(model=names(dplyr::select(thirteen_trends_dinos, 
                                                             trend1:trend3, trend.st1:trend.st10)), intercept=NA, 
                                   interaction=NA, AIC=NA)

# Run loop to fill for.warm (coefficients and p-values)
for (i in warm_fourteen_trends$model) {
  sum <- summary(modelsw4[[i]])
  warm_fourteen_trends[warm_fourteen_trends$model==i, "intercept"]<-
    paste(round(sum$coefficients[1, 1], 3), sep=" ", "±",
          round(sum$coefficients[1, 2], 3),
          ifelse(sum$coefficients[1, 4]<0.001, "***", 
                 ifelse(sum$coefficients[1, 4]<0.01, "**", 
                        ifelse(sum$coefficients[1, 4]<0.05, "*", ""))))
  warm_fourteen_trends[warm_fourteen_trends$model==i, "interaction"]<-
    paste(round(sum$coefficients[2, 1], 3), sep=" ", "±",
          round(sum$coefficients[2, 2], 3),
          ifelse(sum$coefficients[2, 4]<0.001, "***", 
                 ifelse(sum$coefficients[2, 4]<0.01, "**", 
                        ifelse(sum$coefficients[2, 4]<0.05, "*", ""))))
  warm_fourteen_trends[warm_fourteen_trends$model==i, "AIC"]<-as.numeric(round(sum$AICtab[[1]], 1))
}

# Add column weith AIC weights
warm_fourteen_trends$dAIC<- as.numeric(round(aicw(warm_fourteen_trends$AIC)$delta, 1))
warm_fourteen_trends$AIC.weights<- as.numeric(signif(aicw(warm_fourteen_trends$AIC)$w, 3))

# Save output as csv
#write.table(for.warm, "for_warm.csv")

# the  best models has the long term trend calculated over 2 stages
# let's check for overdispersion with the function from markus hall
source("overdisp_fun.R")

overdisp_fun(modelsw4[[which(warm_fourteen_trends$dAIC==0)[1]]])   # not overdispersed, yayyy :D

summary(modelsw4[[which(warm_fourteen_trends$dAIC==0)[1]]])

s <- summary(modelsw4[[which(warm_fourteen_trends$dAIC==0)[1]]])

capture.output(s, file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/Model_summaries/Dinos&Reps_warm.txt")

coefm1 <- coef(modelsw4[[which(warm_fourteen_trends$dAIC==0)[1]]])
coefm1$Genus[1:10,]


#Make a visreg cross-sectional plot by long-term of the best model

warm.p1<-visreg(modelsw4[[which(warm_fourteen_trends$dAIC==0)[1]]], "warming", scale="response", 
                
                by=paste(warm_fourteen_trends[which(warm_fourteen_trends$dAIC==0)[1],1]), 
#                breaks = c(-0.2,0, 0.1),
                
                rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5)) 

warm.p1<-warm.p1 + labs(x="Warming", y="Extinction response") +
  
  theme_classic() + theme(text = element_text(size = 14))


# let's take the values of the visreg and split them into cooling and warming and compare the 
# exctinction response of each by means of boxplots
vis_out_w4 <- visreg(modelsw4[[which(warm_fourteen_trends$dAIC==0)[1]]], "warming", scale="response", 
                     
                     by=paste(warm_fourteen_trends[which(warm_fourteen_trends$dAIC==0)[1],1]), 
 #                    breaks = c(-0.1,0, 0.1),
                     
                     rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5), plot = F) 


vis_out_w4$fit$ctrend <- ifelse(vis_out_w4$fit$trend.st10<0, "Trend = Cooling", "Trend = Warming")

write.table(vis_out_w4$fit, 
            file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/Visreg/Dinos&Reps_warm.csv")


# Basic box plot
ggplot(vis_out_w4$fit, aes(x=ctrend, y=visregFit, fill=ctrend)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#56B4E9", "coral2")) +
  labs(title="Short-Term Warming", x="", y = "P(Extinction)")+
  theme_classic() +
  theme(legend.position = "none", axis.text=element_text(size=13),
        axis.title=element_text(size=14), plot.title = element_text(size=16))



# both are not normally distributed, so we are using a Mann-Whitney U test
wilcox.test(vis_out_w4$fit$visregFit[vis_out_w4$fit$ctrend=="Trend = Warming"],
            vis_out_w4$fit$visregFit[vis_out_w4$fit$ctrend=="Trend = Cooling"],
            paired = F,  conf.int = T)

## Cooling
# Iterate through each long term temperature change
vars = names(dplyr::select(thirteen_trends_dinos, trend1:trend3, trend.st1:trend.st10)) 
modelsc4 = lapply(setNames(vars, vars), function(var) {
  form = paste("extinct~cooling:", var, "+(1|Genus)")
  glmer(form, data=thirteen_trends_dinos, family="binomial")
})

# Make data frame for model output
cool_fourteen_trends <-data.frame(model=names(dplyr::select(thirteen_trends_dinos, trend1:trend3, trend.st1:trend.st10)), intercept=NA, 
                                 interaction=NA, AIC=NA)


# Run loop to fill for.cool (coefficients and p-values)
for (i in cool_fourteen_trends$model) {
  sum <- summary(modelsc4[[i]])
  cool_fourteen_trends[cool_fourteen_trends$model==i, "intercept"]<-
    paste(round(sum$coefficients[1, 1], 3), sep=" ", "±",
          round(sum$coefficients[1, 2], 3),
          ifelse(sum$coefficients[1, 4]<0.001, "***", 
                 ifelse(sum$coefficients[1, 4]<0.01, "**", 
                        ifelse(sum$coefficients[1, 4]<0.05, "*", ""))))
  cool_fourteen_trends[cool_fourteen_trends$model==i, "interaction"]<-
    paste(round(sum$coefficients[2, 1], 3), sep=" ", "±",
          round(sum$coefficients[2, 2], 3),
          ifelse(sum$coefficients[2, 4]<0.001, "***", 
                 ifelse(sum$coefficients[2, 4]<0.01, "**", 
                        ifelse(sum$coefficients[2, 4]<0.05, "*", ""))))
  cool_fourteen_trends[cool_fourteen_trends$model==i, "AIC"]<-as.numeric(round(sum$AICtab[[1]], 1))
}

# Add column weith AIC weights
cool_fourteen_trends$dAIC<- as.numeric(round(aicw(cool_fourteen_trends$AIC)$delta, 1))
cool_fourteen_trends$AIC.weights<- as.numeric(signif(aicw(cool_fourteen_trends$AIC)$w, 3))

# Save output as csv
#write.table(for.cool, "for_cool.csv")

# the  best models has the long term trend calculated over 2 stages
# let's check for overdispersion with the function from markus hall
overdisp_fun(modelsc4[[which(cool_fourteen_trends$dAIC==0)]])   # not overdispersed, yayyy :D

summary(modelsc4[[which(cool_fourteen_trends$dAIC==0)]])

s <- summary(modelsc4[[which(cool_fourteen_trends$dAIC==0)]])

capture.output(s, file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/Model_summaries/Dinos&Reps_cool.txt")

coefm2 <- coef(modelsc4[[which(cool_fourteen_trends$dAIC==0)]])
coefm2$Genus[1:10,]

#for_cool <- read.table("for_cool.csv", sep=" ")
#for_warm <- read.table("for_warm.csv", sep = " ")


#Make a visreg cross-sectional plot by long-term of the best model

cool.p1<-visreg(modelsc4[[which(cool_fourteen_trends$dAIC==0)]], "cooling", scale="response", 
                
                by=paste(cool_fourteen_trends[which(cool_fourteen_trends$dAIC==0),1]), 
#                               breaks = c(-0.2,0.2),
                
                rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5)) 

cool.p1<-cool.p1 + labs(x="Cooling", y="Extinction response") +
  
  theme_classic() + theme(text = element_text(size = 14))

#ggsave("dia_cool_best.pdf", cool.p1)

# let's take the values of the visreg and split them into cooling and warming and compare the 
# exctinction response of each by means of boxplots
vis_out_c4 <- visreg(modelsc4[[which(cool_fourteen_trends$dAIC==0)]], "cooling", scale="response", 
                     
                     by=paste(cool_fourteen_trends[which(cool_fourteen_trends$dAIC==0),1]), 
                                          breaks = c(-0.1,0.1),
                     
                     rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5), plot = F) 



vis_out_c4$fit$ctrend <- ifelse(vis_out_c4$fit$trend.st3<0,"Trend = Cooling", "Trend = Warming")

write.table(vis_out_c4$fit, 
            file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/Visreg/Dinos&Reps_cool.csv")


# Basic box plot
ggplot(vis_out_c4$fit, aes(x=ctrend, y=visregFit, fill=ctrend)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#56B4E9", "coral2")) +
  labs(title="Short-Term Cooling", x="", y = "P(Extinction)")+
  theme_classic() +
  theme(legend.position = "none", axis.text=element_text(size=13),
        axis.title=element_text(size=14), plot.title = element_text(size=16))


# both are not normally distributed, so we are using a Mann-Whitney U test
wilcox.test(vis_out_c4$fit$visregFit[vis_out_c4$fit$ctrend=="Trend = Cooling"],
            vis_out_c4$fit$visregFit[vis_out_c4$fit$ctrend=="Trend = Warming"], 
            paired = F, conf.int = T)

# Save output as csv
#write.table(warm_fourteen_trends, "modell_summary_D&B_w.csv")
#write.table(cool_fourteen_trends, "modell_summary_D&B_c.csv")

# Plot model quality (dAIC vs long-term trend)

# Shows how long the signal of the past influences extinction of focal species. 

#cool.p2<-ggplot(cool_fourteen_trends, aes(x=c(1:length(model)), y=dAIC))

#cool.p2 <- cool.p2 + geom_linerange(aes(ymin=0, ymax=dAIC), 
#                                    colour="#d8ebed", linetype=3) + 
#  geom_line(colour="#65a3a4", size=1.5) +
#  scale_x_continuous("Long-term duration (My)", 1:14, limits=c(1,14)) + 
#  ylab(expression(delta*"AIC")) + 
#  theme_classic() + theme(text = element_text(size = 14))

#ggsave(file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/AIC_plots/dinos_cool.pdf",cool.p2)



# Plot model quality (dAIC vs long-term trend)

# Shows how long the signal of the past influences extinction of focal species. 

#warm.p2<-ggplot(warm_fourteen_trends, aes(x=c(1:length(model)), y= dAIC))

#warm.p2<-warm.p2 + geom_linerange(aes(ymin=0, ymax=dAIC), 
#                                  colour="#d8ebed", linetype=3) + 
#  geom_line(colour="#65a3a4", size=1.5) +
#  scale_x_continuous("Long-term duration (My)", 1:14, limits=c(1,14)) + 
#  ylab(expression(delta*"AIC")) + theme_classic() + theme(text = element_text(size = 14))

#ggsave(file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/AIC_plots/dinos_warm.pdf",warm.p2)

