#### 

### Let's try Jeroen's GLMM with my range through data with information about short temp-change and long
# temp-trends (up to 10 stages, without the respective short-term change)


#remove workspace
rm(list = ls())

# Load range through and temperature data
source("9fourteen_trends.R")
source("12foram_planktonic.R")

# combine them 
fourteen_trends_foram_plank <- fourteen_trends_foram_plank[, -c(14:15)]

fourteen_trends <- rbind(fourteen_trends, fourteen_trends_foram_plank)
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

# replace the NA's with zeros as we have no trend there
fourteen_trends$trend1[is.na(fourteen_trends$trend1)] <- 0
fourteen_trends$trend2[is.na(fourteen_trends$trend2)] <- 0
fourteen_trends$trend3[is.na(fourteen_trends$trend3)] <- 0
fourteen_trends$trend4[is.na(fourteen_trends$trend4)] <- 0

# Split short term temperature change into warming and cooling:
fourteen_trends$cooling<-ifelse(fourteen_trends$change.prev<0, fourteen_trends$change.prev, NA)
fourteen_trends$warming<-ifelse(fourteen_trends$change.prev>0, fourteen_trends$change.prev, NA)

## Warming
# Iterate through each long term temperature change
vars = names(dplyr::select(fourteen_trends, c(trend1:trend4,trend.st1:trend.st10))) 
modelsw = lapply(setNames(vars, vars), function(var) {
  form = paste("extinct~warming:", var, "+(1|Genus)") 
  glmer(form, data=fourteen_trends, family="binomial")
})

# Make data frame for model output
warm_fourteen_trends <- data.frame(model=names(dplyr::select(fourteen_trends, 
                                                             c(trend1:trend4,trend.st1:trend.st10))), intercept=NA, 
                                   interaction=NA, AIC=NA)

# Run loop to fill for.warm (coefficients and p-values)
for (i in warm_fourteen_trends$model) {
  sum <- summary(modelsw[[i]])
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

overdisp_fun(modelsw[[which(warm_fourteen_trends$dAIC==0)]])   # not overdispersed, yayyy :D

summary(modelsw[[which(warm_fourteen_trends$dAIC==0)]])

s <- summary(modelsw[[which(warm_fourteen_trends$dAIC==0)]])

capture.output(s, file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/Model_summaries/Forams_warm.txt")

#coef(modelsw[[which(warm_fourteen_trends$dAIC==0)]])

# create a data frame where we save all effects plus all overdisp
# build an empty data frame to save results for the new data
warm_fourteen_trends$lower_CI <- NA
warm_fourteen_trends$estimate <- NA
warm_fourteen_trends$upper_CI <- NA
warm_fourteen_trends$overdisp <- NA

for (i in 1:14) {
  warm_fourteen_trends$overdisp[i] <- overdisp_fun(modelsw[[i]])[4] 
  
  vis_out_w <- visreg(modelsw[[i]], "warming", scale="response", 
                        
                        by=paste(warm_fourteen_trends[i,1]), 
                        #                breaks = c(-0.5,0.5),
                        
                        rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5), plot = F) 
  
  
  vis_out_w$fit$ctrend <- ifelse(vis_out_w$fit[2]<0, "Trend = Cooling", "Trend = Warming")
  
  wil_test <- wilcox.test(vis_out_w$fit$visregFit[vis_out_w$fit$ctrend=="Trend = Warming"],
              vis_out_w$fit$visregFit[vis_out_w$fit$ctrend=="Trend = Cooling"],
              paired = F,  conf.int = T)
  warm_fourteen_trends$lower_CI[i] <- wil_test$conf.int[1]
  warm_fourteen_trends$upper_CI[i] <- wil_test$conf.int[2]
  warm_fourteen_trends$estimate[i] <- wil_test$estimate
  
}

#Make a visreg cross-sectional plot by long-term of the best model
# if we set scale= response we change the y axis from log odds to probality

warm.p1<-visreg(modelsw[[which(warm_fourteen_trends$dAIC==0)]], "warming", 
                scale="response", 
                #                breaks = c(-0.175,0.014),
                
                by=paste(warm_fourteen_trends[which(warm_fourteen_trends$dAIC==0),1]),
                
                rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5)) 

warm.p1<-warm.p1 + labs(x="Warming", y="P (Extinction)") +
  
  theme_classic() + theme(text = element_text(size = 14))

visreg2d(modelsw[[which(warm_fourteen_trends$dAIC==0)]], "warming", "trend.st4", 
         scale="response", color.palette=colorRampPalette(c("green", "red", "darkred")),
         plot.type="image")

# let's take the values of the visreg and split them into cooling and warming and compare the 
# exctinction response of each by means of boxplots
vis_out_w14 <- visreg(modelsw[[which(warm_fourteen_trends$dAIC==0)]], "warming", scale="response", 
                      
                      by=paste(warm_fourteen_trends[which(warm_fourteen_trends$dAIC==0),1]), 
                      #                breaks = c(-0.5,0.5),
                      
                      rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5), plot = F) 


vis_out_w14$fit$ctrend <- ifelse(vis_out_w14$fit$trend.st4<0, "Trend = Cooling", "Trend = Warming")

write.table(vis_out_w14$fit, 
            file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/Visreg/Forams_warm.csv")

# Basic box plot
windows(7,5)
ggplot(vis_out_w14$fit, aes(x=ctrend, y=visregFit, fill=ctrend)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#56B4E9", "coral2")) +
  labs(title="Short-Term Warming", x="", y = "P(Extinction)")+
  theme_classic() +
  theme(legend.position = "none", axis.text=element_text(size=13),
        axis.title=element_text(size=14), plot.title = element_text(size=16))



# both are not normally distributed, so we are using a Mann-Whitney U test
wilcox.test(vis_out_w14$fit$visregFit[vis_out_w14$fit$ctrend=="Trend = Warming"],
            vis_out_w14$fit$visregFit[vis_out_w14$fit$ctrend=="Trend = Cooling"],
            paired = F,  conf.int = T)

warm.p2<-ggplot(warm_fourteen_trends, aes(x=c(1:length(modelsw)), y=dAIC))

warm.p2 + 
  geom_line(colour="#65a3a4", size=1.5) +
  ylab("dAIC") + 
  geom_line(data = warm_fourteen_trends, aes(y= estimate*75), colour="coral2", size=1.5) +
  my_theme 

## Cooling
# Iterate through each long term temperature change
vars = names(dplyr::select(fourteen_trends, c(trend1:trend4,trend.st1:trend.st10))) 
modelsc = lapply(setNames(vars, vars), function(var) {
  form = paste("extinct~cooling:", var, "+(1|Genus)") 
  # let us try it with crossed random effect with different slopes for Genus, trend and interaction
  glmer(form, data=fourteen_trends, family="binomial")
})

# Make data frame for model output
cool_fourteen_trends <-data.frame(model=names(dplyr::select(fourteen_trends, 
                                                            c(trend1:trend4,trend.st1:trend.st10))), intercept=NA, 
                                  interaction=NA, AIC=NA)

# Run loop to fill for.cool (coefficients and p-values)
for (i in cool_fourteen_trends$model) {
  sum <- summary(modelsc[[i]])
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
# create a data frame where we save all effects plus all overdisp
# build an empty data frame to save results for the new data
cool_fourteen_trends$lower_CI <- NA
cool_fourteen_trends$estimate <- NA
cool_fourteen_trends$upper_CI <- NA
cool_fourteen_trends$overdisp <- NA

for (i in 1:14) {
  cool_fourteen_trends$overdisp[i] <- overdisp_fun(modelsw[[i]])[4] 
  
  vis_out_c <- visreg(modelsc[[i]], "cooling", scale="response", 
                      
                      by=paste(cool_fourteen_trends[i,1]), 
                      #                breaks = c(-0.5,0.5),
                      
                      rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5), plot = F) 
  
  
  vis_out_c$fit$ctrend <- ifelse(vis_out_c$fit[2]<0, "Trend = Cooling", "Trend = Warming")
  
  wil_test <- wilcox.test(vis_out_c$fit$visregFit[vis_out_c$fit$ctrend=="Trend = Cooling"],
                          vis_out_c$fit$visregFit[vis_out_c$fit$ctrend=="Trend = Warming"],
                          paired = F,  conf.int = T)
  cool_fourteen_trends$lower_CI[i] <- wil_test$conf.int[1]
  cool_fourteen_trends$upper_CI[i] <- wil_test$conf.int[2]
  cool_fourteen_trends$estimate[i] <- wil_test$estimate
  
}
# the  best models has the long term trend calculated over 2 stages
# let's check for overdispersion with the function from markus hall
overdisp_fun(modelsc[[which(cool_fourteen_trends$dAIC==0)]]) # again, not overdispersed!

summary(modelsc[[which(cool_fourteen_trends$dAIC==0)]])

s <- summary(modelsc[[which(cool_fourteen_trends$dAIC==0)]])

capture.output(s, file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/Model_summaries/Forams_cool.txt")

#coef(modelsc[[which(cool_fourteen_trends$dAIC==0)]])
#for_cool <- read.table("for_cool.csv", sep=" ")
#for_warm <- read.table("for_warm.csv", sep = " ")


#Make a visreg cross-sectional plot by long-term of the best model

cool.p1<-visreg(modelsc[[which(cool_fourteen_trends$dAIC==0)]], "cooling", 
                scale="response", 
                
                by=paste(cool_fourteen_trends[which(cool_fourteen_trends$dAIC==0),1]),
                
                rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5)) 

cool.p1<-cool.p1 + labs(x="Cooling", y="Extinction response") +
  
  theme_classic() + theme(text = element_text(size = 14))

# let's take the values of the visreg and split them into cooling and warming and compare the 
# exctinction response of each by means of boxplots
vis_out_c14 <- visreg(modelsc[[which(cool_fourteen_trends$dAIC==0)]], "cooling", scale="response", 
                      
                      by=paste(cool_fourteen_trends[which(cool_fourteen_trends$dAIC==0),1]), 
                      #                breaks = c(-0.5,0.5),
                      
                      rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5), plot = F) 


vis_out_c14$fit$ctrend <- ifelse(vis_out_c14$fit$trend.st3<0, "Trend = Cooling", "Trend = Warming")

#write.table(vis_out_c14$fit, 
#            file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/Visreg/Forams_cool.csv")

# Basic box plot
ggplot(vis_out_c14$fit, aes(x=ctrend, y=visregFit, fill=ctrend)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#56B4E9", "coral2")) +
  labs(title="Short-Term Cooling", x="", y = "P(Extinction)")+
  theme_classic() +
  theme(legend.position = "none", axis.text=element_text(size=13),
        axis.title=element_text(size=14), plot.title = element_text(size=16))



# both are not normally distributed, so we are using a Mann-Whitney U test
wilcox.test(vis_out_c14$fit$visregFit[vis_out_c14$fit$ctrend=="Trend = Cooling"],
            vis_out_c14$fit$visregFit[vis_out_c14$fit$ctrend=="Trend = Warming"],
            paired = F,  conf.int = T)


# Save output as csv
#write.table(warm_fourteen_trends, "modell_summary_LBF_w.csv")
#write.table(cool_fourteen_trends, "modell_summary_LBF_c.csv")


# Plot model quality (dAIC vs long-term trend)

# Shows how long the signal of the past influences extinction of focal species. 

#cool.p2<-ggplot(cool_fourteen_trends, aes(x=c(1:length(model)), y=dAIC))

#cool.p2 <- cool.p2 + geom_linerange(aes(ymin=0, ymax=dAIC), 
#                                    colour="#d8ebed", linetype=3) + 
#  geom_line(colour="#65a3a4", size=1.5) +
#  scale_x_continuous("Long-term duration (My)", 1:14, limits=c(1,14)) + 
#  ylab(expression(delta*"AIC")) + 
#  theme_classic() + theme(text = element_text(size = 14))

#ggsave(file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/AIC_plots/forams_cool.pdf",cool.p2)



# Plot model quality (dAIC vs long-term trend)

# Shows how long the signal of the past influences extinction of focal species. 

#warm.p2<-ggplot(warm_fourteen_trends, aes(x=c(1:length(model)), y= dAIC))

#warm.p2<-warm.p2 + geom_linerange(aes(ymin=0, ymax=dAIC), 
#                                  colour="#d8ebed", linetype=3) + 
#  geom_line(colour="#65a3a4", size=1.5) +
#  scale_x_continuous("Long-term duration (My)", 1:14, limits=c(1,14)) + 
#  ylab(expression(delta*"AIC")) + theme_classic() + theme(text = element_text(size = 14))

#ggsave(file = "C:/Users/gmath/Documents/4.Semester/Masterarbeit/AIC_plots/forams_warm.pdf",warm.p2)
