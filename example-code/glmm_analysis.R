# Here we show how to 

# 1 Prepare fossil data for Generalised Mixed Effect Models (GLMM's)
  # 2 Calculate the GLMM
    # 2.1 For warming
      # 2.2 For cooling
        # 3 Analyse the performance of the GLMM
          # 4 Compare palaeoclimate interaction models to traditional ones
            # 5 Extract the results from the GLMM 
              # 5.1 For Warming
                # 5.2 For Cooling




# set working directory to where files are stored using the "rstudioapi" package

# Getting the path of this script
current_path = rstudioapi::getActiveDocumentContext()$path 

# setting it as working directory 
setwd(dirname(current_path ))

# please cross-check that this is the path where you store the other 
# files used for this script
getwd()


# loading packages
library(ggplot2)
library(lme4)
library(visreg)
library(geiger)
library(stargazer)
library(here)


# function needed for overdispersion test

overdisp_fun <- function(model) {
  ## number of variance parameters in an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m) * (nrow(m) + 1)/2
  }
  # The next two lines calculate the residual degrees of freedom
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model)) - model.df
  # extracts the Pearson residuals
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  # Generates a p-value. If less than 0.05, the data are overdispersed.
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}


# load data
# if you have processed other data than bivalves in the "data_preparation.R" file, you can load it in here and 
# use it for the GLMM. But keep in mind that some of these contain an additional taxonmical trend, which needs 
# to be considered in the following calculation as well
thirteen_trends_bivalves <- read.table(here("example-code/data/thirteen_trends_bivalves.csv"))




# 1 Prepare fossil data for Generalised Mixed Effect Models (GLMM' --------


# replace NA's with zeros as we have no trend there
thirteen_trends_bivalves$trend1[is.na(thirteen_trends_bivalves$trend1)] <- 0
thirteen_trends_bivalves$trend2[is.na(thirteen_trends_bivalves$trend2)] <- 0
thirteen_trends_bivalves$trend3[is.na(thirteen_trends_bivalves$trend3)] <- 0

# Split short term temperature change into warming and cooling, to calculate the results for each:
thirteen_trends_bivalves$cooling<-ifelse(thirteen_trends_bivalves$change.prev<0, thirteen_trends_bivalves$change.prev, NA)
thirteen_trends_bivalves$warming<-ifelse(thirteen_trends_bivalves$change.prev>0, thirteen_trends_bivalves$change.prev, NA)



# 2 Calculate the GLMM ----------------------------------------------------



# 2.1  For Warming --------------------------------------------------------

# Iterate through each long term temperature change
vars = names(dplyr::select(thirteen_trends_bivalves, trend1:trend3, trend.st1:trend.st10)) 
modelsw4 = lapply(setNames(vars, vars), function(var) {
  form = paste("extinct~warming:", var, "+(1|Genus)")
  glmer(form, data=thirteen_trends_bivalves, family="binomial")
})


# Make data frame for model output
warm_thirteen_trends <- data.frame(model=names(dplyr::select(thirteen_trends_bivalves,
                                                             trend1:trend3, trend.st1:trend.st10)), 
                                   intercept=NA, interaction=NA, AIC=NA)

# Run loop to fill for.warm (coefficients and p-values)
for (i in warm_thirteen_trends$model) {
  sum <- summary(modelsw4[[i]])
  warm_thirteen_trends[warm_thirteen_trends$model==i, "intercept"]<-
    paste(round(sum$coefficients[1, 1], 3), sep=" ", "±",
          round(sum$coefficients[1, 2], 3),
          ifelse(sum$coefficients[1, 4]<0.001, "***", 
                 ifelse(sum$coefficients[1, 4]<0.01, "**", 
                        ifelse(sum$coefficients[1, 4]<0.05, "*", ""))))
  warm_thirteen_trends[warm_thirteen_trends$model==i, "interaction"]<-
    paste(round(sum$coefficients[2, 1], 3), sep=" ", "±",
          round(sum$coefficients[2, 2], 3),
          ifelse(sum$coefficients[2, 4]<0.001, "***", 
                 ifelse(sum$coefficients[2, 4]<0.01, "**", 
                        ifelse(sum$coefficients[2, 4]<0.05, "*", ""))))
  warm_thirteen_trends[warm_thirteen_trends$model==i, "AIC"]<-as.numeric(round(sum$AICtab[[1]], 1))
}

# Add column weith AIC weights
warm_thirteen_trends$dAIC<- as.numeric(round(aicw(warm_thirteen_trends$AIC)$delta, 1))
warm_thirteen_trends$AIC.weights<- as.numeric(signif(aicw(warm_thirteen_trends$AIC)$w, 3))




# 2.2 For Cooling ---------------------------------------------------------


# Iterate through each long term temperature change
vars = names(dplyr::select(thirteen_trends_bivalves, trend1:trend3, trend.st1:trend.st10)) 
modelsc4 = lapply(setNames(vars, vars), function(var) {
  form = paste("extinct~cooling:", var, "+(1|Genus)")
  glmer(form, data=thirteen_trends_bivalves, family="binomial")
})

# Make data frame for model output
cool_thirteen_trends <-data.frame(model=names(dplyr::select(thirteen_trends_bivalves,
                                                            trend1:trend3, trend.st1:trend.st10)), 
                                  intercept=NA, interaction=NA, AIC=NA)


# Run loop to fill for.cool (coefficients and p-values)
for (i in cool_thirteen_trends$model) {
  sum <- summary(modelsc4[[i]])
  cool_thirteen_trends[cool_thirteen_trends$model==i, "intercept"]<-
    paste(round(sum$coefficients[1, 1], 3), sep=" ", "±",
          round(sum$coefficients[1, 2], 3),
          ifelse(sum$coefficients[1, 4]<0.001, "***", 
                 ifelse(sum$coefficients[1, 4]<0.01, "**", 
                        ifelse(sum$coefficients[1, 4]<0.05, "*", ""))))
  cool_thirteen_trends[cool_thirteen_trends$model==i, "interaction"]<-
    paste(round(sum$coefficients[2, 1], 3), sep=" ", "±",
          round(sum$coefficients[2, 2], 3),
          ifelse(sum$coefficients[2, 4]<0.001, "***", 
                 ifelse(sum$coefficients[2, 4]<0.01, "**", 
                        ifelse(sum$coefficients[2, 4]<0.05, "*", ""))))
  cool_thirteen_trends[cool_thirteen_trends$model==i, "AIC"]<-as.numeric(round(sum$AICtab[[1]], 1))
}

# Add column with AIC weights
cool_thirteen_trends$dAIC<- as.numeric(round(aicw(cool_thirteen_trends$AIC)$delta, 1))
cool_thirteen_trends$AIC.weights<- as.numeric(signif(aicw(cool_thirteen_trends$AIC)$w, 3))




# 3 Analyse the performance of the GLMM -----------------------------------


# let's check for overdispersion with the function from markus huff, calculating the sum of squared 
# Pearson residuals and comparing it to the residual degrees of freedom.
# Model is regarded as overdispersed when p <= 0.05. 
# https://rdrr.io/github/markushuff/PsychHelperFunctions/man/overdisp_fun.html

overdisp_fun(modelsw4[[which(warm_thirteen_trends$dAIC==0)]])   

overdisp_fun(modelsc4[[which(cool_thirteen_trends$dAIC==0)]])   

# summary of the final model
warm_mod <- summary(modelsw4[[which(warm_thirteen_trends$dAIC==0)]])

cool_mod <- summary(modelsc4[[which(cool_thirteen_trends$dAIC==0)]])





# 4 Compare palaeoclimate interaction models to traditional ones ----------


# traditional palao-analysis considers only the change from one bin to the next one, which is the 
# short-term change in our models
warm_mod_trad <- glmer("extinct~warming +(1|Genus)", data=thirteen_trends_bivalves, family="binomial")

cool_mod_trad <- glmer("extinct~cooling +(1|Genus)", data=thirteen_trends_bivalves, family="binomial")

# compare model performance by means of AIC values
summary(warm_mod_trad)$AICtab # traditional models considering short-term change only
warm_mod$AICtab # models including a long-term perspective (short-term and long-term)

summary(cool_mod_trad)$AICtab # short-term change only
cool_mod$AICtab # short-term and long-term




# 5 Extract the results from the GLMM  ------------------------------------




# 5.1 For Warming ---------------------------------------------------------


# extract the results from the GLMM using the "visreg" function
vis_out_w4 <- visreg(modelsw4[[which(warm_thirteen_trends$dAIC==0)]], "warming", scale="response", 
                     by=paste(warm_thirteen_trends[which(warm_thirteen_trends$dAIC==0),1]), 
                     rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5), plot = F) 


# take the values of a visreg and split them into cooling and warming and
vis_out_w4$fit$ctrend <- ifelse(vis_out_w4$fit$trend3<0, "Trend = Cooling", "Trend = Warming")


# compare the exctinction response of warming-cooling and warming-warming by means of boxplots
ggplot(vis_out_w4$fit, aes(x=ctrend, y=visregFit, fill=ctrend)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#56B4E9", "coral2")) +
  labs(title="Short-Term Warming", x="", y = "P(Extinction)")+
  theme_classic() +
  theme(legend.position = "none", axis.text=element_text(size=13),
        axis.title=element_text(size=14), plot.title = element_text(size=16))


# test for normally distributed data
shapiro.test(vis_out_w4$fit$visregFit[vis_out_w4$fit$ctrend=="Trend = Warming"]) # not normally distributed
shapiro.test(vis_out_w4$fit$visregFit[vis_out_w4$fit$ctrend=="Trend = Cooling"]) # not normally distributed

# both are not normally distributed, so we are using a Mann-Whitney U test
wilcox.test(vis_out_w4$fit$visregFit[vis_out_w4$fit$ctrend=="Trend = Warming"],
            vis_out_w4$fit$visregFit[vis_out_w4$fit$ctrend=="Trend = Cooling"],
            paired = F,  conf.int = T)




# 5.2 For Cooling ---------------------------------------------------------


# extract the results from the GLMM using the "visreg" function
vis_out_c4 <- visreg(modelsc4[[which(cool_thirteen_trends$dAIC== 0)]], "cooling", scale="response", 
                     by=paste(cool_thirteen_trends[which(cool_thirteen_trends$dAIC==0),1]), 
                     rug=2, strip.names=F, gg=T, line=list(col="#65a3a4", size=1.5), plot = F) 


# take the values of a visreg and split them into cooling and warming and
vis_out_c4$fit$ctrend <- ifelse(vis_out_c4$fit$trend.st3<0,"Trend = Cooling", "Trend = Warming")


# compare the exctinction response of warming-cooling and warming-warming by means of boxplots
ggplot(vis_out_c4$fit, aes(x=ctrend, y=visregFit, fill=ctrend)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#56B4E9", "coral2")) +
  labs(title="Short-Term Cooling", x="", y = "P(Extinction)")+
  theme_classic() +
  theme(legend.position = "none", axis.text=element_text(size=13),
        axis.title=element_text(size=14), plot.title = element_text(size=16))

# test for normally distributed data
shapiro.test(vis_out_c4$fit$visregFit[vis_out_c4$fit$ctrend=="Trend = Cooling"]) # not normally distributed
shapiro.test(vis_out_c4$fit$visregFit[vis_out_c4$fit$ctrend=="Trend = Warming"]) # not normally distributed


# both are not normally distributed, so we are using a Mann-Whitney U test
wilcox.test(vis_out_c4$fit$visregFit[vis_out_c4$fit$ctrend=="Trend = Cooling"],
            vis_out_c4$fit$visregFit[vis_out_c4$fit$ctrend=="Trend = Warming"], 
            paired = F, conf.int = T)




