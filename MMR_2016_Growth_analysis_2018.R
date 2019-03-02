rm(list = ls())
library(ggplot2)
library(ggthemes)
library(plyr)
library(Hmisc)
library(outliers)
library(datasets)
library(lattice)
library(reshape)
library(reshape2)
library(tidyr)
library(nlme)
library(stats)
library(data.table)
library(lawstat)
library(car)
library(dplyr)
library(base)
library(rminer)
library(tidyverse)
library(lme4)
library(lsmeans)
library(lmerTest)

# standard error function
se <- function(x) sd(x)/sqrt(length(x))

# function to visualize the boxplot, histogram, QQ-plot, and kernel density estimate plot for residuals
normality.plots <- function(x) {
  par(mfrow=c(2,2))
  hist(residuals(x), main = "Histogram", xlab = "Values")
  boxplot(residuals(x), main = "Boxplot", ylab = "Values")
  qqnorm(residuals(x))
  qqline(residuals(x))
  plot(density(residuals(x)), main = "Kernel Density Estimate")
}

# function to visaulize the fitted vs. residuals of model and the absolute residuals vs. fitted
var.plots <- function(x) {
  par(mfrow=c(1,2))
  plot(x = fitted(x), y = residuals(x), xlab="Fitted", ylab="Residuals")
  abline(h=0)
  title("Residuals vs. Fitted")
  plot(x = fitted(x), y = abs(residuals(x)), xlab="Fitted", ylab="Absolute Residuals")
  abline(h=0)
  title("Absolute Residuals vs. Fitted")
}

# completeFun removes NA values in a specific column
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

##############################################################################################################################
##### OPENING DF's
##############################################################################################################################

# opening Weight csv file
Weight <- read.csv("Data/MMR_2016_Weight.csv", encoding = 'UTF-8', na.strings=c("", "NA"))

# opening Surface Area csv file
SurfaceArea <- read.csv("Data/MMR_2016_Surface_area.csv", encoding = 'UTF-8', na.strings=c("", "NA"))

# opening Replicate Treatment csv file
Reps <- read.csv("Data/MMR_2016_Replicate_treatment_assignments.csv", encoding = 'UTF-8', na.strings=c("", "NA"))

##############################################################################################################################
##### RESTRUCTURING  AND MERGING DF'S
##############################################################################################################################

glimpse(Weight)

# converting ID to factor
Weight$ID <- as.factor(Weight$ID)

# converting Diff_weight_g to numeric
Weight$Diff_weight_g <- as.numeric((as.character(Weight$Diff_weight_g)))

# converting Diff_weight_mg to numeric
Weight$Diff_weight_mg <- as.numeric((as.character(Weight$Diff_weight_mg)))

glimpse(Reps)

# converting ID to factor
Reps$ID <- as.factor(Reps$ID)

# converting Temperature to factor
Reps$Temperature <- as.factor(Reps$Temperature)

# converting Tank to factor
Reps$Tank <- as.factor(Reps$Tank)

# reorganiziing Nutrient factor levels
Reps$Nutrient <- ordered(Reps$Nutrient, levels = c("Control", "Ammonium", "Nitrate"))

glimpse(Reps)

# converting Nutrient to factor
Reps$Nutrient <- as.factor(Reps$Nutrient)

# check for NA values
is.na(Weight$Diff_weight_mg) 

# removing NA values
Weight <- completeFun(Weight, "Diff_weight_mg")

# check for NA values
is.na(Weight$Diff_weight_mg) # no NA values

# looking at SurfaceArea df
glimpse(SurfaceArea)

# converting ID to factor
SurfaceArea$ID <- as.factor(SurfaceArea$ID)

# merging Weight and SurfaceArea df's by ID
Weight <- merge.data.frame(Weight, SurfaceArea, by = 'ID')

# merging Weight and Reps df's by ID
Weight <- merge.data.frame(Reps, Weight, by = 'ID')

glimpse(Weight)

# standardizing the difference in weight (mg) by surface area (cm^2)
Weight$Stand_diff <- Weight$Diff_weight_mg / Weight$Surface_area_cm2
# units now mg/cm^2

# standardizing mg/cm^2 by the experiment duration (25 days)
Weight$Stand_diff <- Weight$Stand_diff / 25
# units now mg cm^-2 day^-1

##############################################################################################################################
##### DATA EXPLORATION
##############################################################################################################################

# box plot of growth by tank
ggplot(Weight, aes(x = Tank, y = Stand_diff))+
  geom_boxplot()+
  ylab(expression(paste("Growth (mg ", cm^"-2 ", day^"-1", ")")))+
  theme_few() 
# outlier from tank 1

# box plot of growth by colony
ggplot(Weight, aes(x = Colony, y = Stand_diff))+
  geom_boxplot()+
  ylab(expression(paste("Growth (mg ", cm^"-2 ", day^"-1", ")")))+
  theme_few() 
# looks like there is an outlier in colony 1 - same outlier as in prevous boxplot

# box plot across treatment
ggplot(Weight, aes(x = Treatment, y = Stand_diff))+
  geom_boxplot()+
  ylab(expression(paste("Growth (mg ", cm^"-2 ", day^"-1", ")")))+
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# outlier is in nitrate, 29, intact

##############################################################################################################################
##### PLOTTING MEAN +/- SE WITH OUTLIER
##############################################################################################################################

# summarizing data to plot
WeightSummary <- ddply(Weight, .(Nutrient, Temperature, Scarred), summarise,
                       'mean'=mean(Stand_diff, na.rm = TRUE),
                       'se'=se(Stand_diff),
                       'N'=length(Stand_diff))

dodge<-position_dodge(width=0.6) # this offsets the points so they don't overlap

ggplot(WeightSummary, aes(x = Temperature, y = mean, colour = Scarred, fill = Scarred))+
  facet_wrap(~Nutrient)+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se), 
                width=0, position=position_dodge(width = 0.9))+
  geom_point(data = Weight, aes(x = Temperature, y = Stand_diff, colour = Scarred), shape = 21, alpha = 0.7,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height=0.2, dodge.width=0.9))+
  scale_color_manual(values = c('black', 'black'), guide=FALSE)+
  scale_fill_manual(values = c('gray73', 'white'), labels=c("Intact", "Wounded"))+
  xlab(expression("Temperature ("*degree*C* ")"))+
  ylab(expression(paste("Growth (mg ", cm^"-2 ", day^"-1", ")")))+
  theme_few(base_size = 12)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black"))

##############################################################################################################################
##### DATA ANALYSIS: MIXED EFFECTS MODEL WITH OUTLIER
##############################################################################################################################

# first need to look at spread of the data
hist(Weight$Stand_diff)
qqnorm(Weight$Stand_diff)
# data does not look that normal - likely due to outlier
# moving forward with linear mixed-effects model
ggplot(Weight, aes(x = Temperature, y = Stand_diff, colour=Tank))+
  facet_wrap(Scarred~Nutrient)+
  geom_point()+
  xlab(expression("Temperature ("*degree*C* ")"))+
  ylab("Growth rate (g per cm^2 per day)")+
  theme_classic()

ggplot(Weight, aes(x = Temperature, y = Stand_diff, colour=Tank))+
  facet_wrap(Scarred~Nutrient)+
  xlab(expression("Temperature ("*degree*C* ")"))+
  ylab("Growth rate (g per cm^2 per day)")+
  geom_boxplot()+
  theme_classic()

# mixed effects model to test effects of Temperature and Nutrient treatments on Healing Rate
# Temperature and Nutrient are fixed effects
# Tank is a random effect
# Colony is a random effect

# Treatments:
# Nutrient (Control, Ammonium, Nitrate) - fixed effect
# Temperature (26 C, 29 C) - fixed effect
# need to include Tank as a random effect: + (1 | Tank)
# need to include Colony as a random effect: + (1 | Tank)
library(lme4) # lme4 package allows for random effects

# first looking at distribution of response variable (healing rate)
par(mfrow=c(1,2))
hist(Weight$Stand_diff)
qqnorm(Weight$Stand_diff)
abline(0,1)
shapiro.test(Weight$Stand_diff) 
# weight does not look normal - definitely an extreme outlier
max(Weight$Stand_diff) # this is the outlier value
# no coral grows 7 mg per cm^2 per day - this is computational error
# will perform analysis with mixed effects model with and without outlier

# interaction plot
par(mfrow=c(1,3))
interaction.plot(Weight$Temperature, Weight$Nutrient, Weight$Stand_diff)
interaction.plot(Weight$Temperature, Weight$Scarred, Weight$Stand_diff)
interaction.plot(Weight$Scarred, Weight$Nutrient, Weight$Stand_diff)

# full model with all random effects and the interaction between temperature and nutrient
full.model <- lmer(Stand_diff ~ Temperature * Nutrient * Scarred + (1|Tank) + (1|Colony), data = Weight, REML = TRUE)
tank.model <- lmer(Stand_diff ~ Temperature * Nutrient * Scarred + (1|Tank), data = Weight, REML = TRUE)
colony.model <- lmer(Stand_diff ~ Temperature * Nutrient * Scarred + (1|Colony), data = Weight, REML = TRUE)

library(lattice)
dotplot(ranef(full.model, condVar = TRUE))
# shows the variance of growth across random effects
# doesn't look like there's an effect of colony

# need to check assumptions for random effects:
# 1) residuals are independently distributed
# 2) residuals are from a normal distribution
par(mfrow=c(1,2))
qqnorm(ranef(full.model)$Colony[,1])
qqline(ranef(full.model)$Colony[,1])
qqnorm(ranef(full.model)$Tank[,1])
qqline(ranef(full.model)$Tank[,1])
# the residuals look relatively normal

# first need to determine random effect structure
library(lmerTest)
ranova(full.model) # looking at the LRT for the random effects
# colony is not significant (p = 0.5772)
# tank is not significant (p = 1) 

# corrected AIC values for all models
library(MuMIn)
AICc(tank.model, colony.model, full.model)
# enough evidence that tank can be dropped

# resulting final model:
final.model <- lmer(Stand_diff ~ Nutrient * Temperature * Scarred + (1|Colony), data = Weight, REML = TRUE)
# looking at residuals for normality and constant variance
normality.plots(final.model) # residuals look fairly normally distributed
var.plots(final.model) # residuals look like they're fairly spread out about the line
# assumptions not met - outlier

# f-test for fixed effects
anova(final.model, type=3, ddf = "Kenward-Roger")
# no significant main or interactive effects

##############################################################################################################################
##### DATA ANALYSIS: MIXED EFFECTS MODEL WITHOUT OUTLIER
##############################################################################################################################

# need to remove the row with the outlier
sort(Weight$Stand_diff)
Weight.removed <- Weight %>% subset(Stand_diff < 6)
sort(Weight.removed$Stand_diff) # outlier successfuly removed
range(Weight.removed$Stand_diff)

ggplot(Weight.removed, aes(x = Temperature, y = Stand_diff, colour=Tank))+
  facet_wrap(Scarred~Nutrient)+
  geom_point()+
  xlab(expression("Temperature ("*degree*C* ")"))+
  ylab("Growth rate (g per cm^2 per day)")+
  theme_classic()

ggplot(Weight.removed, aes(x = Temperature, y = Stand_diff, colour=Tank))+
  facet_wrap(Scarred~Nutrient)+
  xlab(expression("Temperature ("*degree*C* ")"))+
  ylab("Growth rate (g per cm^2 per day)")+
  geom_boxplot()+
  theme_classic()

# mixed effects model to test effects of Temperature and Nutrient treatments on Healing Rate
# Temperature and Nutrient are fixed effects
# Tank is a random effect
# Colony is a random effect

# Treatments:
# Nutrient (Control, Ammonium, Nitrate) - fixed effect
# Temperature (26 C, 29 C) - fixed effect
# need to include Tank as a random effect: + (1 | Tank)
# need to include Colony as a random effect: + (1 | Tank)

# first looking at distribution of response variable (healing rate)
par(mfrow=c(1,2))
hist(Weight.removed$Stand_diff)
qqnorm(Weight.removed$Stand_diff)
abline(0,1)
shapiro.test(Weight.removed$Stand_diff) 
# weight does look a lot more normal after removing the outlier

# interaction plot
par(mfrow=c(1,3))
interaction.plot(Weight.removed$Temperature, Weight.removed$Nutrient, Weight.removed$Stand_diff)
interaction.plot(Weight.removed$Temperature, Weight.removed$Scarred, Weight.removed$Stand_diff)
interaction.plot(Weight.removed$Scarred, Weight.removed$Nutrient, Weight.removed$Stand_diff)
# doesn't really look like there's any interactions going on

# full model with all random effects and the interaction between temperature and nutrient
full.model.removed <- lmer(Stand_diff ~ Temperature * Nutrient * Scarred + (1|Tank) + (1|Colony), data = Weight.removed, REML = TRUE)
tank.model.removed <- lmer(Stand_diff ~ Temperature * Nutrient * Scarred + (1|Tank), data = Weight.removed, REML = TRUE)
colony.model.removed <- lmer(Stand_diff ~ Temperature * Nutrient * Scarred + (1|Colony), data = Weight.removed, REML = TRUE)

dotplot(ranef(full.model.removed, condVar = TRUE))

# need to check assumptions for random effects:
# 1) residuals are independently distributed
# 2) residuals are from a normal distribution
par(mfrow=c(1,2))
qqnorm(ranef(full.model.removed)$Colony[,1])
qqline(ranef(full.model.removed)$Colony[,1])
qqnorm(ranef(full.model.removed)$Tank[,1])
qqline(ranef(full.model.removed)$Tank[,1])
# the residuals look relatively normal

# first need to determine random effect structure
ranova(full.model.removed) # looking at the LRT for the random effects
# colony is significant (p = 0.006279)
# tank is not significant (p = 0.8167) 

# corrected AIC values for all models
AICc(tank.model.removed, colony.model.removed, full.model.removed)
# enough evidence that tank can be dropped

# resulting final model:
final.model.removed <- lmer(Stand_diff ~ Temperature * Nutrient * Scarred + (1|Colony), data = Weight.removed, REML = TRUE)
# looking at residuals for normality and constant variance
normality.plots(final.model.removed) # residuals look fairly normally distributed
var.plots(final.model.removed) # residuals look like they're fairly spread out about the line
# assumptions met 

# f-test for fixed effects
anova(final.model.removed, type=3, ddf = "Kenward-Roger")
# no significant main or interactive effects

##############################################################################################################################
##### PLOTTING MEAN +/- SE WITHOUT OUTLIER
##############################################################################################################################

# summarizing data to plot
WeightSummaryRemoved <- ddply(Weight.removed, .(Nutrient, Temperature, Scarred), summarise,
                       'mean'=mean(Stand_diff, na.rm = TRUE),
                       'se'=se(Stand_diff),
                       'N'=length(Stand_diff))

dodge<-position_dodge(width=0.6) # this offsets the points so they don't overlap

names(Weight.removed)

growth <- ggplot(WeightSummaryRemoved, aes(x = Temperature, y = mean, colour = Scarred, fill = Scarred))+
  facet_wrap(~Nutrient)+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se), 
                    width=0, position=position_dodge(width = 0.9))+
  geom_point(data = Weight.removed, aes(x = Temperature, y = Stand_diff, colour = Scarred), shape = 21, alpha = 0.7,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height=0.2, dodge.width=0.9))+
  scale_color_manual(values = c('black', 'black'), guide=FALSE)+
  scale_fill_manual(values = c('gray73', 'white'), labels=c("Intact", "Wounded"))+
  xlab(expression("Temperature ("*degree*C* ")"))+
  ylab(expression(paste("Growth (mg ", cm^"-2 ", day^"-1", ")")))+
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5))+
  theme_few(base_size = 12)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black"))+
  theme(legend.position = c(0.1, 0.89))+
  theme(legend.title = element_blank())
growth

##############################################################################################################################
##### MANUSCRIPT PLOT
##############################################################################################################################

ggsave("growth.png", growth, 
       path = "Output/",
       width = 6, height = 4, units = "in")