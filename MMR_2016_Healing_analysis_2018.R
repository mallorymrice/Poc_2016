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
se <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))

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

##############################################################################################################################
##### OPENING DF's
##############################################################################################################################

# opening up Healing csv files
Healing <- read.csv("Data/MMR_2016_Healing.csv", encoding = 'UTF-8')

# opening Replicate Treatment csv file
Reps <- read.csv("Data/MMR_2016_Replicate_treatment_assignments.csv", encoding = 'UTF-8', na.strings=c("", "NA"))

##############################################################################################################################
##### RESTRUCTURING DF
##############################################################################################################################

glimpse(Healing)

# calculating healing rate
# need to substract final scar area from initial scar area
# and then divide by the time period corals were allowed to heal (24 days)
# healing rate: mm^2 day^-1
Healing$HealingRate <- ( Healing$Initial_scar_area - Healing$Final_scar_area ) / 24
sort(Healing$HealingRate)

# convering ID to factor
Healing$ID <- as.factor(Healing$ID)

glimpse(Healing)

glimpse(Reps)

# convering ID to factor
Reps$ID <- as.factor(Reps$ID)

# converting Temperature to factor
Reps$Temperature <- as.factor(Reps$Temperature)

# converting Tank to factor
Reps$Tank <- as.factor(Reps$Tank)

# reorganizing Nutrient factor levels
Reps$Nutrient <- ordered(Reps$Nutrient, levels = c("Control", "Ammonium", "Nitrate"))

# merging data frames
Healing <- merge.data.frame(Reps, Healing, by = "ID")

glimpse(Healing)

##############################################################################################################################
##### PLOTTING HEALING RATE STANDARDIZED BY DAY (mm^2 day^-1) 
##############################################################################################################################

# summarizing data to plot
HealingSummary <- ddply(Healing, .(Nutrient, Temperature, Scarred), summarise,
                        'mean'=mean(HealingRate, na.rm = TRUE),
                        'se'=se(HealingRate))

dodge<-position_dodge(width=0.6) # this offsets the points so they don't overlap

# plotting Healing Rate by treatment
healing <- ggplot(HealingSummary, aes(x = Temperature, y = mean, colour = Scarred, fill = Scarred))+
  facet_wrap(~Nutrient)+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se), 
                width=0, position=position_dodge(width = 0.9))+
  geom_point(data = Healing, aes(x = Temperature, y = HealingRate, colour = Scarred), shape = 21, alpha = 0.7,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height=0.2, dodge.width=0.9))+
  scale_color_manual(values = c('black'), guide=FALSE)+
  scale_fill_manual(values = c('gray73'), guide = FALSE)+
  xlab(expression("Temperature ("*degree*C* ")"))+
  ylab(bquote('Healing Rate ('*mm^2*' '*day^-1*')'))+
 # scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2.5, 0.5))+
  theme_few(base_size = 12)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black"))
healing

####################################################################################################
##### DATA ANALYSIS: MIXED EFFECTS MODEL
####################################################################################################

# mixed effects model to test effects of Temperature and Nutrient treatments on Healing Rate
# Temperature and Nutrient are fixed effects
# Tank is a random effect
# Colony is a random effect

# Treatments:
# Nutrient (Control, Ammonium, Nitrate) - fixed effect
# Temperature (26 C, 29 C) - fixed effect
# need to include Tank as a random effect: + (1 | Tank)
# need to include Colony as a random effect: + (1 | Tank)

ggplot(Healing, aes(x = Temperature, y = HealingRate, colour=Tank))+
  facet_wrap(~Nutrient)+
  geom_boxplot()+
  theme_classic()

ggplot(Healing, aes(x = Temperature, y = HealingRate, colour=Colony))+
  facet_wrap(~Nutrient)+
  geom_point()+
  theme_classic()

# first looking at distribution of response variable (healing rate)
par(mfrow=c(1,2))
hist(Healing$HealingRate)
qqnorm(Healing$HealingRate)
abline(0,1)
shapiro.test(Healing$HealingRate) # statistically normally distributed
# healing rate looks relatively normal
# can perform a linear mixed effects model on data

# interaction plot
par(mfrow=c(1,1))
interaction.plot(Healing$Temperature, Healing$Nutrient, Healing$HealingRate)
# definitely looks like there is an interaction going on between nutrient and temperature treatments

library(lme4)
# full model with all random effects and the interaction between temperature and nutrient
full.model <- lmer(HealingRate ~ Temperature * Nutrient + (1|Tank) + (1|Colony), data = Healing, REML = TRUE)
tank.model <- lmer(HealingRate ~ Temperature * Nutrient + (1|Tank), data = Healing, REML = TRUE)
colony.model <- lmer(HealingRate ~ Temperature * Nutrient + (1|Colony), data = Healing, REML = TRUE)

library(lattice)
dotplot(ranef(full.model, condVar = TRUE))
# shows the variance of healing rate across random effects

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
# colony is significant (p = 0.00391)
# tank is not significant (p = 0.15168)

# corrected AIC values for all models
library(MuMIn)
AICc(tank.model, colony.model, full.model)
# enough evidence that tank can be dropped

# resulting final model:
final.model <- lmer(HealingRate ~ Nutrient * Temperature + (1|Colony), data = Healing, REML = TRUE)
# looking at residuals for normality and constant variance
normality.plots(full.model) # residuals look fairly normally distributed
var.plots(full.model) # residuals look like they're fairly spread out about the line
# assumptions met

# f-test for fixed effects
anova(final.model, type=3, ddf = "Kenward-Roger")
# significant interaction

# post-hoc tests
library(emmeans)
emmeans(final.model, list(pairwise ~ Temperature * Nutrient), adjust = "tukey")

##############################################################################################################################
##### MANUSCRIPT PLOT
##############################################################################################################################

ggsave("healing.pdf", healing, 
       path = "Output/",
       width = 3, height = 3, units = "in")

##############################################################################################################################
##### CALCULATING PERCENT CHANGES
##############################################################################################################################

HealingSummary
((1.6723313 - 0.5654271) / 1.6723313 ) * 100 # 26 vs. 29 for ambient # 66% decrease
((1.5448958 - 0.5654271) / 1.5448958 ) * 100 # 29 control vs. 29 ammonium
((1.5133920 - 0.5654271) / 1.5133920 ) * 100 # 29 control vs. 29 nitrate
