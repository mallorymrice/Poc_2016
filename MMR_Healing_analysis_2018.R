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

##############################################################################################################################
##### OPENING DF's
##############################################################################################################################

# opening up Healing csv file
Healing <- read.csv(file="MMR_Healing.csv", encoding = 'UTF-8')

##############################################################################################################################
##### RESTRUCTURING DF
##############################################################################################################################

glimpse(Healing)
# converting Temperature to factor
Healing$Temperature <- as.factor(Healing$Temperature)

# converting Tank to factor
Healing$Tank <- as.factor(Healing$Tank)

# reorganizing Nutrient factor levels
Healing$Nutrient <- ordered(Healing$Nutrient, levels = c("Control", "Ammonium", "Nitrate"))

##############################################################################################################################
##### PLOTTING HEALING RATE STANDARDIZED BY DAY (mm^2 day^-1) 
##############################################################################################################################

# summarizing data to plot
HealingSummary <- ddply(Healing, .(Nutrient, Temperature, Corallivory), summarise,
                        'mean'=mean(HealingRate),
                        'se'=se(HealingRate))

dodge<-position_dodge(width=0.6) # this offsets the points so they don't overlap

# plotting Healing Rate by treatment
ggplot(HealingSummary, aes(x = Temperature, y = mean, colour = Corallivory, fill = Corallivory))+
  facet_wrap(~Nutrient)+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se), 
                width=0, position=position_dodge(width = 0.9))+
  geom_point(data = Healing, aes(x = Temperature, y = HealingRate, colour = Corallivory), shape = 21, alpha = 0.7,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height=0.2, dodge.width=0.9))+
  scale_color_manual(values = c('black'), guide=FALSE)+
  scale_fill_manual(values = c('gray73'), guide = FALSE)+
  xlab(expression("Temperature ("*degree*C* ")"))+
  ylab(bquote('Healing Rate ('*mm^2*' '*day^-1*')'))+
 # scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2.5, 0.5))+
  theme_few(base_size = 12)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black"))

##############################################################################################################################
##### INTERACTION PLOT OF HEALING RATE STANDARDIZED BY DAY (mm^2 day^-1) 
##############################################################################################################################

dodge<-position_dodge(width=0.6) # this offsets the points so they don't overlap

# plotting Healing Rate by treatment - BW plot
ggplot(HealingSummary, aes(x = Temperature, y = mean, colour=Nutrient))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width=0, position=position_dodge(width = 0.7))+
  geom_point(position = position_dodge(width = 0.7), size = 4)+
  xlab(expression("Temperature ("*degree*C* ")"))+
  ylab(bquote('Healing Rate ('*mm^2*' '*day^-1*')'))+
  scale_y_continuous(limits = c(0, 2), breaks = seq(0,2,0.5))+
  theme_few(base_size = 12)

##############################################################################################################################
##### PLOTTING HEALING PERCENT (%)
##############################################################################################################################

# summarizing data to plot
PercentHealedSummary <- ddply(Healing, .(Nutrient, Temperature, Corallivory), summarise,
                        'mean'=mean(PercentHealed),
                        'se'=se(PercentHealed))
dodge<-position_dodge(width=0.6) # this offsets the points so they don't overlap

# plotting Healing Rate by treatment - BW plot
ggplot(PercentHealedSummary, aes(x = Temperature, y = mean, fill=Corallivory))+
  facet_wrap(~Nutrient)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width=0, position=position_dodge(width = 0.7))+
  scale_shape_manual(values = c(21))+
  scale_fill_manual(values = c('white'))+
  geom_point(aes(shape=Corallivory), position = position_dodge(width = 0.7), size = 4)+
  xlab(expression("Temperature ("*degree*C* ")"))+
  ylab("Percent Healed (%)")+
  theme_few(base_size = 12)+
  theme(legend.position="none")

####################################################################################################
##### PLOT OF DIFF IN SCAR AREA
####################################################################################################

Healing$Difference <- Healing$Initial_scar_area - Healing$Final_scar_area

# summarizing data to plot
DifferenceSummary <- ddply(Healing, .(Nutrient, Temperature, Corallivory), summarise,
                              'mean'=mean(Difference),
                              'se'=se(Difference))
dodge<-position_dodge(width=0.6) # this offsets the points so they don't overlap

# plotting Healing Rate by treatment - BW plot
ggplot(DifferenceSummary, aes(x = Temperature, y = mean, fill=Corallivory))+
  facet_wrap(~Nutrient)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width=0, position=position_dodge(width = 0.7))+
  scale_shape_manual(values = c(21))+
  scale_fill_manual(values = c('white'))+
  geom_point(aes(shape=Corallivory), position = position_dodge(width = 0.7), size = 4)+
  xlab(expression("Temperature ("*degree*C* ")"))+
  ylab("Difference in scar area (mm^2)")+
  theme_few(base_size = 12)+
  theme(legend.position="none")

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
# determining if the random effect of tank and colony are significant in the model
anova(tank.model, full.model) # looking at the random effect of colony
# colony significant in the model (p = 0.02364)
anova(colony.model, full.model) # looking at hte random effect of tank
# tank not significant in the model (p = 0.5807)
library(lmerTest)
ranova(tank.model) # tank is not significant (p = 0.1835)
ranova(colony.model) # colony is significant (p = 0.03647)
ranova(full.model) # looking at the LRT for the random effects
# colony is significant (p = 0.02159)
# tank is not significant (p = 0.10213) but this p-value is a little iffy to me

# AIC values for all models
AIC(tank.model, colony.model, full.model)
# AIC value for tank model is ~3 values higher than that of colony model
# AIC value for full model (tank and colony as random effects) is similar to that of the colony model
# AIC values suggest that tank is not important to include
library(MuMIn)
AICc(tank.model, colony.model, full.model)
# same pattern for corrected AIC values
# enough evidence that tank can be dropped

# resulting final model:
# Healing Rate ~ Temperature + Nutrient + Temperature*Nutrient + (1|Colony)
final.model <- lmer(HealingRate ~ Nutrient * Temperature + (1|Colony), data = Healing, REML = TRUE)
# looking at residuals for normality and constant variance
normality.plots(final.model) # residuals look fairly normally distributed
var.plots(final.model) # residuals look like they're fairly spread out about the line
# assumptions met

# f-test for fixed effects
anova(final.model, type=3, ddf = "Kenward-Roger")
anova(final.model, type=3, dff = "Satterthwaite")
# significant interaction

# loading lsmeans package to do post-hoc tests
library(lsmeans)
# determing what treatments are driving the significant interaction
LSMeans <- lsmeans(final.model, ~ Nutrient * Temperature, options = list(estName = "HealingRate"))
LSMeans
lsmeansLT(final.model)
difflsmeans(final.model)
pairs(LSMeans)
org.sum <- summary(LSMeans, infer = c(TRUE,TRUE), level = .95)
org.sum

# conditional R^2
library(MuMIn)
r.squaredGLMM(final.model)
# R2m is the R^2 for fixed effects
# R2c is the R^2 for random effects
# here the random effects explain ~17% of the variation in healing rate

# comparing this to the R^2 for the model with tank as random effect
r.squaredGLMM(full.model)
# when added as a random effect tank explains ~13% of the variance in healing rate

####################################################################################################
##### DATA ANALYSIS: SPLIT PLOT DESIGN
####################################################################################################

split.plot <- lmer(HealingRate ~ Nutrient * Temperature + (1|Colony) + (1|Tank/Corallivory), data = Healing, REML = TRUE)
