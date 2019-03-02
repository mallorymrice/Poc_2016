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
library(reshape)
library(rminer)
library(plotrix)
library(tidyverse)

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

# opening Zoox csv file
Zoox <- read.csv("Data/MMR_2016_Symbiodinium_counts.csv", encoding = 'UTF-8', na.strings=c("", "NA"))

# opening Surface Area csv file
SurfaceArea <- read.csv("Data/MMR_2016_Surface_area.csv", encoding = 'UTF-8', na.strings=c("", "NA"))

# opening Replicate Treatment csv file
Reps <- read.csv("Data/MMR_2016_Replicate_treatment_assignments.csv", encoding = 'UTF-8', na.strings=c("", "NA"))

##############################################################################################################################
##### CALCULATING ZOOX DENSITY
##############################################################################################################################

glimpse(Zoox)
# converting ID to factor
Zoox$ID <- as.factor(Zoox$ID)

# need to change # zoox column from integar to numeric
Zoox$No.Zoox
(Zoox$No.Zoox <- as.numeric(Zoox$No.Zoox)) # conversion successful, values didn't change

# need to convert # zoox to million zoox
# steps:
# 1. calculate the average # of zoox counted by dividing the # zoox by # of squares counted (counted 5 squares total)
# 2. calculate the measured cell density (million cells per mL) by multiplying the average # zoox by the dilution volume (10 mL)
# and dividing by the volume of the small square (mL)
# step 1: calculate average # of zoox
Zoox$No.Zoox
(Zoox$average.zoox <- Zoox$No.Zoox / 5) # counted 5 total squares
# step 2: calculate the measured cell density
# smaller square is 0.004 mm^3 so 0.004 mm^3 * 1 mL/ 1000 mm^3 = 0.000004
# hemocytometer used: https://www.labplanet.com/propper-counting-chamber-hemacytometer-090001.html 
# volume of small square:
square.volume <- 0.004 * 1 / 1000
(Zoox$Million.zoox <- ( Zoox$average.zoox * 10 ) / square.volume)
# pipetted 10 uL from 10 mL sample --> reason for 10 being dilution volume
# source for calculations: https://www.hemocytometer.org/hemocytometer-calculation/ 

##############################################################################################################################
##### EXPLORING DATA
##############################################################################################################################

glimpse(Zoox)

# need to summarize the replicate counts
Zoox.summary <- ddply(Zoox, .(ID), summarise,
                      'Million.zoox'=mean(Million.zoox))

Zoox.summary

glimpse(SurfaceArea)

# converting ID to factor
SurfaceArea$ID <- as.factor(SurfaceArea$ID)

# merging Zoox and Surface Area dfs
Zoox.SA <- merge.data.frame(Zoox.summary, SurfaceArea, by = "ID")
glimpse(Zoox.SA)

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

# merging df's
merge <- merge.data.frame(Reps, Zoox.SA, by = "ID")

glimpse(merge)

# need to normalize zoox count by coral surface area
# looking at distribution of data
summary(merge)
# doesn't look like there's any wonky data

# now normalizing zoox counts (million cell) by surface area
# to do this divide zoox counts by surface area
merge$zoox.density <- merge$Million.zoox / merge$Surface_area_cm2
merge$zoox.density
# units: cells * 10^5 per cm^2
glimpse(merge)

range(merge$zoox.density)

##############################################################################################################################
##### PLOTTING DATA
##############################################################################################################################

Summary <- ddply(merge, .(Nutrient, Temperature, Scarred), summarise,
                 'mean'=mean(zoox.density, na.rm = TRUE),
                 'se'=se(zoox.density),
                 'N'=length(zoox.density))
Summary # changing values to scientific notation
max(Summary$mean)
Summary$mean <- Summary$mean / 10^5
Summary$se <- Summary$se / 10^5
merge$zoox.density2 <- merge$zoox.density / 10^5

zoox <- ggplot(Summary, aes(x = Temperature, y = mean, colour = Scarred, fill = Scarred))+
  facet_wrap(~Nutrient)+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width=0, position=position_dodge(width = 0.9))+
  geom_point(data = merge, aes(x = Temperature, y = zoox.density2, colour = Scarred), shape = 21, alpha = 0.7,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height=0.2, dodge.width=0.9))+
  scale_color_manual(values = c('black', 'black'), guide=FALSE)+
  scale_fill_manual(values = c('gray73', 'white'))+
  xlab(expression("Temperature ("*degree*C* ")"))+
  ylab(expression(italic(Symbiodinium)~"density"~{"(cells x 10"}^5*~{cm}^-2*")"~""))+
  theme_few(base_size = 12)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black"))
zoox

##############################################################################################################################
##### ANALYZING DATA: MIXED EFFECTS MODEL
##############################################################################################################################

# looking at distribution of zoox density
hist(merge$zoox.density)
qqnorm(merge$zoox.density)
qqline(merge$zoox.density)
shapiro.test(merge$zoox.density) # data not normally distributed
boxplot(merge$zoox.density ~ merge$Treatment)

glimpse(merge)

library(lme4)
full.model <- lmer(zoox.density ~ Temperature * Nutrient * Scarred + (1|Tank) + (1|Colony), data = merge, REML = TRUE)
tank.model <- lmer(zoox.density ~ Temperature * Nutrient * Scarred + (1|Tank), data = merge, REML = TRUE)
colony.model <- lmer(zoox.density ~ Temperature * Nutrient * Scarred + (1|Colony), data = merge, REML = TRUE)

normality.plots(full.model) # looks mostly normal
var.plots(full.model) # residuals look OK

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
# colony is significant (p = 0.001367)
# tank is not significant (p = 0.135810)

# corrected AIC values for all models
library(MuMIn)
AICc(tank.model, colony.model, full.model)
# enough evidence that tank can be dropped

# resulting final model:
final.model <- lmer(zoox.density2 ~ Temperature * Nutrient * Scarred + (1|Colony), data = merge, REML = TRUE)
# looking at residuals for normality and constant variance
normality.plots(final.model) # residuals look fairly normally distributed
var.plots(final.model) # residuals look like they're fairly spread out about the line
# assumptions met

# f-test for fixed effects
anova(final.model, type=3, ddf = "Kenward-Roger")
# significant effect of temperature (p = 0.005987)
# singificant effect of nutrients (p = 0.013047)

# post-hoc tests
library(emmeans)
emmeans(final.model, list(pairwise ~ Temperature), adjust = "tukey")
emmeans(final.model, list(pairwise ~ Nutrient), adjust = "tukey")
emmeans(final.model, list(pairwise ~ Temperature * Nutrient), adjust = "tukey")

##############################################################################################################################
##### MANUSCRIPT PLOT
##############################################################################################################################

ggsave("zoox.png", zoox, 
       path = "Output/",
       width = 6, height = 4, units = "in")

##############################################################################################################################
##### CALCULATING PERCENT CHANGES
##############################################################################################################################

ddply(merge, .(Temperature), summarise,
                 'mean'=mean(zoox.density2, na.rm = TRUE),
                 'se'=se(zoox.density2),
                 'N'=length(zoox.density2))
(4.624448 - 3.810846) / 3.810846

ddply(merge, .(Nutrient), summarise,
      'mean'=mean(zoox.density2, na.rm = TRUE),
      'se'=se(zoox.density2),
      'N'=length(zoox.density2))
(4.815775 - 3.691771) / 3.691771 # ammonium vs. control
(4.165281 - 3.691771) / 3.691771 # nitrate vs. control
(4.815775 - 4.165281) / 4.815775 # ammonium vs. nitrate

ddply(merge, .(Scarred), summarise,
      'mean'=mean(zoox.density2, na.rm = TRUE),
      'se'=se(zoox.density2),
      'N'=length(zoox.density2))

Summary
(3.677672 - 3.388319) / 3.388319 # intact at 26 vs. intact at 29
