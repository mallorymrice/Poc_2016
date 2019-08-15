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

# opening up nubbin size csv files
size <- read.csv("Data/MMR_2016_Nubbin_initial_size.csv", encoding = 'UTF-8')

# opening Replicate Treatment csv file
Reps <- read.csv("Data/MMR_2016_Replicate_treatment_assignments.csv", encoding = 'UTF-8', na.strings=c("", "NA"))

##############################################################################################################################
##### RESTRUCTURING DF
##############################################################################################################################

glimpse(size)

# convering ID to factor
size$ID <- as.factor(size$ID)
glimpse(size)

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
size <- merge.data.frame(Reps, size, by = "ID")

glimpse(size)

##############################################################################################################################
##### PLOTTING NUBBIN SIZE BY TREATMENT
##############################################################################################################################

# summarizing data to plot
summary <- ddply(size, .(Treatment), summarise,
                        'mean'=mean(Initial_nubbin_height, na.rm = TRUE),
                        'se'=se(Initial_nubbin_height))

dodge<-position_dodge(width=0.6) # this offsets the points so they don't overlap

ggplot(summary, aes(x = Treatment, y = mean))+
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se), 
                width=0, position=position_dodge(width = 0.9))+
  ylab("Initial nubbin height (cm)")+
  theme_few(base_size = 12)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black"))

####################################################################################################
##### DATA ANALYSIS: ONE-WAY ANOVA
####################################################################################################

par(mfrow=c(1,2))
hist(size$Initial_nubbin_height)
qqnorm(size$Initial_nubbin_height)
abline(0,1)
shapiro.test(size$Initial_nubbin_height) # statistically not normally distributed
leveneTest(size$Initial_nubbin_height ~ size$Treatment) # equal variance
# using a nonparametric one-way ANOVA

# boxplots
ggplot(size, aes(x = Treatment, y = Initial_nubbin_height))+
  geom_boxplot()+
  ylab("Initial nubbin height (cm)")+
  theme_few(base_size = 12)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black"))

# Kruskal-Wallis test
kruskal.test(Initial_nubbin_height ~ Treatment, data = size)
# no differences in initial nubbin height across treatments