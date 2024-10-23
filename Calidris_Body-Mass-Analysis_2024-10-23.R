#-----------------------------#
# Calidris Body Mass Analysis #
#     Created 10/23/2024      #
#     Modified 10/23/2024     #
#-----------------------------#

# load packages
library(dplyr)
library(ggplot2)
library(AICcmodavg)
library(RColorBrewer)
library(trtools)
library(tidyverse)

?saveRDS

# ---------------------------------------------------------------------------- #

# Data Processing ####

# Read data
setwd("processed_data")
birds <- read.csv("LM_ShorebirdsALLNeg.csv")

# Subset data for Calidris sp.
calidris <- subset(birds, Species %in% c("LeastSandpiper", 
                                     "PectoralSandpiper", 
                                     "SemipalmatedSandpiper"))

# Remove birds that do not have neonicotinoid concentrations
calidris <- calidris[!is.na(calidris$OverallNeonic), ]

# Remove birds that do not have mass
calidris <- calidris[!is.na(calidris$Mass), ]

# Make neonicotinoid detection column
calidris$Detection <- ifelse(calidris$OverallNeonic > 0, "Detection", "Non-detection")

# Convert categorical variables to factor
calidris$Year <- as.factor(calidris$Year)
calidris$Season <- as.factor(calidris$Season)
calidris$Sex <- as.factor(calidris$Sex)
calidris$Detection <- as.factor(calidris$Detection)

# Check for variable correlation
cor(calidris$Julian, calidris$MigDate)

# ---------------------------------------------------------------------------- #

## Covariate Restrictions ####

# Date into Season and Julian can't be in the same model because of high correlation (0.861)
# Year and Season can't be in the same model because of lack of replication of Season within a year
# Julian and Season can't be in the same model because Season bins Julian
# Year and Species can't be in the same model because all LESA were captured in 2023

# ---------------------------------------------------------------------------- #

# Data Exploration & Linear Regression Model Assumptions ####

## Assumption 1: Zero Expectations of Errors (Residuals) ####

### Residuals Plots 

# Calculate residuals
m <- lm(Mass ~ OverallNeonic, data = calidris)
calidris$yhat <- predict(m)
calidris$residuals <- residuals(m)
calidris$rstudent <- rstudent(m)
head(calidris)

# Simple residual plot
plot(predict(m), rstudent(m))

# Pretty residual plot
ggplot(calidris, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# Bad megaphone pattern... 

# Log transformation of neonics

# Calculate residuals 
m.log <- lm(Mass ~ (log10(OverallNeonic+0.0001)), data = calidris)
calidris$yhat <- predict(m.log)
calidris$residuals <- residuals(m.log)
calidris$rstudent <- rstudent(m.log)
head(calidris)

# Simple residual plot
plot(predict(m.log), rstudent(m.log))

# Pretty residual plot
ggplot(calidris, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# Bad megaphone pattern - violates assumption

# ---------------------------------------------------------------------------- #

## Assumption 2: Equality of error variances (Homoscedasticity) #### 

# No transformation
m <- lm(Mass ~ OverallNeonic, data = calidris)
calidris$yhat <- predict(m)
calidris$residuals <- residuals(m)
calidris$rstudent <- rstudent(m)

# Plot homoscedasticity
ggplot(calidris, aes(x = yhat, y = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual")

# No violations because variability doesn't increase as predicted value increases

# Log transformation
m.log <- lm(Mass ~ (log10(OverallNeonic+0.0001)), data = calidris)
calidris$yhat <- predict(m.log)
calidris$residuals <- residuals(m.log)
calidris$rstudent <- rstudent(m.log)

### PLOT HOMOSCEDASTICITY
ggplot(calidris, aes(x = yhat, y = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual")

# No violations because variability doesn't increase as predicted value increases.
# Log transformation does appear to spread residuals more evenly. 

# ---------------------------------------------------------------------------- #

## Assumption 3: Independence of Errors/Observations ####

# Distribution of one error/observation does not depend on another error/observation

# ---------------------------------------------------------------------------- #

## Assumption 4: Normality ####

m <- lm(Mass ~ OverallNeonic, data = calidris)
hist(rstandard(m))
qqnorm(rstandard(m))

# Does not severely violate normality assumption
# Regression models do not many any assumptions about the distrbution of explanatory 
# variables (e.g., Concentrations), but their distribution does affect inferences. 

# ---------------------------------------------------------------------------- #