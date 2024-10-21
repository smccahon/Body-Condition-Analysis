#------------------------------------#
# Least Sandpiper Body Mass Analysis #
#       Created 10/21/2024           #
#       Modified 10/21/2024          #
#------------------------------------#

# load packages
library(dplyr)
library(ggplot2)
library(AICcmodavg)
library(RColorBrewer)
library(trtools)
library(tidyverse)

# ---------------------------------------------------------------------------- #

# Data Processing ####

# Read data
setwd("processed_data")
birds <- read.csv("LM_ShorebirdsALLNeg.csv")

# Subset data for Least Sandpiper
lesa <- subset(birds, Species %in% c("LeastSandpiper"))

# Remove birds that do not have neonicotinoid concentrations
lesa <- lesa[!is.na(lesa$OverallNeonic), ]

# Make neonicotinoid detection column
lesa$Detection <- ifelse(lesa$OverallNeonic > 0, "Detection", "Non-detection")

# Convert categorical variables to factor
lesa$Season <- as.factor(lesa$Season)
lesa$Detection <- as.factor(lesa$Detection)

# Check for variable correlation
cor(lesa$Julian, lesa$MigDate)

# ---------------------------------------------------------------------------- #

## Covariate Restrictions ####

# Date into Season and Julian can't be in the same model (correlation = 0.90)
# Sex and Year not included because all birds captured were males in 2023
# Julian and Season can't be in the same model because Season bins Julian

# ---------------------------------------------------------------------------- #

# Data Exploration & Linear Regression Model Assumptions ####

## Assumption 1: Zero Expectations of Errors (Residuals) ####

### Residuals Plots 

# Calculate residuals
m <- lm(Mass ~ OverallNeonic, data = lesa)
lesa$yhat <- predict(m, data = lesa)
lesa$residuals <- residuals(m)
lesa$rstudent <- rstudent(m)
head(lesa)

# Simple residual plot
plot(predict(m), rstudent(m))

# Pretty residual plot
ggplot(lesa, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# Dr. Johnson confirmed these residual plots are good even though more 
# variable around 0. The "megaphone" pattern is not concerning as it's mostly 
# a reflection of non-detects (the opposite direction is more concerning). 

# Log transformation of neonics

# Calculate residuals 
m.log <- lm(Mass ~ (log10(OverallNeonic+0.0001)), data = lesa)
lesa$yhat <- predict(m.log)
lesa$residuals <- residuals(m.log)
lesa$rstudent <- rstudent(m.log)
head(lesa)

# Simple residual plot
plot(predict(m.log), rstudent(m.log))

# Pretty residual plot
ggplot(lesa, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# Looks slightly better with logarithmic transformation.

# ---------------------------------------------------------------------------- #

## Assumption 2: Equality of error variances (Homoscedasticity) #### 

# No transformation
m <- lm(Mass ~ OverallNeonic, data = lesa)
lesa$yhat <- predict(m)
lesa$residuals <- residuals(m)
lesa$rstudent <- rstudent(m)

# Plot homoscedasticity
ggplot(lesa, aes(x = yhat, y = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual")

# No violations because variability doesn't increase as predicted value increases

# Log transformation
m.log <- lm(Mass ~ (log10(OverallNeonic+0.0001)), data = lesa)
lesa$yhat <- predict(m.log)
lesa$residuals <- residuals(m.log)
lesa$rstudent <- rstudent(m.log)

### PLOT HOMOSCEDASTICITY
ggplot(lesa, aes(x = yhat, y = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual")

# No violations because variability doesn't increase as predicted value increases.
# Log transformation does appear to spread residuals more evenly. 

# ---------------------------------------------------------------------------- #

## Assumption 3: Independence of Errors/Observations ####

# Distribution of one error/observation does not depend on another error/observation

# ---------------------------------------------------------------------------- #

## Assumption 4: Normality ####

m <- lm(Mass ~ OverallNeonic, data = lesa)
hist(rstandard(m))
qqnorm(rstandard(m))

# Does not severely violate normality assumption
# Regression models do not many any assumptions about the distrbution of explanatory 
# variables (e.g., Concentrations), but their distribution does affect inferences. 

# ---------------------------------------------------------------------------- #

# First Stage AICc Model Selection: Mass ####

### Controlling Variables: Season + Julian + MigDate

### Null and global models
m.null <- lm(Mass ~ 1, data = lesa)
m.global <- lm(Mass ~ Season*MigDate, data = lesa)

### Additive combinations
m1 <- lm(Mass ~ Season, data = lesa)
m2 <- lm(Mass ~ Julian, data = lesa)
m3 <- lm(Mass ~ MigDate, data = lesa)
m4 <- lm(Mass ~ Season + MigDate, data = lesa)

# ---------------------------------------------------------------------------- #

### AICc Model Selection: First Stage ####
models <- list(m.null, m.global, m1, m2, m3, m4)

mod.names <- c('m.null', 'm.global', 'm1', 'm2', 'm3', 'm4')

aictab(models, modnames = mod.names)

# ---------------------------------------------------------------------------- #

### Top Model Summaries: First Stage ####
options(digits = 3)

# Mass ~ Season + MigDate (m4) <- INFORMED NULL MODEL
cbind(summary(m4)$coefficients, confint(m4))

# Mass ~ Season * MigDate (m.global)
cbind(summary(m.global)$coefficients, confint(m.global))

# ---------------------------------------------------------------------------- #

### Conclusion: First Stage ####

# 1. On average, birds in the spring had significantly higher body mass compared to birds in the fall (m4). 
# 2. For each additional day into the migration season, body mass increases significantly (m4).
# 3. The effect of migration date on body mass does not significantly differ between spring and fall. 

# ---------------------------------------------------------------------------- #

# Second Stage AICc Model Selection: Mass ~ Concentration ####

### Null and global models
m.null <- lm(Mass ~ 1, data = lesa)
m.informednull <- lm(Mass ~ Season + MigDate, data = lesa)
m.globalinformednull <- lm(Mass ~ Season*OverallNeonic + MigDate*OverallNeonic, data = lesa)
m.globalseason <- lm(Mass ~ Season + MigDate + OverallNeonic +
                       Season*MigDate + Season*OverallNeonic + MigDate*OverallNeonic, data = lesa)

m.globaljulian <- lm(Mass ~ Julian*OverallNeonic, data = lesa)

### Single combinations of informed and global model covariates
m1 <- lm(Mass ~ OverallNeonic, data = lesa)
m2 <- lm(Mass ~ OverallNeonic + Season, data = lesa)
m3 <- lm(Mass ~ OverallNeonic + Julian, data = lesa)
m4 <- lm(Mass ~ OverallNeonic + MigDate, data = lesa)

### Two additive combinations of informed and global model covariates
m5 <- lm(Mass ~ OverallNeonic + Season + MigDate, data = lesa)

### Two-way interactions (no additive combinations)
m6 <- lm(Mass ~ OverallNeonic*Season, data = lesa)
m7 <- lm(Mass ~ OverallNeonic*Julian, data = lesa)
m8 <- lm(Mass ~ OverallNeonic*MigDate, data = lesa)

### Two-way interactions with a single additive combinations of informed and global model covariates
m9 <- lm(Mass ~ OverallNeonic*Season + MigDate, data = lesa)
m10 <- lm(Mass ~ OverallNeonic*MigDate + Season, data = lesa)
m11 <- lm(Mass ~ OverallNeonic + Season*MigDate, data = lesa)
m12 <- lm(Mass ~ (OverallNeonic * Season) + Season*MigDate, data = lesa)
m13 <- lm(Mass ~ (OverallNeonic * MigDate) + Season*MigDate, data = lesa)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

# Manually list the specially named models
special_models <- list(m.null, m.informednull, m.globalinformednull,
                       m.globalseason)

special_names <- c('m.null', 'm.informednull', 'm.globalinformednull',
                   'm.globalseason')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Add the other models (m1 to m13) to the list
for (i in c(1:13)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

### Top Model Summaries: Second Stage ####

options(digits = 3)

# Mass ~ OverallNeonic + MigDate (m4)
cbind(summary(m4)$coefficients, confint(m4))

# Mass ~ OverallNeonic * MigDate (m8)
cbind(summary(m8)$coefficients, confint(m8))

# Mass ~ OverallNeonic + Season + MigDate (m5)
cbind(summary(m5)$coefficients, confint(m5))

# Mass ~ Season + MigDate (m.informednull)
cbind(summary(m.informednull)$coefficients, confint(m.informednull))

# ----------------------------------------------------------------------------- #

### Conclusions: Second Stage without Log Transformation ####

# 1. There is a positive association between neonicotinoid concentration and body mass (m4).
# 2. Each additional day into migration is associated with an increase in body mass (m4, m5, m.informednull).
# 3. The effect of migration date on mass does not vary with neonicotinoid concentration (m8).
# 4. Birds are heavier in the spring compared to the fall (m.informednull). 

# **Q: Why is there a positive association between neonicotinoid concentrations and body mass? What do these birds with detections and high concentrations have in common? 

# **A: There are only four birds with detections. Two are from spring (semipermanent wetland, Virgil site [77% ag]) and two are from fall (temporary wetland, Glinz [54% ag] and Muske [46% ag]). All are from 2023 and are male. 

# **A: This positive association is likely due to small sample size. 

# ----------------------------------------------------------------------------- #

# Second Stage AICc Model Selection: Mass ~ Log10(Concentration) ####

### Null and global models

m.null <- lm(Mass ~ 1, data = lesa)
m.informednull <- lm(Mass ~ Season + MigDate, data = lesa)
m.globalinformednull <- lm(Mass ~ Season*log10(OverallNeonic + 0.0001) + MigDate*log10(OverallNeonic + 0.0001), data = lesa)
m.globalseason <- lm(Mass ~ Season + MigDate + log10(OverallNeonic + 0.0001) +
                       Season*MigDate + Season*log10(OverallNeonic + 0.0001) + MigDate*log10(OverallNeonic + 0.0001), data = lesa)

m.globaljulian <- lm(Mass ~ Julian*log10(OverallNeonic + 0.0001), data = lesa)

### Single combinations of informed and global model covariates
m1 <- lm(Mass ~ log10(OverallNeonic + 0.0001), data = lesa)
m2 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season, data = lesa)
m3 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Julian, data = lesa)
m4 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + MigDate, data = lesa)

### two additive combinations of informed and global model covariates 
m5 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season + MigDate, data = lesa)

### two-way interactions (no additive combinations)
m6 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season, data = lesa)
m7 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian, data = lesa)
m8 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate, data = lesa)

### two-way interactions with a single additive combinations of informed and global model covariates
m9 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season + MigDate, data = lesa)
m10 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Season, data = lesa)
m11 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season*MigDate, data = lesa)
m12 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Season*MigDate, data = lesa)
m13 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season*MigDate, data = lesa)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

#Manually list the specially named models
special_models <- list(m.null, m.informednull, m.globalinformednull,
                       m.globalseason)

special_names <- c('m.null', 'm.informednull', 'm.globalinformednull',
                   'm.globalseason')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Dynamically add the other models (m1 to m360) to the list
for (i in c(1:13)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

### Top Model Summaries: Second Stage ####

options(digits = 3)

# Mass ~ Log(Neonic) + MigDate (m4)
cbind(summary(m4)$coefficients, confint(m4))

# Mass ~ Log(Neonic) + Season + MigDate (m5)
cbind(summary(m5)$coefficients, confint(m5))

# Mass ~ Log(Neonic) (m1)
cbind(summary(m1)$coefficients, confint(m1))

# ----------------------------------------------------------------------------- #

### Conclusion: Second Stage with Log Transformation ####

# 1. There is a positive association between body mass and log of neonicotinoid concentration (m1, m5, m4). 
# 2. There is a positive association between body mass and date into migration season (m5).

# **The positive associations of neonic and body mass is likely due to small sample size. 

# ----------------------------------------------------------------------------- #

# Second Stage AICc Model Selection: Mass ~ Detection ####

### null and global models
m.null <- lm(Mass ~ 1, data = lesa)
m.informednull <- lm(Mass ~ Season + MigDate, data = lesa)
m.globalinformednull <- lm(Mass ~ Season*Detection + MigDate*Detection, data = lesa)
m.globalseason <- lm(Mass ~ Season + MigDate + Detection +
                       Season*MigDate + Season*Detection + MigDate*Detection, data = lesa)

### single combinations of informed and global model covariates with neonicotinoid concentrations
m1 <- lm(Mass ~ Detection, data = lesa)
m2 <- lm(Mass ~ Detection + Season, data = lesa)
m3 <- lm(Mass ~ Detection + Julian, data = lesa)
m4 <- lm(Mass ~ Detection + MigDate, data = lesa)

### two additive combinations of informed and global model covariates with neonicotinoid concentrations
m5 <- lm(Mass ~ Detection + Season + MigDate, data = lesa)

### two-way interactions (no additive combinations)
m6 <- lm(Mass ~ Detection*Season, data = lesa)
m7 <- lm(Mass ~ Detection*Julian, data = lesa)
m8 <- lm(Mass ~ Detection*MigDate, data = lesa)

### two-way interactions with a single additive combinations of informed and global model covariates with neonicotinoid concentrations
m9 <- lm(Mass ~ Detection*Season + MigDate, data = lesa)
m10 <- lm(Mass ~ Detection*MigDate + Season, data = lesa)
m11 <- lm(Mass ~ Detection + Season*MigDate, data = lesa)
m12 <- lm(Mass ~ (Detection * Season) + Season*MigDate, data = lesa)
m13 <- lm(Mass ~ (Detection * MigDate) + Season*MigDate, data = lesa)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

# Manually list the special models
special_models <- list(m.null, m.informednull, m.globalinformednull,
                       m.globalseason)

special_names <- c('m.null', 'm.informednull', 'm.globalinformednull',
                   'm.globalseason')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Dynamically add the other models (m1 to m360) to the list
for (i in c(1:13)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

## Top Model Summaries: Second Stage ####

options(digits = 3)

# Mass ~ Detection + MigDate (m4)
cbind(summary(m4)$coefficients, confint(m4)) 

# Mass ~ Detection + Season + MigDate (m5)
cbind(summary(m5)$coefficients, confint(m5)) 

# Mass ~ Detection (m1)
cbind(summary(m1)$coefficients, confint(m1)) 

### Conclusion: Second Stage with Detection ####

# ----------------------------------------------------------------------------- #

# 1. Birds with detections had significantly higher body mass (m1, m4).
# 2. Birds with later migration dates had significantly higher body mass (m5).

# **There were only four birds with detections. Positive association is likely due to small sample size.

# ----------------------------------------------------------------------------- #

# Model Comparisons: Concentration, Log(Concentration), Detection ####
m.null <- lm(Mass ~ 1, data = lesa)
m.informednull <- lm(Mass ~ Season + MigDate, data = lesa)
m.c <- lm(Mass ~ OverallNeonic, data = lesa)
m.l <- lm(Mass ~ log10(OverallNeonic + 0.0001), data = lesa)
m.d <- lm(Mass ~ Detection, data = lesa)

models <- list(m.c, m.l, m.d, m.null, m.informednull)

mod.names <- c('m.c', 'm.l', 'm.d', 'm.null', 'm.informednull')

aictab(models, modnames = mod.names)

# ----------------------------------------------------------------------------- #

## Top Model Summaries: Model Comparisons ####

cbind(summary(m.l)$coefficients, confint(m.l))
cbind(summary(m.c)$coefficients, confint(m.c))
cbind(summary(m.d)$coefficients, confint(m.d))

# ----------------------------------------------------------------------------- #

### Conclusion: Model Comparisons ####

# 1. Log(Neonic) and Detection models received similar support (AICc [0.00 - 0.34]). 
# 2. Use model with detections as that makes more sense biologically. 

# ----------------------------------------------------------------------------- #
