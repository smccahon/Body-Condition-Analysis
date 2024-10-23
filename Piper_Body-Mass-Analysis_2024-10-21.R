#---------------------------------#
# Small Piper Body Mass Analysis  #
#       Created 10/21/2024        #
#       Modified 10/23/2024       #
#---------------------------------#

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

# Subset data for pipers
piper <- subset(birds, Species %in% c("SemipalmatedSandpiper", "LeastSandpiper"))

# Remove birds that do not have neonicotinoid concentrations
piper <- piper[!is.na(piper$OverallNeonic), ]

# Make neonicotinoid detection column
piper$Detection <- ifelse(piper$OverallNeonic > 0, "Detection", "Non-detection")

# Convert categorical variables to factor
piper$Season <- as.factor(piper$Season)
piper$Year <- as.factor(piper$Year)
piper$Sex <- as.factor(piper$Sex)
piper$Detection <- as.factor(piper$Detection)

# Check for variable correlation
cor(piper$Julian, piper$MigDate)

# ---------------------------------------------------------------------------- #

## Covariate Restrictions ####

# Date into Season and Julian can't be in the same model because of high correlation (0.924)
# Year and Season can't be in the same model because of lack of replication of Season within a year
# Julian and Season can't be in the same model because Season bins Julian
# Year and Species can't be in the same model because all LESA were captured in 2023

# ---------------------------------------------------------------------------- #

# Data Exploration & Linear Regression Model Assumptions ####

## Assumption 1: Zero Expectations of Errors (Residuals) ####

### Residuals Plots 

# Calculate residuals
m <- lm(Mass ~ OverallNeonic, data = piper)
piper$yhat <- predict(m)
piper$residuals <- residuals(m)
piper$rstudent <- rstudent(m)
head(piper)

# Simple residual plot
plot(predict(m), rstudent(m))

# Pretty residual plot
ggplot(piper, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# Dr. Johnson confirmed these residual plots are good even though more 
# variable around 0. The "megaphone" pattern is not concerning as it's mostly 
# a reflection of non-detects (the opposite direction is more concerning). 

# Log transformation of neonics

# Calculate residuals 
m.log <- lm(Mass ~ (log10(OverallNeonic+0.0001)), data = piper)
piper$yhat <- predict(m.log)
piper$residuals <- residuals(m.log)
piper$rstudent <- rstudent(m.log)
head(piper)

# Simple residual plot
plot(predict(m.log), rstudent(m.log))

# Pretty residual plot
ggplot(piper, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# Looks better with logarithmic transformation.

# ---------------------------------------------------------------------------- #

## Assumption 2: Equality of error variances (Homoscedasticity) #### 

# No transformation
m <- lm(Mass ~ OverallNeonic, data = piper)
piper$yhat <- predict(m)
piper$residuals <- residuals(m)
piper$rstudent <- rstudent(m)

# Plot homoscedasticity
ggplot(piper, aes(x = yhat, y = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual")

# No violations because variability doesn't increase as predicted value increases

# Log transformation
m.log <- lm(Mass ~ (log10(OverallNeonic+0.0001)), data = piper)
piper$yhat <- predict(m.log)
piper$residuals <- residuals(m.log)
piper$rstudent <- rstudent(m.log)

### PLOT HOMOSCEDASTICITY
ggplot(piper, aes(x = yhat, y = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual")

# No violations because variability doesn't increase as predicted value increases.
# Log transformation does appear to spread residuals more evenly. 

# ---------------------------------------------------------------------------- #

## Assumption 3: Independence of Errors/Observations ####

# Distribution of one error/observation does not depend on another error/observation

# ---------------------------------------------------------------------------- #

## Assumption 4: Normality ####

m <- lm(Mass ~ OverallNeonic, data = piper)
hist(rstandard(m))
qqnorm(rstandard(m))

# Does not severely violate normality assumption
# Regression models do not many any assumptions about the distrbution of explanatory 
# variables (e.g., Concentrations), but their distribution does affect inferences. 

# ---------------------------------------------------------------------------- #

# First Stage AICc Model Selection: Mass ####

### Controlling Variables: Year + Sex + Season + Julian + MigDate + Species

### null and global models
m.null <- lm(Mass ~ 1, data= piper)
m.globalyearjulian <- lm(Mass ~ Year + Sex + Julian + 
                           Year*Sex + Year*Julian +
                           Sex*Julian + Sex*Species, data = piper)

m.globalspeciesjulian <-lm(Mass ~ Species + Sex + Julian +
                             Species*Sex + Species*Julian + Sex*Julian +
                             Species, data = piper)

m.globalyearmig <- lm(Mass ~ Year + Sex + MigDate +
                        Year*Sex + Year*MigDate +
                        Sex*MigDate, data = piper)

m.globalseasonmig <- lm(Mass ~ Sex + Season + Species + MigDate +
                          Sex*Season + Sex*Species + Sex*MigDate +
                          Season*Species + Season*MigDate + Species*MigDate, data = piper)

### additive combinations (with restrictions, see above)
m1 <- lm(Mass ~ Year, data = piper)
m2 <- lm(Mass ~ Season, data = piper)
m3 <- lm(Mass ~ Sex, data = piper)
m4 <- lm(Mass ~ Species, data = piper)
m5 <- lm(Mass ~ Julian, data = piper)
m6 <- lm(Mass ~ MigDate, data = piper)

m7 <- lm(Mass ~ Year + Sex, data = piper)
m8 <- lm(Mass ~ Year + Julian, data = piper)
m9  <- lm(Mass ~ Year + MigDate, data = piper)
m10 <- lm(Mass ~ Season + Sex, data = piper)
m11 <- lm(Mass ~ Season + Species, data = piper)
m12 <- lm(Mass ~ Season + MigDate, data = piper)
m13 <- lm(Mass ~ Sex + Species, data = piper)
m14 <- lm(Mass ~ Sex + Julian, data = piper)
m15 <- lm(Mass ~ Sex + MigDate, data = piper)
m16 <- lm(Mass ~ Species + Julian, data = piper)
m17 <- lm(Mass ~ Species + MigDate, data = piper)

m18 <- lm(Mass ~ Year + Sex + Julian, data = piper)
m19 <- lm(Mass ~ Year + Sex + MigDate, data = piper)

m20 <- lm(Mass ~ Sex + Season + Species, data = piper)
m21 <- lm(Mass ~ Sex + Season + MigDate, data = piper)
m22 <- lm(Mass ~ Sex + Julian + Species, data = piper)
m23 <- lm(Mass ~ Sex + Species + MigDate, data = piper)

m24 <- lm(Mass ~ Season + Species + MigDate, data = piper)

m25 <- lm(Mass ~ Sex + Season + Species + MigDate, data = piper)

### two-way interactions (with restrictions, see above)
m26 <- lm(Mass ~ Year * Sex, data = piper)
m27 <- lm(Mass ~ Year * Julian, data = piper)
m28 <- lm(Mass ~ Year * MigDate, data = piper)

m29 <- lm(Mass ~ Sex * Season, data = piper)
m30 <- lm(Mass ~ Sex * Julian, data = piper)
m31 <- lm(Mass ~ Sex * Species, data = piper)
m32 <- lm(Mass ~ Sex * MigDate, data = piper)

m33 <- lm(Mass ~ Season * Species, data = piper)
m34 <- lm(Mass ~ Season * MigDate, data = piper)
m35 <- lm(Mass ~ Species * Julian, data = piper)

### two-way interactions with one additive combination (with restrictions, see above)

#Year
m36 <- lm(Mass ~ (Year * Sex) + Julian, data = piper)
m37 <- lm(Mass ~ (Year * Sex) + MigDate, data = piper)

m38 <- lm(Mass ~ (Year * Julian) + Sex, data = piper)

m39 <- lm(Mass ~ (Year * MigDate) + Sex, data = piper)

#Sex
m40 <- lm(Mass ~ (Sex * Season) + Species, data = piper)
m41 <- lm(Mass ~ (Sex * Season) + MigDate, data = piper)

m42 <- lm(Mass ~ (Sex * Julian) + Year, data = piper)
m43 <- lm(Mass ~ (Sex * Julian) + Species, data = piper)

m44 <- lm(Mass ~ (Sex * Species) + Season, data = piper)
m45 <- lm(Mass ~ (Sex * Species) + Julian, data = piper)
m46 <- lm(Mass ~ (Sex * Species) + MigDate, data = piper)

m47 <- lm(Mass ~ (Sex * MigDate) + Year, data = piper)
m48 <- lm(Mass ~ (Sex * MigDate) + Season, data = piper)
m49 <- lm(Mass ~ (Sex * MigDate) + Species, data = piper)

#Season
m50 <- lm(Mass ~ (Season * Species) + Sex, data = piper)
m51 <- lm(Mass ~ (Season * Species) + MigDate, data = piper)

m52 <- lm(Mass ~ (Season * MigDate) + Sex, data = piper)
m53 <- lm(Mass ~ (Season * MigDate) + Species, data = piper)

#Species
m54 <- lm(Mass ~ (Species * Julian) + Sex, data = piper)

### two-way interactions with two additive combinations (with restrictions, see above)
m55 <- lm(Mass ~ (Sex * Season) + Species + MigDate, data = piper)
m56 <- lm(Mass ~ Sex + (Season * Species) + MigDate, data = piper)
m57 <- lm(Mass ~ Sex + Season + (Species * MigDate), data = piper)

m58 <- lm(Mass ~ (Sex * MigDate) + Species + Season, data = piper)
m59 <- lm(Mass ~ (Sex * Species) + Season + MigDate, data = piper)

### two-way interactions with multiple interactions (with restrictions, see above)

#two-way
m60 <- lm(Mass ~ Year*Sex + Year*Julian, data = piper)
m61 <- lm(Mass ~ Year*Sex + Year*MigDate, data = piper)
m62 <- lm(Mass ~ Year*Sex + Sex*Julian, data = piper)
m64 <- lm(Mass ~ Year*Sex + Sex*MigDate, data = piper)

m65 <- lm(Mass ~ Year*Julian + Year*Sex, data = piper)
m66 <- lm(Mass ~ Year*Julian + Sex*Julian, data = piper)

m67 <- lm(Mass ~ Year*MigDate + Sex*MigDate, data = piper)

m68 <- lm(Mass ~ Sex*Season + Sex*Species, data = piper)
m69 <- lm(Mass ~ Sex*Season + Sex*MigDate, data = piper)
m70 <- lm(Mass ~ Sex*Season + Season*Species, data = piper)

m71 <- lm(Mass ~ Sex*Julian + Sex*Species, data = piper)
m72 <- lm(Mass ~ Sex*Julian + Species*Julian, data = piper)

m73 <- lm(Mass ~ Sex*Species + Sex*MigDate, data = piper)
m74 <- lm(Mass ~ Sex*Species + Season*Species, data = piper)
m75 <- lm(Mass ~ Sex*Species + Season*MigDate, data = piper)
m76 <- lm(Mass ~ Sex*Species + Species*Julian, data = piper)

m77 <- lm(Mass ~ Sex*MigDate + Season*Species, data = piper)
m78 <- lm(Mass ~ Sex*MigDate + Season*MigDate, data = piper)

m79 <- lm(Mass ~ Season*Species + Season*MigDate, data = piper)

#three-way
m80 <- lm(Mass ~ Year*Sex + Year*Julian + Sex*Julian, data = piper)
m81 <- lm(Mass ~ Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m82 <- lm(Mass ~ Year*Julian + Year*Sex + Sex*Julian, data = piper)

m83 <- lm(Mass ~ Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m84 <- lm(Mass ~ Sex*Season + Sex*Species + Season*Species, data = piper)
m85 <- lm(Mass ~ Sex*Season + Sex*Species + Season*MigDate, data = piper)
m86 <- lm(Mass ~ Sex*Season + Sex*Species + Species*Julian, data = piper)
m87 <- lm(Mass ~ Sex*Season + Sex*MigDate + Season*Species, data = piper)
m88 <- lm(Mass ~ Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m89 <- lm(Mass ~ Sex*Season + Season*Species + Season*MigDate, data = piper)

m90 <- lm(Mass ~ Sex*Julian + Sex*Species + Species*Julian, data = piper)

m91 <- lm(Mass ~ Sex*Species + Sex*MigDate + Season*Species, data = piper)
m92 <- lm(Mass ~ Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m93 <- lm(Mass ~ Sex*Species + Season*Species + Season*MigDate, data = piper)

m94 <- lm(Mass ~ Sex*MigDate + Season*Species + Season*MigDate, data = piper)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

#Manually list the special models
special_models <- list(m.null, m.globalspeciesjulian,
                       m.globalyearjulian, m.globalseasonmig, m.globalyearmig)

special_names <- c('m.null', 'm.globalspeciesjulian',
                   'm.globalyearjulian', 'm.globalseasonmig', 'm.globalyearmig')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Dynamically add the other models (m1 to m158) to the list
for (i in c(1:94)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

### Top Model Summaries: Second Stage ####

options(digits = 3)

# Mass ~ Species + MigDate (m17)
cbind(summary(m17)$coefficients, confint(m17))

# Mass ~ Species (m4)
cbind(summary(m4)$coefficients, confint(m4))

# Mass ~ Season * MigDate + Species (m53)
model_summary <- summary(m53)$coefficients
conf_intervals <- confint(m53, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Season + Species + MigDate (m24)
cbind(summary(m24)$coefficients, confint(m24))

# Mass ~ Sex * Season + Species + MigDate (m55)
cbind(summary(m55)$coefficients, confint(m55))

# ----------------------------------------------------------------------------- #

### Conclusion: First Stage ####

# 1. Semipalmated Sandpipers have statistically significant higher body mass. 

# ----------------------------------------------------------------------------- #

# Second Stage AICc Model Selection: Mass ~ Concentration ####

### null and global models
m.null <- lm(Mass ~ 1, data = piper)
m.informednull <- lm(Mass ~ Species + MigDate, data = piper)

m.globalinformednull <- lm(Mass ~ Species + MigDate + OverallNeonic +
                             Species * MigDate + Species * OverallNeonic + 
                             MigDate * OverallNeonic, data = piper)

m.globalspecies <- lm(Mass ~ Sex + Season + Species + MigDate + OverallNeonic + 
                          Sex*Season + Sex*Species + Sex*MigDate + Sex*OverallNeonic +
                          Season*Species + Season*MigDate + Season*OverallNeonic + 
                          Species*MigDate + Species*OverallNeonic + MigDate*OverallNeonic, data = piper)

m.globalyear <-  lm(Mass ~ Sex + Season + Year + MigDate + OverallNeonic + 
                  Sex*Season + Sex*Year + Sex*MigDate + Sex*OverallNeonic +
                  Season*Year + Season*MigDate + Season*OverallNeonic + 
                  Year*MigDate + Year*OverallNeonic + MigDate*OverallNeonic, data = piper)

### single combinations of informed and global model covariates with neonicotinoid concentrations
m1 <- lm(Mass ~ OverallNeonic, data = piper)
m2 <- lm(Mass ~ OverallNeonic + Year, data = piper)
m3 <- lm(Mass ~ OverallNeonic + Sex, data = piper)
m4 <- lm(Mass ~ OverallNeonic + Species, data = piper)
m5 <- lm(Mass ~ OverallNeonic + Season, data = piper)
m6 <- lm(Mass ~ OverallNeonic + Julian, data = piper)
m7 <- lm(Mass ~ OverallNeonic + MigDate, data = piper)

### two additive combinations of informed and global model covariates with neonicotinoid concentrations
m8 <- lm(Mass ~ OverallNeonic + Year + Sex, data = piper)
m9 <- lm(Mass ~ OverallNeonic + Year + Julian, data = piper)
m10 <- lm(Mass ~ OverallNeonic + Year + MigDate, data = piper)

m11 <- lm(Mass ~ OverallNeonic + Sex + Species, data = piper)
m12 <- lm(Mass ~ OverallNeonic + Sex + Season, data = piper)
m13 <- lm(Mass ~ OverallNeonic + Sex + Julian, data = piper)
m14 <- lm(Mass ~ OverallNeonic + Sex + MigDate, data = piper)

m15 <- lm(Mass ~ OverallNeonic + Species + Season, data = piper)
m16 <- lm(Mass ~ OverallNeonic + Species + Julian, data = piper)
m17 <- lm(Mass ~ OverallNeonic + Species + MigDate, data = piper)

m18 <- lm(Mass ~ OverallNeonic + Season + MigDate, data = piper)

### three additive combinations of informed and global model covariates with neonicotinoid concentrations
m19 <- lm(Mass ~ OverallNeonic + Year + Sex + Julian, data = piper)
m20 <- lm(Mass ~ OverallNeonic + Year + Sex + MigDate, data = piper)

m21 <- lm(Mass ~ OverallNeonic + Sex + Species + Season, data = piper)
m22 <- lm(Mass ~ OverallNeonic + Sex + Species + Julian, data = piper)
m23 <- lm(Mass ~ OverallNeonic + Sex + Species + MigDate, data = piper)

m24 <- lm(Mass ~ OverallNeonic + Sex + Season + MigDate, data = piper)

m25 <- lm(Mass ~ OverallNeonic + Species + Season + MigDate, data = piper)

### four additive combinations of informed and global model covariates with neonicotinoid concentrations
m26 <- lm(Mass ~ OverallNeonic + Season + Sex + Species + MigDate, data = piper)

### two-way interactions (no additive combinations)
m27 <- lm(Mass ~ OverallNeonic*Year, data = piper)
m28 <- lm(Mass ~ OverallNeonic*Sex, data = piper)
m29 <- lm(Mass ~ OverallNeonic*Species, data = piper)
m30 <- lm(Mass ~ OverallNeonic*Season, data = piper)
m31 <- lm(Mass ~ OverallNeonic*Julian, data = piper)
m32 <- lm(Mass ~ OverallNeonic*MigDate, data = piper)

### two-way interactions with a single additive combinations of informed and global model covariates with neonicotinoid concentrations
m33 <- lm(Mass ~ OverallNeonic*Year + MigDate, data = piper)
m34 <- lm(Mass ~ OverallNeonic*Year + Julian, data = piper)
m35 <- lm(Mass ~ OverallNeonic*Year + Sex, data = piper)

m36 <- lm(Mass ~ OverallNeonic*Sex + Year, data = piper)
m37 <- lm(Mass ~ OverallNeonic*Sex + Season, data = piper)
m38 <- lm(Mass ~ OverallNeonic*Sex + MigDate, data = piper)
m39 <- lm(Mass ~ OverallNeonic*Sex + Julian, data = piper)
m40 <- lm(Mass ~ OverallNeonic*Sex + Species, data = piper)

m41 <- lm(Mass ~ OverallNeonic*Species + Season, data = piper)
m42 <- lm(Mass ~ OverallNeonic*Species + MigDate, data = piper)
m43 <- lm(Mass ~ OverallNeonic*Species + Julian, data = piper)
m44 <- lm(Mass ~ OverallNeonic*Species + Sex, data = piper)

m45 <- lm(Mass ~ OverallNeonic*Season + MigDate, data = piper)
m46 <- lm(Mass ~ OverallNeonic*Season + Species, data = piper)
m47 <- lm(Mass ~ OverallNeonic*Season + Sex, data = piper)

m48 <- lm(Mass ~ OverallNeonic*Julian + Year, data = piper)
m49 <- lm(Mass ~ OverallNeonic*Julian + Species, data = piper)
m50 <- lm(Mass ~ OverallNeonic*Julian + Sex, data = piper)

m51 <- lm(Mass ~ OverallNeonic*MigDate + Year, data = piper)
m52 <- lm(Mass ~ OverallNeonic*MigDate + Season, data = piper)
m53 <- lm(Mass ~ OverallNeonic*MigDate + Species, data = piper)
m54 <- lm(Mass ~ OverallNeonic*MigDate + Sex, data = piper)

### two-way interactions with two additive combination of informed and global model covariates with neonicotinoid concentrations
#Year
m55 <- lm(Mass ~ OverallNeonic + (Year * Sex) + Julian, data = piper)
m56 <- lm(Mass ~ (OverallNeonic * Year) + Sex + Julian, data = piper)
m57 <- lm(Mass ~ (OverallNeonic * Sex) + Year + Julian, data = piper)
m58 <- lm(Mass ~ (OverallNeonic * Julian) + Sex + Julian, data = piper)

m59 <- lm(Mass ~ OverallNeonic + (Year * Sex) + MigDate, data = piper)
m60 <- lm(Mass ~ (OverallNeonic * Year) + Sex + MigDate, data = piper)
m61 <- lm(Mass ~ (OverallNeonic * Sex) + Year + MigDate, data = piper)
m62 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Year, data = piper)

m63 <- lm(Mass ~ OverallNeonic + (Year * Julian) + Sex, data = piper)
m64 <- lm(Mass ~ (OverallNeonic * Year) + Julian + Sex, data = piper)
m65 <- lm(Mass ~ (OverallNeonic * Julian) + Year + Sex, data = piper)
m66 <- lm(Mass ~ (OverallNeonic * Sex) + Julian + Year, data = piper)

m67 <- lm(Mass ~ OverallNeonic + (Year * MigDate) + Sex, data = piper)
m68 <- lm(Mass ~ (OverallNeonic * Year) + MigDate + Sex, data = piper)
m69 <- lm(Mass ~ (OverallNeonic * MigDate) + Year + Sex, data = piper)
m70 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Year, data = piper)

#Sex
m71 <- lm(Mass ~ OverallNeonic + (Sex * Season) + Species, data = piper)
m72 <- lm(Mass ~ (OverallNeonic * Sex) + Season + Species, data = piper)
m73 <- lm(Mass ~ (OverallNeonic * Season) + Sex + Species, data = piper)
m74 <- lm(Mass ~ (OverallNeonic * Species) + Season + Sex, data = piper)

m75 <- lm(Mass ~ OverallNeonic + (Sex * Season) + MigDate, data = piper)
m76 <- lm(Mass ~ (OverallNeonic * Sex) + Season + MigDate, data = piper)
m77 <- lm(Mass ~ (OverallNeonic * Season) + Sex + MigDate, data = piper)
m78 <- lm(Mass ~ (OverallNeonic * MigDate) + Season + Sex, data = piper)

m79 <- lm(Mass ~ OverallNeonic + (Sex * Julian) + Year, data = piper)
m80 <- lm(Mass ~ (OverallNeonic * Sex) + Julian + Year, data = piper)
m81 <- lm(Mass ~ (OverallNeonic * Julian) + Sex + Year, data = piper)
m82 <- lm(Mass ~ (OverallNeonic * Year) + Julian + Sex, data = piper)

m83 <- lm(Mass ~ OverallNeonic + (Sex * Julian) + Species, data = piper)
m84 <- lm(Mass ~ (OverallNeonic * Sex) + Julian + Species, data = piper)
m85 <- lm(Mass ~ (OverallNeonic * Julian) + Sex + Species, data = piper)
m86 <- lm(Mass ~ (OverallNeonic * Species) + Julian + Sex, data = piper)

m87 <- lm(Mass ~ OverallNeonic + (Sex * Species) + Season, data = piper)
m88 <- lm(Mass ~ (OverallNeonic * Sex) + Species + Season, data = piper)
m89 <- lm(Mass ~ (OverallNeonic * Species) + Sex + Season, data = piper)
m90 <- lm(Mass ~ (OverallNeonic * Season) + Species + Sex, data = piper)

m91 <- lm(Mass ~ OverallNeonic + (Sex * Species) + Julian, data = piper)
m92 <- lm(Mass ~ (OverallNeonic * Sex) + Species + Julian, data = piper)
m93 <- lm(Mass ~ (OverallNeonic * Species) + Sex + Julian, data = piper)
m94 <- lm(Mass ~ (OverallNeonic * Julian) + Species + Sex, data = piper)

m95 <- lm(Mass ~ OverallNeonic + (Sex * Species) + MigDate, data = piper)
m96 <- lm(Mass ~ (OverallNeonic * Sex) + Species + MigDate, data = piper)
m97 <- lm(Mass ~ (OverallNeonic * Species) + Sex + MigDate, data = piper)
m98 <- lm(Mass ~ (OverallNeonic * MigDate) + Species + Sex, data = piper)

m99 <- lm(Mass ~ OverallNeonic + (Sex * MigDate) + Year, data = piper)
m100 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Year, data = piper)
m101 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Year, data = piper)
m102<- lm(Mass ~ (OverallNeonic * Year) + MigDate + Sex, data = piper)

m103 <- lm(Mass ~ OverallNeonic + (Sex * MigDate) + Season, data = piper)
m104 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Season, data = piper)
m105 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Season, data = piper)
m106 <- lm(Mass ~ (OverallNeonic * Season) + MigDate + Sex, data = piper)

m107 <- lm(Mass ~ OverallNeonic + (Sex * MigDate) + Species, data = piper)
m108 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Species, data = piper)
m109 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Species, data = piper)
m110 <- lm(Mass ~ (OverallNeonic * Species) + MigDate + Sex, data = piper)

#Season
m111 <- lm(Mass ~ OverallNeonic + (Season * Species) + Sex, data = piper)
m112 <- lm(Mass ~ (OverallNeonic * Season) + Species + Sex, data = piper)
m113 <- lm(Mass ~ (OverallNeonic * Species) + Season + Sex, data = piper)
m114 <- lm(Mass ~ (OverallNeonic * Sex) + Species + Season, data = piper)

m115 <- lm(Mass ~ OverallNeonic + (Season * Species) + MigDate, data = piper)
m116 <- lm(Mass ~ (OverallNeonic * Season) + Species + MigDate, data = piper)
m117 <- lm(Mass ~ (OverallNeonic * Species) + Season + MigDate, data = piper)
m118 <- lm(Mass ~ (OverallNeonic * MigDate) + Species + Season, data = piper)

m119 <- lm(Mass ~ OverallNeonic + (Season * MigDate) + Sex, data = piper)
m120 <- lm(Mass ~ (OverallNeonic * Season) + MigDate + Sex, data = piper)
m121 <- lm(Mass ~ (OverallNeonic * MigDate) + Season + Sex, data = piper)
m122 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Season, data = piper)

m123 <- lm(Mass ~ OverallNeonic + (Season * MigDate) + Species, data = piper)
m124 <- lm(Mass ~ (OverallNeonic * Season) + MigDate + Species, data = piper)
m125 <- lm(Mass ~ (OverallNeonic * MigDate) + Season + Species, data = piper)
m126 <- lm(Mass ~ (OverallNeonic * Species) + MigDate + Season, data = piper)

#Species)
m127 <- lm(Mass ~ (OverallNeonic * Species) + Julian + Sex, data = piper)
m128 <- lm(Mass ~ (OverallNeonic * Julian) + Species + Sex, data = piper)
m129 <- lm(Mass ~ (OverallNeonic * Sex) + Julian + Species, data = piper)

### two-way interactions with three additive combinations (with restrictions, see above)
# Year + Sex + Species + Season + Julian + MigDate
m130 <- lm(Mass ~ OverallNeonic + (Sex * Season) + Species + MigDate, data = piper)
m131 <- lm(Mass ~ OverallNeonic + Sex + (Season * Species) + MigDate, data = piper)
m132 <- lm(Mass ~ OverallNeonic + Sex + Season + (Species * MigDate), data = piper)
m133 <- lm(Mass ~ (OverallNeonic * Sex) + Season + Species + MigDate, data = piper)
m134 <- lm(Mass ~ (OverallNeonic * Season) + Sex + Species + MigDate, data = piper)
m135 <- lm(Mass ~ (OverallNeonic * Species) + Season + Sex + MigDate, data = piper)
m136 <- lm(Mass ~ (OverallNeonic * MigDate) + Season + Species + Sex, data = piper)

m137 <- lm(Mass ~ OverallNeonic + (Sex * MigDate) + Species + Season, data = piper)
m138 <- lm(Mass ~ OverallNeonic + (Sex * Species) + Season + MigDate, data = piper)

### neonicotinoid concentrations with one interaction (with restrictions, see above)
m139 <- lm(Mass ~ OverallNeonic + Year*Sex, data = piper)
m140 <- lm(Mass ~ OverallNeonic + Year*Season, data = piper)
m141 <- lm(Mass ~ OverallNeonic + Year*Julian, data = piper)
m142 <- lm(Mass ~ OverallNeonic + Year*MigDate, data = piper)

m143 <- lm(Mass ~ OverallNeonic + Sex*Species, data = piper)
m144 <- lm(Mass ~ OverallNeonic + Sex*Season, data = piper)
m145 <- lm(Mass ~ OverallNeonic + Sex*Julian, data = piper)
m146 <- lm(Mass ~ OverallNeonic + Sex*MigDate, data = piper)

m147 <- lm(Mass ~ OverallNeonic + Species*Season, data = piper)
m148 <- lm(Mass ~ OverallNeonic + Species*Julian, data = piper)
m149 <- lm(Mass ~ OverallNeonic + Species*MigDate, data = piper)

m150 <- lm(Mass ~ OverallNeonic + Season*MigDate, data = piper)

### neonicotinoid concentrations with two interactions (with restrictions, see above)
m151 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*Julian, data = piper)
m152 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*MigDate, data = piper)
m153 <- lm(Mass ~ OverallNeonic + Year*Sex + Sex*Julian, data = piper)
m154 <- lm(Mass ~ OverallNeonic + Year*Sex + Sex*MigDate, data = piper)

m155 <- lm(Mass ~ OverallNeonic + Year*Julian + Year*Sex, data = piper)
m156 <- lm(Mass ~ OverallNeonic + Year*Julian + Sex*Julian, data = piper)

m157 <- lm(Mass ~ OverallNeonic + Year*MigDate + Sex*MigDate, data = piper)

m158 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*Species, data = piper)
m159 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*MigDate, data = piper)
m160 <- lm(Mass ~ OverallNeonic + Sex*Season + Season*Species, data = piper)

m161 <- lm(Mass ~ OverallNeonic + Sex*Julian + Sex*Species, data = piper)
m162 <- lm(Mass ~ OverallNeonic + Sex*Julian + Species*Julian, data = piper)

m163 <- lm(Mass ~ OverallNeonic + Sex*Species + Sex*MigDate, data = piper)
m164 <- lm(Mass ~ OverallNeonic + Sex*Species + Season*Species, data = piper)
m165 <- lm(Mass ~ OverallNeonic + Sex*Species + Season*MigDate, data = piper)
m166 <- lm(Mass ~ OverallNeonic + Sex*Species + Species*Julian, data = piper)

m167 <- lm(Mass ~ OverallNeonic + Sex*MigDate + Season*Species, data = piper)
m168 <- lm(Mass ~ OverallNeonic + Sex*MigDate + Season*MigDate, data = piper)
m169 <- lm(Mass ~ OverallNeonic + Season*Species + Season*MigDate, data = piper)

### neonicotinoid concentrations with three interactions (with restrictions, see above)
m170 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*Julian + Sex*Julian, data = piper)
m171 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m172 <- lm(Mass ~ OverallNeonic + Year*Julian + Year*Sex + Sex*Julian, data = piper)
m173 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m174 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*Species + Season*Species, data = piper)
m175 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m176 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*Species + Species*Julian, data = piper)
m177 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m178 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m179 <- lm(Mass ~ OverallNeonic + Sex*Season + Season*Species + Season*MigDate, data = piper)
m180 <- lm(Mass ~ OverallNeonic + Sex*Julian + Sex*Species + Species*Julian, data = piper)
m181 <- lm(Mass ~ OverallNeonic + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m182 <- lm(Mass ~ OverallNeonic + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m183 <- lm(Mass ~ OverallNeonic + Sex*Species + Season*Species + Season*MigDate, data = piper)
m184 <- lm(Mass ~ OverallNeonic + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

### add neonicotinoid interaction with all other interaction combinations
#add single interaction
#Year
m185 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex, data = piper)
m186 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian, data = piper)
m187 <- lm(Mass ~ (OverallNeonic * Year) + Year*MigDate, data = piper)
m188 <- lm(Mass ~ (OverallNeonic * Year) + Sex*Season, data = piper)
m189 <- lm(Mass ~ (OverallNeonic * Year) + Sex*Julian, data = piper)
m190 <- lm(Mass ~ (OverallNeonic * Year) + Sex*MigDate, data = piper)
m191 <- lm(Mass ~ (OverallNeonic * Year) + Season*MigDate, data = piper)

#Sex
m192 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex, data = piper)
m193 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian, data = piper)
m194 <- lm(Mass ~ (OverallNeonic * Sex) + Year*MigDate, data = piper)
m195 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species, data = piper)
m196 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season, data = piper)
m197 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Julian, data = piper)
m198 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*MigDate, data = piper)
m199 <- lm(Mass ~ (OverallNeonic * Sex) + Species*Season, data = piper)
m200 <- lm(Mass ~ (OverallNeonic * Sex) + Species*Julian, data = piper)
m201 <- lm(Mass ~ (OverallNeonic * Sex) + Species*MigDate, data = piper)
m202 <- lm(Mass ~ (OverallNeonic * Sex) + Season*MigDate, data = piper)

#Species
m203 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species, data = piper)
m204 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season, data = piper)
m205 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Julian, data = piper)
m206 <- lm(Mass ~ (OverallNeonic * Species) + Sex*MigDate, data = piper)
m207 <- lm(Mass ~ (OverallNeonic * Species) + Species*Season, data = piper)
m208 <- lm(Mass ~ (OverallNeonic * Species) + Species*Julian, data = piper)
m209 <- lm(Mass ~ (OverallNeonic * Species) + Species*MigDate, data = piper)
m210 <- lm(Mass ~ (OverallNeonic * Species) + Season*MigDate, data = piper)

#Season
m211 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species, data = piper)
m212 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season, data = piper)
m213 <- lm(Mass ~ (OverallNeonic * Season) + Sex*MigDate, data = piper)
m214 <- lm(Mass ~ (OverallNeonic * Season) + Species*Season, data = piper)
m215 <- lm(Mass ~ (OverallNeonic * Season) + Species*MigDate, data = piper)
m216 <- lm(Mass ~ (OverallNeonic * Season) + Season*MigDate, data = piper)

#Julian
m217 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex, data = piper)
m218 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian, data = piper)
m219 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Species, data = piper)
m220 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Julian, data = piper)
m221 <- lm(Mass ~ (OverallNeonic * Julian) + Species*Julian, data = piper)

#MigDate
m222 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex, data = piper)
m223 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*MigDate, data = piper)
m224 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species, data = piper)
m225 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season, data = piper)
m226 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*MigDate, data = piper)
m227 <- lm(Mass ~ (OverallNeonic * MigDate) + Species*Season, data = piper)
m228 <- lm(Mass ~ (OverallNeonic * MigDate) + Species*MigDate, data = piper)
m229 <- lm(Mass ~ (OverallNeonic * MigDate) + Season*MigDate, data = piper)

#add two interactions
### neonicotinoid concentrations with two interactions (with restrictions, see above)
#Year
m230 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*Julian, data = piper)
m231 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*MigDate, data = piper)
m232 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Sex*Julian, data = piper)
m233 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Sex*MigDate, data = piper)
m234 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian + Year*Sex, data = piper)
m235 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian + Sex*Julian, data = piper)
m236 <- lm(Mass ~ (OverallNeonic * Year) + Year*MigDate + Sex*MigDate, data = piper)
m237 <- lm(Mass ~ (OverallNeonic * Year) + Sex*Julian + Sex*Species, data = piper)
m238 <- lm(Mass ~ (OverallNeonic * Year) + Sex*Species + Sex*MigDate, data = piper)

#Sex
m239 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*Julian, data = piper)
m240 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*MigDate, data = piper)
m241 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Sex*Julian, data = piper)
m242 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Sex*MigDate, data = piper)
m243 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Year*Sex, data = piper)
m244 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Sex*Julian, data = piper)
m245 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Sex*Species, data = piper)
m246 <- lm(Mass ~ (OverallNeonic * Sex) + Year*MigDate + Sex*MigDate, data = piper)
m247 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*Species, data = piper)
m248 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*MigDate, data = piper)
m249 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Season*Species, data = piper)
m250 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Julian + Sex*Species, data = piper)
m251 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Julian + Species*Julian, data = piper)
m252 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Sex*MigDate, data = piper)
m253 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Season*Species, data = piper)
m254 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Season*MigDate, data = piper)
m255 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Species*Julian, data = piper)
m256 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*MigDate + Season*Species, data = piper)
m257 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*MigDate + Season*MigDate, data = piper)
m258 <- lm(Mass ~ (OverallNeonic * Sex) + Season*Species + Season*MigDate, data = piper)

#Species
m259 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*Species, data = piper)
m260 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*MigDate, data = piper)
m261 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Season*Species, data = piper)
m262 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Julian + Sex*Species, data = piper)
m263 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Julian + Species*Julian, data = piper)
m264 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Sex*MigDate, data = piper)
m265 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Season*Species, data = piper)
m266 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Season*MigDate, data = piper)
m267 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Species*Julian, data = piper)
m268 <- lm(Mass ~ (OverallNeonic * Species) + Sex*MigDate + Season*Species, data = piper)
m269 <- lm(Mass ~ (OverallNeonic * Species) + Sex*MigDate + Season*MigDate, data = piper)
m270 <- lm(Mass ~ (OverallNeonic * Species) + Season*Species + Season*MigDate, data = piper)

#Season
m271 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*Species, data = piper)
m272 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*MigDate, data = piper)
m273 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Season*Species, data = piper)
m274 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species + Sex*MigDate, data = piper)
m275 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species + Season*Species, data = piper)
m276 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species + Season*MigDate, data = piper)
m277 <- lm(Mass ~ (OverallNeonic * Season) + Sex*MigDate + Season*Species, data = piper)
m278 <- lm(Mass ~ (OverallNeonic * Season) + Sex*MigDate + Season*MigDate, data = piper)
m279 <- lm(Mass ~ (OverallNeonic * Season) + Season*Species + Season*MigDate, data = piper)

#Julian
m280 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex + Year*Julian, data = piper)
m281 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex + Sex*Julian, data = piper)
m282 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian + Year*Sex, data = piper)
m283 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian + Sex*Julian, data = piper)
m284 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Julian + Sex*Species, data = piper)
m285 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Julian + Species*Julian, data = piper)
m286 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Species + Season*Species, data = piper)
m287 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Species + Species*Julian, data = piper)

#MigDate
m288 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex + Year*MigDate, data = piper)
m289 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex + Sex*MigDate, data = piper)
m290 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*MigDate + Sex*MigDate, data = piper)
m291 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*Species, data = piper)
m292 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*MigDate, data = piper)
m293 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Season*Species, data = piper)
m294 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species + Sex*MigDate, data = piper)
m295 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species + Season*Species, data = piper)
m296 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species + Season*MigDate, data = piper)
m297 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*MigDate + Season*Species, data = piper)
m298 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*MigDate + Season*MigDate, data = piper)
m299 <- lm(Mass ~ (OverallNeonic * MigDate) + Season*Species + Season*MigDate, data = piper)

#add three interactions
#Year
m300 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*Julian + Sex*Julian, data = piper)
m301 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m302 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian + Year*Sex + Sex*Julian, data = piper)

#Sex
m303 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*Julian + Sex*Julian, data = piper)
m304 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m305 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Year*Sex + Sex*Julian, data = piper)
m306 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m307 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*Species + Season*Species, data = piper)
m308 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m309 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*Species + Species*Julian, data = piper)
m310 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m311 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m312 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m313 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Julian + Sex*Species + Species*Julian, data = piper)
m314 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m315 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m316 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m317 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

#Species
m318 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m319 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*Species + Season*Species, data = piper)
m320 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m321 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*Species + Species*Julian, data = piper)
m322 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m323 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m324 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m325 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Julian + Sex*Species + Species*Julian, data = piper)
m326 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m327 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m328 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m329 <- lm(Mass ~ (OverallNeonic * Species) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

#Season
m330 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m331 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*Species + Season*Species, data = piper)
m332 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m333 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m334 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m335 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m336 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m337 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m338 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m339 <- lm(Mass ~ (OverallNeonic * Season) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

#Julian
m340 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex + Year*Julian + Sex*Julian, data = piper)
m341 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian + Year*Sex + Sex*Julian, data = piper)
m342 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Julian + Sex*Species + Species*Julian, data = piper)

#MigDate
m343 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m344 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*MigDate + Sex*Species + Sex*MigDate, data = piper)
m345 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m346 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*Species + Season*Species, data = piper)
m347 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m348 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m349 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m350 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m351 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m352 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m353 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m354 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

#Manually list the special models
special_models <- list(m.null, m.globalspecies, m.globalinformednull, m.globalyear,
                       m.informednull)

special_names <- c('m.null', 'm.globalspecies', 'm.globalinformednull', 'm.globalyear',
                   'm.informednull')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Dynamically add the other models (m1 to m360) to the list
for (i in c(1:354)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

## Top Model Summaries: Second Stage ####

# Mass ~ Neonic * MigDate + Sex*Species + Sex*MigDate (m294)
model_summary <- summary(m294)$coefficients
conf_intervals <- confint(m294, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Neonic * MigDate + Species (m53)
model_summary <- summary(m53)$coefficients
conf_intervals <- confint(m53, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Neonic * Species + Sex * MigDate (m206)
model_summary <- summary(m206)$coefficients
conf_intervals <- confint(m206, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Neonic * Species + Sex * Species + Sex * MigDate (m264)
model_summary <- summary(m264)$coefficients
conf_intervals <- confint(m264, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Neonic * Species + MigDate (m42)
model_summary <- summary(m42)$coefficients
conf_intervals <- confint(m42, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# ----------------------------------------------------------------------------- #

### Conclusions: Second Stage ####

# 1. Semipalmated sandpipers have significantly higher body mass (m294, m53, m206, m264, m42).
# 2. Males have significantly higher body mass (m294, m206, m264).
# 3. Positive association between date into migration season and body mass (m294, m53, m206, m264, m42).
# 4. Neonicotinoid concentrations are positively associated with body mass (m206, m264, m42).
# 5. Neonicotinoid concentrations vary with date into migration season (m294). The positive coefficient suggests that as migration date increases, the effect of neonicotinoid exposure on body mass also increases. In other words, the impact of later migration dates on body mass is greater at higher levels of neonic exposure. 

# ----------------------------------------------------------------------------- #

# Second Stage AICc Model Selection: Mass ~ Log(Concentration) ####

### null and global models
m.null <- lm(Mass ~ 1, data = piper)
m.informednull <- lm(Mass ~ Species + MigDate, data = piper)

m.globalinformednull <- lm(Mass ~ Species + MigDate + log10(OverallNeonic + 0.0001) +
                             Species * MigDate + Species * log10(OverallNeonic + 0.0001) + 
                             MigDate * log10(OverallNeonic + 0.0001), data = piper)

m.globalspecies <- lm(Mass ~ Sex + Season + Species + MigDate + log10(OverallNeonic + 0.0001) + 
                        Sex*Season + Sex*Species + Sex*MigDate + Sex*log10(OverallNeonic + 0.0001) +
                        Season*Species + Season*MigDate + Season*log10(OverallNeonic + 0.0001) + 
                        Species*MigDate + Species*log10(OverallNeonic + 0.0001) + MigDate*log10(OverallNeonic + 0.0001), data = piper)

m.globalyear <-  lm(Mass ~ Sex + Season + Year + MigDate + log10(OverallNeonic + 0.0001) + 
                      Sex*Season + Sex*Year + Sex*MigDate + Sex*log10(OverallNeonic + 0.0001) +
                      Season*Year + Season*MigDate + Season*log10(OverallNeonic + 0.0001) + 
                      Year*MigDate + Year*log10(OverallNeonic + 0.0001) + MigDate*log10(OverallNeonic + 0.0001), data = piper)

### single combinations of informed and global model covariates with neonicotinoid concentrations
m1 <- lm(Mass ~ log10(OverallNeonic + 0.0001), data = piper)
m2 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year, data = piper)
m3 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex, data = piper)
m4 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species, data = piper)
m5 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season, data = piper)
m6 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Julian, data = piper)
m7 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + MigDate, data = piper)

### two additive combinations of informed and global model covariates with neonicotinoid concentrations
m8 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Sex, data = piper)
m9 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Julian, data = piper)
m10 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + MigDate, data = piper)

m11 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Species, data = piper)
m12 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Season, data = piper)
m13 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Julian, data = piper)
m14 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + MigDate, data = piper)

m15 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species + Season, data = piper)
m16 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species + Julian, data = piper)
m17 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species + MigDate, data = piper)

m18 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season + MigDate, data = piper)

### three additive combinations of informed and global model covariates with neonicotinoid concentrations
m19 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Sex + Julian, data = piper)
m20 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Sex + MigDate, data = piper)

m21 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Species + Season, data = piper)
m22 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Species + Julian, data = piper)
m23 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Species + MigDate, data = piper)

m24 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Season + MigDate, data = piper)

m25 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species + Season + MigDate, data = piper)

### four additive combinations of informed and global model covariates with neonicotinoid concentrations
m26 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season + Sex + Species + MigDate, data = piper)

### two-way interactions (no additive combinations)
m27 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year, data = piper)
m28 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex, data = piper)
m29 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Species, data = piper)
m30 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season, data = piper)
m31 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian, data = piper)
m32 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate, data = piper)

### two-way interactions with a single additive combinations of informed and global model covariates with neonicotinoid concentrations
m33 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year + MigDate, data = piper)
m34 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year + Julian, data = piper)
m35 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year + Sex, data = piper)

m36 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Year, data = piper)
m37 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Season, data = piper)
m38 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + MigDate, data = piper)
m39 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Julian, data = piper)
m40 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Species, data = piper)

m41 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Species + Season, data = piper)
m42 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Species + MigDate, data = piper)
m43 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Species + Julian, data = piper)
m44 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Species + Sex, data = piper)

m45 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season + MigDate, data = piper)
m46 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season + Species, data = piper)
m47 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season + Sex, data = piper)

m48 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian + Year, data = piper)
m49 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian + Species, data = piper)
m50 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian + Sex, data = piper)

m51 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Year, data = piper)
m52 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Season, data = piper)
m53 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Species, data = piper)
m54 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Sex, data = piper)

### two-way interactions with two additive combination of informed and global model covariates with neonicotinoid concentrations
#Year
m55 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * Sex) + Julian, data = piper)
m56 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex + Julian, data = piper)
m57 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year + Julian, data = piper)
m58 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex + Julian, data = piper)

m59 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * Sex) + MigDate, data = piper)
m60 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex + MigDate, data = piper)
m61 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year + MigDate, data = piper)
m62 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Year, data = piper)

m63 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * Julian) + Sex, data = piper)
m64 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Julian + Sex, data = piper)
m65 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year + Sex, data = piper)
m66 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Julian + Year, data = piper)

m67 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * MigDate) + Sex, data = piper)
m68 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + MigDate + Sex, data = piper)
m69 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year + Sex, data = piper)
m70 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Year, data = piper)

#Sex
m71 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Season) + Species, data = piper)
m72 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season + Species, data = piper)
m73 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex + Species, data = piper)
m74 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Season + Sex, data = piper)

m75 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Season) + MigDate, data = piper)
m76 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season + MigDate, data = piper)
m77 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex + MigDate, data = piper)
m78 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season + Sex, data = piper)

m79 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Julian) + Year, data = piper)
m80 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Julian + Year, data = piper)
m81 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex + Year, data = piper)
m82 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Julian + Sex, data = piper)

m83 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Julian) + Species, data = piper)
m84 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Julian + Species, data = piper)
m85 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex + Species, data = piper)
m86 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Julian + Sex, data = piper)

m87 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Species) + Season, data = piper)
m88 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species + Season, data = piper)
m89 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex + Season, data = piper)
m90 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Species + Sex, data = piper)

m91 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Species) + Julian, data = piper)
m92 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species + Julian, data = piper)
m93 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex + Julian, data = piper)
m94 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Species + Sex, data = piper)

m95 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Species) + MigDate, data = piper)
m96 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species + MigDate, data = piper)
m97 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex + MigDate, data = piper)
m98 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Species + Sex, data = piper)

m99 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * MigDate) + Year, data = piper)
m100 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Year, data = piper)
m101 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Year, data = piper)
m102<- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + MigDate + Sex, data = piper)

m103 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * MigDate) + Season, data = piper)
m104 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Season, data = piper)
m105 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Season, data = piper)
m106 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + MigDate + Sex, data = piper)

m107 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * MigDate) + Species, data = piper)
m108 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Species, data = piper)
m109 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Species, data = piper)
m110 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + MigDate + Sex, data = piper)

#Season
m111 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Season * Species) + Sex, data = piper)
m112 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Species + Sex, data = piper)
m113 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Season + Sex, data = piper)
m114 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species + Season, data = piper)

m115 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Season * Species) + MigDate, data = piper)
m116 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Species + MigDate, data = piper)
m117 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Season + MigDate, data = piper)
m118 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Species + Season, data = piper)

m119 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Season * MigDate) + Sex, data = piper)
m120 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + MigDate + Sex, data = piper)
m121 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season + Sex, data = piper)
m122 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Season, data = piper)

m123 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Season * MigDate) + Species, data = piper)
m124 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + MigDate + Species, data = piper)
m125 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season + Species, data = piper)
m126 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + MigDate + Season, data = piper)

#Species)
m127 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Julian + Sex, data = piper)
m128 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Species + Sex, data = piper)
m129 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Julian + Species, data = piper)

### two-way interactions with three additive combinations (with restrictions, see above)
# Year + Sex + Species + Season + Julian + MigDate
m130 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Season) + Species + MigDate, data = piper)
m131 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + (Season * Species) + MigDate, data = piper)
m132 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Season + (Species * MigDate), data = piper)
m133 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season + Species + MigDate, data = piper)
m134 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex + Species + MigDate, data = piper)
m135 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Season + Sex + MigDate, data = piper)
m136 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season + Species + Sex, data = piper)

m137 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * MigDate) + Species + Season, data = piper)
m138 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Species) + Season + MigDate, data = piper)

### neonicotinoid concentrations with one interaction (with restrictions, see above)
m139 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex, data = piper)
m140 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Season, data = piper)
m141 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian, data = piper)
m142 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*MigDate, data = piper)

m143 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species, data = piper)
m144 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season, data = piper)
m145 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Julian, data = piper)
m146 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*MigDate, data = piper)

m147 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species*Season, data = piper)
m148 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species*Julian, data = piper)
m149 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species*MigDate, data = piper)

m150 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season*MigDate, data = piper)

### neonicotinoid concentrations with two interactions (with restrictions, see above)
m151 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*Julian, data = piper)
m152 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*MigDate, data = piper)
m153 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Sex*Julian, data = piper)
m154 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Sex*MigDate, data = piper)

m155 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian + Year*Sex, data = piper)
m156 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian + Sex*Julian, data = piper)

m157 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*MigDate + Sex*MigDate, data = piper)

m158 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*Species, data = piper)
m159 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*MigDate, data = piper)
m160 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Season*Species, data = piper)

m161 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Julian + Sex*Species, data = piper)
m162 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Julian + Species*Julian, data = piper)

m163 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Sex*MigDate, data = piper)
m164 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Season*Species, data = piper)
m165 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Season*MigDate, data = piper)
m166 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Species*Julian, data = piper)

m167 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*MigDate + Season*Species, data = piper)
m168 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*MigDate + Season*MigDate, data = piper)
m169 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season*Species + Season*MigDate, data = piper)

### neonicotinoid concentrations with three interactions (with restrictions, see above)
m170 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*Julian + Sex*Julian, data = piper)
m171 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m172 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian + Year*Sex + Sex*Julian, data = piper)
m173 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m174 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*Species + Season*Species, data = piper)
m175 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m176 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*Species + Species*Julian, data = piper)
m177 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m178 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m179 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m180 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Julian + Sex*Species + Species*Julian, data = piper)
m181 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m182 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m183 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m184 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

### add neonicotinoid interaction with all other interaction combinations
#add single interaction
#Year
m185 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex, data = piper)
m186 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian, data = piper)
m187 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*MigDate, data = piper)
m188 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*Season, data = piper)
m189 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*Julian, data = piper)
m190 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*MigDate, data = piper)
m191 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Season*MigDate, data = piper)

#Sex
m192 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex, data = piper)
m193 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian, data = piper)
m194 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*MigDate, data = piper)
m195 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species, data = piper)
m196 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season, data = piper)
m197 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Julian, data = piper)
m198 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*MigDate, data = piper)
m199 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species*Season, data = piper)
m200 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species*Julian, data = piper)
m201 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species*MigDate, data = piper)
m202 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season*MigDate, data = piper)

#Species
m203 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species, data = piper)
m204 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season, data = piper)
m205 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Julian, data = piper)
m206 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*MigDate, data = piper)
m207 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Species*Season, data = piper)
m208 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Species*Julian, data = piper)
m209 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Species*MigDate, data = piper)
m210 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Season*MigDate, data = piper)

#Season
m211 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species, data = piper)
m212 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season, data = piper)
m213 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*MigDate, data = piper)
m214 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Species*Season, data = piper)
m215 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Species*MigDate, data = piper)
m216 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Season*MigDate, data = piper)

#Julian
m217 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex, data = piper)
m218 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian, data = piper)
m219 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Species, data = piper)
m220 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Julian, data = piper)
m221 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Species*Julian, data = piper)

#MigDate
m222 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex, data = piper)
m223 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*MigDate, data = piper)
m224 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species, data = piper)
m225 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season, data = piper)
m226 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*MigDate, data = piper)
m227 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Species*Season, data = piper)
m228 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Species*MigDate, data = piper)
m229 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season*MigDate, data = piper)

#add two interactions
### neonicotinoid concentrations with two interactions (with restrictions, see above)
#Year
m230 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*Julian, data = piper)
m231 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*MigDate, data = piper)
m232 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Sex*Julian, data = piper)
m233 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Sex*MigDate, data = piper)
m234 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian + Year*Sex, data = piper)
m235 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian + Sex*Julian, data = piper)
m236 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*MigDate + Sex*MigDate, data = piper)
m237 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*Julian + Sex*Species, data = piper)
m238 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*Species + Sex*MigDate, data = piper)

#Sex
m239 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*Julian, data = piper)
m240 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*MigDate, data = piper)
m241 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Sex*Julian, data = piper)
m242 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Sex*MigDate, data = piper)
m243 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Year*Sex, data = piper)
m244 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Sex*Julian, data = piper)
m245 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Sex*Species, data = piper)
m246 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*MigDate + Sex*MigDate, data = piper)
m247 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*Species, data = piper)
m248 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*MigDate, data = piper)
m249 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Season*Species, data = piper)
m250 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Julian + Sex*Species, data = piper)
m251 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Julian + Species*Julian, data = piper)
m252 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Sex*MigDate, data = piper)
m253 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Season*Species, data = piper)
m254 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Season*MigDate, data = piper)
m255 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Species*Julian, data = piper)
m256 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*MigDate + Season*Species, data = piper)
m257 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*MigDate + Season*MigDate, data = piper)
m258 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season*Species + Season*MigDate, data = piper)

#Species
m259 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*Species, data = piper)
m260 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*MigDate, data = piper)
m261 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Season*Species, data = piper)
m262 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Julian + Sex*Species, data = piper)
m263 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Julian + Species*Julian, data = piper)
m264 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Sex*MigDate, data = piper)
m265 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Season*Species, data = piper)
m266 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Season*MigDate, data = piper)
m267 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Species*Julian, data = piper)
m268 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*MigDate + Season*Species, data = piper)
m269 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*MigDate + Season*MigDate, data = piper)
m270 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Season*Species + Season*MigDate, data = piper)

#Season
m271 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*Species, data = piper)
m272 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*MigDate, data = piper)
m273 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Season*Species, data = piper)
m274 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species + Sex*MigDate, data = piper)
m275 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species + Season*Species, data = piper)
m276 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species + Season*MigDate, data = piper)
m277 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*MigDate + Season*Species, data = piper)
m278 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*MigDate + Season*MigDate, data = piper)
m279 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Season*Species + Season*MigDate, data = piper)

#Julian
m280 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex + Year*Julian, data = piper)
m281 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex + Sex*Julian, data = piper)
m282 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian + Year*Sex, data = piper)
m283 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian + Sex*Julian, data = piper)
m284 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Julian + Sex*Species, data = piper)
m285 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Julian + Species*Julian, data = piper)
m286 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Species + Season*Species, data = piper)
m287 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Species + Species*Julian, data = piper)

#MigDate
m288 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex + Year*MigDate, data = piper)
m289 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex + Sex*MigDate, data = piper)
m290 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*MigDate + Sex*MigDate, data = piper)
m291 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*Species, data = piper)
m292 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*MigDate, data = piper)
m293 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Season*Species, data = piper)
m294 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species + Sex*MigDate, data = piper)
m295 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species + Season*Species, data = piper)
m296 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species + Season*MigDate, data = piper)
m297 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*MigDate + Season*Species, data = piper)
m298 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*MigDate + Season*MigDate, data = piper)
m299 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season*Species + Season*MigDate, data = piper)

#add three interactions
#Year
m300 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*Julian + Sex*Julian, data = piper)
m301 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m302 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian + Year*Sex + Sex*Julian, data = piper)

#Sex
m303 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*Julian + Sex*Julian, data = piper)
m304 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m305 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Year*Sex + Sex*Julian, data = piper)
m306 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m307 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*Species + Season*Species, data = piper)
m308 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m309 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*Species + Species*Julian, data = piper)
m310 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m311 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m312 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m313 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Julian + Sex*Species + Species*Julian, data = piper)
m314 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m315 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m316 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m317 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

#Species
m318 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m319 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*Species + Season*Species, data = piper)
m320 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m321 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*Species + Species*Julian, data = piper)
m322 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m323 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m324 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m325 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Julian + Sex*Species + Species*Julian, data = piper)
m326 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m327 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m328 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m329 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

#Season
m330 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m331 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*Species + Season*Species, data = piper)
m332 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m333 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m334 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m335 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m336 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m337 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m338 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m339 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

#Julian
m340 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex + Year*Julian + Sex*Julian, data = piper)
m341 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian + Year*Sex + Sex*Julian, data = piper)
m342 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Julian + Sex*Species + Species*Julian, data = piper)

#MigDate
m343 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m344 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*MigDate + Sex*Species + Sex*MigDate, data = piper)
m345 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m346 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*Species + Season*Species, data = piper)
m347 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m348 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m349 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m350 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m351 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m352 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m353 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m354 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

#Manually list the special models
special_models <- list(m.null, m.globalspecies, m.globalinformednull, m.globalyear,
                       m.informednull)

special_names <- c('m.null', 'm.globalspecies', 'm.globalinformednull', 'm.globalyear',
                   'm.informednull')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Dynamically add the other models (m1 to m360) to the list
for (i in c(1:354)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

### Top Model Summaries: Second Stage ####

options(digits = 3)

# Mass ~ Log(Neonic) + Sex * MigDate + Species (m107)
cbind(summary(m107)$coefficients, confint(m107))

# Mass ~ Log(Neonic) + Sex * Species + Sex * MigDate (m163)
model_summary <- summary(m163)$coefficients
conf_intervals <- confint(m163, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Log(Neonic) + Sex * Season + Species + MigDate (m130)
cbind(summary(m130)$coefficients, confint(m130))

# Mass ~ Log(Neonic) + Sex * Season + Sex * Species + Season * MigDate (m175)
cbind(summary(m175)$coefficients, confint(m175))
model_summary <- summary(m175)$coefficients
conf_intervals <- confint(m175, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Log(Neonic) + Sex * Julian + Species (m83)
model_summary <- summary(m83)$coefficients
conf_intervals <- confint(m83, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Log(Neonic) + Sex * Julian + Sex * Species (m161)
model_summary <- summary(m161)$coefficients
conf_intervals <- confint(m161, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Log(Neonic) + Species + MigDate (m17)
model_summary <- summary(m17)$coefficients
conf_intervals <- confint(m17, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Log(Neonic) * MigDate + Sex * Species + Sex * MigDate (m294)
model_summary <- summary(m294)$coefficients
conf_intervals <- confint(m294, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# ----------------------------------------------------------------------------- #

### Conclusion: Second Stage with Log Transformation ####

# Main conclusion: body mass is positively associated with log(neonicotinoid).
# Notes: There are only 11 birds with detections...

# ----------------------------------------------------------------------------- #

# Second Stage AICc Model Selection: Mass ~ Detection ####

### null and global models
m.null <- lm(Mass ~ 1, data = piper)
m.informednull <- lm(Mass ~ Species + MigDate, data = piper)

m.globalinformednull <- lm(Mass ~ Species + MigDate + Detection +
                             Species * MigDate + Species * Detection + 
                             MigDate * Detection, data = piper)

m.globalspecies <- lm(Mass ~ Sex + Season + Species + MigDate + Detection + 
                        Sex*Season + Sex*Species + Sex*MigDate + Sex*Detection +
                        Season*Species + Season*MigDate + Season*Detection + 
                        Species*MigDate + Species*Detection + MigDate*Detection, data = piper)

m.globalyear <-  lm(Mass ~ Sex + Season + Year + MigDate + Detection + 
                      Sex*Season + Sex*Year + Sex*MigDate + Sex*Detection +
                      Season*Year + Season*MigDate + Season*Detection + 
                      Year*MigDate + Year*Detection + MigDate*Detection, data = piper)

### single combinations of informed and global model covariates with neonicotinoid concentrations
m1 <- lm(Mass ~ Detection, data = piper)
m2 <- lm(Mass ~ Detection + Year, data = piper)
m3 <- lm(Mass ~ Detection + Sex, data = piper)
m4 <- lm(Mass ~ Detection + Species, data = piper)
m5 <- lm(Mass ~ Detection + Season, data = piper)
m6 <- lm(Mass ~ Detection + Julian, data = piper)
m7 <- lm(Mass ~ Detection + MigDate, data = piper)

### two additive combinations of informed and global model covariates with neonicotinoid concentrations
m8 <- lm(Mass ~ Detection + Year + Sex, data = piper)
m9 <- lm(Mass ~ Detection + Year + Julian, data = piper)
m10 <- lm(Mass ~ Detection + Year + MigDate, data = piper)

m11 <- lm(Mass ~ Detection + Sex + Species, data = piper)
m12 <- lm(Mass ~ Detection + Sex + Season, data = piper)
m13 <- lm(Mass ~ Detection + Sex + Julian, data = piper)
m14 <- lm(Mass ~ Detection + Sex + MigDate, data = piper)

m15 <- lm(Mass ~ Detection + Species + Season, data = piper)
m16 <- lm(Mass ~ Detection + Species + Julian, data = piper)
m17 <- lm(Mass ~ Detection + Species + MigDate, data = piper)

m18 <- lm(Mass ~ Detection + Season + MigDate, data = piper)

### three additive combinations of informed and global model covariates with neonicotinoid concentrations
m19 <- lm(Mass ~ Detection + Year + Sex + Julian, data = piper)
m20 <- lm(Mass ~ Detection + Year + Sex + MigDate, data = piper)

m21 <- lm(Mass ~ Detection + Sex + Species + Season, data = piper)
m22 <- lm(Mass ~ Detection + Sex + Species + Julian, data = piper)
m23 <- lm(Mass ~ Detection + Sex + Species + MigDate, data = piper)

m24 <- lm(Mass ~ Detection + Sex + Season + MigDate, data = piper)

m25 <- lm(Mass ~ Detection + Species + Season + MigDate, data = piper)

### four additive combinations of informed and global model covariates with neonicotinoid concentrations
m26 <- lm(Mass ~ Detection + Season + Sex + Species + MigDate, data = piper)

### two-way interactions (no additive combinations)
m27 <- lm(Mass ~ Detection*Year, data = piper)
m28 <- lm(Mass ~ Detection*Sex, data = piper)
m29 <- lm(Mass ~ Detection*Species, data = piper)
m30 <- lm(Mass ~ Detection*Season, data = piper)
m31 <- lm(Mass ~ Detection*Julian, data = piper)
m32 <- lm(Mass ~ Detection*MigDate, data = piper)

### two-way interactions with a single additive combinations of informed and global model covariates with neonicotinoid concentrations
m33 <- lm(Mass ~ Detection*Year + MigDate, data = piper)
m34 <- lm(Mass ~ Detection*Year + Julian, data = piper)
m35 <- lm(Mass ~ Detection*Year + Sex, data = piper)

m36 <- lm(Mass ~ Detection*Sex + Year, data = piper)
m37 <- lm(Mass ~ Detection*Sex + Season, data = piper)
m38 <- lm(Mass ~ Detection*Sex + MigDate, data = piper)
m39 <- lm(Mass ~ Detection*Sex + Julian, data = piper)
m40 <- lm(Mass ~ Detection*Sex + Species, data = piper)

m41 <- lm(Mass ~ Detection*Species + Season, data = piper)
m42 <- lm(Mass ~ Detection*Species + MigDate, data = piper)
m43 <- lm(Mass ~ Detection*Species + Julian, data = piper)
m44 <- lm(Mass ~ Detection*Species + Sex, data = piper)

m45 <- lm(Mass ~ Detection*Season + MigDate, data = piper)
m46 <- lm(Mass ~ Detection*Season + Species, data = piper)
m47 <- lm(Mass ~ Detection*Season + Sex, data = piper)

m48 <- lm(Mass ~ Detection*Julian + Year, data = piper)
m49 <- lm(Mass ~ Detection*Julian + Species, data = piper)
m50 <- lm(Mass ~ Detection*Julian + Sex, data = piper)

m51 <- lm(Mass ~ Detection*MigDate + Year, data = piper)
m52 <- lm(Mass ~ Detection*MigDate + Season, data = piper)
m53 <- lm(Mass ~ Detection*MigDate + Species, data = piper)
m54 <- lm(Mass ~ Detection*MigDate + Sex, data = piper)

### two-way interactions with two additive combination of informed and global model covariates with neonicotinoid concentrations
#Year
m55 <- lm(Mass ~ Detection + (Year * Sex) + Julian, data = piper)
m56 <- lm(Mass ~ (Detection * Year) + Sex + Julian, data = piper)
m57 <- lm(Mass ~ (Detection * Sex) + Year + Julian, data = piper)
m58 <- lm(Mass ~ (Detection * Julian) + Sex + Julian, data = piper)

m59 <- lm(Mass ~ Detection + (Year * Sex) + MigDate, data = piper)
m60 <- lm(Mass ~ (Detection * Year) + Sex + MigDate, data = piper)
m61 <- lm(Mass ~ (Detection * Sex) + Year + MigDate, data = piper)
m62 <- lm(Mass ~ (Detection * MigDate) + Sex + Year, data = piper)

m63 <- lm(Mass ~ Detection + (Year * Julian) + Sex, data = piper)
m64 <- lm(Mass ~ (Detection * Year) + Julian + Sex, data = piper)
m65 <- lm(Mass ~ (Detection * Julian) + Year + Sex, data = piper)
m66 <- lm(Mass ~ (Detection * Sex) + Julian + Year, data = piper)

m67 <- lm(Mass ~ Detection + (Year * MigDate) + Sex, data = piper)
m68 <- lm(Mass ~ (Detection * Year) + MigDate + Sex, data = piper)
m69 <- lm(Mass ~ (Detection * MigDate) + Year + Sex, data = piper)
m70 <- lm(Mass ~ (Detection * Sex) + MigDate + Year, data = piper)

#Sex
m71 <- lm(Mass ~ Detection + (Sex * Season) + Species, data = piper)
m72 <- lm(Mass ~ (Detection * Sex) + Season + Species, data = piper)
m73 <- lm(Mass ~ (Detection * Season) + Sex + Species, data = piper)
m74 <- lm(Mass ~ (Detection * Species) + Season + Sex, data = piper)

m75 <- lm(Mass ~ Detection + (Sex * Season) + MigDate, data = piper)
m76 <- lm(Mass ~ (Detection * Sex) + Season + MigDate, data = piper)
m77 <- lm(Mass ~ (Detection * Season) + Sex + MigDate, data = piper)
m78 <- lm(Mass ~ (Detection * MigDate) + Season + Sex, data = piper)

m79 <- lm(Mass ~ Detection + (Sex * Julian) + Year, data = piper)
m80 <- lm(Mass ~ (Detection * Sex) + Julian + Year, data = piper)
m81 <- lm(Mass ~ (Detection * Julian) + Sex + Year, data = piper)
m82 <- lm(Mass ~ (Detection * Year) + Julian + Sex, data = piper)

m83 <- lm(Mass ~ Detection + (Sex * Julian) + Species, data = piper)
m84 <- lm(Mass ~ (Detection * Sex) + Julian + Species, data = piper)
m85 <- lm(Mass ~ (Detection * Julian) + Sex + Species, data = piper)
m86 <- lm(Mass ~ (Detection * Species) + Julian + Sex, data = piper)

m87 <- lm(Mass ~ Detection + (Sex * Species) + Season, data = piper)
m88 <- lm(Mass ~ (Detection * Sex) + Species + Season, data = piper)
m89 <- lm(Mass ~ (Detection * Species) + Sex + Season, data = piper)
m90 <- lm(Mass ~ (Detection * Season) + Species + Sex, data = piper)

m91 <- lm(Mass ~ Detection + (Sex * Species) + Julian, data = piper)
m92 <- lm(Mass ~ (Detection * Sex) + Species + Julian, data = piper)
m93 <- lm(Mass ~ (Detection * Species) + Sex + Julian, data = piper)
m94 <- lm(Mass ~ (Detection * Julian) + Species + Sex, data = piper)

m95 <- lm(Mass ~ Detection + (Sex * Species) + MigDate, data = piper)
m96 <- lm(Mass ~ (Detection * Sex) + Species + MigDate, data = piper)
m97 <- lm(Mass ~ (Detection * Species) + Sex + MigDate, data = piper)
m98 <- lm(Mass ~ (Detection * MigDate) + Species + Sex, data = piper)

m99 <- lm(Mass ~ Detection + (Sex * MigDate) + Year, data = piper)
m100 <- lm(Mass ~ (Detection * Sex) + MigDate + Year, data = piper)
m101 <- lm(Mass ~ (Detection * MigDate) + Sex + Year, data = piper)
m102<- lm(Mass ~ (Detection * Year) + MigDate + Sex, data = piper)

m103 <- lm(Mass ~ Detection + (Sex * MigDate) + Season, data = piper)
m104 <- lm(Mass ~ (Detection * Sex) + MigDate + Season, data = piper)
m105 <- lm(Mass ~ (Detection * MigDate) + Sex + Season, data = piper)
m106 <- lm(Mass ~ (Detection * Season) + MigDate + Sex, data = piper)

m107 <- lm(Mass ~ Detection + (Sex * MigDate) + Species, data = piper)
m108 <- lm(Mass ~ (Detection * Sex) + MigDate + Species, data = piper)
m109 <- lm(Mass ~ (Detection * MigDate) + Sex + Species, data = piper)
m110 <- lm(Mass ~ (Detection * Species) + MigDate + Sex, data = piper)

#Season
m111 <- lm(Mass ~ Detection + (Season * Species) + Sex, data = piper)
m112 <- lm(Mass ~ (Detection * Season) + Species + Sex, data = piper)
m113 <- lm(Mass ~ (Detection * Species) + Season + Sex, data = piper)
m114 <- lm(Mass ~ (Detection * Sex) + Species + Season, data = piper)

m115 <- lm(Mass ~ Detection + (Season * Species) + MigDate, data = piper)
m116 <- lm(Mass ~ (Detection * Season) + Species + MigDate, data = piper)
m117 <- lm(Mass ~ (Detection * Species) + Season + MigDate, data = piper)
m118 <- lm(Mass ~ (Detection * MigDate) + Species + Season, data = piper)

m119 <- lm(Mass ~ Detection + (Season * MigDate) + Sex, data = piper)
m120 <- lm(Mass ~ (Detection * Season) + MigDate + Sex, data = piper)
m121 <- lm(Mass ~ (Detection * MigDate) + Season + Sex, data = piper)
m122 <- lm(Mass ~ (Detection * Sex) + MigDate + Season, data = piper)

m123 <- lm(Mass ~ Detection + (Season * MigDate) + Species, data = piper)
m124 <- lm(Mass ~ (Detection * Season) + MigDate + Species, data = piper)
m125 <- lm(Mass ~ (Detection * MigDate) + Season + Species, data = piper)
m126 <- lm(Mass ~ (Detection * Species) + MigDate + Season, data = piper)

#Species)
m127 <- lm(Mass ~ (Detection * Species) + Julian + Sex, data = piper)
m128 <- lm(Mass ~ (Detection * Julian) + Species + Sex, data = piper)
m129 <- lm(Mass ~ (Detection * Sex) + Julian + Species, data = piper)

### two-way interactions with three additive combinations (with restrictions, see above)
# Year + Sex + Species + Season + Julian + MigDate
m130 <- lm(Mass ~ Detection + (Sex * Season) + Species + MigDate, data = piper)
m131 <- lm(Mass ~ Detection + Sex + (Season * Species) + MigDate, data = piper)
m132 <- lm(Mass ~ Detection + Sex + Season + (Species * MigDate), data = piper)
m133 <- lm(Mass ~ (Detection * Sex) + Season + Species + MigDate, data = piper)
m134 <- lm(Mass ~ (Detection * Season) + Sex + Species + MigDate, data = piper)
m135 <- lm(Mass ~ (Detection * Species) + Season + Sex + MigDate, data = piper)
m136 <- lm(Mass ~ (Detection * MigDate) + Season + Species + Sex, data = piper)

m137 <- lm(Mass ~ Detection + (Sex * MigDate) + Species + Season, data = piper)
m138 <- lm(Mass ~ Detection + (Sex * Species) + Season + MigDate, data = piper)

### neonicotinoid concentrations with one interaction (with restrictions, see above)
m139 <- lm(Mass ~ Detection + Year*Sex, data = piper)
m140 <- lm(Mass ~ Detection + Year*Season, data = piper)
m141 <- lm(Mass ~ Detection + Year*Julian, data = piper)
m142 <- lm(Mass ~ Detection + Year*MigDate, data = piper)

m143 <- lm(Mass ~ Detection + Sex*Species, data = piper)
m144 <- lm(Mass ~ Detection + Sex*Season, data = piper)
m145 <- lm(Mass ~ Detection + Sex*Julian, data = piper)
m146 <- lm(Mass ~ Detection + Sex*MigDate, data = piper)

m147 <- lm(Mass ~ Detection + Species*Season, data = piper)
m148 <- lm(Mass ~ Detection + Species*Julian, data = piper)
m149 <- lm(Mass ~ Detection + Species*MigDate, data = piper)

m150 <- lm(Mass ~ Detection + Season*MigDate, data = piper)

### neonicotinoid concentrations with two interactions (with restrictions, see above)
m151 <- lm(Mass ~ Detection + Year*Sex + Year*Julian, data = piper)
m152 <- lm(Mass ~ Detection + Year*Sex + Year*MigDate, data = piper)
m153 <- lm(Mass ~ Detection + Year*Sex + Sex*Julian, data = piper)
m154 <- lm(Mass ~ Detection + Year*Sex + Sex*MigDate, data = piper)

m155 <- lm(Mass ~ Detection + Year*Julian + Year*Sex, data = piper)
m156 <- lm(Mass ~ Detection + Year*Julian + Sex*Julian, data = piper)

m157 <- lm(Mass ~ Detection + Year*MigDate + Sex*MigDate, data = piper)

m158 <- lm(Mass ~ Detection + Sex*Season + Sex*Species, data = piper)
m159 <- lm(Mass ~ Detection + Sex*Season + Sex*MigDate, data = piper)
m160 <- lm(Mass ~ Detection + Sex*Season + Season*Species, data = piper)

m161 <- lm(Mass ~ Detection + Sex*Julian + Sex*Species, data = piper)
m162 <- lm(Mass ~ Detection + Sex*Julian + Species*Julian, data = piper)

m163 <- lm(Mass ~ Detection + Sex*Species + Sex*MigDate, data = piper)
m164 <- lm(Mass ~ Detection + Sex*Species + Season*Species, data = piper)
m165 <- lm(Mass ~ Detection + Sex*Species + Season*MigDate, data = piper)
m166 <- lm(Mass ~ Detection + Sex*Species + Species*Julian, data = piper)

m167 <- lm(Mass ~ Detection + Sex*MigDate + Season*Species, data = piper)
m168 <- lm(Mass ~ Detection + Sex*MigDate + Season*MigDate, data = piper)
m169 <- lm(Mass ~ Detection + Season*Species + Season*MigDate, data = piper)

### neonicotinoid concentrations with three interactions (with restrictions, see above)
m170 <- lm(Mass ~ Detection + Year*Sex + Year*Julian + Sex*Julian, data = piper)
m171 <- lm(Mass ~ Detection + Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m172 <- lm(Mass ~ Detection + Year*Julian + Year*Sex + Sex*Julian, data = piper)
m173 <- lm(Mass ~ Detection + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m174 <- lm(Mass ~ Detection + Sex*Season + Sex*Species + Season*Species, data = piper)
m175 <- lm(Mass ~ Detection + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m176 <- lm(Mass ~ Detection + Sex*Season + Sex*Species + Species*Julian, data = piper)
m177 <- lm(Mass ~ Detection + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m178 <- lm(Mass ~ Detection + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m179 <- lm(Mass ~ Detection + Sex*Season + Season*Species + Season*MigDate, data = piper)
m180 <- lm(Mass ~ Detection + Sex*Julian + Sex*Species + Species*Julian, data = piper)
m181 <- lm(Mass ~ Detection + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m182 <- lm(Mass ~ Detection + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m183 <- lm(Mass ~ Detection + Sex*Species + Season*Species + Season*MigDate, data = piper)
m184 <- lm(Mass ~ Detection + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

### add neonicotinoid interaction with all other interaction combinations
#add single interaction
#Year
m185 <- lm(Mass ~ (Detection * Year) + Year*Sex, data = piper)
m186 <- lm(Mass ~ (Detection * Year) + Year*Julian, data = piper)
m187 <- lm(Mass ~ (Detection * Year) + Year*MigDate, data = piper)
m188 <- lm(Mass ~ (Detection * Year) + Sex*Season, data = piper)
m189 <- lm(Mass ~ (Detection * Year) + Sex*Julian, data = piper)
m190 <- lm(Mass ~ (Detection * Year) + Sex*MigDate, data = piper)
m191 <- lm(Mass ~ (Detection * Year) + Season*MigDate, data = piper)

#Sex
m192 <- lm(Mass ~ (Detection * Sex) + Year*Sex, data = piper)
m193 <- lm(Mass ~ (Detection * Sex) + Year*Julian, data = piper)
m194 <- lm(Mass ~ (Detection * Sex) + Year*MigDate, data = piper)
m195 <- lm(Mass ~ (Detection * Sex) + Sex*Species, data = piper)
m196 <- lm(Mass ~ (Detection * Sex) + Sex*Season, data = piper)
m197 <- lm(Mass ~ (Detection * Sex) + Sex*Julian, data = piper)
m198 <- lm(Mass ~ (Detection * Sex) + Sex*MigDate, data = piper)
m199 <- lm(Mass ~ (Detection * Sex) + Species*Season, data = piper)
m200 <- lm(Mass ~ (Detection * Sex) + Species*Julian, data = piper)
m201 <- lm(Mass ~ (Detection * Sex) + Species*MigDate, data = piper)
m202 <- lm(Mass ~ (Detection * Sex) + Season*MigDate, data = piper)

#Species
m203 <- lm(Mass ~ (Detection * Species) + Sex*Species, data = piper)
m204 <- lm(Mass ~ (Detection * Species) + Sex*Season, data = piper)
m205 <- lm(Mass ~ (Detection * Species) + Sex*Julian, data = piper)
m206 <- lm(Mass ~ (Detection * Species) + Sex*MigDate, data = piper)
m207 <- lm(Mass ~ (Detection * Species) + Species*Season, data = piper)
m208 <- lm(Mass ~ (Detection * Species) + Species*Julian, data = piper)
m209 <- lm(Mass ~ (Detection * Species) + Species*MigDate, data = piper)
m210 <- lm(Mass ~ (Detection * Species) + Season*MigDate, data = piper)

#Season
m211 <- lm(Mass ~ (Detection * Season) + Sex*Species, data = piper)
m212 <- lm(Mass ~ (Detection * Season) + Sex*Season, data = piper)
m213 <- lm(Mass ~ (Detection * Season) + Sex*MigDate, data = piper)
m214 <- lm(Mass ~ (Detection * Season) + Species*Season, data = piper)
m215 <- lm(Mass ~ (Detection * Season) + Species*MigDate, data = piper)
m216 <- lm(Mass ~ (Detection * Season) + Season*MigDate, data = piper)

#Julian
m217 <- lm(Mass ~ (Detection * Julian) + Year*Sex, data = piper)
m218 <- lm(Mass ~ (Detection * Julian) + Year*Julian, data = piper)
m219 <- lm(Mass ~ (Detection * Julian) + Sex*Species, data = piper)
m220 <- lm(Mass ~ (Detection * Julian) + Sex*Julian, data = piper)
m221 <- lm(Mass ~ (Detection * Julian) + Species*Julian, data = piper)

#MigDate
m222 <- lm(Mass ~ (Detection * MigDate) + Year*Sex, data = piper)
m223 <- lm(Mass ~ (Detection * MigDate) + Year*MigDate, data = piper)
m224 <- lm(Mass ~ (Detection * MigDate) + Sex*Species, data = piper)
m225 <- lm(Mass ~ (Detection * MigDate) + Sex*Season, data = piper)
m226 <- lm(Mass ~ (Detection * MigDate) + Sex*MigDate, data = piper)
m227 <- lm(Mass ~ (Detection * MigDate) + Species*Season, data = piper)
m228 <- lm(Mass ~ (Detection * MigDate) + Species*MigDate, data = piper)
m229 <- lm(Mass ~ (Detection * MigDate) + Season*MigDate, data = piper)

#add two interactions
### neonicotinoid concentrations with two interactions (with restrictions, see above)
#Year
m230 <- lm(Mass ~ (Detection * Year) + Year*Sex + Year*Julian, data = piper)
m231 <- lm(Mass ~ (Detection * Year) + Year*Sex + Year*MigDate, data = piper)
m232 <- lm(Mass ~ (Detection * Year) + Year*Sex + Sex*Julian, data = piper)
m233 <- lm(Mass ~ (Detection * Year) + Year*Sex + Sex*MigDate, data = piper)
m234 <- lm(Mass ~ (Detection * Year) + Year*Julian + Year*Sex, data = piper)
m235 <- lm(Mass ~ (Detection * Year) + Year*Julian + Sex*Julian, data = piper)
m236 <- lm(Mass ~ (Detection * Year) + Year*MigDate + Sex*MigDate, data = piper)
m237 <- lm(Mass ~ (Detection * Year) + Sex*Julian + Sex*Species, data = piper)
m238 <- lm(Mass ~ (Detection * Year) + Sex*Species + Sex*MigDate, data = piper)

#Sex
m239 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Year*Julian, data = piper)
m240 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Year*MigDate, data = piper)
m241 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Sex*Julian, data = piper)
m242 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Sex*MigDate, data = piper)
m243 <- lm(Mass ~ (Detection * Sex) + Year*Julian + Year*Sex, data = piper)
m244 <- lm(Mass ~ (Detection * Sex) + Year*Julian + Sex*Julian, data = piper)
m245 <- lm(Mass ~ (Detection * Sex) + Year*Julian + Sex*Species, data = piper)
m246 <- lm(Mass ~ (Detection * Sex) + Year*MigDate + Sex*MigDate, data = piper)
m247 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Sex*Species, data = piper)
m248 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Sex*MigDate, data = piper)
m249 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Season*Species, data = piper)
m250 <- lm(Mass ~ (Detection * Sex) + Sex*Julian + Sex*Species, data = piper)
m251 <- lm(Mass ~ (Detection * Sex) + Sex*Julian + Species*Julian, data = piper)
m252 <- lm(Mass ~ (Detection * Sex) + Sex*Species + Sex*MigDate, data = piper)
m253 <- lm(Mass ~ (Detection * Sex) + Sex*Species + Season*Species, data = piper)
m254 <- lm(Mass ~ (Detection * Sex) + Sex*Species + Season*MigDate, data = piper)
m255 <- lm(Mass ~ (Detection * Sex) + Sex*Species + Species*Julian, data = piper)
m256 <- lm(Mass ~ (Detection * Sex) + Sex*MigDate + Season*Species, data = piper)
m257 <- lm(Mass ~ (Detection * Sex) + Sex*MigDate + Season*MigDate, data = piper)
m258 <- lm(Mass ~ (Detection * Sex) + Season*Species + Season*MigDate, data = piper)

#Species
m259 <- lm(Mass ~ (Detection * Species) + Sex*Season + Sex*Species, data = piper)
m260 <- lm(Mass ~ (Detection * Species) + Sex*Season + Sex*MigDate, data = piper)
m261 <- lm(Mass ~ (Detection * Species) + Sex*Season + Season*Species, data = piper)
m262 <- lm(Mass ~ (Detection * Species) + Sex*Julian + Sex*Species, data = piper)
m263 <- lm(Mass ~ (Detection * Species) + Sex*Julian + Species*Julian, data = piper)
m264 <- lm(Mass ~ (Detection * Species) + Sex*Species + Sex*MigDate, data = piper)
m265 <- lm(Mass ~ (Detection * Species) + Sex*Species + Season*Species, data = piper)
m266 <- lm(Mass ~ (Detection * Species) + Sex*Species + Season*MigDate, data = piper)
m267 <- lm(Mass ~ (Detection * Species) + Sex*Species + Species*Julian, data = piper)
m268 <- lm(Mass ~ (Detection * Species) + Sex*MigDate + Season*Species, data = piper)
m269 <- lm(Mass ~ (Detection * Species) + Sex*MigDate + Season*MigDate, data = piper)
m270 <- lm(Mass ~ (Detection * Species) + Season*Species + Season*MigDate, data = piper)

#Season
m271 <- lm(Mass ~ (Detection * Season) + Sex*Season + Sex*Species, data = piper)
m272 <- lm(Mass ~ (Detection * Season) + Sex*Season + Sex*MigDate, data = piper)
m273 <- lm(Mass ~ (Detection * Season) + Sex*Season + Season*Species, data = piper)
m274 <- lm(Mass ~ (Detection * Season) + Sex*Species + Sex*MigDate, data = piper)
m275 <- lm(Mass ~ (Detection * Season) + Sex*Species + Season*Species, data = piper)
m276 <- lm(Mass ~ (Detection * Season) + Sex*Species + Season*MigDate, data = piper)
m277 <- lm(Mass ~ (Detection * Season) + Sex*MigDate + Season*Species, data = piper)
m278 <- lm(Mass ~ (Detection * Season) + Sex*MigDate + Season*MigDate, data = piper)
m279 <- lm(Mass ~ (Detection * Season) + Season*Species + Season*MigDate, data = piper)

#Julian
m280 <- lm(Mass ~ (Detection * Julian) + Year*Sex + Year*Julian, data = piper)
m281 <- lm(Mass ~ (Detection * Julian) + Year*Sex + Sex*Julian, data = piper)
m282 <- lm(Mass ~ (Detection * Julian) + Year*Julian + Year*Sex, data = piper)
m283 <- lm(Mass ~ (Detection * Julian) + Year*Julian + Sex*Julian, data = piper)
m284 <- lm(Mass ~ (Detection * Julian) + Sex*Julian + Sex*Species, data = piper)
m285 <- lm(Mass ~ (Detection * Julian) + Sex*Julian + Species*Julian, data = piper)
m286 <- lm(Mass ~ (Detection * Julian) + Sex*Species + Season*Species, data = piper)
m287 <- lm(Mass ~ (Detection * Julian) + Sex*Species + Species*Julian, data = piper)

#MigDate
m288 <- lm(Mass ~ (Detection * MigDate) + Year*Sex + Year*MigDate, data = piper)
m289 <- lm(Mass ~ (Detection * MigDate) + Year*Sex + Sex*MigDate, data = piper)
m290 <- lm(Mass ~ (Detection * MigDate) + Year*MigDate + Sex*MigDate, data = piper)
m291 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Sex*Species, data = piper)
m292 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Sex*MigDate, data = piper)
m293 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Season*Species, data = piper)
m294 <- lm(Mass ~ (Detection * MigDate) + Sex*Species + Sex*MigDate, data = piper)
m295 <- lm(Mass ~ (Detection * MigDate) + Sex*Species + Season*Species, data = piper)
m296 <- lm(Mass ~ (Detection * MigDate) + Sex*Species + Season*MigDate, data = piper)
m297 <- lm(Mass ~ (Detection * MigDate) + Sex*MigDate + Season*Species, data = piper)
m298 <- lm(Mass ~ (Detection * MigDate) + Sex*MigDate + Season*MigDate, data = piper)
m299 <- lm(Mass ~ (Detection * MigDate) + Season*Species + Season*MigDate, data = piper)

#add three interactions
#Year
m300 <- lm(Mass ~ (Detection * Year) + Year*Sex + Year*Julian + Sex*Julian, data = piper)
m301 <- lm(Mass ~ (Detection * Year) + Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m302 <- lm(Mass ~ (Detection * Year) + Year*Julian + Year*Sex + Sex*Julian, data = piper)

#Sex
m303 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Year*Julian + Sex*Julian, data = piper)
m304 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m305 <- lm(Mass ~ (Detection * Sex) + Year*Julian + Year*Sex + Sex*Julian, data = piper)
m306 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m307 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Sex*Species + Season*Species, data = piper)
m308 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m309 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Sex*Species + Species*Julian, data = piper)
m310 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m311 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m312 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m313 <- lm(Mass ~ (Detection * Sex) + Sex*Julian + Sex*Species + Species*Julian, data = piper)
m314 <- lm(Mass ~ (Detection * Sex) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m315 <- lm(Mass ~ (Detection * Sex) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m316 <- lm(Mass ~ (Detection * Sex) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m317 <- lm(Mass ~ (Detection * Sex) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

#Species
m318 <- lm(Mass ~ (Detection * Species) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m319 <- lm(Mass ~ (Detection * Species) + Sex*Season + Sex*Species + Season*Species, data = piper)
m320 <- lm(Mass ~ (Detection * Species) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m321 <- lm(Mass ~ (Detection * Species) + Sex*Season + Sex*Species + Species*Julian, data = piper)
m322 <- lm(Mass ~ (Detection * Species) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m323 <- lm(Mass ~ (Detection * Species) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m324 <- lm(Mass ~ (Detection * Species) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m325 <- lm(Mass ~ (Detection * Species) + Sex*Julian + Sex*Species + Species*Julian, data = piper)
m326 <- lm(Mass ~ (Detection * Species) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m327 <- lm(Mass ~ (Detection * Species) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m328 <- lm(Mass ~ (Detection * Species) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m329 <- lm(Mass ~ (Detection * Species) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

#Season
m330 <- lm(Mass ~ (Detection * Season) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m331 <- lm(Mass ~ (Detection * Season) + Sex*Season + Sex*Species + Season*Species, data = piper)
m332 <- lm(Mass ~ (Detection * Season) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m333 <- lm(Mass ~ (Detection * Season) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m334 <- lm(Mass ~ (Detection * Season) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m335 <- lm(Mass ~ (Detection * Season) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m336 <- lm(Mass ~ (Detection * Season) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m337 <- lm(Mass ~ (Detection * Season) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m338 <- lm(Mass ~ (Detection * Season) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m339 <- lm(Mass ~ (Detection * Season) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

#Julian
m340 <- lm(Mass ~ (Detection * Julian) + Year*Sex + Year*Julian + Sex*Julian, data = piper)
m341 <- lm(Mass ~ (Detection * Julian) + Year*Julian + Year*Sex + Sex*Julian, data = piper)
m342 <- lm(Mass ~ (Detection * Julian) + Sex*Julian + Sex*Species + Species*Julian, data = piper)

#MigDate
m343 <- lm(Mass ~ (Detection * MigDate) + Year*Sex + Year*MigDate + Sex*MigDate, data = piper)
m344 <- lm(Mass ~ (Detection * MigDate) + Year*MigDate + Sex*Species + Sex*MigDate, data = piper)
m345 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Sex*Species + Sex*MigDate, data = piper)
m346 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Sex*Species + Season*Species, data = piper)
m347 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Sex*Species + Season*MigDate, data = piper)
m348 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Sex*MigDate + Season*Species, data = piper)
m349 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Sex*MigDate + Season*MigDate, data = piper)
m350 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Season*Species + Season*MigDate, data = piper)
m351 <- lm(Mass ~ (Detection * MigDate) + Sex*Species + Sex*MigDate + Season*Species, data = piper)
m352 <- lm(Mass ~ (Detection * MigDate) + Sex*Species + Sex*MigDate + Season*MigDate, data = piper)
m353 <- lm(Mass ~ (Detection * MigDate) + Sex*Species + Season*Species + Season*MigDate, data = piper)
m354 <- lm(Mass ~ (Detection * MigDate) + Sex*MigDate + Season*Species + Season*MigDate, data = piper)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

#Manually list the special models
special_models <- list(m.null, m.globalspecies, m.globalinformednull, m.globalyear,
                       m.informednull)

special_names <- c('m.null', 'm.globalspecies', 'm.globalinformednull', 'm.globalyear',
                   'm.informednull')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Dynamically add the other models (m1 to m360) to the list
for (i in c(1:354)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

### Top Model Summaries: Second Stage ####

options(digits = 3)

# Mass ~ Detection + Sex * Season + Sex * Species + Season * MigDate (m175)
model_summary <- summary(m175)$coefficients
conf_intervals <- confint(m175, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Detection + (Sex * Season) + Species + MigDate (m130)
model_summary <- summary(m130)$coefficients
conf_intervals <- confint(m130, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Detection + Sex * MigDate + Species (m107)
model_summary <- summary(m107)$coefficients
conf_intervals <- confint(m107, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Detection + Sex * Species + Sex * MigDate (m163)
model_summary <- summary(m163)$coefficients
conf_intervals <- confint(m163, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# ----------------------------------------------------------------------------- #

# Model comparisons: Concentration, Log(Concentration), Detection
m.c <- lm(Mass ~ OverallNeonic, data = piper)
m.l <- lm(Mass ~ log10(OverallNeonic + 0.0001), data = piper)
m.d <- lm(Mass ~ Detection, data = piper)
m.null <- lm(Mass ~ 1, data = piper)

special_models <- list(m.null, m.c, m.l, m.d)

special_names <- c('m.null', 'm.c', 'm.l', 'm.d')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Calculate AIC table
aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

# Plotting ####

## Added-Variable Plots ####

library(car)
avPlots(m175)

?avPlots

## Effect Plots ####
library(effects)
plot(allEffects(m175))

## ggplot Visuals ####
library(ggplot2)
ggplot(piper, aes(x = Detection, y = Mass)) +
  geom_boxplot() + 
  theme_classic() + 
  xlab("Neonicotinoid Detection Status") + 
  ylab("Mass")

# ----------------------------------------------------------------------------- #

# Why is there a positive relationship between neonics and mass? ####

## Are birds with detections earlier migrants? ####
# Yes! But there's a lot of variability. When are birds most exposed to neonics? Look at water and invert data.
ggplot(piper, aes(x = Detection, y = MigDate)) +
  geom_boxplot() + 
  theme_classic() + 
  xlab("Neonicotinoid Detection Status") + 
  ylab("Date into Stopover Season")

# Assess normality assumption (t-test vs. Mann-Whitney)
hist(piper$MigDate)
shapiro.test(piper$MigDate)
wilcox.test(MigDate ~ Detection, data = piper)


