#---------------------------------#
# Calidris sp. Body Mass Analysis #
#       Created 10/21/2024        #
#       Modified 10/21/2024       #
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

# Subset data for Least Sandcalidris
calidris <- subset(birds, Species %in% c("SemipalmatedSandpiper", "LeastSandpiper"))

# Remove birds that do not have neonicotinoid concentrations
calidris <- calidris[!is.na(calidris$OverallNeonic), ]

# Make neonicotinoid detection column
calidris$Detection <- ifelse(calidris$OverallNeonic > 0, "Detection", "Non-detection")

# Convert categorical variables to factor
calidris$Season <- as.factor(calidris$Season)
calidris$Year <- as.factor(calidris$Year)
calidris$Sex <- as.factor(calidris$Sex)
calidris$Detection <- as.factor(calidris$Detection)

# Check for variable correlation
cor(calidris$Julian, calidris$MigDate)

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

# Dr. Johnson confirmed these residual plots are good even though more 
# variable around 0. The "megaphone" pattern is not concerning as it's mostly 
# a reflection of non-detects (the opposite direction is more concerning). 

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

# Looks better with logarithmic transformation.

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

# First Stage AICc Model Selection: Mass ####

### Controlling Variables: Year + Sex + Season + Julian + MigDate + Species

### null and global models
m.null <- lm(Mass ~ 1, data= calidris)
m.globalyearjulian <- lm(Mass ~ Year + Sex + Julian + 
                           Year*Sex + Year*Julian +
                           Sex*Julian + Sex*Species, data = calidris)

m.globalspeciesjulian <-lm(Mass ~ Species + Sex + Julian +
                             Species*Sex + Species*Julian + Sex*Julian +
                             Species, data = calidris)

m.globalyearmig <- lm(Mass ~ Year + Sex + MigDate +
                        Year*Sex + Year*MigDate +
                        Sex*MigDate, data = calidris)

m.globalseasonmig <- lm(Mass ~ Sex + Season + Species + MigDate +
                          Sex*Season + Sex*Species + Sex*MigDate +
                          Season*Species + Season*MigDate + Species*MigDate, data = calidris)

### additive combinations (with restrictions, see above)
m1 <- lm(Mass ~ Year, data = calidris)
m2 <- lm(Mass ~ Season, data = calidris)
m3 <- lm(Mass ~ Sex, data = calidris)
m4 <- lm(Mass ~ Species, data = calidris)
m5 <- lm(Mass ~ Julian, data = calidris)
m6 <- lm(Mass ~ MigDate, data = calidris)

m7 <- lm(Mass ~ Year + Sex, data = calidris)
m8 <- lm(Mass ~ Year + Julian, data = calidris)
m9  <- lm(Mass ~ Year + MigDate, data = calidris)
m10 <- lm(Mass ~ Season + Sex, data = calidris)
m11 <- lm(Mass ~ Season + Species, data = calidris)
m12 <- lm(Mass ~ Season + MigDate, data = calidris)
m13 <- lm(Mass ~ Sex + Species, data = calidris)
m14 <- lm(Mass ~ Sex + Julian, data = calidris)
m15 <- lm(Mass ~ Sex + MigDate, data = calidris)
m16 <- lm(Mass ~ Species + Julian, data = calidris)
m17 <- lm(Mass ~ Species + MigDate, data = calidris)

m18 <- lm(Mass ~ Year + Sex + Julian, data = calidris)
m19 <- lm(Mass ~ Year + Sex + MigDate, data = calidris)

m20 <- lm(Mass ~ Sex + Season + Species, data = calidris)
m21 <- lm(Mass ~ Sex + Season + MigDate, data = calidris)
m22 <- lm(Mass ~ Sex + Julian + Species, data = calidris)
m23 <- lm(Mass ~ Sex + Species + MigDate, data = calidris)

m24 <- lm(Mass ~ Season + Species + MigDate, data = calidris)

m25 <- lm(Mass ~ Sex + Season + Species + MigDate, data = calidris)

### two-way interactions (with restrictions, see above)
m26 <- lm(Mass ~ Year * Sex, data = calidris)
m27 <- lm(Mass ~ Year * Julian, data = calidris)
m28 <- lm(Mass ~ Year * MigDate, data = calidris)

m29 <- lm(Mass ~ Sex * Season, data = calidris)
m30 <- lm(Mass ~ Sex * Julian, data = calidris)
m31 <- lm(Mass ~ Sex * Species, data = calidris)
m32 <- lm(Mass ~ Sex * MigDate, data = calidris)

m33 <- lm(Mass ~ Season * Species, data = calidris)
m34 <- lm(Mass ~ Season * MigDate, data = calidris)
m35 <- lm(Mass ~ Species * Julian, data = calidris)

### two-way interactions with one additive combination (with restrictions, see above)

#Year
m36 <- lm(Mass ~ (Year * Sex) + Julian, data = calidris)
m37 <- lm(Mass ~ (Year * Sex) + MigDate, data = calidris)

m38 <- lm(Mass ~ (Year * Julian) + Sex, data = calidris)

m39 <- lm(Mass ~ (Year * MigDate) + Sex, data = calidris)

#Sex
m40 <- lm(Mass ~ (Sex * Season) + Species, data = calidris)
m41 <- lm(Mass ~ (Sex * Season) + MigDate, data = calidris)

m42 <- lm(Mass ~ (Sex * Julian) + Year, data = calidris)
m43 <- lm(Mass ~ (Sex * Julian) + Species, data = calidris)

m44 <- lm(Mass ~ (Sex * Species) + Season, data = calidris)
m45 <- lm(Mass ~ (Sex * Species) + Julian, data = calidris)
m46 <- lm(Mass ~ (Sex * Species) + MigDate, data = calidris)

m47 <- lm(Mass ~ (Sex * MigDate) + Year, data = calidris)
m48 <- lm(Mass ~ (Sex * MigDate) + Season, data = calidris)
m49 <- lm(Mass ~ (Sex * MigDate) + Species, data = calidris)

#Season
m50 <- lm(Mass ~ (Season * Species) + Sex, data = calidris)
m51 <- lm(Mass ~ (Season * Species) + MigDate, data = calidris)

m52 <- lm(Mass ~ (Season * MigDate) + Sex, data = calidris)
m53 <- lm(Mass ~ (Season * MigDate) + Species, data = calidris)

#Species
m54 <- lm(Mass ~ (Species * Julian) + Sex, data = calidris)

### two-way interactions with two additive combinations (with restrictions, see above)
m55 <- lm(Mass ~ (Sex * Season) + Species + MigDate, data = calidris)
m56 <- lm(Mass ~ Sex + (Season * Species) + MigDate, data = calidris)
m57 <- lm(Mass ~ Sex + Season + (Species * MigDate), data = calidris)

m58 <- lm(Mass ~ (Sex * MigDate) + Species + Season, data = calidris)
m59 <- lm(Mass ~ (Sex * Species) + Season + MigDate, data = calidris)

### two-way interactions with multiple interactions (with restrictions, see above)

#two-way
m60 <- lm(Mass ~ Year*Sex + Year*Julian, data = calidris)
m61 <- lm(Mass ~ Year*Sex + Year*MigDate, data = calidris)
m62 <- lm(Mass ~ Year*Sex + Sex*Julian, data = calidris)
m64 <- lm(Mass ~ Year*Sex + Sex*MigDate, data = calidris)

m65 <- lm(Mass ~ Year*Julian + Year*Sex, data = calidris)
m66 <- lm(Mass ~ Year*Julian + Sex*Julian, data = calidris)

m67 <- lm(Mass ~ Year*MigDate + Sex*MigDate, data = calidris)

m68 <- lm(Mass ~ Sex*Season + Sex*Species, data = calidris)
m69 <- lm(Mass ~ Sex*Season + Sex*MigDate, data = calidris)
m70 <- lm(Mass ~ Sex*Season + Season*Species, data = calidris)

m71 <- lm(Mass ~ Sex*Julian + Sex*Species, data = calidris)
m72 <- lm(Mass ~ Sex*Julian + Species*Julian, data = calidris)

m73 <- lm(Mass ~ Sex*Species + Sex*MigDate, data = calidris)
m74 <- lm(Mass ~ Sex*Species + Season*Species, data = calidris)
m75 <- lm(Mass ~ Sex*Species + Season*MigDate, data = calidris)
m76 <- lm(Mass ~ Sex*Species + Species*Julian, data = calidris)

m77 <- lm(Mass ~ Sex*MigDate + Season*Species, data = calidris)
m78 <- lm(Mass ~ Sex*MigDate + Season*MigDate, data = calidris)

m79 <- lm(Mass ~ Season*Species + Season*MigDate, data = calidris)

#three-way
m80 <- lm(Mass ~ Year*Sex + Year*Julian + Sex*Julian, data = calidris)
m81 <- lm(Mass ~ Year*Sex + Year*MigDate + Sex*MigDate, data = calidris)
m82 <- lm(Mass ~ Year*Julian + Year*Sex + Sex*Julian, data = calidris)

m83 <- lm(Mass ~ Sex*Season + Sex*Species + Sex*MigDate, data = calidris)
m84 <- lm(Mass ~ Sex*Season + Sex*Species + Season*Species, data = calidris)
m85 <- lm(Mass ~ Sex*Season + Sex*Species + Season*MigDate, data = calidris)
m86 <- lm(Mass ~ Sex*Season + Sex*Species + Species*Julian, data = calidris)
m87 <- lm(Mass ~ Sex*Season + Sex*MigDate + Season*Species, data = calidris)
m88 <- lm(Mass ~ Sex*Season + Sex*MigDate + Season*MigDate, data = calidris)
m89 <- lm(Mass ~ Sex*Season + Season*Species + Season*MigDate, data = calidris)

m90 <- lm(Mass ~ Sex*Julian + Sex*Species + Species*Julian, data = calidris)

m91 <- lm(Mass ~ Sex*Species + Sex*MigDate + Season*Species, data = calidris)
m92 <- lm(Mass ~ Sex*Species + Sex*MigDate + Season*MigDate, data = calidris)
m93 <- lm(Mass ~ Sex*Species + Season*Species + Season*MigDate, data = calidris)

m94 <- lm(Mass ~ Sex*MigDate + Season*Species + Season*MigDate, data = calidris)

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
m.null <- lm(Mass ~ 1, data = calidris)
m.informednull <- lm(Mass ~ Species + MigDate, data = calidris)

m.globalinformednull <- lm(Mass ~ Species + MigDate + OverallNeonic +
                             Species * MigDate + Species * OverallNeonic + 
                             MigDate * OverallNeonic, data = calidris)

m.globalspecies <- lm(Mass ~ Sex + Season + Species + MigDate + OverallNeonic + 
                          Sex*Season + Sex*Species + Sex*MigDate + Sex*OverallNeonic +
                          Season*Species + Season*MigDate + Season*OverallNeonic + 
                          Species*MigDate + Species*OverallNeonic + MigDate*OverallNeonic, data = calidris)

m.globalyear <-  lm(Mass ~ Sex + Season + Year + MigDate + OverallNeonic + 
                  Sex*Season + Sex*Year + Sex*MigDate + Sex*OverallNeonic +
                  Season*Year + Season*MigDate + Season*OverallNeonic + 
                  Year*MigDate + Year*OverallNeonic + MigDate*OverallNeonic, data = calidris)

### single combinations of informed and global model covariates with neonicotinoid concentrations
m1 <- lm(Mass ~ OverallNeonic, data = calidris)
m2 <- lm(Mass ~ OverallNeonic + Year, data = calidris)
m3 <- lm(Mass ~ OverallNeonic + Sex, data = calidris)
m4 <- lm(Mass ~ OverallNeonic + Species, data = calidris)
m5 <- lm(Mass ~ OverallNeonic + Season, data = calidris)
m6 <- lm(Mass ~ OverallNeonic + Julian, data = calidris)
m7 <- lm(Mass ~ OverallNeonic + MigDate, data = calidris)

### two additive combinations of informed and global model covariates with neonicotinoid concentrations
m8 <- lm(Mass ~ OverallNeonic + Year + Sex, data = calidris)
m9 <- lm(Mass ~ OverallNeonic + Year + Julian, data = calidris)
m10 <- lm(Mass ~ OverallNeonic + Year + MigDate, data = calidris)

m11 <- lm(Mass ~ OverallNeonic + Sex + Species, data = calidris)
m12 <- lm(Mass ~ OverallNeonic + Sex + Season, data = calidris)
m13 <- lm(Mass ~ OverallNeonic + Sex + Julian, data = calidris)
m14 <- lm(Mass ~ OverallNeonic + Sex + MigDate, data = calidris)

m15 <- lm(Mass ~ OverallNeonic + Species + Season, data = calidris)
m16 <- lm(Mass ~ OverallNeonic + Species + Julian, data = calidris)
m17 <- lm(Mass ~ OverallNeonic + Species + MigDate, data = calidris)

m18 <- lm(Mass ~ OverallNeonic + Season + MigDate, data = calidris)

### three additive combinations of informed and global model covariates with neonicotinoid concentrations
m19 <- lm(Mass ~ OverallNeonic + Year + Sex + Julian, data = calidris)
m20 <- lm(Mass ~ OverallNeonic + Year + Sex + MigDate, data = calidris)

m21 <- lm(Mass ~ OverallNeonic + Sex + Species + Season, data = calidris)
m22 <- lm(Mass ~ OverallNeonic + Sex + Species + Julian, data = calidris)
m23 <- lm(Mass ~ OverallNeonic + Sex + Species + MigDate, data = calidris)

m24 <- lm(Mass ~ OverallNeonic + Sex + Season + MigDate, data = calidris)

m25 <- lm(Mass ~ OverallNeonic + Species + Season + MigDate, data = calidris)

### four additive combinations of informed and global model covariates with neonicotinoid concentrations
m26 <- lm(Mass ~ OverallNeonic + Season + Sex + Species + MigDate, data = calidris)

### two-way interactions (no additive combinations)
m27 <- lm(Mass ~ OverallNeonic*Year, data = calidris)
m28 <- lm(Mass ~ OverallNeonic*Sex, data = calidris)
m29 <- lm(Mass ~ OverallNeonic*Species, data = calidris)
m30 <- lm(Mass ~ OverallNeonic*Season, data = calidris)
m31 <- lm(Mass ~ OverallNeonic*Julian, data = calidris)
m32 <- lm(Mass ~ OverallNeonic*MigDate, data = calidris)

### two-way interactions with a single additive combinations of informed and global model covariates with neonicotinoid concentrations
m33 <- lm(Mass ~ OverallNeonic*Year + MigDate, data = calidris)
m34 <- lm(Mass ~ OverallNeonic*Year + Julian, data = calidris)
m35 <- lm(Mass ~ OverallNeonic*Year + Sex, data = calidris)

m36 <- lm(Mass ~ OverallNeonic*Sex + Year, data = calidris)
m37 <- lm(Mass ~ OverallNeonic*Sex + Season, data = calidris)
m38 <- lm(Mass ~ OverallNeonic*Sex + MigDate, data = calidris)
m39 <- lm(Mass ~ OverallNeonic*Sex + Julian, data = calidris)
m40 <- lm(Mass ~ OverallNeonic*Sex + Species, data = calidris)

m41 <- lm(Mass ~ OverallNeonic*Species + Season, data = calidris)
m42 <- lm(Mass ~ OverallNeonic*Species + MigDate, data = calidris)
m43 <- lm(Mass ~ OverallNeonic*Species + Julian, data = calidris)
m44 <- lm(Mass ~ OverallNeonic*Species + Sex, data = calidris)

m45 <- lm(Mass ~ OverallNeonic*Season + MigDate, data = calidris)
m46 <- lm(Mass ~ OverallNeonic*Season + Species, data = calidris)
m47 <- lm(Mass ~ OverallNeonic*Season + Sex, data = calidris)

m48 <- lm(Mass ~ OverallNeonic*Julian + Year, data = calidris)
m49 <- lm(Mass ~ OverallNeonic*Julian + Species, data = calidris)
m50 <- lm(Mass ~ OverallNeonic*Julian + Sex, data = calidris)

m51 <- lm(Mass ~ OverallNeonic*MigDate + Year, data = calidris)
m52 <- lm(Mass ~ OverallNeonic*MigDate + Season, data = calidris)
m53 <- lm(Mass ~ OverallNeonic*MigDate + Species, data = calidris)
m54 <- lm(Mass ~ OverallNeonic*MigDate + Sex, data = calidris)

### two-way interactions with two additive combination of informed and global model covariates with neonicotinoid concentrations
#Year
m55 <- lm(Mass ~ OverallNeonic + (Year * Sex) + Julian, data = calidris)
m56 <- lm(Mass ~ (OverallNeonic * Year) + Sex + Julian, data = calidris)
m57 <- lm(Mass ~ (OverallNeonic * Sex) + Year + Julian, data = calidris)
m58 <- lm(Mass ~ (OverallNeonic * Julian) + Sex + Julian, data = calidris)

m59 <- lm(Mass ~ OverallNeonic + (Year * Sex) + MigDate, data = calidris)
m60 <- lm(Mass ~ (OverallNeonic * Year) + Sex + MigDate, data = calidris)
m61 <- lm(Mass ~ (OverallNeonic * Sex) + Year + MigDate, data = calidris)
m62 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Year, data = calidris)

m63 <- lm(Mass ~ OverallNeonic + (Year * Julian) + Sex, data = calidris)
m64 <- lm(Mass ~ (OverallNeonic * Year) + Julian + Sex, data = calidris)
m65 <- lm(Mass ~ (OverallNeonic * Julian) + Year + Sex, data = calidris)
m66 <- lm(Mass ~ (OverallNeonic * Sex) + Julian + Year, data = calidris)

m67 <- lm(Mass ~ OverallNeonic + (Year * MigDate) + Sex, data = calidris)
m68 <- lm(Mass ~ (OverallNeonic * Year) + MigDate + Sex, data = calidris)
m69 <- lm(Mass ~ (OverallNeonic * MigDate) + Year + Sex, data = calidris)
m70 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Year, data = calidris)

#Sex
m71 <- lm(Mass ~ OverallNeonic + (Sex * Season) + Species, data = calidris)
m72 <- lm(Mass ~ (OverallNeonic * Sex) + Season + Species, data = calidris)
m73 <- lm(Mass ~ (OverallNeonic * Season) + Sex + Species, data = calidris)
m74 <- lm(Mass ~ (OverallNeonic * Species) + Season + Sex, data = calidris)

m75 <- lm(Mass ~ OverallNeonic + (Sex * Season) + MigDate, data = calidris)
m76 <- lm(Mass ~ (OverallNeonic * Sex) + Season + MigDate, data = calidris)
m77 <- lm(Mass ~ (OverallNeonic * Season) + Sex + MigDate, data = calidris)
m78 <- lm(Mass ~ (OverallNeonic * MigDate) + Season + Sex, data = calidris)

m79 <- lm(Mass ~ OverallNeonic + (Sex * Julian) + Year, data = calidris)
m80 <- lm(Mass ~ (OverallNeonic * Sex) + Julian + Year, data = calidris)
m81 <- lm(Mass ~ (OverallNeonic * Julian) + Sex + Year, data = calidris)
m82 <- lm(Mass ~ (OverallNeonic * Year) + Julian + Sex, data = calidris)

m83 <- lm(Mass ~ OverallNeonic + (Sex * Julian) + Species, data = calidris)
m84 <- lm(Mass ~ (OverallNeonic * Sex) + Julian + Species, data = calidris)
m85 <- lm(Mass ~ (OverallNeonic * Julian) + Sex + Species, data = calidris)
m86 <- lm(Mass ~ (OverallNeonic * Species) + Julian + Sex, data = calidris)

m87 <- lm(Mass ~ OverallNeonic + (Sex * Species) + Season, data = calidris)
m88 <- lm(Mass ~ (OverallNeonic * Sex) + Species + Season, data = calidris)
m89 <- lm(Mass ~ (OverallNeonic * Species) + Sex + Season, data = calidris)
m90 <- lm(Mass ~ (OverallNeonic * Season) + Species + Sex, data = calidris)

m91 <- lm(Mass ~ OverallNeonic + (Sex * Species) + Julian, data = calidris)
m92 <- lm(Mass ~ (OverallNeonic * Sex) + Species + Julian, data = calidris)
m93 <- lm(Mass ~ (OverallNeonic * Species) + Sex + Julian, data = calidris)
m94 <- lm(Mass ~ (OverallNeonic * Julian) + Species + Sex, data = calidris)

m95 <- lm(Mass ~ OverallNeonic + (Sex * Species) + MigDate, data = calidris)
m96 <- lm(Mass ~ (OverallNeonic * Sex) + Species + MigDate, data = calidris)
m97 <- lm(Mass ~ (OverallNeonic * Species) + Sex + MigDate, data = calidris)
m98 <- lm(Mass ~ (OverallNeonic * MigDate) + Species + Sex, data = calidris)

m99 <- lm(Mass ~ OverallNeonic + (Sex * MigDate) + Year, data = calidris)
m100 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Year, data = calidris)
m101 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Year, data = calidris)
m102<- lm(Mass ~ (OverallNeonic * Year) + MigDate + Sex, data = calidris)

m103 <- lm(Mass ~ OverallNeonic + (Sex * MigDate) + Season, data = calidris)
m104 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Season, data = calidris)
m105 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Season, data = calidris)
m106 <- lm(Mass ~ (OverallNeonic * Season) + MigDate + Sex, data = calidris)

m107 <- lm(Mass ~ OverallNeonic + (Sex * MigDate) + Species, data = calidris)
m108 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Species, data = calidris)
m109 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Species, data = calidris)
m110 <- lm(Mass ~ (OverallNeonic * Species) + MigDate + Sex, data = calidris)

#Season
m111 <- lm(Mass ~ OverallNeonic + (Season * Species) + Sex, data = calidris)
m112 <- lm(Mass ~ (OverallNeonic * Season) + Species + Sex, data = calidris)
m113 <- lm(Mass ~ (OverallNeonic * Species) + Season + Sex, data = calidris)
m114 <- lm(Mass ~ (OverallNeonic * Sex) + Species + Season, data = calidris)

m115 <- lm(Mass ~ OverallNeonic + (Season * Species) + MigDate, data = calidris)
m116 <- lm(Mass ~ (OverallNeonic * Season) + Species + MigDate, data = calidris)
m117 <- lm(Mass ~ (OverallNeonic * Species) + Season + MigDate, data = calidris)
m118 <- lm(Mass ~ (OverallNeonic * MigDate) + Species + Season, data = calidris)

m119 <- lm(Mass ~ OverallNeonic + (Season * MigDate) + Sex, data = calidris)
m120 <- lm(Mass ~ (OverallNeonic * Season) + MigDate + Sex, data = calidris)
m121 <- lm(Mass ~ (OverallNeonic * MigDate) + Season + Sex, data = calidris)
m122 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Season, data = calidris)

m123 <- lm(Mass ~ OverallNeonic + (Season * MigDate) + Species, data = calidris)
m124 <- lm(Mass ~ (OverallNeonic * Season) + MigDate + Species, data = calidris)
m125 <- lm(Mass ~ (OverallNeonic * MigDate) + Season + Species, data = calidris)
m126 <- lm(Mass ~ (OverallNeonic * Species) + MigDate + Season, data = calidris)

#Species)
m127 <- lm(Mass ~ (OverallNeonic * Species) + Julian + Sex, data = calidris)
m128 <- lm(Mass ~ (OverallNeonic * Julian) + Species + Sex, data = calidris)
m129 <- lm(Mass ~ (OverallNeonic * Sex) + Julian + Species, data = calidris)

### two-way interactions with three additive combinations (with restrictions, see above)
# Year + Sex + Species + Season + Julian + MigDate
m130 <- lm(Mass ~ OverallNeonic + (Sex * Season) + Species + MigDate, data = calidris)
m131 <- lm(Mass ~ OverallNeonic + Sex + (Season * Species) + MigDate, data = calidris)
m132 <- lm(Mass ~ OverallNeonic + Sex + Season + (Species * MigDate), data = calidris)
m133 <- lm(Mass ~ (OverallNeonic * Sex) + Season + Species + MigDate, data = calidris)
m134 <- lm(Mass ~ (OverallNeonic * Season) + Sex + Species + MigDate, data = calidris)
m135 <- lm(Mass ~ (OverallNeonic * Species) + Season + Sex + MigDate, data = calidris)
m136 <- lm(Mass ~ (OverallNeonic * MigDate) + Season + Species + Sex, data = calidris)

m137 <- lm(Mass ~ OverallNeonic + (Sex * MigDate) + Species + Season, data = calidris)
m138 <- lm(Mass ~ OverallNeonic + (Sex * Species) + Season + MigDate, data = calidris)

### neonicotinoid concentrations with one interaction (with restrictions, see above)
m139 <- lm(Mass ~ OverallNeonic + Year*Sex, data = calidris)
m140 <- lm(Mass ~ OverallNeonic + Year*Season, data = calidris)
m141 <- lm(Mass ~ OverallNeonic + Year*Julian, data = calidris)
m142 <- lm(Mass ~ OverallNeonic + Year*MigDate, data = calidris)

m143 <- lm(Mass ~ OverallNeonic + Sex*Species, data = calidris)
m144 <- lm(Mass ~ OverallNeonic + Sex*Season, data = calidris)
m145 <- lm(Mass ~ OverallNeonic + Sex*Julian, data = calidris)
m146 <- lm(Mass ~ OverallNeonic + Sex*MigDate, data = calidris)

m147 <- lm(Mass ~ OverallNeonic + Species*Season, data = calidris)
m148 <- lm(Mass ~ OverallNeonic + Species*Julian, data = calidris)
m149 <- lm(Mass ~ OverallNeonic + Species*MigDate, data = calidris)

m150 <- lm(Mass ~ OverallNeonic + Season*MigDate, data = calidris)

### neonicotinoid concentrations with two interactions (with restrictions, see above)
m151 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*Julian, data = calidris)
m152 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*MigDate, data = calidris)
m153 <- lm(Mass ~ OverallNeonic + Year*Sex + Sex*Julian, data = calidris)
m154 <- lm(Mass ~ OverallNeonic + Year*Sex + Sex*MigDate, data = calidris)

m155 <- lm(Mass ~ OverallNeonic + Year*Julian + Year*Sex, data = calidris)
m156 <- lm(Mass ~ OverallNeonic + Year*Julian + Sex*Julian, data = calidris)

m157 <- lm(Mass ~ OverallNeonic + Year*MigDate + Sex*MigDate, data = calidris)

m158 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*Species, data = calidris)
m159 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*MigDate, data = calidris)
m160 <- lm(Mass ~ OverallNeonic + Sex*Season + Season*Species, data = calidris)

m161 <- lm(Mass ~ OverallNeonic + Sex*Julian + Sex*Species, data = calidris)
m162 <- lm(Mass ~ OverallNeonic + Sex*Julian + Species*Julian, data = calidris)

m163 <- lm(Mass ~ OverallNeonic + Sex*Species + Sex*MigDate, data = calidris)
m164 <- lm(Mass ~ OverallNeonic + Sex*Species + Season*Species, data = calidris)
m165 <- lm(Mass ~ OverallNeonic + Sex*Species + Season*MigDate, data = calidris)
m166 <- lm(Mass ~ OverallNeonic + Sex*Species + Species*Julian, data = calidris)

m167 <- lm(Mass ~ OverallNeonic + Sex*MigDate + Season*Species, data = calidris)
m168 <- lm(Mass ~ OverallNeonic + Sex*MigDate + Season*MigDate, data = calidris)
m169 <- lm(Mass ~ OverallNeonic + Season*Species + Season*MigDate, data = calidris)

### neonicotinoid concentrations with three interactions (with restrictions, see above)
m170 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*Julian + Sex*Julian, data = calidris)
m171 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*MigDate + Sex*MigDate, data = calidris)
m172 <- lm(Mass ~ OverallNeonic + Year*Julian + Year*Sex + Sex*Julian, data = calidris)
m173 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*Species + Sex*MigDate, data = calidris)
m174 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*Species + Season*Species, data = calidris)
m175 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*Species + Season*MigDate, data = calidris)
m176 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*Species + Species*Julian, data = calidris)
m177 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*MigDate + Season*Species, data = calidris)
m178 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*MigDate + Season*MigDate, data = calidris)
m179 <- lm(Mass ~ OverallNeonic + Sex*Season + Season*Species + Season*MigDate, data = calidris)
m180 <- lm(Mass ~ OverallNeonic + Sex*Julian + Sex*Species + Species*Julian, data = calidris)
m181 <- lm(Mass ~ OverallNeonic + Sex*Species + Sex*MigDate + Season*Species, data = calidris)
m182 <- lm(Mass ~ OverallNeonic + Sex*Species + Sex*MigDate + Season*MigDate, data = calidris)
m183 <- lm(Mass ~ OverallNeonic + Sex*Species + Season*Species + Season*MigDate, data = calidris)
m184 <- lm(Mass ~ OverallNeonic + Sex*MigDate + Season*Species + Season*MigDate, data = calidris)

### add neonicotinoid interaction with all other interaction combinations
#add single interaction
#Year
m185 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex, data = calidris)
m186 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian, data = calidris)
m187 <- lm(Mass ~ (OverallNeonic * Year) + Year*MigDate, data = calidris)
m188 <- lm(Mass ~ (OverallNeonic * Year) + Sex*Season, data = calidris)
m189 <- lm(Mass ~ (OverallNeonic * Year) + Sex*Julian, data = calidris)
m190 <- lm(Mass ~ (OverallNeonic * Year) + Sex*MigDate, data = calidris)
m191 <- lm(Mass ~ (OverallNeonic * Year) + Season*MigDate, data = calidris)

#Sex
m192 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex, data = calidris)
m193 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian, data = calidris)
m194 <- lm(Mass ~ (OverallNeonic * Sex) + Year*MigDate, data = calidris)
m195 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species, data = calidris)
m196 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season, data = calidris)
m197 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Julian, data = calidris)
m198 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*MigDate, data = calidris)
m199 <- lm(Mass ~ (OverallNeonic * Sex) + Species*Season, data = calidris)
m200 <- lm(Mass ~ (OverallNeonic * Sex) + Species*Julian, data = calidris)
m201 <- lm(Mass ~ (OverallNeonic * Sex) + Species*MigDate, data = calidris)
m202 <- lm(Mass ~ (OverallNeonic * Sex) + Season*MigDate, data = calidris)

#Species
m203 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species, data = calidris)
m204 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season, data = calidris)
m205 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Julian, data = calidris)
m206 <- lm(Mass ~ (OverallNeonic * Species) + Sex*MigDate, data = calidris)
m207 <- lm(Mass ~ (OverallNeonic * Species) + Species*Season, data = calidris)
m208 <- lm(Mass ~ (OverallNeonic * Species) + Species*Julian, data = calidris)
m209 <- lm(Mass ~ (OverallNeonic * Species) + Species*MigDate, data = calidris)
m210 <- lm(Mass ~ (OverallNeonic * Species) + Season*MigDate, data = calidris)

#Season
m211 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species, data = calidris)
m212 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season, data = calidris)
m213 <- lm(Mass ~ (OverallNeonic * Season) + Sex*MigDate, data = calidris)
m214 <- lm(Mass ~ (OverallNeonic * Season) + Species*Season, data = calidris)
m215 <- lm(Mass ~ (OverallNeonic * Season) + Species*MigDate, data = calidris)
m216 <- lm(Mass ~ (OverallNeonic * Season) + Season*MigDate, data = calidris)

#Julian
m217 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex, data = calidris)
m218 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian, data = calidris)
m219 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Species, data = calidris)
m220 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Julian, data = calidris)
m221 <- lm(Mass ~ (OverallNeonic * Julian) + Species*Julian, data = calidris)

#MigDate
m222 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex, data = calidris)
m223 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*MigDate, data = calidris)
m224 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species, data = calidris)
m225 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season, data = calidris)
m226 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*MigDate, data = calidris)
m227 <- lm(Mass ~ (OverallNeonic * MigDate) + Species*Season, data = calidris)
m228 <- lm(Mass ~ (OverallNeonic * MigDate) + Species*MigDate, data = calidris)
m229 <- lm(Mass ~ (OverallNeonic * MigDate) + Season*MigDate, data = calidris)

#add two interactions
### neonicotinoid concentrations with two interactions (with restrictions, see above)
#Year
m230 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*Julian, data = calidris)
m231 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*MigDate, data = calidris)
m232 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Sex*Julian, data = calidris)
m233 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Sex*MigDate, data = calidris)
m234 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian + Year*Sex, data = calidris)
m235 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian + Sex*Julian, data = calidris)
m236 <- lm(Mass ~ (OverallNeonic * Year) + Year*MigDate + Sex*MigDate, data = calidris)
m237 <- lm(Mass ~ (OverallNeonic * Year) + Sex*Julian + Sex*Species, data = calidris)
m238 <- lm(Mass ~ (OverallNeonic * Year) + Sex*Species + Sex*MigDate, data = calidris)

#Sex
m239 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*Julian, data = calidris)
m240 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*MigDate, data = calidris)
m241 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Sex*Julian, data = calidris)
m242 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Sex*MigDate, data = calidris)
m243 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Year*Sex, data = calidris)
m244 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Sex*Julian, data = calidris)
m245 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Sex*Species, data = calidris)
m246 <- lm(Mass ~ (OverallNeonic * Sex) + Year*MigDate + Sex*MigDate, data = calidris)
m247 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*Species, data = calidris)
m248 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*MigDate, data = calidris)
m249 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Season*Species, data = calidris)
m250 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Julian + Sex*Species, data = calidris)
m251 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Julian + Species*Julian, data = calidris)
m252 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Sex*MigDate, data = calidris)
m253 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Season*Species, data = calidris)
m254 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Season*MigDate, data = calidris)
m255 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Species*Julian, data = calidris)
m256 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*MigDate + Season*Species, data = calidris)
m257 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*MigDate + Season*MigDate, data = calidris)
m258 <- lm(Mass ~ (OverallNeonic * Sex) + Season*Species + Season*MigDate, data = calidris)

#Species
m259 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*Species, data = calidris)
m260 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*MigDate, data = calidris)
m261 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Season*Species, data = calidris)
m262 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Julian + Sex*Species, data = calidris)
m263 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Julian + Species*Julian, data = calidris)
m264 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Sex*MigDate, data = calidris)
m265 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Season*Species, data = calidris)
m266 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Season*MigDate, data = calidris)
m267 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Species*Julian, data = calidris)
m268 <- lm(Mass ~ (OverallNeonic * Species) + Sex*MigDate + Season*Species, data = calidris)
m269 <- lm(Mass ~ (OverallNeonic * Species) + Sex*MigDate + Season*MigDate, data = calidris)
m270 <- lm(Mass ~ (OverallNeonic * Species) + Season*Species + Season*MigDate, data = calidris)

#Season
m271 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*Species, data = calidris)
m272 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*MigDate, data = calidris)
m273 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Season*Species, data = calidris)
m274 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species + Sex*MigDate, data = calidris)
m275 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species + Season*Species, data = calidris)
m276 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species + Season*MigDate, data = calidris)
m277 <- lm(Mass ~ (OverallNeonic * Season) + Sex*MigDate + Season*Species, data = calidris)
m278 <- lm(Mass ~ (OverallNeonic * Season) + Sex*MigDate + Season*MigDate, data = calidris)
m279 <- lm(Mass ~ (OverallNeonic * Season) + Season*Species + Season*MigDate, data = calidris)

#Julian
m280 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex + Year*Julian, data = calidris)
m281 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex + Sex*Julian, data = calidris)
m282 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian + Year*Sex, data = calidris)
m283 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian + Sex*Julian, data = calidris)
m284 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Julian + Sex*Species, data = calidris)
m285 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Julian + Species*Julian, data = calidris)
m286 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Species + Season*Species, data = calidris)
m287 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Species + Species*Julian, data = calidris)

#MigDate
m288 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex + Year*MigDate, data = calidris)
m289 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex + Sex*MigDate, data = calidris)
m290 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*MigDate + Sex*MigDate, data = calidris)
m291 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*Species, data = calidris)
m292 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*MigDate, data = calidris)
m293 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Season*Species, data = calidris)
m294 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species + Sex*MigDate, data = calidris)
m295 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species + Season*Species, data = calidris)
m296 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species + Season*MigDate, data = calidris)
m297 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*MigDate + Season*Species, data = calidris)
m298 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*MigDate + Season*MigDate, data = calidris)
m299 <- lm(Mass ~ (OverallNeonic * MigDate) + Season*Species + Season*MigDate, data = calidris)

#add three interactions
#Year
m300 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*Julian + Sex*Julian, data = calidris)
m301 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*MigDate + Sex*MigDate, data = calidris)
m302 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian + Year*Sex + Sex*Julian, data = calidris)

#Sex
m303 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*Julian + Sex*Julian, data = calidris)
m304 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*MigDate + Sex*MigDate, data = calidris)
m305 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Year*Sex + Sex*Julian, data = calidris)
m306 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*Species + Sex*MigDate, data = calidris)
m307 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*Species + Season*Species, data = calidris)
m308 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*Species + Season*MigDate, data = calidris)
m309 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*Species + Species*Julian, data = calidris)
m310 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*MigDate + Season*Species, data = calidris)
m311 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*MigDate + Season*MigDate, data = calidris)
m312 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Season*Species + Season*MigDate, data = calidris)
m313 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Julian + Sex*Species + Species*Julian, data = calidris)
m314 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Sex*MigDate + Season*Species, data = calidris)
m315 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Sex*MigDate + Season*MigDate, data = calidris)
m316 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Species + Season*Species + Season*MigDate, data = calidris)
m317 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*MigDate + Season*Species + Season*MigDate, data = calidris)

#Species
m318 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*Species + Sex*MigDate, data = calidris)
m319 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*Species + Season*Species, data = calidris)
m320 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*Species + Season*MigDate, data = calidris)
m321 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*Species + Species*Julian, data = calidris)
m322 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*MigDate + Season*Species, data = calidris)
m323 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Sex*MigDate + Season*MigDate, data = calidris)
m324 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Season + Season*Species + Season*MigDate, data = calidris)
m325 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Julian + Sex*Species + Species*Julian, data = calidris)
m326 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Sex*MigDate + Season*Species, data = calidris)
m327 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Sex*MigDate + Season*MigDate, data = calidris)
m328 <- lm(Mass ~ (OverallNeonic * Species) + Sex*Species + Season*Species + Season*MigDate, data = calidris)
m329 <- lm(Mass ~ (OverallNeonic * Species) + Sex*MigDate + Season*Species + Season*MigDate, data = calidris)

#Season
m330 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*Species + Sex*MigDate, data = calidris)
m331 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*Species + Season*Species, data = calidris)
m332 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*Species + Season*MigDate, data = calidris)
m333 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*MigDate + Season*Species, data = calidris)
m334 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*MigDate + Season*MigDate, data = calidris)
m335 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Season*Species + Season*MigDate, data = calidris)
m336 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species + Sex*MigDate + Season*Species, data = calidris)
m337 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species + Sex*MigDate + Season*MigDate, data = calidris)
m338 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Species + Season*Species + Season*MigDate, data = calidris)
m339 <- lm(Mass ~ (OverallNeonic * Season) + Sex*MigDate + Season*Species + Season*MigDate, data = calidris)

#Julian
m340 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex + Year*Julian + Sex*Julian, data = calidris)
m341 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian + Year*Sex + Sex*Julian, data = calidris)
m342 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Julian + Sex*Species + Species*Julian, data = calidris)

#MigDate
m343 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex + Year*MigDate + Sex*MigDate, data = calidris)
m344 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*MigDate + Sex*Species + Sex*MigDate, data = calidris)
m345 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*Species + Sex*MigDate, data = calidris)
m346 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*Species + Season*Species, data = calidris)
m347 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*Species + Season*MigDate, data = calidris)
m348 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*MigDate + Season*Species, data = calidris)
m349 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*MigDate + Season*MigDate, data = calidris)
m350 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Season*Species + Season*MigDate, data = calidris)
m351 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species + Sex*MigDate + Season*Species, data = calidris)
m352 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species + Sex*MigDate + Season*MigDate, data = calidris)
m353 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Species + Season*Species + Season*MigDate, data = calidris)
m354 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*MigDate + Season*Species + Season*MigDate, data = calidris)

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
m.null <- lm(Mass ~ 1, data = calidris)
m.informednull <- lm(Mass ~ Species + MigDate, data = calidris)

m.globalinformednull <- lm(Mass ~ Species + MigDate + log10(OverallNeonic + 0.0001) +
                             Species * MigDate + Species * log10(OverallNeonic + 0.0001) + 
                             MigDate * log10(OverallNeonic + 0.0001), data = calidris)

m.globalspecies <- lm(Mass ~ Sex + Season + Species + MigDate + log10(OverallNeonic + 0.0001) + 
                        Sex*Season + Sex*Species + Sex*MigDate + Sex*log10(OverallNeonic + 0.0001) +
                        Season*Species + Season*MigDate + Season*log10(OverallNeonic + 0.0001) + 
                        Species*MigDate + Species*log10(OverallNeonic + 0.0001) + MigDate*log10(OverallNeonic + 0.0001), data = calidris)

m.globalyear <-  lm(Mass ~ Sex + Season + Year + MigDate + log10(OverallNeonic + 0.0001) + 
                      Sex*Season + Sex*Year + Sex*MigDate + Sex*log10(OverallNeonic + 0.0001) +
                      Season*Year + Season*MigDate + Season*log10(OverallNeonic + 0.0001) + 
                      Year*MigDate + Year*log10(OverallNeonic + 0.0001) + MigDate*log10(OverallNeonic + 0.0001), data = calidris)

### single combinations of informed and global model covariates with neonicotinoid concentrations
m1 <- lm(Mass ~ log10(OverallNeonic + 0.0001), data = calidris)
m2 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year, data = calidris)
m3 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex, data = calidris)
m4 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species, data = calidris)
m5 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season, data = calidris)
m6 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Julian, data = calidris)
m7 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + MigDate, data = calidris)

### two additive combinations of informed and global model covariates with neonicotinoid concentrations
m8 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Sex, data = calidris)
m9 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Julian, data = calidris)
m10 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + MigDate, data = calidris)

m11 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Species, data = calidris)
m12 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Season, data = calidris)
m13 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Julian, data = calidris)
m14 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + MigDate, data = calidris)

m15 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species + Season, data = calidris)
m16 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species + Julian, data = calidris)
m17 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species + MigDate, data = calidris)

m18 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season + MigDate, data = calidris)

### three additive combinations of informed and global model covariates with neonicotinoid concentrations
m19 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Sex + Julian, data = calidris)
m20 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Sex + MigDate, data = calidris)

m21 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Species + Season, data = calidris)
m22 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Species + Julian, data = calidris)
m23 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Species + MigDate, data = calidris)

m24 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Season + MigDate, data = calidris)

m25 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species + Season + MigDate, data = calidris)

### four additive combinations of informed and global model covariates with neonicotinoid concentrations
m26 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season + Sex + Species + MigDate, data = calidris)

### two-way interactions (no additive combinations)
m27 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year, data = calidris)
m28 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex, data = calidris)
m29 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Species, data = calidris)
m30 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season, data = calidris)
m31 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian, data = calidris)
m32 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate, data = calidris)

### two-way interactions with a single additive combinations of informed and global model covariates with neonicotinoid concentrations
m33 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year + MigDate, data = calidris)
m34 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year + Julian, data = calidris)
m35 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year + Sex, data = calidris)

m36 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Year, data = calidris)
m37 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Season, data = calidris)
m38 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + MigDate, data = calidris)
m39 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Julian, data = calidris)
m40 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Species, data = calidris)

m41 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Species + Season, data = calidris)
m42 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Species + MigDate, data = calidris)
m43 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Species + Julian, data = calidris)
m44 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Species + Sex, data = calidris)

m45 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season + MigDate, data = calidris)
m46 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season + Species, data = calidris)
m47 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season + Sex, data = calidris)

m48 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian + Year, data = calidris)
m49 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian + Species, data = calidris)
m50 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian + Sex, data = calidris)

m51 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Year, data = calidris)
m52 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Season, data = calidris)
m53 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Species, data = calidris)
m54 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Sex, data = calidris)

### two-way interactions with two additive combination of informed and global model covariates with neonicotinoid concentrations
#Year
m55 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * Sex) + Julian, data = calidris)
m56 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex + Julian, data = calidris)
m57 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year + Julian, data = calidris)
m58 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex + Julian, data = calidris)

m59 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * Sex) + MigDate, data = calidris)
m60 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex + MigDate, data = calidris)
m61 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year + MigDate, data = calidris)
m62 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Year, data = calidris)

m63 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * Julian) + Sex, data = calidris)
m64 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Julian + Sex, data = calidris)
m65 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year + Sex, data = calidris)
m66 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Julian + Year, data = calidris)

m67 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * MigDate) + Sex, data = calidris)
m68 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + MigDate + Sex, data = calidris)
m69 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year + Sex, data = calidris)
m70 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Year, data = calidris)

#Sex
m71 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Season) + Species, data = calidris)
m72 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season + Species, data = calidris)
m73 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex + Species, data = calidris)
m74 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Season + Sex, data = calidris)

m75 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Season) + MigDate, data = calidris)
m76 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season + MigDate, data = calidris)
m77 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex + MigDate, data = calidris)
m78 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season + Sex, data = calidris)

m79 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Julian) + Year, data = calidris)
m80 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Julian + Year, data = calidris)
m81 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex + Year, data = calidris)
m82 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Julian + Sex, data = calidris)

m83 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Julian) + Species, data = calidris)
m84 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Julian + Species, data = calidris)
m85 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex + Species, data = calidris)
m86 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Julian + Sex, data = calidris)

m87 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Species) + Season, data = calidris)
m88 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species + Season, data = calidris)
m89 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex + Season, data = calidris)
m90 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Species + Sex, data = calidris)

m91 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Species) + Julian, data = calidris)
m92 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species + Julian, data = calidris)
m93 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex + Julian, data = calidris)
m94 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Species + Sex, data = calidris)

m95 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Species) + MigDate, data = calidris)
m96 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species + MigDate, data = calidris)
m97 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex + MigDate, data = calidris)
m98 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Species + Sex, data = calidris)

m99 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * MigDate) + Year, data = calidris)
m100 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Year, data = calidris)
m101 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Year, data = calidris)
m102<- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + MigDate + Sex, data = calidris)

m103 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * MigDate) + Season, data = calidris)
m104 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Season, data = calidris)
m105 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Season, data = calidris)
m106 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + MigDate + Sex, data = calidris)

m107 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * MigDate) + Species, data = calidris)
m108 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Species, data = calidris)
m109 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Species, data = calidris)
m110 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + MigDate + Sex, data = calidris)

#Season
m111 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Season * Species) + Sex, data = calidris)
m112 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Species + Sex, data = calidris)
m113 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Season + Sex, data = calidris)
m114 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species + Season, data = calidris)

m115 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Season * Species) + MigDate, data = calidris)
m116 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Species + MigDate, data = calidris)
m117 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Season + MigDate, data = calidris)
m118 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Species + Season, data = calidris)

m119 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Season * MigDate) + Sex, data = calidris)
m120 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + MigDate + Sex, data = calidris)
m121 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season + Sex, data = calidris)
m122 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Season, data = calidris)

m123 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Season * MigDate) + Species, data = calidris)
m124 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + MigDate + Species, data = calidris)
m125 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season + Species, data = calidris)
m126 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + MigDate + Season, data = calidris)

#Species)
m127 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Julian + Sex, data = calidris)
m128 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Species + Sex, data = calidris)
m129 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Julian + Species, data = calidris)

### two-way interactions with three additive combinations (with restrictions, see above)
# Year + Sex + Species + Season + Julian + MigDate
m130 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Season) + Species + MigDate, data = calidris)
m131 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + (Season * Species) + MigDate, data = calidris)
m132 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Season + (Species * MigDate), data = calidris)
m133 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season + Species + MigDate, data = calidris)
m134 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex + Species + MigDate, data = calidris)
m135 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Season + Sex + MigDate, data = calidris)
m136 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season + Species + Sex, data = calidris)

m137 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * MigDate) + Species + Season, data = calidris)
m138 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Species) + Season + MigDate, data = calidris)

### neonicotinoid concentrations with one interaction (with restrictions, see above)
m139 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex, data = calidris)
m140 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Season, data = calidris)
m141 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian, data = calidris)
m142 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*MigDate, data = calidris)

m143 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species, data = calidris)
m144 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season, data = calidris)
m145 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Julian, data = calidris)
m146 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*MigDate, data = calidris)

m147 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species*Season, data = calidris)
m148 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species*Julian, data = calidris)
m149 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Species*MigDate, data = calidris)

m150 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season*MigDate, data = calidris)

### neonicotinoid concentrations with two interactions (with restrictions, see above)
m151 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*Julian, data = calidris)
m152 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*MigDate, data = calidris)
m153 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Sex*Julian, data = calidris)
m154 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Sex*MigDate, data = calidris)

m155 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian + Year*Sex, data = calidris)
m156 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian + Sex*Julian, data = calidris)

m157 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*MigDate + Sex*MigDate, data = calidris)

m158 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*Species, data = calidris)
m159 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*MigDate, data = calidris)
m160 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Season*Species, data = calidris)

m161 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Julian + Sex*Species, data = calidris)
m162 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Julian + Species*Julian, data = calidris)

m163 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Sex*MigDate, data = calidris)
m164 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Season*Species, data = calidris)
m165 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Season*MigDate, data = calidris)
m166 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Species*Julian, data = calidris)

m167 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*MigDate + Season*Species, data = calidris)
m168 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*MigDate + Season*MigDate, data = calidris)
m169 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season*Species + Season*MigDate, data = calidris)

### neonicotinoid concentrations with three interactions (with restrictions, see above)
m170 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*Julian + Sex*Julian, data = calidris)
m171 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*MigDate + Sex*MigDate, data = calidris)
m172 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian + Year*Sex + Sex*Julian, data = calidris)
m173 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*Species + Sex*MigDate, data = calidris)
m174 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*Species + Season*Species, data = calidris)
m175 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*Species + Season*MigDate, data = calidris)
m176 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*Species + Species*Julian, data = calidris)
m177 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*MigDate + Season*Species, data = calidris)
m178 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*MigDate + Season*MigDate, data = calidris)
m179 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Season*Species + Season*MigDate, data = calidris)
m180 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Julian + Sex*Species + Species*Julian, data = calidris)
m181 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Sex*MigDate + Season*Species, data = calidris)
m182 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Sex*MigDate + Season*MigDate, data = calidris)
m183 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Species + Season*Species + Season*MigDate, data = calidris)
m184 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*MigDate + Season*Species + Season*MigDate, data = calidris)

### add neonicotinoid interaction with all other interaction combinations
#add single interaction
#Year
m185 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex, data = calidris)
m186 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian, data = calidris)
m187 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*MigDate, data = calidris)
m188 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*Season, data = calidris)
m189 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*Julian, data = calidris)
m190 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*MigDate, data = calidris)
m191 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Season*MigDate, data = calidris)

#Sex
m192 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex, data = calidris)
m193 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian, data = calidris)
m194 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*MigDate, data = calidris)
m195 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species, data = calidris)
m196 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season, data = calidris)
m197 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Julian, data = calidris)
m198 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*MigDate, data = calidris)
m199 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species*Season, data = calidris)
m200 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species*Julian, data = calidris)
m201 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Species*MigDate, data = calidris)
m202 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season*MigDate, data = calidris)

#Species
m203 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species, data = calidris)
m204 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season, data = calidris)
m205 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Julian, data = calidris)
m206 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*MigDate, data = calidris)
m207 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Species*Season, data = calidris)
m208 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Species*Julian, data = calidris)
m209 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Species*MigDate, data = calidris)
m210 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Season*MigDate, data = calidris)

#Season
m211 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species, data = calidris)
m212 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season, data = calidris)
m213 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*MigDate, data = calidris)
m214 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Species*Season, data = calidris)
m215 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Species*MigDate, data = calidris)
m216 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Season*MigDate, data = calidris)

#Julian
m217 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex, data = calidris)
m218 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian, data = calidris)
m219 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Species, data = calidris)
m220 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Julian, data = calidris)
m221 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Species*Julian, data = calidris)

#MigDate
m222 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex, data = calidris)
m223 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*MigDate, data = calidris)
m224 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species, data = calidris)
m225 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season, data = calidris)
m226 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*MigDate, data = calidris)
m227 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Species*Season, data = calidris)
m228 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Species*MigDate, data = calidris)
m229 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season*MigDate, data = calidris)

#add two interactions
### neonicotinoid concentrations with two interactions (with restrictions, see above)
#Year
m230 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*Julian, data = calidris)
m231 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*MigDate, data = calidris)
m232 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Sex*Julian, data = calidris)
m233 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Sex*MigDate, data = calidris)
m234 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian + Year*Sex, data = calidris)
m235 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian + Sex*Julian, data = calidris)
m236 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*MigDate + Sex*MigDate, data = calidris)
m237 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*Julian + Sex*Species, data = calidris)
m238 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*Species + Sex*MigDate, data = calidris)

#Sex
m239 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*Julian, data = calidris)
m240 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*MigDate, data = calidris)
m241 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Sex*Julian, data = calidris)
m242 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Sex*MigDate, data = calidris)
m243 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Year*Sex, data = calidris)
m244 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Sex*Julian, data = calidris)
m245 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Sex*Species, data = calidris)
m246 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*MigDate + Sex*MigDate, data = calidris)
m247 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*Species, data = calidris)
m248 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*MigDate, data = calidris)
m249 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Season*Species, data = calidris)
m250 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Julian + Sex*Species, data = calidris)
m251 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Julian + Species*Julian, data = calidris)
m252 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Sex*MigDate, data = calidris)
m253 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Season*Species, data = calidris)
m254 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Season*MigDate, data = calidris)
m255 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Species*Julian, data = calidris)
m256 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*MigDate + Season*Species, data = calidris)
m257 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*MigDate + Season*MigDate, data = calidris)
m258 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season*Species + Season*MigDate, data = calidris)

#Species
m259 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*Species, data = calidris)
m260 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*MigDate, data = calidris)
m261 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Season*Species, data = calidris)
m262 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Julian + Sex*Species, data = calidris)
m263 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Julian + Species*Julian, data = calidris)
m264 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Sex*MigDate, data = calidris)
m265 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Season*Species, data = calidris)
m266 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Season*MigDate, data = calidris)
m267 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Species*Julian, data = calidris)
m268 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*MigDate + Season*Species, data = calidris)
m269 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*MigDate + Season*MigDate, data = calidris)
m270 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Season*Species + Season*MigDate, data = calidris)

#Season
m271 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*Species, data = calidris)
m272 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*MigDate, data = calidris)
m273 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Season*Species, data = calidris)
m274 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species + Sex*MigDate, data = calidris)
m275 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species + Season*Species, data = calidris)
m276 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species + Season*MigDate, data = calidris)
m277 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*MigDate + Season*Species, data = calidris)
m278 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*MigDate + Season*MigDate, data = calidris)
m279 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Season*Species + Season*MigDate, data = calidris)

#Julian
m280 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex + Year*Julian, data = calidris)
m281 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex + Sex*Julian, data = calidris)
m282 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian + Year*Sex, data = calidris)
m283 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian + Sex*Julian, data = calidris)
m284 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Julian + Sex*Species, data = calidris)
m285 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Julian + Species*Julian, data = calidris)
m286 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Species + Season*Species, data = calidris)
m287 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Species + Species*Julian, data = calidris)

#MigDate
m288 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex + Year*MigDate, data = calidris)
m289 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex + Sex*MigDate, data = calidris)
m290 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*MigDate + Sex*MigDate, data = calidris)
m291 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*Species, data = calidris)
m292 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*MigDate, data = calidris)
m293 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Season*Species, data = calidris)
m294 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species + Sex*MigDate, data = calidris)
m295 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species + Season*Species, data = calidris)
m296 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species + Season*MigDate, data = calidris)
m297 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*MigDate + Season*Species, data = calidris)
m298 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*MigDate + Season*MigDate, data = calidris)
m299 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season*Species + Season*MigDate, data = calidris)

#add three interactions
#Year
m300 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*Julian + Sex*Julian, data = calidris)
m301 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*MigDate + Sex*MigDate, data = calidris)
m302 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian + Year*Sex + Sex*Julian, data = calidris)

#Sex
m303 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*Julian + Sex*Julian, data = calidris)
m304 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*MigDate + Sex*MigDate, data = calidris)
m305 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Year*Sex + Sex*Julian, data = calidris)
m306 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*Species + Sex*MigDate, data = calidris)
m307 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*Species + Season*Species, data = calidris)
m308 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*Species + Season*MigDate, data = calidris)
m309 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*Species + Species*Julian, data = calidris)
m310 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*MigDate + Season*Species, data = calidris)
m311 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*MigDate + Season*MigDate, data = calidris)
m312 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Season*Species + Season*MigDate, data = calidris)
m313 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Julian + Sex*Species + Species*Julian, data = calidris)
m314 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Sex*MigDate + Season*Species, data = calidris)
m315 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Sex*MigDate + Season*MigDate, data = calidris)
m316 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Species + Season*Species + Season*MigDate, data = calidris)
m317 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*MigDate + Season*Species + Season*MigDate, data = calidris)

#Species
m318 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*Species + Sex*MigDate, data = calidris)
m319 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*Species + Season*Species, data = calidris)
m320 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*Species + Season*MigDate, data = calidris)
m321 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*Species + Species*Julian, data = calidris)
m322 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*MigDate + Season*Species, data = calidris)
m323 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Sex*MigDate + Season*MigDate, data = calidris)
m324 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Season + Season*Species + Season*MigDate, data = calidris)
m325 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Julian + Sex*Species + Species*Julian, data = calidris)
m326 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Sex*MigDate + Season*Species, data = calidris)
m327 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Sex*MigDate + Season*MigDate, data = calidris)
m328 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*Species + Season*Species + Season*MigDate, data = calidris)
m329 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Species) + Sex*MigDate + Season*Species + Season*MigDate, data = calidris)

#Season
m330 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*Species + Sex*MigDate, data = calidris)
m331 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*Species + Season*Species, data = calidris)
m332 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*Species + Season*MigDate, data = calidris)
m333 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*MigDate + Season*Species, data = calidris)
m334 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*MigDate + Season*MigDate, data = calidris)
m335 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Season*Species + Season*MigDate, data = calidris)
m336 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species + Sex*MigDate + Season*Species, data = calidris)
m337 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species + Sex*MigDate + Season*MigDate, data = calidris)
m338 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Species + Season*Species + Season*MigDate, data = calidris)
m339 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*MigDate + Season*Species + Season*MigDate, data = calidris)

#Julian
m340 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex + Year*Julian + Sex*Julian, data = calidris)
m341 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian + Year*Sex + Sex*Julian, data = calidris)
m342 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Julian + Sex*Species + Species*Julian, data = calidris)

#MigDate
m343 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex + Year*MigDate + Sex*MigDate, data = calidris)
m344 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*MigDate + Sex*Species + Sex*MigDate, data = calidris)
m345 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*Species + Sex*MigDate, data = calidris)
m346 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*Species + Season*Species, data = calidris)
m347 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*Species + Season*MigDate, data = calidris)
m348 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*MigDate + Season*Species, data = calidris)
m349 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*MigDate + Season*MigDate, data = calidris)
m350 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Season*Species + Season*MigDate, data = calidris)
m351 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species + Sex*MigDate + Season*Species, data = calidris)
m352 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species + Sex*MigDate + Season*MigDate, data = calidris)
m353 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Species + Season*Species + Season*MigDate, data = calidris)
m354 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*MigDate + Season*Species + Season*MigDate, data = calidris)

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
model_summary <- summary(m107)$coefficients
conf_intervals <- confint(m107, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Log(Neonic) + Sex * Species + Sex * MigDate (m163)
model_summary <- summary(m163)$coefficients
conf_intervals <- confint(m163, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Log(Neonic) + Sex * Season + Species + MigDate (m130)
model_summary <- summary(m130)$coefficients
conf_intervals <- confint(m130, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# Mass ~ Log(Neonic) + Sex * Season + Sex * Species + Season * MigDate (m175)
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



