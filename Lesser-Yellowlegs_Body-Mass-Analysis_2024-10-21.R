#--------------------------------------#
# Lesser Yellowlegs Body Mass Analysis #
#         Created 10/21/2024           #
#         Modified 10/21/2024          #
#--------------------------------------#

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

# Subset data for Lesser Yellowlegs
leye <- subset(birds, Species %in% c("LesserYellowlegs"))

# Make neonicotinoid detection column
leye$Detection <- ifelse(leye$OverallNeonic > 0, "Detection", "Non-detection")

# Convert categorical variables to factor
leye$Year <- as.factor(leye$Year)
leye$Season <- as.factor(leye$Season)
leye$Sex <- as.factor(leye$Sex)
leye$Detection <- as.factor(leye$Detection)

# Check for variable correlation
cor(leye$Julian, leye$MigDate)

# ---------------------------------------------------------------------------- #

## Covariate Restrictions ####

# Date into Season and Julian can't be in the same model (correlation = 0.89)
# Year and Season can't be in the same model because of lack of replication of Season within a year
# Julian and Season can't be in the same model because Season bins Julian

# ---------------------------------------------------------------------------- #

# Data Exploration & Linear Regression Model Assumptions ####

## Assumption 1: Zero Expectations of Errors (Residuals) ####

### Residuals Plots 

# Calculate residuals
m <- lm(Mass ~ OverallNeonic, data = leye)
leye$yhat <- predict(m)
leye$residuals <- residuals(m)
leye$rstudent <- rstudent(m)
head(leye)

# Simple residual plot
plot(predict(m), rstudent(m))

# Pretty residual plot
ggplot(leye, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# Dr. Johnson confirmed these residual plots are good even though more 
# variable around 0. The "megaphone" pattern is not concerning as it's mostly 
# a reflection of non-detects (the opposite direction is more concerning). 

# Log transformation of neonics

# Calculate residuals 
m.log <- lm(Mass ~ (log10(OverallNeonic+0.0001)), data = leye)
leye$yhat <- predict(m.log)
leye$residuals <- residuals(m.log)
leye$rstudent <- rstudent(m.log)
head(leye)

# Simple residual plot
plot(predict(m.log), rstudent(m.log))

# Pretty residual plot
ggplot(leye, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# Looks better with logarithmic transformation.

# ---------------------------------------------------------------------------- #

## Assumption 2: Equality of error variances (Homoscedasticity) #### 

# No transformation
m <- lm(Mass ~ OverallNeonic, data = leye)
leye$yhat <- predict(m)
leye$residuals <- residuals(m)
leye$rstudent <- rstudent(m)

# Plot homoscedasticity
ggplot(leye, aes(x = yhat, y = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual")

# No violations because variability doesn't increase as predicted value increases

# Log transformation
m.log <- lm(Mass ~ (log10(OverallNeonic+0.0001)), data = leye)
leye$yhat <- predict(m.log)
leye$residuals <- residuals(m.log)
leye$rstudent <- rstudent(m.log)

### PLOT HOMOSCEDASTICITY
ggplot(leye, aes(x = yhat, y = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual")

# No violations because variability doesn't increase as predicted value increases.
# Log transformation does appear to spread residuals more evenly. 

# ---------------------------------------------------------------------------- #

## Assumption 3: Independence of Errors/Observations ####

# Distribution of one error/observation does not depend on another error/observation

# ---------------------------------------------------------------------------- #

## Assumption 4: Normality ####

m <- lm(Mass ~ OverallNeonic, data = leye)
hist(rstandard(m))
qqnorm(rstandard(m))

# Does not severely violate normality assumption
# Regression models do not many any assumptions about the distrbution of explanatory 
# variables (e.g., Concentrations), but their distribution does affect inferences. 

# ---------------------------------------------------------------------------- #

# First Stage AICc Model Selection: Mass ####

### Controlling Variables: Year + Sex + Season + Julian + MigDate
 
### Null and global models 
m.null <- lm(Mass ~ 1, data= leye)
m.globalyearjulian <- lm(Mass ~ Year + Sex + Julian + 
                           Year*Sex + Year*Julian +  
                           Sex*Julian, data = leye)

m.globalyearmig <- lm(Mass ~ Year + Sex + MigDate +
                        Year*Sex + Year*MigDate +
                        Sex*MigDate, data = leye)

m.globalseasonmig <- lm(Mass ~ Sex + Season + MigDate +
                          Sex*Season + Sex*MigDate +
                          Season*MigDate, data = leye)

### Additive models 
m1 <- lm(Mass ~ Year, data = leye)
m2 <- lm(Mass ~ Season, data = leye)
m3 <- lm(Mass ~ Sex, data = leye)
m4 <- lm(Mass ~ Julian, data = leye)
m5 <- lm(Mass ~ MigDate, data = leye)

m6 <- lm(Mass ~ Year + Sex, data = leye)
m7 <- lm(Mass ~ Year + Julian, data = leye)
m8  <- lm(Mass ~ Year + MigDate, data = leye)
m9 <- lm(Mass ~ Season + Sex, data = leye)
m10 <- lm(Mass ~ Season + MigDate, data = leye)
m11 <- lm(Mass ~ Sex + Julian, data = leye)
m12 <- lm(Mass ~ Sex + MigDate, data = leye)

m13 <- lm(Mass ~ Year + Sex + Julian, data = leye)
m14 <- lm(Mass ~ Year + Sex + MigDate, data = leye)

m15 <- lm(Mass ~ Sex + Season + MigDate, data = leye)

### Two-way interactions 
m16 <- lm(Mass ~ Year * Sex, data = leye)
m17 <- lm(Mass ~ Year * Julian, data = leye)
m18 <- lm(Mass ~ Year * MigDate, data = leye)

m19 <- lm(Mass ~ Sex * Season, data = leye)
m20 <- lm(Mass ~ Sex * Julian, data = leye)
m21 <- lm(Mass ~ Sex * MigDate, data = leye)

m22 <- lm(Mass ~ Season * MigDate, data = leye)

### Two-way interactions with one additive combination 

# Year
m23 <- lm(Mass ~ (Year * Sex) + Julian, data = leye)
m24 <- lm(Mass ~ (Year * Sex) + MigDate, data = leye)
m25 <- lm(Mass ~ (Year * Julian) + Sex, data = leye)
m26 <- lm(Mass ~ (Year * MigDate) + Sex, data = leye)

# Sex
m27 <- lm(Mass ~ (Sex * Season) + MigDate, data = leye)
m28 <- lm(Mass ~ (Sex * Julian) + Year, data = leye)
m29 <- lm(Mass ~ (Sex * MigDate) + Year, data = leye)
m30 <- lm(Mass ~ (Sex * MigDate) + Season, data = leye)

# Season
m31 <- lm(Mass ~ (Season * MigDate) + Sex, data = leye)

### Two-way interactions with two additive combinations
m32 <- lm(Mass ~ (Year * Sex) + Julian, data = leye)
m33 <- lm(Mass ~ Year + (Sex * Julian), data = leye)

m34 <- lm(Mass ~ (Year * Julian) + Sex, data = leye)

m35 <- lm(Mass ~ (Year * Sex) + MigDate, data = leye)

m36 <- lm(Mass ~ (Year * MigDate) + Sex, data = leye)
m37 <- lm(Mass ~ Year + (MigDate * Sex), data = leye)

m38 <- lm(Mass ~ (Sex * Season) + MigDate, data = leye)
m39 <- lm(Mass ~ Sex + (Season * MigDate), data = leye)

m40 <- lm(Mass ~ (Sex * MigDate) + Season, data = leye)

### Two-way interactions with multiple interactions

# Two-way
m41 <- lm(Mass ~ Year*Sex + Year*Julian, data = leye)
m42 <- lm(Mass ~ Year*Sex + Year*MigDate, data = leye)
m43 <- lm(Mass ~ Year*Sex + Sex*Julian, data = leye)
m44 <- lm(Mass ~ Year*Sex + Sex*MigDate, data = leye)
m45 <- lm(Mass ~ Year*Julian + Year*Sex, data = leye)
m46 <- lm(Mass ~ Year*Julian + Sex*Julian, data = leye)
m47 <- lm(Mass ~ Year*MigDate + Sex*MigDate, data = leye)
m48 <- lm(Mass ~ Sex*Season + Sex*MigDate, data = leye)
m49 <- lm(Mass ~ Sex*MigDate + Season*MigDate, data = leye)

# Three-way
m50 <- lm(Mass ~ Year*Sex + Year*Julian + Sex*Julian, data = leye)
m51 <- lm(Mass ~ Year*Julian + Year*Sex + Sex*Julian, data = leye)
m52 <- lm(Mass ~ Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

# ---------------------------------------------------------------------------- #

## AICc Model Selection: First Stage ####
models <- list(m.null, m.globalyearjulian, m.globalseasonmig, m.globalyearmig, 
               m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, 
               m16, m17, m18, m19, m20, m21, m22, m23, m24, m25, m26, m27, m28, 
               m29, m30, m31, m32, m33, m34, m35, m36, m37, m38, m39, m40, m41,
               m42, m43, m44, m45, m46, m47, m48, m49, m50, m51, m52)

mod.names <- c('m.null', 'm.globalyearjulian', 'm.globalseasonmig', 'm.globalyearmig', 
               'm1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 
               'm11', 'm12', 'm13', 'm14', 'm15', 'm16', 'm17', 'm18', 'm19', 
               'm20', 'm21', 'm22', 'm23', 'm24', 'm25', 'm26', 'm27', 'm28', 
               'm29', 'm30', 'm31', 'm32', 'm33', 'm34', 'm35', 'm36', 'm37', 
               'm38', 'm39', 'm40', 'm41', 'm42', 'm43', 'm44', 'm45', 'm46', 
               'm47', 'm48', 'm49', 'm50', 'm51', 'm52')

aictab(models, modnames = mod.names)

# Shortened AICc Model Selection With Top Models
models <- list(m1, m6, m18, m8, m7, m.globalyearmig, m.globalyearjulian, 
               m.null, m.globalseasonmig)

mod.names <- c('m1', 'm6', 'm18', 'm8', 'm7', 'm.globalyearmig', 
               'm.globalyearjulian', 'm.null', 'm.globalseasonmig')

aictab(models, modnames = mod.names)

# ---------------------------------------------------------------------------- #

### Top Model Summaries: First Stage ####
options(digits = 3)

# Mass ~ Year (m1) <- INFORMED NULL MODEL
cbind(summary(m1)$coefficients, confint(m1))

# ---------------------------------------------------------------------------- #

### Conclusion: First Stage ####
# Mass ~ Year is informed null model. All other models have Delta AICc > 2.

# ---------------------------------------------------------------------------- #

# Second Stage AICc Model Selection: Mass ~ Concentration ####

### Null and global models
m.null <- lm(Mass ~ 1, data = leye)
m.informednull <- lm(Mass ~ Year, data = leye)
m.globalinformednull <- lm(Mass ~ Year * OverallNeonic, data = leye)

m.globalyearjulian <- lm(Mass ~ Year + Sex + Julian + OverallNeonic + 
                           Year*Sex + Year*Julian + Year*OverallNeonic + 
                           Sex*Julian + Sex*OverallNeonic + 
                           Julian*OverallNeonic, data = leye)

m.globalyearmig <- lm(Mass ~ Year + Sex + MigDate + OverallNeonic +
                        Year*Sex + Year*MigDate + Year*OverallNeonic + 
                        Sex*MigDate + Sex*OverallNeonic + 
                        MigDate*OverallNeonic, data = leye)

m.globalseasonmig <- lm(Mass ~ Sex + Season + MigDate + OverallNeonic + 
                          Sex*Season + Sex*MigDate + Sex*OverallNeonic +
                          Season*MigDate + Season*OverallNeonic + 
                          MigDate*OverallNeonic, data = leye)

### Single combinations of informed and global model covariates
m1 <- lm(Mass ~ OverallNeonic, data = leye)
m2 <- lm(Mass ~ OverallNeonic + Year, data = leye)
m3 <- lm(Mass ~ OverallNeonic + Sex, data = leye)
m4 <- lm(Mass ~ OverallNeonic + Season, data = leye)
m5 <- lm(Mass ~ OverallNeonic + Julian, data = leye)
m6 <- lm(Mass ~ OverallNeonic + MigDate, data = leye)

### Two additive combinations of informed and global model covariates
m7 <- lm(Mass ~ OverallNeonic + Year + Sex, data = leye)
m8 <- lm(Mass ~ OverallNeonic + Year + Julian, data = leye)
m9 <- lm(Mass ~ OverallNeonic + Year + MigDate, data = leye)

m10 <- lm(Mass ~ OverallNeonic + Sex + Season, data = leye)
m11 <- lm(Mass ~ OverallNeonic + Sex + Julian, data = leye)
m12 <- lm(Mass ~ OverallNeonic + Sex + MigDate, data = leye)

m13 <- lm(Mass ~ OverallNeonic + Season + MigDate, data = leye)

### Three additive combinations of informed and global model covariates
m14 <- lm(Mass ~ OverallNeonic + Year + Sex + Julian, data = leye)
m15 <- lm(Mass ~ OverallNeonic + Year + Sex + MigDate, data = leye)
m16 <- lm(Mass ~ OverallNeonic + Sex + Season + MigDate, data = leye)

### Two-way interactions (no additive combinations)
m17 <- lm(Mass ~ OverallNeonic*Year, data = leye)
m18 <- lm(Mass ~ OverallNeonic*Sex, data = leye)
m19 <- lm(Mass ~ OverallNeonic*Season, data = leye)
m20 <- lm(Mass ~ OverallNeonic*Julian, data = leye)
m21 <- lm(Mass ~ OverallNeonic*MigDate, data = leye)

### Two-way interactions with a single additive combinations of informed and global model covariates 
m22 <- lm(Mass ~ OverallNeonic*Year + MigDate, data = leye)
m23 <- lm(Mass ~ OverallNeonic*Year + Julian, data = leye)
m24 <- lm(Mass ~ OverallNeonic*Year + Sex, data = leye)

m25 <- lm(Mass ~ OverallNeonic*Sex + Year, data = leye)
m26 <- lm(Mass ~ OverallNeonic*Sex + Season, data = leye)
m27 <- lm(Mass ~ OverallNeonic*Sex + MigDate, data = leye)
m28 <- lm(Mass ~ OverallNeonic*Sex + Julian, data = leye)

m29 <- lm(Mass ~ OverallNeonic*Season + MigDate, data = leye)
m30 <- lm(Mass ~ OverallNeonic*Season + Sex, data = leye)

m31 <- lm(Mass ~ OverallNeonic*Julian + Year, data = leye)
m32 <- lm(Mass ~ OverallNeonic*Julian + Sex, data = leye)

m33 <- lm(Mass ~ OverallNeonic*MigDate + Year, data = leye)
m34 <- lm(Mass ~ OverallNeonic*MigDate + Season, data = leye)
m35 <- lm(Mass ~ OverallNeonic*MigDate + Sex, data = leye)

### Two-way interactions with two additive combination of informed and global model covariates
# Year
m36 <- lm(Mass ~ OverallNeonic + (Year * Sex) + Julian, data = leye)
m37 <- lm(Mass ~ (OverallNeonic * Year) + Sex + Julian, data = leye)
m38 <- lm(Mass ~ (OverallNeonic * Sex) + Year + Julian, data = leye)
m39 <- lm(Mass ~ (OverallNeonic * Julian) + Sex + Julian, data = leye)

m40 <- lm(Mass ~ OverallNeonic + (Year * Sex) + MigDate, data = leye)
m41 <- lm(Mass ~ (OverallNeonic * Year) + Sex + MigDate, data = leye)
m42 <- lm(Mass ~ (OverallNeonic * Sex) + Year + MigDate, data = leye)
m43 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Year, data = leye)

m44 <- lm(Mass ~ OverallNeonic + (Year * Julian) + Sex, data = leye)
m45 <- lm(Mass ~ (OverallNeonic * Year) + Julian + Sex, data = leye)
m46 <- lm(Mass ~ (OverallNeonic * Julian) + Year + Sex, data = leye)
m47 <- lm(Mass ~ (OverallNeonic * Sex) + Julian + Year, data = leye)

# Sex
m48 <- lm(Mass ~ OverallNeonic + (Sex * Season) + MigDate, data = leye)
m49 <- lm(Mass ~ (OverallNeonic * Sex) + Season + MigDate, data = leye)
m50 <- lm(Mass ~ (OverallNeonic * Season) + Sex + MigDate, data = leye)
m51 <- lm(Mass ~ (OverallNeonic * MigDate) + Season + Sex, data = leye)

m52 <- lm(Mass ~ OverallNeonic + (Sex * Julian) + Year, data = leye)
m53 <- lm(Mass ~ (OverallNeonic * Sex) + Julian + Year, data = leye)
m54 <- lm(Mass ~ (OverallNeonic * Julian) + Sex + Year, data = leye)
m55 <- lm(Mass ~ (OverallNeonic * Year) + Julian + Sex, data = leye)

m56 <- lm(Mass ~ OverallNeonic + (Sex * MigDate) + Year, data = leye)
m57 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Year, data = leye)
m58 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Year, data = leye)
m59 <- lm(Mass ~ (OverallNeonic * Year) + MigDate + Sex, data = leye)

m60 <- lm(Mass ~ OverallNeonic + (Sex * MigDate) + Season, data = leye)
m61 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Season, data = leye)
m62 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Season, data = leye)
m63 <- lm(Mass ~ (OverallNeonic * Season) + MigDate + Sex, data = leye)

### Neonicotinoid concentrations with one interaction
m64 <- lm(Mass ~ OverallNeonic + Year*Sex, data = leye)
m65 <- lm(Mass ~ OverallNeonic + Year*Season, data = leye)
m66 <- lm(Mass ~ OverallNeonic + Year*Julian, data = leye)
m67 <- lm(Mass ~ OverallNeonic + Year*MigDate, data = leye)
m68 <- lm(Mass ~ OverallNeonic + Sex*Season, data = leye)
m69 <- lm(Mass ~ OverallNeonic + Sex*Julian, data = leye)
m70 <- lm(Mass ~ OverallNeonic + Sex*MigDate, data = leye)
m71 <- lm(Mass ~ OverallNeonic + Season*MigDate, data = leye)

### Neonicotinoid concentrations with two interactions
m72 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*Julian, data = leye)
m73 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*MigDate, data = leye)
m74 <- lm(Mass ~ OverallNeonic + Year*Sex + Sex*Julian, data = leye)
m75 <- lm(Mass ~ OverallNeonic + Year*Sex + Sex*MigDate, data = leye)
m76 <- lm(Mass ~ OverallNeonic + Year*Julian + Year*Sex, data = leye)
m77 <- lm(Mass ~ OverallNeonic + Year*Julian + Sex*Julian, data = leye)
m78 <- lm(Mass ~ OverallNeonic + Year*MigDate + Sex*MigDate, data = leye)
m79 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*MigDate, data = leye)
m80 <- lm(Mass ~ OverallNeonic + Sex*MigDate + Season*MigDate, data = leye)

### Neonicotinoid concentrations with three interactions
m81 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*Julian + Sex*Julian, data = leye)
m82 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*MigDate + Sex*MigDate, data = leye)
m83 <- lm(Mass ~ OverallNeonic + Year*Julian + Year*Sex + Sex*Julian, data = leye)
m84 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

### Add neonicotinoid interaction with all other interaction combinations
# Add single interaction
#Year
m85 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex, data = leye)
m86 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian, data = leye)
m87 <- lm(Mass ~ (OverallNeonic * Year) + Year*MigDate, data = leye)
m88 <- lm(Mass ~ (OverallNeonic * Year) + Sex*Season, data = leye)
m89 <- lm(Mass ~ (OverallNeonic * Year) + Sex*Julian, data = leye)
m90 <- lm(Mass ~ (OverallNeonic * Year) + Sex*MigDate, data = leye)
m91 <- lm(Mass ~ (OverallNeonic * Year) + Season*MigDate, data = leye)

#Sex
m92 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex, data = leye)
m93 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian, data = leye)
m94 <- lm(Mass ~ (OverallNeonic * Sex) + Year*MigDate, data = leye)
m95 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season, data = leye)
m96 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Julian, data = leye)
m97 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*MigDate, data = leye)
m98 <- lm(Mass ~ (OverallNeonic * Sex) + Season*MigDate, data = leye)

#Season
m99 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season, data = leye)
m100 <- lm(Mass ~ (OverallNeonic * Season) + Sex*MigDate, data = leye)
m101 <- lm(Mass ~ (OverallNeonic * Season) + Season*MigDate, data = leye)

#Julian
m102 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex, data = leye)
m103 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian, data = leye)
m104 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Julian, data = leye)

#MigDate
m105 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex, data = leye)
m106 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*MigDate, data = leye)
m107 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season, data = leye)
m108 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*MigDate, data = leye)
m109 <- lm(Mass ~ (OverallNeonic * MigDate) + Season*MigDate, data = leye)

# Add two interactions
### Neonicotinoid concentrations with two interactions
# Year
m110 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*Julian, data = leye)
m111 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*MigDate, data = leye)
m112 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Sex*Julian, data = leye)
m113 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Sex*MigDate, data = leye)
m114 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian + Year*Sex, data = leye)
m115 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian + Sex*Julian, data = leye)
m116 <- lm(Mass ~ (OverallNeonic * Year) + Year*MigDate + Sex*MigDate, data = leye)

# Sex
m117 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*Julian, data = leye)
m118 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*MigDate, data = leye)
m119 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Sex*Julian, data = leye)
m120 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Sex*MigDate, data = leye)
m121 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Year*Sex, data = leye)
m122 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Sex*Julian, data = leye)
m123 <- lm(Mass ~ (OverallNeonic * Sex) + Year*MigDate + Sex*MigDate, data = leye)
m124 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*MigDate, data = leye)
m125 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*MigDate + Season*MigDate, data = leye)

# Season
m126 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*MigDate, data = leye)
m127 <- lm(Mass ~ (OverallNeonic * Season) + Sex*MigDate + Season*MigDate, data = leye)

# Julian
m128 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex + Year*Julian, data = leye)
m129 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex + Sex*Julian, data = leye)
m130 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian + Year*Sex, data = leye)
m131 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian + Sex*Julian, data = leye)

# MigDate
m132 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex + Year*MigDate, data = leye)
m133 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex + Sex*MigDate, data = leye)
m134 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*MigDate + Sex*MigDate, data = leye)
m135 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*MigDate, data = leye)
m136 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*MigDate + Season*MigDate, data = leye)

# Add three interactions
# Year
m137 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*Julian + Sex*Julian, data = leye)
m138 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*MigDate + Sex*MigDate, data = leye)

# Sex
m139 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*Julian + Sex*Julian, data = leye)
m140 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*MigDate + Sex*MigDate, data = leye)
m141 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Year*Sex + Sex*Julian, data = leye)
m142 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

# Season
m143 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

# Julian
m144 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex + Year*Julian + Sex*Julian, data = leye)
m145 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian + Year*Sex + Sex*Julian, data = leye)

# MigDate
m146 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex + Year*MigDate + Sex*MigDate, data = leye)
m147 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

# Manually list the specially named models
special_models <- list(m.null, m.informednull, m.globalinformednull,
                       m.globalyearjulian, m.globalseasonmig, m.globalyearmig)

special_names <- c('m.null', 'm.informednull', 'm.globalinformednull',
                   'm.globalyearjulian', 'm.globalseasonmig', 'm.globalyearmig')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Add the other models (m1 to m147) to the list
for (i in c(1:147)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aic_table <- aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

### Top Model Summaries: Second Stage ####
cbind(summary(m.informednull)$coefficients, confint(m.informednull))

options(digits = 3)

# Mass ~ Year (m1) <- INFORMED NULL MODEL
cbind(summary(m1)$coefficients, confint(m1))

# Mass ~ OverallNeonic + Year (m2)
cbind(summary(m2)$coefficients, confint(m2)) # confidence intervals overlap 0

# ----------------------------------------------------------------------------- #

### Conclusions: Second Stage without Log Transformation ####

# Informed null model (Mass ~ Year) came out as the top model. Neonicotinoid 
# concentrations do not further explain any variation in mass beyond year. 

# ----------------------------------------------------------------------------- #

# Second Stage AICc Model Selection: Mass ~ Log10(Concentration) ####

### Null and global models
m.null <- lm(Mass ~ 1, data = leye)
m.informednull <- lm(Mass ~ Year, data = leye)
m.globalinformednull <- lm(Mass ~ Year * log10(OverallNeonic + 0.0001), data = leye)

m.globalyearjulian <- lm(Mass ~ Year + Sex + Julian + log10(OverallNeonic + 0.0001) + 
                           Year*Sex + Year*Julian + Year*log10(OverallNeonic + 0.0001) + 
                           Sex*Julian + Sex*log10(OverallNeonic + 0.0001) + 
                           Julian*log10(OverallNeonic + 0.0001), data = leye)

m.globalyearmig <- lm(Mass ~ Year + Sex + MigDate + log10(OverallNeonic + 0.0001) +
                        Year*Sex + Year*MigDate + Year*log10(OverallNeonic + 0.0001) + 
                        Sex*MigDate + Sex*log10(OverallNeonic + 0.0001) + 
                        MigDate*log10(OverallNeonic + 0.0001), data = leye)

m.globalseasonmig <- lm(Mass ~ Sex + Season + MigDate + log10(OverallNeonic + 0.0001) + 
                          Sex*Season + Sex*MigDate + Sex*log10(OverallNeonic + 0.0001) +
                          Season*MigDate + Season*log10(OverallNeonic + 0.0001) + 
                          MigDate*log10(OverallNeonic + 0.0001), data = leye)

### single combinations of informed and global model covariates with neonicotinoid concentrations
m1 <- lm(Mass ~ log10(OverallNeonic + 0.0001), data = leye)
m2 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year, data = leye)
m3 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex, data = leye)
m4 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season, data = leye)
m5 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Julian, data = leye)
m6 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + MigDate, data = leye)

### two additive combinations of informed and global model covariates with neonicotinoid concentrations
m7 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Sex, data = leye)
m8 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Julian, data = leye)
m9 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + MigDate, data = leye)

m10 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Season, data = leye)
m11 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Julian, data = leye)
m12 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + MigDate, data = leye)

m13 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season + MigDate, data = leye)

### three additive combinations of informed and global model covariates with neonicotinoid concentrations
m14 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Sex + Julian, data = leye)
m15 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Sex + MigDate, data = leye)
m16 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Season + MigDate, data = leye)

### two-way interactions (no additive combinations)
m17 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year, data = leye)
m18 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex, data = leye)
m19 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season, data = leye)
m20 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian, data = leye)
m21 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate, data = leye)

### two-way interactions with a single additive combinations of informed and global model covariates with neonicotinoid concentrations
m22 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year + MigDate, data = leye)
m23 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year + Julian, data = leye)
m24 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year + Sex, data = leye)

m25 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Year, data = leye)
m26 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Season, data = leye)
m27 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + MigDate, data = leye)
m28 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Julian, data = leye)

m29 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season + MigDate, data = leye)
m30 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season + Sex, data = leye)

m31 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian + Year, data = leye)
m32 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian + Sex, data = leye)

m33 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Year, data = leye)
m34 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Season, data = leye)
m35 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Sex, data = leye)

### two-way interactions with two additive combination of informed and global model covariates with neonicotinoid concentrations
#Year
m36 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * Sex) + Julian, data = leye)
m37 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex + Julian, data = leye)
m38 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year + Julian, data = leye)
m39 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex + Julian, data = leye)

m40 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * Sex) + MigDate, data = leye)
m41 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex + MigDate, data = leye)
m42 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year + MigDate, data = leye)
m43 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Year, data = leye)

m44 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * Julian) + Sex, data = leye)
m45 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Julian + Sex, data = leye)
m46 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year + Sex, data = leye)
m47 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Julian + Year, data = leye)

#Sex
m48 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Season) + MigDate, data = leye)
m49 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season + MigDate, data = leye)
m50 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex + MigDate, data = leye)
m51 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season + Sex, data = leye)

m52 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Julian) + Year, data = leye)
m53 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Julian + Year, data = leye)
m54 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex + Year, data = leye)
m55 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Julian + Sex, data = leye)

m56 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * MigDate) + Year, data = leye)
m57 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Year, data = leye)
m58 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Year, data = leye)
m59 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + MigDate + Sex, data = leye)

m60 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * MigDate) + Season, data = leye)
m61 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Season, data = leye)
m62 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Season, data = leye)
m63 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + MigDate + Sex, data = leye)

### neonicotinoid concentrations with one interaction (with restrictions, see above)
m64 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex, data = leye)
m65 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Season, data = leye)
m66 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian, data = leye)
m67 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*MigDate, data = leye)
m68 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season, data = leye)
m69 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Julian, data = leye)
m70 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*MigDate, data = leye)
m71 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season*MigDate, data = leye)

### neonicotinoid concentrations with two interactions (with restrictions, see above)
m72 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*Julian, data = leye)
m73 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*MigDate, data = leye)
m74 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Sex*Julian, data = leye)
m75 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Sex*MigDate, data = leye)
m76 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian + Year*Sex, data = leye)
m77 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian + Sex*Julian, data = leye)
m78 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*MigDate + Sex*MigDate, data = leye)
m79 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*MigDate, data = leye)
m80 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*MigDate + Season*MigDate, data = leye)

### neonicotinoid concentrations with three interactions (with restrictions, see above)
m81 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*Julian + Sex*Julian, data = leye)
m82 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*MigDate + Sex*MigDate, data = leye)
m83 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian + Year*Sex + Sex*Julian, data = leye)
m84 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

### add neonicotinoid interaction with all other interaction combinations
#add single interaction
#Year
m85 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex, data = leye)
m86 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian, data = leye)
m87 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*MigDate, data = leye)
m88 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*Season, data = leye)
m89 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*Julian, data = leye)
m90 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*MigDate, data = leye)
m91 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Season*MigDate, data = leye)

#Sex
m92 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex, data = leye)
m93 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian, data = leye)
m94 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*MigDate, data = leye)
m95 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season, data = leye)
m96 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Julian, data = leye)
m97 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*MigDate, data = leye)
m98 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season*MigDate, data = leye)

#Season
m99 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season, data = leye)
m100 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*MigDate, data = leye)
m101 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Season*MigDate, data = leye)

#Julian
m102 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex, data = leye)
m103 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian, data = leye)
m104 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Julian, data = leye)

#MigDate
m105 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex, data = leye)
m106 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*MigDate, data = leye)
m107 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season, data = leye)
m108 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*MigDate, data = leye)
m109 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season*MigDate, data = leye)

#add two interactions
### neonicotinoid concentrations with two interactions (with restrictions, see above)
#Year
m110 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*Julian, data = leye)
m111 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*MigDate, data = leye)
m112 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Sex*Julian, data = leye)
m113 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Sex*MigDate, data = leye)
m114 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian + Year*Sex, data = leye)
m115 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian + Sex*Julian, data = leye)
m116 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*MigDate + Sex*MigDate, data = leye)

#Sex
m117 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*Julian, data = leye)
m118 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*MigDate, data = leye)
m119 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Sex*Julian, data = leye)
m120 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Sex*MigDate, data = leye)
m121 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Year*Sex, data = leye)
m122 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Sex*Julian, data = leye)
m123 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*MigDate + Sex*MigDate, data = leye)
m124 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*MigDate, data = leye)
m125 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*MigDate + Season*MigDate, data = leye)

#Season
m126 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*MigDate, data = leye)
m127 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*MigDate + Season*MigDate, data = leye)

#Julian
m128 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex + Year*Julian, data = leye)
m129 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex + Sex*Julian, data = leye)
m130 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian + Year*Sex, data = leye)
m131 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian + Sex*Julian, data = leye)

#MigDate
m132 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex + Year*MigDate, data = leye)
m133 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex + Sex*MigDate, data = leye)
m134 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*MigDate + Sex*MigDate, data = leye)
m135 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*MigDate, data = leye)
m136 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*MigDate + Season*MigDate, data = leye)

#add three interactions
#Year
m137 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*Julian + Sex*Julian, data = leye)
m138 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*MigDate + Sex*MigDate, data = leye)

#Sex
m139 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*Julian + Sex*Julian, data = leye)
m140 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*MigDate + Sex*MigDate, data = leye)
m141 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Year*Sex + Sex*Julian, data = leye)
m142 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

#Season
m143 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

#Julian
m144 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex + Year*Julian + Sex*Julian, data = leye)
m145 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian + Year*Sex + Sex*Julian, data = leye)

#MigDate
m146 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex + Year*MigDate + Sex*MigDate, data = leye)
m147 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

# Manually list the specially named models
special_models <- list(m.null, m.informednull, m.globalinformednull,
                       m.globalyearjulian, m.globalseasonmig, m.globalyearmig)

special_names <- c('m.null', 'm.informednull', 'm.globalinformednull',
                   'm.globalyearjulian', 'm.globalseasonmig', 'm.globalyearmig')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Add the other models (m1 to m360) to the list
for (i in c(1:147)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aic_table <- aictab(cand.set = models, modnames = mod_names)

# ------------------------------------------------------------------------------ #

## Top Model Summaries: Second Stage ####

# Mass ~ Year (m1) <- INFORMED NULL MODEL
cbind(summary(m1)$coefficients, confint(m1))

# Mass ~ Year + Log10(Neonic) (m2)
cbind(summary(m2)$coefficients, confint(m2)) # confidence intervals overlap 0.

# ------------------------------------------------------------------------------ #

### Conclusion: Second Stage with Log Transformation ####
# Informed null model (Mass ~ Year) came out as the top model. Neonicotinoid 
# concentrations with log transformation do not further explain any variation in mass beyond year. 

# ------------------------------------------------------------------------------ #

# Second Stage AICc Model Selection: Mass ~ Detection ####

### Null and global models
m.null <- lm(Mass ~ 1, data = leye)
m.informednull <- lm(Mass ~ Year, data = leye)
m.globalinformednull <- lm(Mass ~ Year * Detection, data = leye)

m.globalyearjulian <- lm(Mass ~ Year + Sex + Julian + Detection + 
                           Year*Sex + Year*Julian + Year*Detection + 
                           Sex*Julian + Sex*Detection + 
                           Julian*Detection, data = leye)

m.globalyearmig <- lm(Mass ~ Year + Sex + MigDate + Detection +
                        Year*Sex + Year*MigDate + Year*Detection + 
                        Sex*MigDate + Sex*Detection + 
                        MigDate*Detection, data = leye)

m.globalseasonmig <- lm(Mass ~ Sex + Season + MigDate + Detection + 
                          Sex*Season + Sex*MigDate + Sex*Detection +
                          Season*MigDate + Season*Detection + 
                          MigDate*Detection, data = leye)

### Single combinations of informed and global model covariates
m1 <- lm(Mass ~ Detection, data = leye)
m2 <- lm(Mass ~ Detection + Year, data = leye)
m3 <- lm(Mass ~ Detection + Sex, data = leye)
m4 <- lm(Mass ~ Detection + Season, data = leye)
m5 <- lm(Mass ~ Detection + Julian, data = leye)
m6 <- lm(Mass ~ Detection + MigDate, data = leye)

### Two additive combinations of informed and global model covariates
m7 <- lm(Mass ~ Detection + Year + Sex, data = leye)
m8 <- lm(Mass ~ Detection + Year + Julian, data = leye)
m9 <- lm(Mass ~ Detection + Year + MigDate, data = leye)

m10 <- lm(Mass ~ Detection + Sex + Season, data = leye)
m11 <- lm(Mass ~ Detection + Sex + Julian, data = leye)
m12 <- lm(Mass ~ Detection + Sex + MigDate, data = leye)

m13 <- lm(Mass ~ Detection + Season + MigDate, data = leye)

### Three additive combinations of informed and global model covariates
m14 <- lm(Mass ~ Detection + Year + Sex + Julian, data = leye)
m15 <- lm(Mass ~ Detection + Year + Sex + MigDate, data = leye)
m16 <- lm(Mass ~ Detection + Sex + Season + MigDate, data = leye)

### Two-way interactions (no additive combinations)
m17 <- lm(Mass ~ Detection*Year, data = leye)
m18 <- lm(Mass ~ Detection*Sex, data = leye)
m19 <- lm(Mass ~ Detection*Season, data = leye)
m20 <- lm(Mass ~ Detection*Julian, data = leye)
m21 <- lm(Mass ~ Detection*MigDate, data = leye)

### Two-way interactions with a single additive combinations of informed and global model covariates
m22 <- lm(Mass ~ Detection*Year + MigDate, data = leye)
m23 <- lm(Mass ~ Detection*Year + Julian, data = leye)
m24 <- lm(Mass ~ Detection*Year + Sex, data = leye)

m25 <- lm(Mass ~ Detection*Sex + Year, data = leye)
m26 <- lm(Mass ~ Detection*Sex + Season, data = leye)
m27 <- lm(Mass ~ Detection*Sex + MigDate, data = leye)
m28 <- lm(Mass ~ Detection*Sex + Julian, data = leye)

m29 <- lm(Mass ~ Detection*Season + MigDate, data = leye)
m30 <- lm(Mass ~ Detection*Season + Sex, data = leye)

m31 <- lm(Mass ~ Detection*Julian + Year, data = leye)
m32 <- lm(Mass ~ Detection*Julian + Sex, data = leye)

m33 <- lm(Mass ~ Detection*MigDate + Year, data = leye)
m34 <- lm(Mass ~ Detection*MigDate + Season, data = leye)
m35 <- lm(Mass ~ Detection*MigDate + Sex, data = leye)

### Two-way interactions with two additive combination of informed and global model covariates
#Year
m36 <- lm(Mass ~ Detection + (Year * Sex) + Julian, data = leye)
m37 <- lm(Mass ~ (Detection * Year) + Sex + Julian, data = leye)
m38 <- lm(Mass ~ (Detection * Sex) + Year + Julian, data = leye)
m39 <- lm(Mass ~ (Detection * Julian) + Sex + Julian, data = leye)

m40 <- lm(Mass ~ Detection + (Year * Sex) + MigDate, data = leye)
m41 <- lm(Mass ~ (Detection * Year) + Sex + MigDate, data = leye)
m42 <- lm(Mass ~ (Detection * Sex) + Year + MigDate, data = leye)
m43 <- lm(Mass ~ (Detection * MigDate) + Sex + Year, data = leye)

m44 <- lm(Mass ~ Detection + (Year * Julian) + Sex, data = leye)
m45 <- lm(Mass ~ (Detection * Year) + Julian + Sex, data = leye)
m46 <- lm(Mass ~ (Detection * Julian) + Year + Sex, data = leye)
m47 <- lm(Mass ~ (Detection * Sex) + Julian + Year, data = leye)

#Sex
m48 <- lm(Mass ~ Detection + (Sex * Season) + MigDate, data = leye)
m49 <- lm(Mass ~ (Detection * Sex) + Season + MigDate, data = leye)
m50 <- lm(Mass ~ (Detection * Season) + Sex + MigDate, data = leye)
m51 <- lm(Mass ~ (Detection * MigDate) + Season + Sex, data = leye)

m52 <- lm(Mass ~ Detection + (Sex * Julian) + Year, data = leye)
m53 <- lm(Mass ~ (Detection * Sex) + Julian + Year, data = leye)
m54 <- lm(Mass ~ (Detection * Julian) + Sex + Year, data = leye)
m55 <- lm(Mass ~ (Detection * Year) + Julian + Sex, data = leye)

m56 <- lm(Mass ~ Detection + (Sex * MigDate) + Year, data = leye)
m57 <- lm(Mass ~ (Detection * Sex) + MigDate + Year, data = leye)
m58 <- lm(Mass ~ (Detection * MigDate) + Sex + Year, data = leye)
m59 <- lm(Mass ~ (Detection * Year) + MigDate + Sex, data = leye)

m60 <- lm(Mass ~ Detection + (Sex * MigDate) + Season, data = leye)
m61 <- lm(Mass ~ (Detection * Sex) + MigDate + Season, data = leye)
m62 <- lm(Mass ~ (Detection * MigDate) + Sex + Season, data = leye)
m63 <- lm(Mass ~ (Detection * Season) + MigDate + Sex, data = leye)

### Neonicotinoid detections with one interaction 
m64 <- lm(Mass ~ Detection + Year*Sex, data = leye)
m65 <- lm(Mass ~ Detection + Year*Season, data = leye)
m66 <- lm(Mass ~ Detection + Year*Julian, data = leye)
m67 <- lm(Mass ~ Detection + Year*MigDate, data = leye)
m68 <- lm(Mass ~ Detection + Sex*Season, data = leye)
m69 <- lm(Mass ~ Detection + Sex*Julian, data = leye)
m70 <- lm(Mass ~ Detection + Sex*MigDate, data = leye)
m71 <- lm(Mass ~ Detection + Season*MigDate, data = leye)

### Neonicotinoid detections with two interactions 
m72 <- lm(Mass ~ Detection + Year*Sex + Year*Julian, data = leye)
m73 <- lm(Mass ~ Detection + Year*Sex + Year*MigDate, data = leye)
m74 <- lm(Mass ~ Detection + Year*Sex + Sex*Julian, data = leye)
m75 <- lm(Mass ~ Detection + Year*Sex + Sex*MigDate, data = leye)
m76 <- lm(Mass ~ Detection + Year*Julian + Year*Sex, data = leye)
m77 <- lm(Mass ~ Detection + Year*Julian + Sex*Julian, data = leye)
m78 <- lm(Mass ~ Detection + Year*MigDate + Sex*MigDate, data = leye)
m79 <- lm(Mass ~ Detection + Sex*Season + Sex*MigDate, data = leye)
m80 <- lm(Mass ~ Detection + Sex*MigDate + Season*MigDate, data = leye)

### Neonicotinoid detections with three interactions (with restrictions, see above)
m81 <- lm(Mass ~ Detection + Year*Sex + Year*Julian + Sex*Julian, data = leye)
m82 <- lm(Mass ~ Detection + Year*Sex + Year*MigDate + Sex*MigDate, data = leye)
m83 <- lm(Mass ~ Detection + Year*Julian + Year*Sex + Sex*Julian, data = leye)
m84 <- lm(Mass ~ Detection + Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

### Add neonicotinoid interaction with all other interaction combinations
# Add single interaction
#Year
m85 <- lm(Mass ~ (Detection * Year) + Year*Sex, data = leye)
m86 <- lm(Mass ~ (Detection * Year) + Year*Julian, data = leye)
m87 <- lm(Mass ~ (Detection * Year) + Year*MigDate, data = leye)
m88 <- lm(Mass ~ (Detection * Year) + Sex*Season, data = leye)
m89 <- lm(Mass ~ (Detection * Year) + Sex*Julian, data = leye)
m90 <- lm(Mass ~ (Detection * Year) + Sex*MigDate, data = leye)
m91 <- lm(Mass ~ (Detection * Year) + Season*MigDate, data = leye)

#Sex
m92 <- lm(Mass ~ (Detection * Sex) + Year*Sex, data = leye)
m93 <- lm(Mass ~ (Detection * Sex) + Year*Julian, data = leye)
m94 <- lm(Mass ~ (Detection * Sex) + Year*MigDate, data = leye)
m95 <- lm(Mass ~ (Detection * Sex) + Sex*Season, data = leye)
m96 <- lm(Mass ~ (Detection * Sex) + Sex*Julian, data = leye)
m97 <- lm(Mass ~ (Detection * Sex) + Sex*MigDate, data = leye)
m98 <- lm(Mass ~ (Detection * Sex) + Season*MigDate, data = leye)

#Season
m99 <- lm(Mass ~ (Detection * Season) + Sex*Season, data = leye)
m100 <- lm(Mass ~ (Detection * Season) + Sex*MigDate, data = leye)
m101 <- lm(Mass ~ (Detection * Season) + Season*MigDate, data = leye)

#Julian
m102 <- lm(Mass ~ (Detection * Julian) + Year*Sex, data = leye)
m103 <- lm(Mass ~ (Detection * Julian) + Year*Julian, data = leye)
m104 <- lm(Mass ~ (Detection * Julian) + Sex*Julian, data = leye)

#MigDate
m105 <- lm(Mass ~ (Detection * MigDate) + Year*Sex, data = leye)
m106 <- lm(Mass ~ (Detection * MigDate) + Year*MigDate, data = leye)
m107 <- lm(Mass ~ (Detection * MigDate) + Sex*Season, data = leye)
m108 <- lm(Mass ~ (Detection * MigDate) + Sex*MigDate, data = leye)
m109 <- lm(Mass ~ (Detection * MigDate) + Season*MigDate, data = leye)

# Add two interactions
### Neonicotinoid detection with two interactions 
#Year
m110 <- lm(Mass ~ (Detection * Year) + Year*Sex + Year*Julian, data = leye)
m111 <- lm(Mass ~ (Detection * Year) + Year*Sex + Year*MigDate, data = leye)
m112 <- lm(Mass ~ (Detection * Year) + Year*Sex + Sex*Julian, data = leye)
m113 <- lm(Mass ~ (Detection * Year) + Year*Sex + Sex*MigDate, data = leye)
m114 <- lm(Mass ~ (Detection * Year) + Year*Julian + Year*Sex, data = leye)
m115 <- lm(Mass ~ (Detection * Year) + Year*Julian + Sex*Julian, data = leye)
m116 <- lm(Mass ~ (Detection * Year) + Year*MigDate + Sex*MigDate, data = leye)

#Sex
m117 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Year*Julian, data = leye)
m118 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Year*MigDate, data = leye)
m119 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Sex*Julian, data = leye)
m120 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Sex*MigDate, data = leye)
m121 <- lm(Mass ~ (Detection * Sex) + Year*Julian + Year*Sex, data = leye)
m122 <- lm(Mass ~ (Detection * Sex) + Year*Julian + Sex*Julian, data = leye)
m123 <- lm(Mass ~ (Detection * Sex) + Year*MigDate + Sex*MigDate, data = leye)
m124 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Sex*MigDate, data = leye)
m125 <- lm(Mass ~ (Detection * Sex) + Sex*MigDate + Season*MigDate, data = leye)

#Season
m126 <- lm(Mass ~ (Detection * Season) + Sex*Season + Sex*MigDate, data = leye)
m127 <- lm(Mass ~ (Detection * Season) + Sex*MigDate + Season*MigDate, data = leye)

#Julian
m128 <- lm(Mass ~ (Detection * Julian) + Year*Sex + Year*Julian, data = leye)
m129 <- lm(Mass ~ (Detection * Julian) + Year*Sex + Sex*Julian, data = leye)
m130 <- lm(Mass ~ (Detection * Julian) + Year*Julian + Year*Sex, data = leye)
m131 <- lm(Mass ~ (Detection * Julian) + Year*Julian + Sex*Julian, data = leye)

#MigDate
m132 <- lm(Mass ~ (Detection * MigDate) + Year*Sex + Year*MigDate, data = leye)
m133 <- lm(Mass ~ (Detection * MigDate) + Year*Sex + Sex*MigDate, data = leye)
m134 <- lm(Mass ~ (Detection * MigDate) + Year*MigDate + Sex*MigDate, data = leye)
m135 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Sex*MigDate, data = leye)
m136 <- lm(Mass ~ (Detection * MigDate) + Sex*MigDate + Season*MigDate, data = leye)

#add three interactions
#Year
m137 <- lm(Mass ~ (Detection * Year) + Year*Sex + Year*Julian + Sex*Julian, data = leye)
m138 <- lm(Mass ~ (Detection * Year) + Year*Sex + Year*MigDate + Sex*MigDate, data = leye)

#Sex
m139 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Year*Julian + Sex*Julian, data = leye)
m140 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Year*MigDate + Sex*MigDate, data = leye)
m141 <- lm(Mass ~ (Detection * Sex) + Year*Julian + Year*Sex + Sex*Julian, data = leye)
m142 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

#Season
m143 <- lm(Mass ~ (Detection * Season) + Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

#Julian
m144 <- lm(Mass ~ (Detection * Julian) + Year*Sex + Year*Julian + Sex*Julian, data = leye)
m145 <- lm(Mass ~ (Detection * Julian) + Year*Julian + Year*Sex + Sex*Julian, data = leye)

#MigDate
m146 <- lm(Mass ~ (Detection * MigDate) + Year*Sex + Year*MigDate + Sex*MigDate, data = leye)
m147 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Sex*MigDate + Season*MigDate, data = leye)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

# Manually list the specially named models
special_models <- list(m.null, m.informednull, m.globalinformednull,
                       m.globalyearjulian, m.globalseasonmig, m.globalyearmig)

special_names <- c('m.null', 'm.informednull', 'm.globalinformednull',
                   'm.globalyearjulian', 'm.globalseasonmig', 'm.globalyearmig')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Dynamically add the other models (m1 to m147) to the list
for (i in c(1:147)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

## Top Model Summaries: Second Stage ####

# Mass ~ Detection + Year (m2)
cbind(summary(m2)$coefficients, confint(m2)) # confidence intervals overlap 0 for detections.

# ----------------------------------------------------------------------------- #

### Conclusion: Second Stage with Detection ####

# Informed null model (Mass ~ Year) came out as the top model. Detection of neonicotinoids 
# does not further explain any variation in mass beyond year. Birds in 2023 had
# significantly lower body mass but this is not related to neonicotinoid concentrations
# or detections.

# ----------------------------------------------------------------------------- #

# Model Comparisons (m2): Concentration, Log(Concentration), Detection ####
m.null <- lm(Mass ~ 1, data = leye)
m.informednull <- lm(Mass ~ Year, data = leye)
m.c <- lm(Mass ~ OverallNeonic + Year, data = leye)
m.l <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year, data = leye)
m.d <- lm(Mass ~ Detection + Year, data = leye)

models <- list(m.c, m.l, m.d, m.null, m.informednull)

mod.names <- c('m.c', 'm.l', 'm.d', 'm.null', 'm.informednull')

aictab(models, modnames = mod.names)

# ----------------------------------------------------------------------------- #

## Top Model Summaries: Model Comparisons of m2 ####

cbind(summary(m.l)$coefficients, confint(m.l))
cbind(summary(m.c)$coefficients, confint(m.c))
cbind(summary(m.d)$coefficients, confint(m.d))

# ----------------------------------------------------------------------------- #

### Conclusion: Model Comparisons of m2 ####

# All three models performed similarily and received similar support 
# (AICc [1.82 - 2.43]). 
# Concentrations, log(concentrations), and detections do not significantly 
# impact mass beyond year.

# ----------------------------------------------------------------------------- #
