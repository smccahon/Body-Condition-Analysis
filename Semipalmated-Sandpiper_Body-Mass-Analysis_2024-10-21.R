#-------------------------------------------#
# Semipalmated Sandpiper Body Mass Analysis #
#           Created 10/21/2024              #
#          Modified 10/21/2024              #
#-------------------------------------------#

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
sesa <- subset(birds, Species %in% c("SemipalmatedSandpiper"))

# Remove birds that do not have neonicotinoid concentrations
sesa <- sesa[!is.na(sesa$OverallNeonic), ]

# Make neonicotinoid detection column
sesa$Detection <- ifelse(sesa$OverallNeonic > 0, "Detection", "Non-detection")

# Convert categorical variables to factor
sesa$Season <- as.factor(sesa$Season)
sesa$Year <- as.factor(sesa$Year)
sesa$Sex <- as.factor(sesa$Sex)
sesa$Detection <- as.factor(sesa$Detection)

# Check for variable correlation
cor(sesa$Julian, sesa$MigDate)

# ---------------------------------------------------------------------------- #

## Covariate Restrictions ####

# Date into Season and Julian can't be in the same model (correlation = 0.93)
# Year and Season can't be in the same model because of lack of replication of Season within a year
# Julian and Season can't be in the same model because Season bins Julian

# ---------------------------------------------------------------------------- #

# Data Exploration & Linear Regression Model Assumptions ####

## Assumption 1: Zero Expectations of Errors (Residuals) ####

### Residuals Plots 

# Calculate residuals
m <- lm(Mass ~ OverallNeonic, data = sesa)
sesa$yhat <- predict(m)
sesa$residuals <- residuals(m)
sesa$rstudent <- rstudent(m)
head(sesa)

# Simple residual plot
plot(predict(m), rstudent(m))

# Pretty residual plot
ggplot(sesa, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# Dr. Johnson confirmed these residual plots are good even though more 
# variable around 0. The "megaphone" pattern is not concerning as it's mostly 
# a reflection of non-detects (the opposite direction is more concerning). 

# Log transformation of neonics

# Calculate residuals 
m.log <- lm(Mass ~ (log10(OverallNeonic+0.0001)), data = sesa)
sesa$yhat <- predict(m.log)
sesa$residuals <- residuals(m.log)
sesa$rstudent <- rstudent(m.log)
head(sesa)

# Simple residual plot
plot(predict(m.log), rstudent(m.log))

# Pretty residual plot
ggplot(sesa, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# Looks better with logarithmic transformation.

# ---------------------------------------------------------------------------- #

## Assumption 2: Equality of error variances (Homoscedasticity) #### 

# No transformation
m <- lm(Mass ~ OverallNeonic, data = sesa)
sesa$yhat <- predict(m)
sesa$residuals <- residuals(m)
sesa$rstudent <- rstudent(m)

# Plot homoscedasticity
ggplot(sesa, aes(x = yhat, y = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual")

# No violations because variability doesn't increase as predicted value increases

# Log transformation
m.log <- lm(Mass ~ (log10(OverallNeonic+0.0001)), data = sesa)
sesa$yhat <- predict(m.log)
sesa$residuals <- residuals(m.log)
sesa$rstudent <- rstudent(m.log)

### PLOT HOMOSCEDASTICITY
ggplot(sesa, aes(x = yhat, y = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual")

# No violations because variability doesn't increase as predicted value increases.
# Log transformation does appear to spread residuals more evenly. 

# ---------------------------------------------------------------------------- #

## Assumption 3: Independence of Errors/Observations ####

# Distribution of one error/observation does not depend on another error/observation

# ---------------------------------------------------------------------------- #

## Assumption 4: Normality ####

m <- lm(Mass ~ OverallNeonic, data = sesa)
hist(rstandard(m))
qqnorm(rstandard(m))

# Does not severely violate normality assumption
# Regression models do not many any assumptions about the distrbution of explanatory 
# variables (e.g., Concentrations), but their distribution does affect inferences. 

# ---------------------------------------------------------------------------- #

# First Stage AICc Model Selection: Mass ####

### Controlling Variables: Year + Sex + Season + Julian + MigDate

### Null and global models 
m.null <- lm(Mass ~ 1, data= sesa)
m.globalyearjulian <- lm(Mass ~ Year + Sex + Julian + 
                           Year*Sex + Year*Julian +  
                           Sex*Julian, data = sesa)

m.globalyearmig <- lm(Mass ~ Year + Sex + MigDate +
                        Year*Sex + Year*MigDate +
                        Sex*MigDate, data = sesa)

m.globalseasonmig <- lm(Mass ~ Sex + Season + MigDate +
                          Sex*Season + Sex*MigDate +
                          Season*MigDate, data = sesa)

### additive combinations (with restrictions, see above)
m1 <- lm(Mass ~ Year, data = sesa)
m2 <- lm(Mass ~ Season, data = sesa)
m3 <- lm(Mass ~ Sex, data = sesa)
m4 <- lm(Mass ~ Julian, data = sesa)
m5 <- lm(Mass ~ MigDate, data = sesa)

m6 <- lm(Mass ~ Year + Sex, data = sesa)
m7 <- lm(Mass ~ Year + Julian, data = sesa)
m8  <- lm(Mass ~ Year + MigDate, data = sesa)
m9 <- lm(Mass ~ Season + Sex, data = sesa)
m10 <- lm(Mass ~ Season + MigDate, data = sesa)
m11 <- lm(Mass ~ Sex + Julian, data = sesa)
m12 <- lm(Mass ~ Sex + MigDate, data = sesa)

m13 <- lm(Mass ~ Year + Sex + Julian, data = sesa)
m14 <- lm(Mass ~ Year + Sex + MigDate, data = sesa)

m15 <- lm(Mass ~ Sex + Season + MigDate, data = sesa)

### two-way interactions (with restrictions, see above)
m16 <- lm(Mass ~ Year * Sex, data = sesa)
m17 <- lm(Mass ~ Year * Julian, data = sesa)
m18 <- lm(Mass ~ Year * MigDate, data = sesa)

m19 <- lm(Mass ~ Sex * Season, data = sesa)
m20 <- lm(Mass ~ Sex * Julian, data = sesa)
m21 <- lm(Mass ~ Sex * MigDate, data = sesa)

m22 <- lm(Mass ~ Season * MigDate, data = sesa)

### two-way interactions with one additive combination (with restrictions, see above)

#Year
m23 <- lm(Mass ~ (Year * Sex) + Julian, data = sesa)
m24 <- lm(Mass ~ (Year * Sex) + MigDate, data = sesa)
m25 <- lm(Mass ~ (Year * Julian) + Sex, data = sesa)
m26 <- lm(Mass ~ (Year * MigDate) + Sex, data = sesa)

#Sex
m27 <- lm(Mass ~ (Sex * Season) + MigDate, data = sesa)
m28 <- lm(Mass ~ (Sex * Julian) + Year, data = sesa)
m29 <- lm(Mass ~ (Sex * MigDate) + Year, data = sesa)
m30 <- lm(Mass ~ (Sex * MigDate) + Season, data = sesa)

#Season
m31 <- lm(Mass ~ (Season * MigDate) + Sex, data = sesa)

### two-way interactions with two additive combinations (with restrictions, see above)
m32 <- lm(Mass ~ (Year * Sex) + Julian, data = sesa)
m33 <- lm(Mass ~ Year + (Sex * Julian), data = sesa)

m34 <- lm(Mass ~ (Year * Julian) + Sex, data = sesa)

m35 <- lm(Mass ~ (Year * Sex) + MigDate, data = sesa)

m36 <- lm(Mass ~ (Year * MigDate) + Sex, data = sesa)
m37 <- lm(Mass ~ Year + (MigDate * Sex), data = sesa)

m38 <- lm(Mass ~ (Sex * Season) + MigDate, data = sesa)
m39 <- lm(Mass ~ Sex + (Season * MigDate), data = sesa)

m40 <- lm(Mass ~ (Sex * MigDate) + Season, data = sesa)

### two-way interactions with multiple interactions (with restrictions, see above)

#two-way
m41 <- lm(Mass ~ Year*Sex + Year*Julian, data = sesa)
m42 <- lm(Mass ~ Year*Sex + Year*MigDate, data = sesa)
m43 <- lm(Mass ~ Year*Sex + Sex*Julian, data = sesa)
m44 <- lm(Mass ~ Year*Sex + Sex*MigDate, data = sesa)
m45 <- lm(Mass ~ Year*Julian + Year*Sex, data = sesa)
m46 <- lm(Mass ~ Year*Julian + Sex*Julian, data = sesa)
m47 <- lm(Mass ~ Year*MigDate + Sex*MigDate, data = sesa)
m48 <- lm(Mass ~ Sex*Season + Sex*MigDate, data = sesa)
m49 <- lm(Mass ~ Sex*MigDate + Season*MigDate, data = sesa)

#three-way
m50 <- lm(Mass ~ Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m51 <- lm(Mass ~ Year*Julian + Year*Sex + Sex*Julian, data = sesa)
m52 <- lm(Mass ~ Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

# ---------------------------------------------------------------------------- #

## AICc Model Selection: First Stage ####

models <- list(m.null, m.globalyearjulian, m.globalseasonmig, m.globalyearmig, 
               m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, 
               m16, m17, m18, m19, m20, m21, m22, m23, m24, m25, m26, m27, m28, 
               m29, m30, m31, m32, m33, m34, m35, m36, m37, m38, m39, m40, m41,
               m42, m43, m44, m45, m46, m47, m48, m49, m50, m51, m52)

mod.names <- c('m.null', 'm.globalyearjulian', 'm.globalseasonmig', 'm.globalyearmig', 
               'm1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 
               'm12', 'm13', 'm14', 'm15', 'm16', 'm17', 'm18', 'm19', 'm20', 
               'm21', 'm22', 'm23', 'm24', 'm25', 'm26', 'm27', 'm28', 'm29', 
               'm30', 'm31', 'm32', 'm33', 'm34', 'm35', 'm36', 'm37', 'm38', 'm39', 'm40', 'm41',
               'm42', 'm43', 'm44', 'm45', 'm46', 'm47', 'm48', 'm49', 'm50', 'm51', 'm52')

aictab(models, modnames = mod.names)

# ---------------------------------------------------------------------------- #

### Top Model Summaries: First Stage ####
options(digits = 3)

# Mass ~ Sex 
cbind(summary(m3)$coefficients, confint(m3))

# ---------------------------------------------------------------------------- #

### Conclusion: First Stage ####
# No variables have any effect on body mass.

# ---------------------------------------------------------------------------- #

# Second Stage AICc Model Selection: Mass ~ Concentration ####

### null and global models
m.null <- lm(Mass ~ 1, data = sesa)

m.globalyearjulian <- lm(Mass ~ Year + Sex + Julian + OverallNeonic + 
                           Year*Sex + Year*Julian + Year*OverallNeonic + 
                           Sex*Julian + Sex*OverallNeonic + 
                           Julian*OverallNeonic, data = sesa)

m.globalyearmig <- lm(Mass ~ Year + Sex + MigDate + OverallNeonic +
                        Year*Sex + Year*MigDate + Year*OverallNeonic + 
                        Sex*MigDate + Sex*OverallNeonic + 
                        MigDate*OverallNeonic, data = sesa)

m.globalseasonmig <- lm(Mass ~ Sex + Season + MigDate + OverallNeonic + 
                          Sex*Season + Sex*MigDate + Sex*OverallNeonic +
                          Season*MigDate + Season*OverallNeonic + 
                          MigDate*OverallNeonic, data = sesa)

### single combinations of informed and global model covariates with neonicotinoid concentrations
m1 <- lm(Mass ~ OverallNeonic, data = sesa)
m2 <- lm(Mass ~ OverallNeonic + Year, data = sesa)
m3 <- lm(Mass ~ OverallNeonic + Sex, data = sesa)
m4 <- lm(Mass ~ OverallNeonic + Season, data = sesa)
m5 <- lm(Mass ~ OverallNeonic + Julian, data = sesa)
m6 <- lm(Mass ~ OverallNeonic + MigDate, data = sesa)

### two additive combinations of informed and global model covariates with neonicotinoid concentrations
m7 <- lm(Mass ~ OverallNeonic + Year + Sex, data = sesa)
m8 <- lm(Mass ~ OverallNeonic + Year + Julian, data = sesa)
m9 <- lm(Mass ~ OverallNeonic + Year + MigDate, data = sesa)

m10 <- lm(Mass ~ OverallNeonic + Sex + Season, data = sesa)
m11 <- lm(Mass ~ OverallNeonic + Sex + Julian, data = sesa)
m12 <- lm(Mass ~ OverallNeonic + Sex + MigDate, data = sesa)

m13 <- lm(Mass ~ OverallNeonic + Season + MigDate, data = sesa)

### three additive combinations of informed and global model covariates with neonicotinoid concentrations
m14 <- lm(Mass ~ OverallNeonic + Year + Sex + Julian, data = sesa)
m15 <- lm(Mass ~ OverallNeonic + Year + Sex + MigDate, data = sesa)
m16 <- lm(Mass ~ OverallNeonic + Sex + Season + MigDate, data = sesa)

### two-way interactions (no additive combinations)
m17 <- lm(Mass ~ OverallNeonic*Year, data = sesa)
m18 <- lm(Mass ~ OverallNeonic*Sex, data = sesa)
m19 <- lm(Mass ~ OverallNeonic*Season, data = sesa)
m20 <- lm(Mass ~ OverallNeonic*Julian, data = sesa)
m21 <- lm(Mass ~ OverallNeonic*MigDate, data = sesa)

### two-way interactions with a single additive combinations of informed and global model covariates with neonicotinoid concentrations
m22 <- lm(Mass ~ OverallNeonic*Year + MigDate, data = sesa)
m23 <- lm(Mass ~ OverallNeonic*Year + Julian, data = sesa)
m24 <- lm(Mass ~ OverallNeonic*Year + Sex, data = sesa)

m25 <- lm(Mass ~ OverallNeonic*Sex + Year, data = sesa)
m26 <- lm(Mass ~ OverallNeonic*Sex + Season, data = sesa)
m27 <- lm(Mass ~ OverallNeonic*Sex + MigDate, data = sesa)
m28 <- lm(Mass ~ OverallNeonic*Sex + Julian, data = sesa)

m29 <- lm(Mass ~ OverallNeonic*Season + MigDate, data = sesa)
m30 <- lm(Mass ~ OverallNeonic*Season + Sex, data = sesa)

m31 <- lm(Mass ~ OverallNeonic*Julian + Year, data = sesa)
m32 <- lm(Mass ~ OverallNeonic*Julian + Sex, data = sesa)

m33 <- lm(Mass ~ OverallNeonic*MigDate + Year, data = sesa)
m34 <- lm(Mass ~ OverallNeonic*MigDate + Season, data = sesa)
m35 <- lm(Mass ~ OverallNeonic*MigDate + Sex, data = sesa)

### two-way interactions with two additive combination of informed and global model covariates with neonicotinoid concentrations
#Year
m36 <- lm(Mass ~ OverallNeonic + (Year * Sex) + Julian, data = sesa)
m37 <- lm(Mass ~ (OverallNeonic * Year) + Sex + Julian, data = sesa)
m38 <- lm(Mass ~ (OverallNeonic * Sex) + Year + Julian, data = sesa)
m39 <- lm(Mass ~ (OverallNeonic * Julian) + Sex + Julian, data = sesa)

m40 <- lm(Mass ~ OverallNeonic + (Year * Sex) + MigDate, data = sesa)
m41 <- lm(Mass ~ (OverallNeonic * Year) + Sex + MigDate, data = sesa)
m42 <- lm(Mass ~ (OverallNeonic * Sex) + Year + MigDate, data = sesa)
m43 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Year, data = sesa)

m44 <- lm(Mass ~ OverallNeonic + (Year * Julian) + Sex, data = sesa)
m45 <- lm(Mass ~ (OverallNeonic * Year) + Julian + Sex, data = sesa)
m46 <- lm(Mass ~ (OverallNeonic * Julian) + Year + Sex, data = sesa)
m47 <- lm(Mass ~ (OverallNeonic * Sex) + Julian + Year, data = sesa)

#Sex
m48 <- lm(Mass ~ OverallNeonic + (Sex * Season) + MigDate, data = sesa)
m49 <- lm(Mass ~ (OverallNeonic * Sex) + Season + MigDate, data = sesa)
m50 <- lm(Mass ~ (OverallNeonic * Season) + Sex + MigDate, data = sesa)
m51 <- lm(Mass ~ (OverallNeonic * MigDate) + Season + Sex, data = sesa)

m52 <- lm(Mass ~ OverallNeonic + (Sex * Julian) + Year, data = sesa)
m53 <- lm(Mass ~ (OverallNeonic * Sex) + Julian + Year, data = sesa)
m54 <- lm(Mass ~ (OverallNeonic * Julian) + Sex + Year, data = sesa)
m55 <- lm(Mass ~ (OverallNeonic * Year) + Julian + Sex, data = sesa)

m56 <- lm(Mass ~ OverallNeonic + (Sex * MigDate) + Year, data = sesa)
m57 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Year, data = sesa)
m58 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Year, data = sesa)
m59 <- lm(Mass ~ (OverallNeonic * Year) + MigDate + Sex, data = sesa)

m60 <- lm(Mass ~ OverallNeonic + (Sex * MigDate) + Season, data = sesa)
m61 <- lm(Mass ~ (OverallNeonic * Sex) + MigDate + Season, data = sesa)
m62 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex + Season, data = sesa)
m63 <- lm(Mass ~ (OverallNeonic * Season) + MigDate + Sex, data = sesa)

### neonicotinoid concentrations with one interaction (with restrictions, see above)
m64 <- lm(Mass ~ OverallNeonic + Year*Sex, data = sesa)
m65 <- lm(Mass ~ OverallNeonic + Year*Season, data = sesa)
m66 <- lm(Mass ~ OverallNeonic + Year*Julian, data = sesa)
m67 <- lm(Mass ~ OverallNeonic + Year*MigDate, data = sesa)
m68 <- lm(Mass ~ OverallNeonic + Sex*Season, data = sesa)
m69 <- lm(Mass ~ OverallNeonic + Sex*Julian, data = sesa)
m70 <- lm(Mass ~ OverallNeonic + Sex*MigDate, data = sesa)
m71 <- lm(Mass ~ OverallNeonic + Season*MigDate, data = sesa)

### neonicotinoid concentrations with two interactions (with restrictions, see above)
m72 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*Julian, data = sesa)
m73 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*MigDate, data = sesa)
m74 <- lm(Mass ~ OverallNeonic + Year*Sex + Sex*Julian, data = sesa)
m75 <- lm(Mass ~ OverallNeonic + Year*Sex + Sex*MigDate, data = sesa)
m76 <- lm(Mass ~ OverallNeonic + Year*Julian + Year*Sex, data = sesa)
m77 <- lm(Mass ~ OverallNeonic + Year*Julian + Sex*Julian, data = sesa)
m78 <- lm(Mass ~ OverallNeonic + Year*MigDate + Sex*MigDate, data = sesa)
m79 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*MigDate, data = sesa)
m80 <- lm(Mass ~ OverallNeonic + Sex*MigDate + Season*MigDate, data = sesa)

### neonicotinoid concentrations with three interactions (with restrictions, see above)
m81 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m82 <- lm(Mass ~ OverallNeonic + Year*Sex + Year*MigDate + Sex*MigDate, data = sesa)
m83 <- lm(Mass ~ OverallNeonic + Year*Julian + Year*Sex + Sex*Julian, data = sesa)
m84 <- lm(Mass ~ OverallNeonic + Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

### add neonicotinoid interaction with all other interaction combinations
#add single interaction
#Year
m85 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex, data = sesa)
m86 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian, data = sesa)
m87 <- lm(Mass ~ (OverallNeonic * Year) + Year*MigDate, data = sesa)
m88 <- lm(Mass ~ (OverallNeonic * Year) + Sex*Season, data = sesa)
m89 <- lm(Mass ~ (OverallNeonic * Year) + Sex*Julian, data = sesa)
m90 <- lm(Mass ~ (OverallNeonic * Year) + Sex*MigDate, data = sesa)
m91 <- lm(Mass ~ (OverallNeonic * Year) + Season*MigDate, data = sesa)

#Sex
m92 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex, data = sesa)
m93 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian, data = sesa)
m94 <- lm(Mass ~ (OverallNeonic * Sex) + Year*MigDate, data = sesa)
m95 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season, data = sesa)
m96 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Julian, data = sesa)
m97 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*MigDate, data = sesa)
m98 <- lm(Mass ~ (OverallNeonic * Sex) + Season*MigDate, data = sesa)

#Season
m99 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season, data = sesa)
m100 <- lm(Mass ~ (OverallNeonic * Season) + Sex*MigDate, data = sesa)
m101 <- lm(Mass ~ (OverallNeonic * Season) + Season*MigDate, data = sesa)

#Julian
m102 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex, data = sesa)
m103 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian, data = sesa)
m104 <- lm(Mass ~ (OverallNeonic * Julian) + Sex*Julian, data = sesa)

#MigDate
m105 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex, data = sesa)
m106 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*MigDate, data = sesa)
m107 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season, data = sesa)
m108 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*MigDate, data = sesa)
m109 <- lm(Mass ~ (OverallNeonic * MigDate) + Season*MigDate, data = sesa)

#add two interactions
### neonicotinoid concentrations with two interactions (with restrictions, see above)
#Year
m110 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*Julian, data = sesa)
m111 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*MigDate, data = sesa)
m112 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Sex*Julian, data = sesa)
m113 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Sex*MigDate, data = sesa)
m114 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian + Year*Sex, data = sesa)
m115 <- lm(Mass ~ (OverallNeonic * Year) + Year*Julian + Sex*Julian, data = sesa)
m116 <- lm(Mass ~ (OverallNeonic * Year) + Year*MigDate + Sex*MigDate, data = sesa)

#Sex
m117 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*Julian, data = sesa)
m118 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*MigDate, data = sesa)
m119 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Sex*Julian, data = sesa)
m120 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Sex*MigDate, data = sesa)
m121 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Year*Sex, data = sesa)
m122 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Sex*Julian, data = sesa)
m123 <- lm(Mass ~ (OverallNeonic * Sex) + Year*MigDate + Sex*MigDate, data = sesa)
m124 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*MigDate, data = sesa)
m125 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*MigDate + Season*MigDate, data = sesa)

#Season
m126 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*MigDate, data = sesa)
m127 <- lm(Mass ~ (OverallNeonic * Season) + Sex*MigDate + Season*MigDate, data = sesa)

#Julian
m128 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex + Year*Julian, data = sesa)
m129 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex + Sex*Julian, data = sesa)
m130 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian + Year*Sex, data = sesa)
m131 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian + Sex*Julian, data = sesa)

#MigDate
m132 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex + Year*MigDate, data = sesa)
m133 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex + Sex*MigDate, data = sesa)
m134 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*MigDate + Sex*MigDate, data = sesa)
m135 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*MigDate, data = sesa)
m136 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*MigDate + Season*MigDate, data = sesa)

#add three interactions
#Year
m137 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m138 <- lm(Mass ~ (OverallNeonic * Year) + Year*Sex + Year*MigDate + Sex*MigDate, data = sesa)

#Sex
m139 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m140 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Sex + Year*MigDate + Sex*MigDate, data = sesa)
m141 <- lm(Mass ~ (OverallNeonic * Sex) + Year*Julian + Year*Sex + Sex*Julian, data = sesa)
m142 <- lm(Mass ~ (OverallNeonic * Sex) + Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

#Season
m143 <- lm(Mass ~ (OverallNeonic * Season) + Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

#Julian
m144 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m145 <- lm(Mass ~ (OverallNeonic * Julian) + Year*Julian + Year*Sex + Sex*Julian, data = sesa)

#MigDate
m146 <- lm(Mass ~ (OverallNeonic * MigDate) + Year*Sex + Year*MigDate + Sex*MigDate, data = sesa)
m147 <- lm(Mass ~ (OverallNeonic * MigDate) + Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

#Manually list the special models
special_models <- list(m.null,
                       m.globalyearjulian, m.globalseasonmig, m.globalyearmig)

special_names <- c('m.null', 
                   'm.globalyearjulian', 'm.globalseasonmig', 'm.globalyearmig')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Dynamically add the other models (m1 to m360) to the list
for (i in c(1:147)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

### Top Model Summaries: Second Stage ####

options(digits = 3)

# Mass ~ OverallNeonic
cbind(summary(m1)$coefficients, confint(m1))

# ----------------------------------------------------------------------------- #

### Conclusion: First Stage ####
# No variables have any effect on body mass.

# ---------------------------------------------------------------------------- #

# Second Stage AICc Model Selection: Mass ~ Log(Concentration) ####

### null and global models
m.null <- lm(Mass ~ 1, data = sesa)
m.globalyearjulian <- lm(Mass ~ Year + Sex + Julian + log10(OverallNeonic + 0.0001) + 
                           Year*Sex + Year*Julian + Year*log10(OverallNeonic + 0.0001) + 
                           Sex*Julian + Sex*log10(OverallNeonic + 0.0001) + 
                           Julian*log10(OverallNeonic + 0.0001), data = sesa)

m.globalyearmig <- lm(Mass ~ Year + Sex + MigDate + log10(OverallNeonic + 0.0001) +
                        Year*Sex + Year*MigDate + Year*log10(OverallNeonic + 0.0001) + 
                        Sex*MigDate + Sex*log10(OverallNeonic + 0.0001) + 
                        MigDate*log10(OverallNeonic + 0.0001), data = sesa)

m.globalseasonmig <- lm(Mass ~ Sex + Season + MigDate + log10(OverallNeonic + 0.0001) + 
                          Sex*Season + Sex*MigDate + Sex*log10(OverallNeonic + 0.0001) +
                          Season*MigDate + Season*log10(OverallNeonic + 0.0001) + 
                          MigDate*log10(OverallNeonic + 0.0001), data = sesa)

### single combinations of informed and global model covariates with neonicotinoid concentrations
m1 <- lm(Mass ~ log10(OverallNeonic + 0.0001), data = sesa)
m2 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year, data = sesa)
m3 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex, data = sesa)
m4 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season, data = sesa)
m5 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Julian, data = sesa)
m6 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + MigDate, data = sesa)

### two additive combinations of informed and global model covariates with neonicotinoid concentrations
m7 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Sex, data = sesa)
m8 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Julian, data = sesa)
m9 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + MigDate, data = sesa)

m10 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Season, data = sesa)
m11 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Julian, data = sesa)
m12 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + MigDate, data = sesa)

m13 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season + MigDate, data = sesa)

### three additive combinations of informed and global model covariates with neonicotinoid concentrations
m14 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Sex + Julian, data = sesa)
m15 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year + Sex + MigDate, data = sesa)
m16 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex + Season + MigDate, data = sesa)

### two-way interactions (no additive combinations)
m17 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year, data = sesa)
m18 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex, data = sesa)
m19 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season, data = sesa)
m20 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian, data = sesa)
m21 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate, data = sesa)

### two-way interactions with a single additive combinations of informed and global model covariates with neonicotinoid concentrations
m22 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year + MigDate, data = sesa)
m23 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year + Julian, data = sesa)
m24 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Year + Sex, data = sesa)

m25 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Year, data = sesa)
m26 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Season, data = sesa)
m27 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + MigDate, data = sesa)
m28 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Sex + Julian, data = sesa)

m29 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season + MigDate, data = sesa)
m30 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Season + Sex, data = sesa)

m31 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian + Year, data = sesa)
m32 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*Julian + Sex, data = sesa)

m33 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Year, data = sesa)
m34 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Season, data = sesa)
m35 <- lm(Mass ~ log10(OverallNeonic + 0.0001)*MigDate + Sex, data = sesa)

### two-way interactions with two additive combination of informed and global model covariates with neonicotinoid concentrations
#Year
m36 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * Sex) + Julian, data = sesa)
m37 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex + Julian, data = sesa)
m38 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year + Julian, data = sesa)
m39 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex + Julian, data = sesa)

m40 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * Sex) + MigDate, data = sesa)
m41 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex + MigDate, data = sesa)
m42 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year + MigDate, data = sesa)
m43 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Year, data = sesa)

m44 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Year * Julian) + Sex, data = sesa)
m45 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Julian + Sex, data = sesa)
m46 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year + Sex, data = sesa)
m47 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Julian + Year, data = sesa)

#Sex
m48 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Season) + MigDate, data = sesa)
m49 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season + MigDate, data = sesa)
m50 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex + MigDate, data = sesa)
m51 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season + Sex, data = sesa)

m52 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * Julian) + Year, data = sesa)
m53 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Julian + Year, data = sesa)
m54 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex + Year, data = sesa)
m55 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Julian + Sex, data = sesa)

m56 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * MigDate) + Year, data = sesa)
m57 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Year, data = sesa)
m58 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Year, data = sesa)
m59 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + MigDate + Sex, data = sesa)

m60 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + (Sex * MigDate) + Season, data = sesa)
m61 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + MigDate + Season, data = sesa)
m62 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex + Season, data = sesa)
m63 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + MigDate + Sex, data = sesa)

### neonicotinoid concentrations with one interaction (with restrictions, see above)
m64 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex, data = sesa)
m65 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Season, data = sesa)
m66 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian, data = sesa)
m67 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*MigDate, data = sesa)
m68 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season, data = sesa)
m69 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Julian, data = sesa)
m70 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*MigDate, data = sesa)
m71 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Season*MigDate, data = sesa)

### neonicotinoid concentrations with two interactions (with restrictions, see above)
m72 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*Julian, data = sesa)
m73 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*MigDate, data = sesa)
m74 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Sex*Julian, data = sesa)
m75 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Sex*MigDate, data = sesa)
m76 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian + Year*Sex, data = sesa)
m77 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian + Sex*Julian, data = sesa)
m78 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*MigDate + Sex*MigDate, data = sesa)
m79 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*MigDate, data = sesa)
m80 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*MigDate + Season*MigDate, data = sesa)

### neonicotinoid concentrations with three interactions (with restrictions, see above)
m81 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m82 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Sex + Year*MigDate + Sex*MigDate, data = sesa)
m83 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Year*Julian + Year*Sex + Sex*Julian, data = sesa)
m84 <- lm(Mass ~ log10(OverallNeonic + 0.0001) + Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

### add neonicotinoid interaction with all other interaction combinations
#add single interaction
#Year
m85 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex, data = sesa)
m86 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian, data = sesa)
m87 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*MigDate, data = sesa)
m88 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*Season, data = sesa)
m89 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*Julian, data = sesa)
m90 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Sex*MigDate, data = sesa)
m91 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Season*MigDate, data = sesa)

#Sex
m92 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex, data = sesa)
m93 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian, data = sesa)
m94 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*MigDate, data = sesa)
m95 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season, data = sesa)
m96 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Julian, data = sesa)
m97 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*MigDate, data = sesa)
m98 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Season*MigDate, data = sesa)

#Season
m99 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season, data = sesa)
m100 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*MigDate, data = sesa)
m101 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Season*MigDate, data = sesa)

#Julian
m102 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex, data = sesa)
m103 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian, data = sesa)
m104 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Sex*Julian, data = sesa)

#MigDate
m105 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex, data = sesa)
m106 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*MigDate, data = sesa)
m107 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season, data = sesa)
m108 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*MigDate, data = sesa)
m109 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Season*MigDate, data = sesa)

#add two interactions
### neonicotinoid concentrations with two interactions (with restrictions, see above)
#Year
m110 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*Julian, data = sesa)
m111 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*MigDate, data = sesa)
m112 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Sex*Julian, data = sesa)
m113 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Sex*MigDate, data = sesa)
m114 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian + Year*Sex, data = sesa)
m115 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Julian + Sex*Julian, data = sesa)
m116 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*MigDate + Sex*MigDate, data = sesa)

#Sex
m117 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*Julian, data = sesa)
m118 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*MigDate, data = sesa)
m119 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Sex*Julian, data = sesa)
m120 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Sex*MigDate, data = sesa)
m121 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Year*Sex, data = sesa)
m122 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Sex*Julian, data = sesa)
m123 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*MigDate + Sex*MigDate, data = sesa)
m124 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*MigDate, data = sesa)
m125 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*MigDate + Season*MigDate, data = sesa)

#Season
m126 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*MigDate, data = sesa)
m127 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*MigDate + Season*MigDate, data = sesa)

#Julian
m128 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex + Year*Julian, data = sesa)
m129 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex + Sex*Julian, data = sesa)
m130 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian + Year*Sex, data = sesa)
m131 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian + Sex*Julian, data = sesa)

#MigDate
m132 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex + Year*MigDate, data = sesa)
m133 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex + Sex*MigDate, data = sesa)
m134 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*MigDate + Sex*MigDate, data = sesa)
m135 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*MigDate, data = sesa)
m136 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*MigDate + Season*MigDate, data = sesa)

#add three interactions
#Year
m137 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m138 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Year) + Year*Sex + Year*MigDate + Sex*MigDate, data = sesa)

#Sex
m139 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m140 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Sex + Year*MigDate + Sex*MigDate, data = sesa)
m141 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Year*Julian + Year*Sex + Sex*Julian, data = sesa)
m142 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Sex) + Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

#Season
m143 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Season) + Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

#Julian
m144 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m145 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * Julian) + Year*Julian + Year*Sex + Sex*Julian, data = sesa)

#MigDate
m146 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Year*Sex + Year*MigDate + Sex*MigDate, data = sesa)
m147 <- lm(Mass ~ (log10(OverallNeonic + 0.0001) * MigDate) + Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

#Manually list the special models
special_models <- list(m.null,
                       m.globalyearjulian, m.globalseasonmig, m.globalyearmig)

special_names <- c('m.null',
                   'm.globalyearjulian', 'm.globalseasonmig', 'm.globalyearmig')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Dynamically add the other models (m1 to m360) to the list
for (i in c(1:147)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

### Conclusion: Second Stage with Log Transformation ####

# No variable has any effect on body mass, including log neonic concentration.

# ----------------------------------------------------------------------------- #

# Second Stage AICc Model Selection: Mass ~ Detection ####

### null and global models
m.null <- lm(Mass ~ 1, data = sesa)
m.globalyearjulian <- lm(Mass ~ Year + Sex + Julian + Detection + 
                           Year*Sex + Year*Julian + Year*Detection + 
                           Sex*Julian + Sex*Detection + 
                           Julian*Detection, data = sesa)

m.globalyearmig <- lm(Mass ~ Year + Sex + MigDate + Detection +
                        Year*Sex + Year*MigDate + Year*Detection + 
                        Sex*MigDate + Sex*Detection + 
                        MigDate*Detection, data = sesa)

m.globalseasonmig <- lm(Mass ~ Sex + Season + MigDate + Detection + 
                          Sex*Season + Sex*MigDate + Sex*Detection +
                          Season*MigDate + Season*Detection + 
                          MigDate*Detection, data = sesa)

### single combinations of informed and global model covariates with neonicotinoid concentrations
m1 <- lm(Mass ~ Detection, data = sesa)
m2 <- lm(Mass ~ Detection + Year, data = sesa)
m3 <- lm(Mass ~ Detection + Sex, data = sesa)
m4 <- lm(Mass ~ Detection + Season, data = sesa)
m5 <- lm(Mass ~ Detection + Julian, data = sesa)
m6 <- lm(Mass ~ Detection + MigDate, data = sesa)

### two additive combinations of informed and global model covariates with neonicotinoid concentrations
m7 <- lm(Mass ~ Detection + Year + Sex, data = sesa)
m8 <- lm(Mass ~ Detection + Year + Julian, data = sesa)
m9 <- lm(Mass ~ Detection + Year + MigDate, data = sesa)

m10 <- lm(Mass ~ Detection + Sex + Season, data = sesa)
m11 <- lm(Mass ~ Detection + Sex + Julian, data = sesa)
m12 <- lm(Mass ~ Detection + Sex + MigDate, data = sesa)

m13 <- lm(Mass ~ Detection + Season + MigDate, data = sesa)

### three additive combinations of informed and global model covariates with neonicotinoid concentrations
m14 <- lm(Mass ~ Detection + Year + Sex + Julian, data = sesa)
m15 <- lm(Mass ~ Detection + Year + Sex + MigDate, data = sesa)
m16 <- lm(Mass ~ Detection + Sex + Season + MigDate, data = sesa)

### two-way interactions (no additive combinations)
m17 <- lm(Mass ~ Detection*Year, data = sesa)
m18 <- lm(Mass ~ Detection*Sex, data = sesa)
m19 <- lm(Mass ~ Detection*Season, data = sesa)
m20 <- lm(Mass ~ Detection*Julian, data = sesa)
m21 <- lm(Mass ~ Detection*MigDate, data = sesa)

### two-way interactions with a single additive combinations of informed and global model covariates with neonicotinoid concentrations
m22 <- lm(Mass ~ Detection*Year + MigDate, data = sesa)
m23 <- lm(Mass ~ Detection*Year + Julian, data = sesa)
m24 <- lm(Mass ~ Detection*Year + Sex, data = sesa)

m25 <- lm(Mass ~ Detection*Sex + Year, data = sesa)
m26 <- lm(Mass ~ Detection*Sex + Season, data = sesa)
m27 <- lm(Mass ~ Detection*Sex + MigDate, data = sesa)
m28 <- lm(Mass ~ Detection*Sex + Julian, data = sesa)

m29 <- lm(Mass ~ Detection*Season + MigDate, data = sesa)
m30 <- lm(Mass ~ Detection*Season + Sex, data = sesa)

m31 <- lm(Mass ~ Detection*Julian + Year, data = sesa)
m32 <- lm(Mass ~ Detection*Julian + Sex, data = sesa)

m33 <- lm(Mass ~ Detection*MigDate + Year, data = sesa)
m34 <- lm(Mass ~ Detection*MigDate + Season, data = sesa)
m35 <- lm(Mass ~ Detection*MigDate + Sex, data = sesa)

### two-way interactions with two additive combination of informed and global model covariates with neonicotinoid concentrations
#Year
m36 <- lm(Mass ~ Detection + (Year * Sex) + Julian, data = sesa)
m37 <- lm(Mass ~ (Detection * Year) + Sex + Julian, data = sesa)
m38 <- lm(Mass ~ (Detection * Sex) + Year + Julian, data = sesa)
m39 <- lm(Mass ~ (Detection * Julian) + Sex + Julian, data = sesa)

m40 <- lm(Mass ~ Detection + (Year * Sex) + MigDate, data = sesa)
m41 <- lm(Mass ~ (Detection * Year) + Sex + MigDate, data = sesa)
m42 <- lm(Mass ~ (Detection * Sex) + Year + MigDate, data = sesa)
m43 <- lm(Mass ~ (Detection * MigDate) + Sex + Year, data = sesa)

m44 <- lm(Mass ~ Detection + (Year * Julian) + Sex, data = sesa)
m45 <- lm(Mass ~ (Detection * Year) + Julian + Sex, data = sesa)
m46 <- lm(Mass ~ (Detection * Julian) + Year + Sex, data = sesa)
m47 <- lm(Mass ~ (Detection * Sex) + Julian + Year, data = sesa)

#Sex
m48 <- lm(Mass ~ Detection + (Sex * Season) + MigDate, data = sesa)
m49 <- lm(Mass ~ (Detection * Sex) + Season + MigDate, data = sesa)
m50 <- lm(Mass ~ (Detection * Season) + Sex + MigDate, data = sesa)
m51 <- lm(Mass ~ (Detection * MigDate) + Season + Sex, data = sesa)

m52 <- lm(Mass ~ Detection + (Sex * Julian) + Year, data = sesa)
m53 <- lm(Mass ~ (Detection * Sex) + Julian + Year, data = sesa)
m54 <- lm(Mass ~ (Detection * Julian) + Sex + Year, data = sesa)
m55 <- lm(Mass ~ (Detection * Year) + Julian + Sex, data = sesa)

m56 <- lm(Mass ~ Detection + (Sex * MigDate) + Year, data = sesa)
m57 <- lm(Mass ~ (Detection * Sex) + MigDate + Year, data = sesa)
m58 <- lm(Mass ~ (Detection * MigDate) + Sex + Year, data = sesa)
m59 <- lm(Mass ~ (Detection * Year) + MigDate + Sex, data = sesa)

m60 <- lm(Mass ~ Detection + (Sex * MigDate) + Season, data = sesa)
m61 <- lm(Mass ~ (Detection * Sex) + MigDate + Season, data = sesa)
m62 <- lm(Mass ~ (Detection * MigDate) + Sex + Season, data = sesa)
m63 <- lm(Mass ~ (Detection * Season) + MigDate + Sex, data = sesa)

### neonicotinoid concentrations with one interaction (with restrictions, see above)
m64 <- lm(Mass ~ Detection + Year*Sex, data = sesa)
m65 <- lm(Mass ~ Detection + Year*Season, data = sesa)
m66 <- lm(Mass ~ Detection + Year*Julian, data = sesa)
m67 <- lm(Mass ~ Detection + Year*MigDate, data = sesa)
m68 <- lm(Mass ~ Detection + Sex*Season, data = sesa)
m69 <- lm(Mass ~ Detection + Sex*Julian, data = sesa)
m70 <- lm(Mass ~ Detection + Sex*MigDate, data = sesa)
m71 <- lm(Mass ~ Detection + Season*MigDate, data = sesa)

### neonicotinoid concentrations with two interactions (with restrictions, see above)
m72 <- lm(Mass ~ Detection + Year*Sex + Year*Julian, data = sesa)
m73 <- lm(Mass ~ Detection + Year*Sex + Year*MigDate, data = sesa)
m74 <- lm(Mass ~ Detection + Year*Sex + Sex*Julian, data = sesa)
m75 <- lm(Mass ~ Detection + Year*Sex + Sex*MigDate, data = sesa)
m76 <- lm(Mass ~ Detection + Year*Julian + Year*Sex, data = sesa)
m77 <- lm(Mass ~ Detection + Year*Julian + Sex*Julian, data = sesa)
m78 <- lm(Mass ~ Detection + Year*MigDate + Sex*MigDate, data = sesa)
m79 <- lm(Mass ~ Detection + Sex*Season + Sex*MigDate, data = sesa)
m80 <- lm(Mass ~ Detection + Sex*MigDate + Season*MigDate, data = sesa)

### neonicotinoid concentrations with three interactions (with restrictions, see above)
m81 <- lm(Mass ~ Detection + Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m82 <- lm(Mass ~ Detection + Year*Sex + Year*MigDate + Sex*MigDate, data = sesa)
m83 <- lm(Mass ~ Detection + Year*Julian + Year*Sex + Sex*Julian, data = sesa)
m84 <- lm(Mass ~ Detection + Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

### add neonicotinoid interaction with all other interaction combinations
#add single interaction
#Year
m85 <- lm(Mass ~ (Detection * Year) + Year*Sex, data = sesa)
m86 <- lm(Mass ~ (Detection * Year) + Year*Julian, data = sesa)
m87 <- lm(Mass ~ (Detection * Year) + Year*MigDate, data = sesa)
m88 <- lm(Mass ~ (Detection * Year) + Sex*Season, data = sesa)
m89 <- lm(Mass ~ (Detection * Year) + Sex*Julian, data = sesa)
m90 <- lm(Mass ~ (Detection * Year) + Sex*MigDate, data = sesa)
m91 <- lm(Mass ~ (Detection * Year) + Season*MigDate, data = sesa)

#Sex
m92 <- lm(Mass ~ (Detection * Sex) + Year*Sex, data = sesa)
m93 <- lm(Mass ~ (Detection * Sex) + Year*Julian, data = sesa)
m94 <- lm(Mass ~ (Detection * Sex) + Year*MigDate, data = sesa)
m95 <- lm(Mass ~ (Detection * Sex) + Sex*Season, data = sesa)
m96 <- lm(Mass ~ (Detection * Sex) + Sex*Julian, data = sesa)
m97 <- lm(Mass ~ (Detection * Sex) + Sex*MigDate, data = sesa)
m98 <- lm(Mass ~ (Detection * Sex) + Season*MigDate, data = sesa)

#Season
m99 <- lm(Mass ~ (Detection * Season) + Sex*Season, data = sesa)
m100 <- lm(Mass ~ (Detection * Season) + Sex*MigDate, data = sesa)
m101 <- lm(Mass ~ (Detection * Season) + Season*MigDate, data = sesa)

#Julian
m102 <- lm(Mass ~ (Detection * Julian) + Year*Sex, data = sesa)
m103 <- lm(Mass ~ (Detection * Julian) + Year*Julian, data = sesa)
m104 <- lm(Mass ~ (Detection * Julian) + Sex*Julian, data = sesa)

#MigDate
m105 <- lm(Mass ~ (Detection * MigDate) + Year*Sex, data = sesa)
m106 <- lm(Mass ~ (Detection * MigDate) + Year*MigDate, data = sesa)
m107 <- lm(Mass ~ (Detection * MigDate) + Sex*Season, data = sesa)
m108 <- lm(Mass ~ (Detection * MigDate) + Sex*MigDate, data = sesa)
m109 <- lm(Mass ~ (Detection * MigDate) + Season*MigDate, data = sesa)

#add two interactions
### neonicotinoid concentrations with two interactions (with restrictions, see above)
#Year
m110 <- lm(Mass ~ (Detection * Year) + Year*Sex + Year*Julian, data = sesa)
m111 <- lm(Mass ~ (Detection * Year) + Year*Sex + Year*MigDate, data = sesa)
m112 <- lm(Mass ~ (Detection * Year) + Year*Sex + Sex*Julian, data = sesa)
m113 <- lm(Mass ~ (Detection * Year) + Year*Sex + Sex*MigDate, data = sesa)
m114 <- lm(Mass ~ (Detection * Year) + Year*Julian + Year*Sex, data = sesa)
m115 <- lm(Mass ~ (Detection * Year) + Year*Julian + Sex*Julian, data = sesa)
m116 <- lm(Mass ~ (Detection * Year) + Year*MigDate + Sex*MigDate, data = sesa)

#Sex
m117 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Year*Julian, data = sesa)
m118 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Year*MigDate, data = sesa)
m119 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Sex*Julian, data = sesa)
m120 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Sex*MigDate, data = sesa)
m121 <- lm(Mass ~ (Detection * Sex) + Year*Julian + Year*Sex, data = sesa)
m122 <- lm(Mass ~ (Detection * Sex) + Year*Julian + Sex*Julian, data = sesa)
m123 <- lm(Mass ~ (Detection * Sex) + Year*MigDate + Sex*MigDate, data = sesa)
m124 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Sex*MigDate, data = sesa)
m125 <- lm(Mass ~ (Detection * Sex) + Sex*MigDate + Season*MigDate, data = sesa)

#Season
m126 <- lm(Mass ~ (Detection * Season) + Sex*Season + Sex*MigDate, data = sesa)
m127 <- lm(Mass ~ (Detection * Season) + Sex*MigDate + Season*MigDate, data = sesa)

#Julian
m128 <- lm(Mass ~ (Detection * Julian) + Year*Sex + Year*Julian, data = sesa)
m129 <- lm(Mass ~ (Detection * Julian) + Year*Sex + Sex*Julian, data = sesa)
m130 <- lm(Mass ~ (Detection * Julian) + Year*Julian + Year*Sex, data = sesa)
m131 <- lm(Mass ~ (Detection * Julian) + Year*Julian + Sex*Julian, data = sesa)

#MigDate
m132 <- lm(Mass ~ (Detection * MigDate) + Year*Sex + Year*MigDate, data = sesa)
m133 <- lm(Mass ~ (Detection * MigDate) + Year*Sex + Sex*MigDate, data = sesa)
m134 <- lm(Mass ~ (Detection * MigDate) + Year*MigDate + Sex*MigDate, data = sesa)
m135 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Sex*MigDate, data = sesa)
m136 <- lm(Mass ~ (Detection * MigDate) + Sex*MigDate + Season*MigDate, data = sesa)

#add three interactions
#Year
m137 <- lm(Mass ~ (Detection * Year) + Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m138 <- lm(Mass ~ (Detection * Year) + Year*Sex + Year*MigDate + Sex*MigDate, data = sesa)

#Sex
m139 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m140 <- lm(Mass ~ (Detection * Sex) + Year*Sex + Year*MigDate + Sex*MigDate, data = sesa)
m141 <- lm(Mass ~ (Detection * Sex) + Year*Julian + Year*Sex + Sex*Julian, data = sesa)
m142 <- lm(Mass ~ (Detection * Sex) + Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

#Season
m143 <- lm(Mass ~ (Detection * Season) + Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

#Julian
m144 <- lm(Mass ~ (Detection * Julian) + Year*Sex + Year*Julian + Sex*Julian, data = sesa)
m145 <- lm(Mass ~ (Detection * Julian) + Year*Julian + Year*Sex + Sex*Julian, data = sesa)

#MigDate
m146 <- lm(Mass ~ (Detection * MigDate) + Year*Sex + Year*MigDate + Sex*MigDate, data = sesa)
m147 <- lm(Mass ~ (Detection * MigDate) + Sex*Season + Sex*MigDate + Season*MigDate, data = sesa)

# ----------------------------------------------------------------------------- #

## AICc Model Selection: Second Stage ####

#Manually list the special models
special_models <- list(m.null,
                       m.globalyearjulian, m.globalseasonmig, m.globalyearmig)

special_names <- c('m.null',
                   'm.globalyearjulian', 'm.globalseasonmig', 'm.globalyearmig')

# Create a list to store all models, starting with the special models
models <- special_models
mod_names <- special_names

# Dynamically add the other models (m1 to m360) to the list
for (i in c(1:147)) {
  model_name <- paste0("m", i)
  models[[length(models) + 1]] <- get(model_name)
  mod_names <- c(mod_names, model_name)
}

# Calculate AIC table
aictab(cand.set = models, modnames = mod_names)

# ----------------------------------------------------------------------------- #

### Conclusion: Second Stage with Detection ####

# No variable has any effect on body mass, including log neonic concentration.

# ----------------------------------------------------------------------------- #