#-----------------------------#
# Killdeer Body Mass Analysis #
#     Created 10/23/2024      #
#     Modified 10/21/2024     #
#-----------------------------#

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

# Subset data for Lesser Yellowlegs
kill <- subset(birds, Species %in% c("Killdeer"))

# Make neonicotinoid detection column
kill$Detection <- ifelse(kill$OverallNeonic > 0, "Detection", "Non-detection")

# Convert categorical variables to factor
kill$Year <- as.factor(kill$Year)
kill$Season <- as.factor(kill$Season)
kill$Sex <- as.factor(kill$Sex)
kill$Detection <- as.factor(kill$Detection)

# Remove birds that do not have neonicotinoid concentrations
kill <- kill[!is.na(kill$OverallNeonic), ]

# ---------------------------------------------------------------------------- #

## Covariate Restrictions ####

# Date into Season NA because Killdeer are a resident species
# Year and Season can't be in the same model because of lack of replication of Season within a year
# Julian and Season can't be in the same model because Season bins Julian

# ---------------------------------------------------------------------------- #

# Data Exploration & Linear Regression Model Assumptions ####

## Assumption 1: Zero Expectations of Errors (Residuals) ####

### Residuals Plots 

# Calculate residuals
m <- lm(Mass ~ OverallNeonic, data = kill)
kill$yhat <- predict(m)
kill$residuals <- residuals(m)
kill$rstudent <- rstudent(m)
head(kill)

# Simple residual plot
plot(predict(m), rstudent(m))

# Pretty residual plot
ggplot(kill, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# clear megaphone pattern: violates assumption

# Calculate residuals 
m.log <- lm(Mass ~ (log10(OverallNeonic+0.0001)), data = kill)
kill$yhat <- predict(m.log)
kill$residuals <- residuals(m.log)
kill$rstudent <- rstudent(m.log)
head(kill)

# Simple residual plot
plot(predict(m.log), rstudent(m.log))

# Pretty residual plot
ggplot(kill, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# clear megaphone pattern: violates assumption

# Calculate residuals
m <- lm(Mass ~ Detection, data = kill)
kill$yhat <- predict(m)
kill$residuals <- residuals(m)
kill$rstudent <- rstudent(m)
head(kill)

# Simple residual plot
plot(predict(m), rstudent(m))

# Pretty residual plot
ggplot(kill, aes(x = yhat, y = rstudent)) + geom_segment(aes(x = yhat, xend = yhat, y = 0, yend = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual") + 
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

# ---------------------------------------------------------------------------- #

## Assumption 2: Equality of error variances (Homoscedasticity) #### 

# No transformation
m <- lm(Mass ~ OverallNeonic, data = kill)
kill$yhat <- predict(m)
kill$residuals <- residuals(m)
kill$rstudent <- rstudent(m)

# Plot homoscedasticity
ggplot(kill, aes(x = yhat, y = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual")

# Log transformation
m.log <- lm(Mass ~ (log10(OverallNeonic+0.0001)), data = kill)
kill$yhat <- predict(m.log)
kill$residuals <- residuals(m.log)
kill$rstudent <- rstudent(m.log)

### PLOT HOMOSCEDASTICITY
ggplot(kill, aes(x = yhat, y = rstudent)) + 
  geom_point() + theme_classic() + labs(x = "Predicted Value", y = "Studentized Residual")

# extreme violation due to outlier; assumption violated

# ---------------------------------------------------------------------------- #

## Assumption 3: Independence of Errors/Observations ####

# Distribution of one error/observation does not depend on another error/observation

# ---------------------------------------------------------------------------- #

## Assumption 4: Normality ####

m <- lm(Mass ~ Detection, data = kill)
hist(rstandard(m))
qqnorm(rstandard(m))

# Does not severely violate normality assumption
# Regression models do not have any assumptions about the distribution of explanatory 
# variables (e.g., Concentrations), but their distribution does affect inferences. 

# ---------------------------------------------------------------------------- #

# Conclusion: Given all the violated assumptions, convert concentrations into detections and non-detections

# ---------------------------------------------------------------------------- #

# First Stage AICc Model Selection: Mass ####

### Controlling Variables: Year + Sex + Season + Julian

### null and global models
m.null <- lm(Mass ~ 1, data= kill)
m.globalyearjulian <- lm(Mass ~ Year*Sex + Year*Julian + Sex*Julian, data = kill)

### additive combinations (with restrictions, see above)
m1 <- lm(Mass ~ Year, data = kill)
m2 <- lm(Mass ~ Season, data = kill)
m3 <- lm(Mass ~ Sex, data = kill)
m4 <- lm(Mass ~ Julian, data = kill)
m5 <- lm(Mass ~ Year + Sex, data = kill)
m6 <- lm(Mass ~ Year + Julian, data = kill)
m7 <- lm(Mass ~ Season + Sex, data = kill)
m8 <- lm(Mass ~ Sex + Julian, data = kill)
m9 <- lm(Mass ~ Year + Sex + Julian, data = kill)

### two-way interactions (with restrictions, see above)
m10 <- lm(Mass ~ Year * Sex, data = kill)
m11 <- lm(Mass ~ Year * Julian, data = kill)
m12 <- lm(Mass ~ Sex * Season, data = kill)
m13 <- lm(Mass ~ Sex * Julian, data = kill)

### two-way interactions with one additive combination (with restrictions, see above)
m14 <- lm(Mass ~ (Year * Sex) + Julian, data = kill)
m15 <- lm(Mass ~ (Sex * Julian) + Year, data = kill)

### two-way interactions with multiple interactions (with restrictions, see above)
#two-way
m16 <- lm(Mass ~ Year*Sex + Year*Julian, data = kill)
m17 <- lm(Mass ~ Year*Sex + Sex*Julian, data = kill)
m18 <- lm(Mass ~ Year*Julian + Year*Sex, data = kill)
m19 <- lm(Mass ~ Year*Julian + Sex*Julian, data = kill)

#three-way
m20 <- lm(Mass ~ Year*Sex + Year*Julian + Sex*Julian, data = kill)
m21 <- lm(Mass ~ Year*Julian + Year*Sex + Sex*Julian, data = kill)

### AIC model selection
models <- list(m.null, m.globalyearjulian, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16,
               m17, m18, m19, m20, m21)

mod.names <- c('m.null', 'm.globalyearjulian', 'm1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13', 'm14', 'm15', 'm16',
               'm17', 'm18', 'm19', 'm20', 'm21')

aictab(models, modnames = mod.names)

# ---------------------------------------------------------------------------- #

## Top Model Summaries: First Stage ####

options(digits = 3)

# Mass ~ Year + Sex + Julian (m9)
cbind(summary(m9)$coefficients, confint(m9))

# Mass ~ Year * Sex + Julian (m14)
model_summary <- summary(m14)$coefficients
conf_intervals <- confint(m14, na.rm = TRUE)
common_terms <- intersect(rownames(model_summary), rownames(conf_intervals))
model_summary <- model_summary[common_terms, , drop = FALSE]
conf_intervals <- conf_intervals[common_terms, , drop = FALSE]
cbind(model_summary, conf_intervals)

# ---------------------------------------------------------------------------- #

# Second Stage AICc Model Selection: Mass ~ Detection ####

### null and global models
m.null <- lm(Mass ~ 1, data= kill)
m.informednull <- lm(Mass ~ Year + Sex + Julian, data = kill)
m.globalinformednull <- lm(Mass ~ Year + Sex + Julian + Detection + 
                             Year * Sex + Year * Julian + 
                             Sex * Julian + Julian * Detection, 
                           data = kill)
m.global <- lm(Mass ~ Year*Sex + Year*Julian + Sex*Julian, data = kill)

### additive combinations
m1 <- lm(Mass ~ Year + Detection, data = kill)
m2 <- lm(Mass ~ Season + Detection, data = kill)
m3 <- lm(Mass ~ Sex + Detection, data = kill)
m4 <- lm(Mass ~ Julian + Detection, data = kill)
m5 <- lm(Mass ~ Year + Sex + Detection, data = kill)
m6 <- lm(Mass ~ Year + Julian + Detection, data = kill)
m7 <- lm(Mass ~ Season + Sex + Julian, data = kill)
m8 <- lm(Mass ~ Year + Sex + Julian + Detection, data = kill)

### two-way interactions 
m9 <- lm(Mass ~ Detection * Julian, data = kill)

### two-way interactions with one additive combination 
m10 <- lm(Mass ~ (Detection * Julian) + Year, data = kill)
m11 <- lm(Mass ~ (Detection * Julian) + Season, data = kill)
m12 <- lm(Mass ~ (Detection * Julian) + Sex, data = kill)

### two-way interactions with multiple interactions 
#two-way
m13 <- lm(Mass ~ Detection + Year*Sex + Year*Julian, data = kill)
m14 <- lm(Mass ~ Detection + Year*Sex + Sex*Julian, data = kill)
m15 <- lm(Mass ~ Detection + Year*Julian + Year*Sex, data = kill)
m16 <- lm(Mass ~ Detection + Year*Julian + Sex*Julian, data = kill)

m17 <- lm(Mass ~ Detection*Julian + Year*Sex + Year*Julian, data = kill)
m18 <- lm(Mass ~ Detection*Julian + Year*Sex + Sex*Julian, data = kill)
m19 <- lm(Mass ~ Detection*Julian + Year*Julian + Year*Sex, data = kill)
m20 <- lm(Mass ~ Detection*Julian + Year*Julian + Sex*Julian, data = kill)

#three-way
m21 <- lm(Mass ~ Detection + Year*Sex + Year*Julian + Sex*Julian, data = kill)
m22 <- lm(Mass ~ Detection + Year*Julian + Year*Sex + Sex*Julian, data = kill)

m23 <- lm(Mass ~ Detection*Julian + Year*Sex + Year*Julian + Sex*Julian, data = kill)
m24 <- lm(Mass ~ Detection*Julian + Year*Julian + Year*Sex + Sex*Julian, data = kill)

### AIC model selection
models <- list(m.null, m.informednull, m.global, m.globalinformednull, m1, m2, 
               m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16,
               m17, m18, m19, m20, m21, m22, m23, m24)

mod.names <- c('m.null', 'm.informednull', 'm.global', 'm.globalinformednull', 
               'm1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 
               'm11', 'm12', 'm13', 'm14', 'm15', 'm16',
               'm17', 'm18', 'm19', 'm20', 'm21', 'm22', 'm23', 'm24')

aictab(models, modnames = mod.names)
