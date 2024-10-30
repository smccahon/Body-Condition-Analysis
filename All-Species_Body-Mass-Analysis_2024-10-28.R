#----------------------------------#
#  All Species Body Mass Analysis  #
#        Created 10/28/2024        #          
#       Modified 10/29/2024        #
#----------------------------------#

# load packages
library(dplyr)
library(ggplot2)
library(AICcmodavg)
library(RColorBrewer)
library(trtools)
library(tidyverse)
library(lme4)
library(car)

# ---------------------------------------------------------------------------- #

# Data Processing & Manipulation ####

# Read data
setwd("processed_data")
birds <- read.csv("LM_ShorebirdsALLNeg.csv")

# Make neonicotinoid detection column (Detection/Non-detection)
birds$Detection <- ifelse(birds$OverallNeonic > 0, "Detection", "Non-detection")

# Make sampling event column
birds <- birds %>% 
  mutate(Event = paste(Season, Year, sep = "_"))

birds <- birds %>% 
  mutate(Event = ifelse(Season %in% c("Spring", "Fall"),
                        paste(Season, Year, sep = "_"),
                        NA))

# Logarithmic transformation of neonics
birds <- birds %>% 
  mutate(LogNeonic = log10(OverallNeonic + 0.0001))

# Reduce species to those with 3 observations
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# Filter out birds that do no have measured neonic concentrations
birds <- subset(birds, !is.na(OverallNeonic))

# Convert categorical variables to factor
birds$Sex <- as.factor(birds$Sex)
birds$Detection <- as.factor(birds$Detection)
birds$Event <- as.factor(birds$Event)
birds$Species <- as.factor(birds$Species)

# Count number of individuals per species with 
# 1. Measured neonic concentrations and 2) at least 3 individuals
table(birds$Species)

# ---------------------------------------------------------------------------- #

# All Species Modeling ####

options(digits = 3)

## Reduced Null Model ####
m.reduced <- lmer(Mass ~ Sex + Event + (1 | Species) + Julian, data = birds)
summary(m.reduced)
confint(m.reduced)

# Sampling event is important: Birds in fall 2023 had significantly lower body mass (likely due to age)

#-------#

## Neonicotinoid Detection Model ####
m.detection <- lmer(Mass ~ Sex + Event + (1 | Species) + Julian + Detection, data = birds)
summary(m.detection)
 confint(m.detection)

# Conclusion: The effect of neonicotinoid detection on body mass is not significant
# Sampling event is the only significant predictor (birds in fall 2023 had significantly lower body mass)

#-------#

### Compare full and reduced models: Detection ####

anova(m.reduced, m.detection)

# Neonicotinoid detection does not further explain any variation in body mass

#-------#

## Neonicotinoid Concentration Model ####
m.conc <- lmer(Mass ~ Sex + Event + (1 | Species) + Julian + OverallNeonic, data = birds)
summary(m.conc)
confint(m.conc)

# Conclusion: The effect of neonicotinoid concentration on body mass is not significant
# Sampling event is the only significant predictor (birds in fall 2023 had significantly lower body mass)

#-------#

### Compare full and reduced models: Concentration ####

anova(m.reduced, m.conc)

# Neonicotinoid concentration does not further explain any variation in body mass

#-------#

## Neonicotinoid Log(Concentration) Model ####
m.log <- lmer(Mass ~ Sex + Event + (1 | Species) + Julian + LogNeonic, data = birds)
summary(m.log)
confint(m.log)

# Conclusion: The effect of log(concentration) on body mass is not significant
# Sampling event is the only significant predictor (birds in fall 2023 had significantly lower body mass)

#-------#

### Compare full and reduced models: Log(Concentration) ####

anova(m.reduced, m.log)

# Log(concentration) does not further explain any variation in body mass

# ---------------------------------------------------------------------------- #

## Agricultural Intensity Model ####

m.fullnullag <- lmer(Mass ~ Sex + Event + Julian + (1 | Species), data = birds)
m.bestnullag <- lmer(Mass ~ Sex + Event + (1 | Species), data = birds)

anova(m.fullnullag, m.bestnullag) # excluding Julian improves model fit

#------#

### Compare full and reduced models: Agricultural Intenstiy ####

m.fullag <- lmer(Mass ~ Sex + Event + (1 | Species) + PercentAg, data = birds)

summary(m.fullag)
confint(m.fullag)

# The effect of surrounding agricultural intensity on body mass is not significant

anova(m.fullag, m.bestnullag)

# Including percent agriculture around capture site does not improve model fit. 

# ---------------------------------------------------------------------------- #

## Identifying best reduced model with species as random effect ####

m.sex <- lmer(Mass ~ Sex + (1 | Species), data = birds, REML = FALSE)
m.julian <- lmer(Mass ~ Julian + (1 | Species), data = birds, REML = FALSE)
m.event <- lmer(Mass ~ Event + (1 | Species), data = birds, REML = FALSE)
m.sexjulian <- lmer(Mass ~ Sex + Julian + (1 | Species), data = birds, REML = FALSE)
m.sexevent <- lmer(Mass ~ Sex + Event + (1 | Species), data = birds, REML = FALSE)
m.julianevent <- lmer(Mass ~ Julian + Event + (1 | Species), data = birds, REML = FALSE)
m.sexjulianevent <- lmer(Mass ~ Julian + Sex + Event + (1 | Species), data = birds, REML = FALSE)

### AICc Model Selection: All Species ####
models <- list(m.sex, m.julian, m.event, m.sexjulian, m.sexevent, m.julianevent,
               m.sexjulianevent)

mod.names <- c('m.sex', 'm.julian', 'm.event', 'm.sexjulian', 'm.sexevent', 'm.julianevent',
               'm.sexjulianevent')

aictab(models, modnames = mod.names)

### Top model summary: Mass ~ Sex + Event + (1 | Species) ####

summary(m.sexevent)
confint(m.sexevent)

#-------#

## Identifying best neonicotinoid model (with species as a random effect) ####

m.detection <- lmer(Mass ~ Detection + (1| Species), data = birds, REML = FALSE)
m.conc <- lmer(Mass ~ OverallNeonic + (1 | Species), data = birds, REML = FALSE)
m.log <- lmer(Mass ~ LogNeonic + (1 | Species), data = birds, REML = FALSE)

### AICc Model Selection: All Species ####
models <- list(m.detection, m.conc, m.log)

mod.names <- c('m.detection', 'm.conc', 'm.log')

aictab(models, modnames = mod.names)

# Conclusion: All models perform similarly (within 2 Delta AICc)

# ---------------------------------------------------------------------------- #


















# ---------------------------------------------------------------------------- #

## Assess interactions: Neonic * Julian ####
m.detectionjulian <- lmer(Mass ~ Detection*Julian + Sex + Event + (1 | Species), data = birds)
summary(m.detectionjulian)
confint(m.detectionjulian)

# Interaction between detection and julian is significant. The positive effect of 
# detection on mass becomes less pronounced as Julian day progresses
# Birds in fall 2023 had significantly lower body mass

m.concjulian <- lmer(Mass ~ OverallNeonic*Julian + Sex + Event + (1 | Species), data = birds)


m.interact <- lm(Mass ~ LogNeonic*MigDate + Sex + Event + MigDate, data = leye)
cbind(summary(m.interact)$coefficients, confint(m.interact))

# Neonicotinoids do not vary with migration date

#-------#

## Assess interactions: Neonic * Sex: LEYE ####
m.interact <- lm(Mass ~ Detection*Sex + Sex + Event + MigDate, data = leye)
cbind(summary(m.interact)$coefficients, confint(m.interact))

m.interact <- lm(Mass ~ OverallNeonic*Sex + Sex + Event + MigDate, data = leye)
cbind(summary(m.interact)$coefficients, confint(m.interact))

m.interact <- lm(Mass ~ LogNeonic*Sex + Sex + Event + MigDate, data = leye)
cbind(summary(m.interact)$coefficients, confint(m.interact))

# Neonicotinoids do not vary with sex

#-------#

## Assess interactions: Neonic * Event: LEYE ####
m.interact <- lm(Mass ~ Detection*Event + Sex + Event + MigDate, data = leye)
summary(m.interact)$coefficients
confint(m.interact)

m.interact <- lm(Mass ~ OverallNeonic*Event + Sex + Event + MigDate, data = leye)
summary(m.interact)$coefficients
confint(m.interact)

m.interact <- lm(Mass ~ LogNeonic*Event + Sex + Event + MigDate, data = leye)
summary(m.interact)$coefficients
confint(m.interact)

# Neonicotinoids do not vary with sampling event

#-------#

# Assess Model Assumptions ####

