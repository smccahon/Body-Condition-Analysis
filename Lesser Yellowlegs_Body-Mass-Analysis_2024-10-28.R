#----------------------#
#  Body Mass Analysis  #
#   Created 10/28/2024 #          
#  Modified 10/28/2024 #
#----------------------#

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



# Convert categorical variables to factor
birds$Sex <- as.factor(birds$Sex)
birds$Detection <- as.factor(birds$Detection)
birds$Event <- as.factor(birds$Event)

# Subset data for Lesser Yellowlegs
leye <- subset(birds, Species %in% c("LesserYellowlegs"))

# ---------------------------------------------------------------------------- #

# Lesser Yellowlegs Modelling ####

## Reduced Null Model: LEYE ####
options(digits = 3)
# Make neonicotinoid detection column (High, Low, Non-detection)
m.null <- lm(Mass ~ Sex + Event + MigDate, data = leye)
cbind(summary(m.null)$coefficients, confint(m.null))

# Lesser Yellowlegs had signficantly lower body mass in Fall 2023
# Age likely has an important effect on body mass

#-------#

## Neonicotinoid Detection Model: LEYE ####
m.detect <- lm(Mass ~ Sex + Event + MigDate + Detection, data = leye)
cbind(summary(m.detect)$coefficients, confint(m.detect)) 

# No effect of neonicotinoid detection on body mass (CI overlap 0)

#-------#

### Comparing full and reduced models: Detection (LEYE) ####
anova(m.null, m.detect)

# Full model does NOT significantly improve fit 

#-------#

## Neonicotinoid Concentration Model: LEYE ####
m.conc <- lm(Mass ~ Sex + Event + MigDate + OverallNeonic, data = leye)
cbind(summary(m.conc)$coefficients, confint(m.conc)) 

# No effect of neonicotinoid concentration on body mass (CI overlap 0)

#-------#

### Comparing full and reduced models: Concentration (LEYE) ####
anova(m.null, m.conc)

# Full model does NOT significantly improve fit 

#-------#

## Log(Neonicotinoid Concentration) Model: LEYE ####
m.log <- lm(Mass ~ Sex + Event + MigDate + LogNeonic, data = leye)
cbind(summary(m.log)$coefficients, confint(m.log)) 

# No effect of neonicotinoid concentration on body mass (CI overlap 0)

#-------#

### Comparing full and reduced models: Log(Concentration) (LEYE) ####
anova(m.null, m.log)

# Full model does NOT significantly improve fit 

#-------#

### Assess interactions: Neonic * MigDate: LEYE ####
m.interact <- lm(Mass ~ Detection*MigDate + Sex + Event + MigDate, data = leye)
cbind(summary(m.interact)$coefficients, confint(m.interact))

m.interact <- lm(Mass ~ OverallNeonic*MigDate + Sex + Event + MigDate, data = leye)
cbind(summary(m.interact)$coefficients, confint(m.interact))

m.interact <- lm(Mass ~ LogNeonic*MigDate + Sex + Event + MigDate, data = leye)
cbind(summary(m.interact)$coefficients, confint(m.interact))

# Neonicotinoids do not vary with migration date

#-------#

### Assess interactions: Neonic * Sex: LEYE ####
m.interact <- lm(Mass ~ Detection*Sex + Sex + Event + MigDate, data = leye)
cbind(summary(m.interact)$coefficients, confint(m.interact))

m.interact <- lm(Mass ~ OverallNeonic*Sex + Sex + Event + MigDate, data = leye)
cbind(summary(m.interact)$coefficients, confint(m.interact))

m.interact <- lm(Mass ~ LogNeonic*Sex + Sex + Event + MigDate, data = leye)
cbind(summary(m.interact)$coefficients, confint(m.interact))

# Neonicotinoids do not vary with sex

#-------#

### Assess interactions: Neonic * Event: LEYE ####
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

# Agricultural Intensity Model ####
m.ag <- lm(Mass ~ Sex + Event + MigDate + PercentAg, data = leye)
cbind(summary(m.ag)$coefficients, confint(m.ag))

m.ag <- lm(Mass ~ Sex + Event + MigDate + AgCategory, data = leye)
cbind(summary(m.ag)$coefficients, confint(m.ag))

#-------#

### Comparing full and reduced models: Log(Concentration) (LEYE) ####
anova(m.null, m.ag)

# ---------------------------------------------------------------------------- #

# Identifying best null model: LEYE ####
m.null <- lm(Mass ~ 1, data = leye)
m.sex <- lm(Mass ~ Sex, data = leye)
m.mig <- lm(Mass ~ MigDate, data = leye)
m.event <- lm(Mass ~ Event, data = leye) # BEST NULL MODEL
m.sexmig <- lm(Mass ~ Sex + MigDate, data = leye)
m.sexevent <- lm(Mass ~ Sex + Event, data = leye)
m.migevent <- lm(Mass ~ MigDate + Event, data = leye)
m.sexmigevent <- lm(Mass ~ MigDate + Sex + Event, data = leye)

# Model summaries
summary(m.event)
summary(m.log)

## AICc Model Selection: LEYE ####
models <- list(m.null, m.sex, m.mig, m.event, m.sexmig, m.sexevent, m.migevent,
               m.sexmigevent)

mod.names <- c('m.null', 'm.sex', 'm.mig', 'm.event', 'm.sexmig', 'm.sexevent', 'm.migevent',
               'm.sexmigevent')

aictab(models, modnames = mod.names)

# Identify best neonic model: LEYE ####
m.conc <- lm(Mass ~ OverallNeonic, data = leye)
m.log <- lm(Mass ~ LogNeonic, data = leye)
m.detect <- lm(Mass ~ Detection, data = leye)

## AICc Model Selection: LEYE ####
models <- list(m.null, m.conc, m.log, m.detect)

mod.names <- c('m.null', 'm.conc', 'm.log', 'm.detect')

aictab(models, modnames = mod.names) # DETECTION AND LOG TRANSFORMATION ARE BEST

## Likelihood Ratio Test: LEYE ####
m.logevent <- lm(Mass ~ LogNeonic + Event, data = leye)
m.detectevent <- lm(Mass ~ Detection + Event, data = leye)
anova(m.logevent, m.event)
anova(m.detectevent, m.event)

m.logevent <- lm(Mass ~ LogNeonic + Event, data = leye)
models <- list(m.log, m.event, m.logevent)
mod.names <- c('m.log', 'm.event', 'm.logevent')
aictab(models, modnames = mod.names)
anova(m.event, m.log)

# Conclusion: Neonicotinoids do not further inform body mass beyond sampling event.
# Evidence: AICc Model Selection, beta coefficients, and null-hypothesis significance testing

# ---------------------------------------------------------------------------- #

# Check Model Assumptions ####

# Diagnostic plot of null model
plot(m.null)

# Conclusion: Assumptions reasonably met

#--------#

# Diagnostic plot of concentration model
plot(m.conc)

# Conclusion: Some deviation from normality; one extreme outlier
# Solution: Transform neonicotinoid variable

#--------#

# Diagnostic plot of log(concentration) model
plot(m.log)

# Conclusion: Assumptions reasonably met

#--------#

# Diagnostic plot of detection model
plot(m.detect)

# Conclusion: Assumptions reasonably met

#--------#

# Diagnostic plot of ag intensification model (continuous variable)
plot(m.ag)

# Conclusion: Assumptions reasonably met

#--------#

# Assess multicollinearity ####
vif(m.null)
vif(m.detect)
vif(m.conc)
vif(m.log)
vif(m.ag)

# Conclusion: Low to moderate multicollinearity

# ---------------------------------------------------------------------------- #