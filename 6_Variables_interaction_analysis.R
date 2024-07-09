library(tidyverse)
library(jtools)
library(interactions)

###Investigating interactions between outcome predicting varialbes

#1- Loading dataframe
data <- read.csv("outcomes&eigengenes.csv")

#2- Interaction analysis
tested.outcome <- data[]
fiti <- lm(tested.outcome ~ Eigengene.PC1 * Eigengene.PC2 * Treatment, data = data)
summ(fiti)

#3- Ploting interaction
ip <- interact_plot(fiti, pred = Eigengene, modx = Treatment)

