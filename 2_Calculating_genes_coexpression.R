library(tidyverse)
library(WGCNA)

###Calculating genes co-expression

#1- Computing optimal soft-power for a scale free topology model fit
sft <- pickSoftThreshold(nano.counts, RsquaredCut = 0.8, powerVector = c(1:30), networkType = "unsigned", verbose = 5)
sft_power <- sft$powerEstimate

#2- Calculating gene co-expression
co_expression.df <- as.data.frame(as.table(adjacency(t(nano.counts), power = sft_power)))

