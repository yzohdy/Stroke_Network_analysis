library(tidyverse)
library(brainGraph)
library(igraph)
library(tmod)

### Correlation between Rich-club genes expression and outcomes

#1- Loading outcomes dataframe
outcomes <- read.csv("Outcomes.csv")

#2- Extracting Rich-club and non-Rich-club genes expression
RC.genes <- as.data.frame(V(RC.g))
data.RC <- filter(data, rownames(nano.counts) %in% paste(rownames(RC.genes)))
data.non.RC.df <- filter(data, !rownames(data) %in% paste(rownames(RC.genes)))

#3- Calculating Eigengene for Rich-club versus randomly selected genes from the DEGs network

#3.1- Calculating Rich-club genes eigengenes
RC.groups <- data.frame(ID=c("a","b"), Title=c("RC", "All"))
RC.genes.list <- list(a =c(paste(rownames(RC.genes))), b = c (paste(unique(DEGs))))
RC.tmod <- makeTmodGS(gs2gene= RC.genes.list, gs= RC.groups)
RC.tmod.set <- RC.tmod[grep("RC", RC.tmod$gs$Title)]
RC.eigv1 <- eigengene(data, g= rownames(data), mset= RC.tmod.set, k=1)
RC.eigv2 <- eigengene(data, g= rownames(data), mset= RC.tmod.set, k=2)
RC.eigv3 <- eigengene(data, g= rownames(data), mset= RC.tmod.set, k=3)
RC.egenes <- data.frame(a = t(RC.eigv1), b= t(RC.eigv2), c= t(RC.eigv3))
RC.egenes.outcomes <- merge(RC.egenes, outcomes, by= 'row.names')


#3.2- Calculating non-Rich-club genes eigengenes
genes.not.RC <- unique(DEGs[!(DEGs %in% rownames(RC.genes))])

rand.EG <- list(outcome)
for (i in 1:1000) {
  rand.groups <- data.frame(ID = c("a", "b"), Title = c("Rand", "All"))
  rand.genes <- sample(genes.not.RC, size = nrow(RC.genes))
  rand.gene.list <- list(a = c(paste(rand.genes)), b = c(paste(unique(DEGs))))
  rand.tmod <- makeTmodGS(gs2gene = rand.gene.list, gs = rand.groups)
  rand.tmod.set <- rand.tmod[grep("Rand", rand.tmod$gs$Title)]
  rand.eigv1 <- eigengene(data, g = rownames(data), mset = rand.tmod.set, k = 1)
  rand.eigv2 <- eigengene(data, g = rownames(data), mset = rand.tmod.set, k = 2)
  rand.eigv3 <- eigengene(data, g = rownames(data), mset = rand.tmod.set, k = 3)
  rand.egenes <- data.frame(a = t(rand.eigv1), b = t(rand.eigv2), c = t(rand.eigv3))
  rand.EG[[i]] <- merge(rand.egenes, outcomes, by= 'row.names')}

#4- Compare the correlation between Rich-club genes expression versus non-Rich-club genes with outcomes
tested.outcome <- outcomes[]

#4.1- Correlating Rich-club genes' eigengenes with outcomes
RC.eigenR2 <- summary(lm(tested.outcome ~ a + a.1, data = RC.egenes.outcomes)) 

#4.2- Correlating random non-Rich-club gene sets eigengenes with outcomes
rand.EG.lm <- list()
rand.EG.lm.R2 <- list()
for(i in 1:length(rand.EG)){
  rand.EG.lm[[i]] <- summary(lm(tested.outcome ~ a + a.1, data = rand.EG[[i]])) 
  rand.EG.lm.R2[[i]] <- rand.EG.lm[[i]]$r.squared} 
rand.EG.lm.R2.df <- data.frame(R2 = t(as.data.frame(rand.EG.lm.R2)), row.names = 1:length(rand.EG))

#4.3 - T.test
t.test(rand.EG.lm.R2.df[,1], mu= RC.eigenR2$r.squared, alternative = "two.sided", paired = F, var.equal = F, conf.level = 0.95)


