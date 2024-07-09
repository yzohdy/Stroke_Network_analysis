library(tidyverse)
library(org.Mm.eg.db)
library(RCy3)
library(brainGraph)
library(igraph)
library(clusterProfiler)

###Graph theory analysis

#1- Plotting network degrees distribution
g.RC <- rich_club_all(DEGs.g)
g.RC <- g.RC[!is.nan(phi), ]
ggplot(g.RC, aes(x = k, y= Nk)) +
  theme_bw() +
  geom_smooth(color= "black", alpha=0) +
  geom_point(color= "black", alpha=0.4) +
  labs(y = "Nodes with degree > k", x= "Degree(k)", title = "Degree(k) Distribution" ) +
  guides(fill = T)

#2- Comparing DEGs network with random networks of similiar degree distribution
#2.1- Construction of random networks
rand.n <- 1000
RR <- list()
for(i in 1:rand.n) {
  RR[[i]] <- sample_degseq(degree(DEGs.g), method = c("vl"))
  E(RR[[i]])$weight <- sample(E(DEGs.g)$weight)}

#2.2- Calculate cluster-coefficent for random graph and compare with DEGs network
clu.phi <- vector()
for (i in 1:rand.n) {clu.phi[i] <- transitivity(RR[[i]])}
clu.phi.t.test <- t.test(clu.phi, mu= transitivity(DEGs.g))

#3- Rich Club analysis
#3.1- Rich-club analysis for random networks
r.RC <- list()
for (i in 1:rand.n) {
  r.RC[[i]] <- rich_club_all(RR[[i]], weighted = T)}

#3.2- Average Rich-club coefficient for random networks 
rand.phi <- matrix(0, nrow = nrow(r.RC[[1]]) , ncol = 1)
for (row in 1:nrow(r.RC[[1]])) {
  row_values <- sapply(r.RC, function(x) x[row,]$phi)
  rand.phi[row, ] <- mean(row_values)}

#3.3- Rich-club analysis for DEGs network
DEGs.g.RC <- rich_club_all(DEGs.g, weighted = T)

#3.4- Normalized Rich-club coefficient
norm.RC.phi <- (DEGs.g.RC$phi/rand.phi)
RC.df <- data.frame(degree = DEGs.g.RC$k, Network.phi = DEGs.g.RC$phi, Rand.phi = rand.phi, Norm.phi = norm.RC.phi)
RC.df <- RC.df[!is.nan(RC.df$Network.phi), ]

#3.5- Idenfitying Hubs nodes based on Rich-club coefficient
DEGs.g.RC.coef <- rich_club_coeff(g, k= 80, weighted = F, A= as_adjacency_matrix(DEGs.g))

#3.6- Constructing Rich-club genes network
RC.g <- DEGs.g.RC.coef$graph
V(RC.g)$degree <- degree(RC.g)
V(RC.g)$node.betweenness <- betweenness(RC.g)
V(RC.g)$strength <- strength(RC.g)

#3.7- Visualizing Rich-club genes network
createNetworkFromIgraph(g.RC, title = "DEGs networks Rich-club genes")

#4- Comparing the degrees between Rich-club and non-Rich-club genes
degree.g.df <- as.data.frame(degree(g))
degree.RC <- filter(degree.g.df, rownames(degree.g.df) %in% paste(V(RC.g)))
degree.not.RC <- filter(degree.g.df, !rownames(degree.g.df) %in% paste(V(RC.g)))

#5- Functional annotation of Rich-club genes
GO.results.RC <- enrichGO(gene = V(RC.g), OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = c("BP","CC"))
