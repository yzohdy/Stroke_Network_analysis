library(STRINGdb)
library(brainGraph)
library(igraph)
library(RCy3)

###Constructing and visualzing interaction network

#1- Identifing gene-gene interactions using STRING-database
string_db <- STRINGdb$new(version="12.0", species=9606,score_threshold=200, network_type="full", input_directory="")
DEGs.id <- string_db$map(DEGs, "geneID", removeUnmappedRows = TRUE)
interaction.network <- string_db$get_interactions(DEGs.id$STRING_id)

#2- Construct network
DEGs.g <-graph.data.frame(interaction.network[,c(1,2)])
E(DEGs.g)$weight <- co_expression.df
V(DEGs.g)$degree <- degree(DEGs.g)
V(DEGs.g)$node.betweenness <- betweenness(DEGs.g)
V(DEGs.g)$node.corness <- coreness(DEGs.g)
V(DEGs.g)$strength <- strength(DEGs.g)
V(DEGs.g)$transitivity <- transitivity(DEGs.g)

#3- Visualize on Cytoscape
createNetworkFromIgraph(DEGs.g, title = "DEGs Interaction Network")
setEdgeLineWidthMapping("weight", widths=c(0,1))
setEdgeColorDefault("grey")
setNodeColorMapping("degree", colors = paletteColorBrewerGreens)
setNodeSizeMapping("strength")


