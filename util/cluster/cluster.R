# script reference http://www.sthda.com/english/wiki/cluster-analysis-in-r-unsupervised-machine-learning#basic-clustering-methods
# clustering after pca
# Install factoextra
##install.packages("factoextra")
# Install cluster package
##install.packages("cluster")
library(factoextra)
library(cluster)

res.dist <- get_dist(USArrests, stand = TRUE, method = "pearson")
fviz_dist(res.dist,gradient = list(low = "#00AFBB",mid = "white" ,high = "#FC4E07"))
