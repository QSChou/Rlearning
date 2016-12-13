# use the factoextra packages for pca plot
library(factoextra)
iris_pca <- princomp(iris[,-5])
fviz_pca(iris_pca)
fviz_pca_ind(iris_pca)
fviz_pca_contrib(iris_pca)
fviz_pca_var(iris_pca)
fviz_pca_biplot(iris_pca)
