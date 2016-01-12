# Below example use the iris data and knn to do classification
# Author: QS Chou Version 0.1 Date: 2016/01/12
# Although normaliztion is not required in the iris data set, below code still includes the normalization step for learning purpose
attach(iris)
normalize <- function(x) {
  num <- x-min(x)
  denum <- max(x) - min(x)
  return (num/denum) 
}
# test normalize; case when divided 0, it will return Null automatically.
n1 <- normalize(c(100,100,100))
# test normalize; usual case
n2 <- normalize(c(200,100,300,29,39,19,10))
# test normalize; negative value case
n3 <- normalize(c(-200,100,300,-29,39,-19,10))
iris_norm <- as.data.frame(lapply(iris[1:4], normalize))
iris_n <- cbind(iris_norm,iris[5])

# sampling data set into training & test
set.seed(1234)
ind <- sample(2, nrow(iris),replace = TRUE, prob = c(0.7,0.3))

# Original data set
iris_training <- iris[ind==1,1:4]
iris_training_labels <- iris[ind==1,5]
iris_test <- iris[ind==2, 1:4]
iris_test_lables <- iris[ind==2,5]

# Normalized data set
iris_n_training <- iris_n[ind==1,1:4]
iris_n_training_labels <- iris_n[ind==1,5]
iris_n_test <- iris_n[ind==2, 1:4]
iris_n_test_lables <- iris_n[ind==2,5]

# Apply knn for learning; it will need to install the 'class' packages

## Try with k= 3 & use the original data set
library(class)
iris_pred <- knn(train = iris_training, test = iris_test, cl = iris_training_labels, k=3)
iris_test_pred <- cbind(iris_test,iris_test_lables,iris_pred)
summary(iris_test_pred)
table(iris_test_pred$iris_test_lables,iris_test_pred$iris_pred)

## Try with k= 5 & use the original data set
library(class)
iris_pred <- knn(train = iris_training, test = iris_test, cl = iris_training_labels, k=5)
iris_test_pred <- cbind(iris_test,iris_test_lables,iris_pred)
summary(iris_test_pred)
table(iris_test_pred$iris_test_lables,iris_test_pred$iris_pred)


## Try with k= 3 & use the normalized data set
library(class)
iris_pred <- knn(train = iris_n_training, test = iris_n_test, cl = iris_n_training_labels, k=3)
iris_test_pred <- cbind(iris_n_test,iris_n_test_lables,iris_pred)
summary(iris_test_pred)
table(iris_test_pred$iris_n_test_lables,iris_test_pred$iris_pred)

## Try with k= 5 & use the normalized data set
library(class)
iris_pred <- knn(train = iris_n_training, test = iris_n_test, cl = iris_n_training_labels, k=11)
iris_test_pred <- cbind(iris_test,iris_n_test_lables,iris_pred)
summary(iris_test_pred)
table(iris_test_pred$iris_n_test_lables,iris_test_pred$iris_pred)
