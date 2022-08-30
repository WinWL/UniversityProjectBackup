library(ggplot2)
library(GGally)
library(plyr)
library(class)
library(rpart) # for the tree
library(rattle) # for the nice picture of the tree
library(gridExtra)
library(e1071)
library(mixture)
library(xtable)

## Importing Data
# rm(list=ls())
setwd("D:/Dropbox/Dropbox/School/SS_4M03_Multivariate_Analysis")
glass_data <- read.csv("glass_dataset/glass_data.csv") 

# Renaming columns
glass_data_column_names <- c("Id", "Ri", "Na", "Mg",
                             "Al", "Si", "K", "Ca",
                             "Ba", "Fe", "Type")
names(glass_data) <- glass_data_column_names
glass_data$Type <- as.factor(glass_data$Type)


# Removing Id column
gd <- glass_data[,-1]
summary(gd)
str(gd)

# Number of each type
length(which(gd$Type == 1)) 
length(which(gd$Type == 2)) 
length(which(gd$Type == 3)) 
length(which(gd$Type == 4)) 
length(which(gd$Type == 5)) 
length(which(gd$Type == 6)) 
length(which(gd$Type == 7)) 


# Only 3 types
gd_3_only <- gd
gd_3_only$Type <- mapvalues(glass_data$Type,
                            from = c("1", "2", "3", "5", "6", "7"),
                            to = c("1", "1", "1", "2", "2", "3"))
# Visualizations
head(gd)
print(xtable(head(gd), type = "latex"), file = "head.txt")

ggpairs(gd, aes(colour = Type, alpha = 0.4))
# ggpairs(gd_3_only, aes(colour = Type,  alpha = 0.4))
# ggpairs(gd[, -c(1,2,5,6,7)], aes(colour = Type, alpha = 0.4))
# ggpairs(gd_3_only[, -c(1,2,5,6,7)], aes(colour = Type,  alpha = 0.4))

# Box plots
a <- ggplot(gd, aes("Boxplot for all", Mg)) +
  xlab("")  + geom_boxplot() +
  scale_x_discrete(breaks=NULL) 
a <- ggplot(gd, aes(y = Ri, x = Type, fill = Type)) + geom_boxplot() +  xlab("Type") + ggtitle("Refractive Index") + theme(legend.position = "none")
b <- ggplot(gd, aes(y = Na, x = Type, fill = Type)) + geom_boxplot() +  xlab("Type") + ggtitle("Sodium") + theme(legend.position = "none")
c <- ggplot(gd, aes(y = Mg, x = Type, fill = Type)) + geom_boxplot() +  xlab("Type") + ggtitle("Magnesium") + theme(legend.position = "none")
d <- ggplot(gd, aes(y = Al, x = Type, fill = Type)) + geom_boxplot() +  xlab("Type") + ggtitle("Aluminum") + theme(legend.position = "none")
e <- ggplot(gd, aes(y = Si, x = Type, fill = Type)) + geom_boxplot() +  xlab("Type") + ggtitle("Silicon") + theme(legend.position = "none")
f <- ggplot(gd, aes(y = K, x = Type, fill = Type)) + geom_boxplot() +  xlab("Type") + ggtitle("Potassium") + theme(legend.position = "none")
g <- ggplot(gd, aes(y = Ca, x = Type, fill = Type)) + geom_boxplot() +  xlab("Type") + ggtitle("Calcium") + theme(legend.position = "none")
h <- ggplot(gd, aes(y = Ba, x = Type, fill = Type)) + geom_boxplot() +  xlab("Type") + ggtitle("Barium") + theme(legend.position = "none")
i <- ggplot(gd, aes(y = Fe, x = Type, fill = Type)) + geom_boxplot() +  xlab("Type") + ggtitle("Iron") + theme(legend.position = "none")
grid.arrange(a, b, c,d,e,f, g,h,i, nrow = 3, widths = c(3,3,3))

ggparcoord(gd, columns=1:9, alphaLines=0.5,groupColumn="Type",  scale="center") + xlab("") + ylab("")
ggparcoord(gd_3_only, columns=1:9, alphaLines=0.5,groupColumn="Type",  scale="center") + xlab("") + ylab("")

# Samples for straified sampling 
#   1 -  69 Type 1 Glass
#  70 - 145 Type 2 Glass
# 146 - 162 Type 3 Glass
# 163 - 175 Type 5 Glass
# 176 - 184 Type 6 Glass
# 185 - 213 Type 7 Glass 

strat <- function(){
  strat_T1 <- sample(c(  1: 69), floor(0.75*( 69 - 0  )))
  strat_T2 <- sample(c( 70:145), floor(0.75*(145 -  69)))
  strat_T3 <- sample(c(146:162), floor(0.75*(162 - 145)))
  strat_T5 <- sample(c(163:175), floor(0.75*(175 - 162)))
  strat_T6 <- sample(c(176:184), floor(0.5*(184 - 175)))
  strat_T7 <- sample(c(185:213), floor(0.75*(213 - 184)))
  train_row <- c(strat_T1, strat_T2, strat_T3, strat_T5, strat_T6, strat_T7)
  return(train_row)
}

train_row <- strat()

# This function calculates the KNN of a given test, train and label vectors
# It checks for 1 to max_n neighbours and prints out the associated misclassification
# rates and the ARI.
# Also prints out the confusion table for the knn that gives max ARI. 
KNN_Analysis <- function(train, test, true_train, true_test, max_n = 10) {
  
  # KNN with 1 neighbour. Get the confusion table and ARI at the same time
  knn_gd <- knn(train, test, true_train, k = 1)
  table_knn <- table(true_test, knn_gd)
  ARI_knn <- classAgreement(table_knn)
  
  # Construct a data frame and list to store values
  knn_df <- data.frame(1, ARI_knn[[1]], ARI_knn[[4]])
  names(knn_df) <- c("Neighbours","Misclassification Rate", "ARI")
  knn_all <- list(list())
  knn_all <- c(knn_all, list(list(knn_gd, table_knn, ARI_knn)))
  
  for (n in c(2:max_n)) {
    # Repeat KNN with more than 1 neighbour
    knn_gd <- knn(train, test, true_train, k = n)
    table_knn <- table(true_test, knn_gd)
    ARI_knn <- classAgreement(table_knn)
    
    # Store results
    df <- data.frame(n, ARI_knn[[1]], ARI_knn[[4]])
    names(df) <- c("Neighbours","Misclassification Rate", "ARI")
    knn_df <- rbind(knn_df, df)
    knn_all <- c(knn_all, list(list(knn_gd, table_knn, ARI_knn)))
  }
  
  # Prints out the relevant information
  # Finds the ARI of highest value and prints out associated ARI and confusion table 
  max_ind <- which(knn_df$ARI == max(knn_df$ARI))
  if (length(max_ind) > 1) { max_ind <- max_ind[1]}
  print(knn_df)
  print(knn_all[[max_ind+1]])
  paste("Max ARI is with neighbours =", max_ind)
}
citation(package = "e1071")
citation(package = "class")


print(xtable(knn_df, type = "latex"), include.rownames=FALSE, file = "knn_base.txt")
print(xtable(knn_all[[max_ind+1]][[2]], type = "latex"), file = "knn_conf.txt")


### Full Type Variable
# Making training and test sets
# And separating out the labels for both
train <- gd[train_row,]; test <- gd[-train_row,] 
true_train <- train$Type; true_test <- test$Type
train <- train[,-10]; test <- test[,-10] 
KNN_Analysis(train = train, test = test, true_train =  true_train, true_test = true_test)

# Removing some variables
x_data <- gd
x_data <- x_data[, -c(1,2,5,6,7)]
train <- x_data[train_row,]; test <- x_data[-train_row,] 
true_train <- train$Type; true_test <- test$Type
train <- train[,-10]; test <- test[,-10] 
KNN_Analysis(train = train, test = test, true_train =  true_train, true_test = true_test)

### Collapsed Type Variable
train <- gd_3_only[train_row,]; test <- gd_3_only[-train_row,] 
true_train <- train$Type; true_test <- test$Type
train <- train[,-10]; test <- test[,-10] 
KNN_Analysis(train = train, test = test, true_train =  true_train, true_test = true_test)

# Removing some variables
x_data <- gd_3_only
x_data <- x_data[, -c(1,2,5,6,7)]
train <- x_data[train_row,]; test <- x_data[-train_row,] 
true_train <- train$Type; true_test <- test$Type
train <- train[,-10]; test <- test[,-10] 
KNN_Analysis(train = train, test = test, true_train =  true_train, true_test = true_test)

# Constructs a classification tree out of the given data
# Assumes the given data is the glass data frame
# Assumes there is a train_row 
# Intended usage is for x to be the glass data set. 
# Intention is for x to change by changing the number of factors in the type variable 
classificationTree <- function(x) {
  
  # Construct trees
  gd_tree <- rpart(Type ~ ., data=x, method="class", parms = list(split = 'gini'))
  gd_tree_2 <- rpart(Type ~., data = x, subset = train_row, parms = list(split = 'gini'))
  
  # Predict on test set
  table_tree <- table(predict(gd_tree_2, test, type = "class"), x[-train_row, "Type"])
  
  # Draw figures
  print("Full Tree")
  fancyRpartPlot(gd_tree, main="Classification Tree for the Glass Data, Full Data")
  printcp(gd_tree)
  
  print("Test/Train Tree")
  fancyRpartPlot(gd_tree_2, main="Classification Tree for the Glass Data, Test Data")
  printcp(gd_tree_2)
  
  # Classification Rate 
  print("Confusion Table + ARI")
  print(table_tree)
  classAgreement(table_tree)
}
citation(package = "rpart")

?rpart
x <- gd
print(xtable(gd_tree_2$cptable, type = "latex"), include.rownames=FALSE, file = "cla_cp.txt")
print(xtable(table_tree, type = "latex"), file = "cla_mat.txt")

train_row <- train_row <- strat()

train <- gd[train_row,]; test <- gd[-train_row,] 
true_train <- train$Type; true_test <- test$Type
train <- train[,-10]; test <- test[,-10] 

# Testing the classification tree with all 6 present glass types and then after collapsing to 3 types
classificationTree(gd)

train <- gd_3_only[train_row,]; test <- gd_3_only[-train_row,] 
true_train <- train$Type; true_test <- test$Type
train <- train[,-10]; test <- test[,-10] 
classificationTree(gd_3_only)

test <- test[, -c(1,2,5,6,7)]
classificationTree(gd[, -c(1,2,5,6,7)])
classificationTree(gd_3_only[, -c(1,2,5,6,7)])