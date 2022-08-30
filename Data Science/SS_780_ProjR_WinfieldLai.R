


# Winfield Lai
# 001327375
# Stats 780 Final Project


# rm(list=ls())
# setwd("D:/Dropbox/Dropbox/School/SS_780_DataScience/Project")
# setwd("E:/Dropbox/Dropbox/School/SS_780_DataScience/Project")

library(dplyr)
library(xtable)
library(ggplot2)
library(reshape)
library(ggpubr)
library(gridExtra)
library(arules)
library(e1071)
library(nnet)
library(mclust)
library(gbm)
library(randomForest)
library(tree)
library(GGally)
library(ggthemes)
library(reshape2)
library(foreach)
library(doParallel)
library(bibtex)
source("std_lift.R")
write.bib(c('arules', 
            'e1071', 
            'randomForest',
            'nnet',
            'tree'), 
          file='RPackages')


#########################################
# Importing in the data
df.raw <- read.csv("expanded.csv")
colnames(df.raw) <- c("Class", 
                   "CapShape",
                   "CapSurface",
                   "CapColor",
                   "Bruises", 
                   "Odor", 
                   "GillAttachment",
                   "GillSpacing",
                   "GillSize",
                   "GillColor",
                   "StalkShape",
                   "StalkRoot",
                   "StalkSurfaceAboveRing",
                   "StalkSurfaceBelowRing",
                   "StalkColorAboveRing",
                   "StalkColorBelowRing",
                   "VeilType",
                   "VeilColor",
                   "RingNumber",
                   "RingColor",
                   "SporePrintColor",
                   "Population",
                   "Habitat")
# Remove Missing data
df.rNA1 <- df.raw %>% filter_all(all_vars(. != "?")) %>% na.omit()
df.rNA <- df.rNA1 %>% 
  apply(MARGIN=2, function(x) as.factor(x)) %>%
  as.data.frame()

# names(df.rNA)
# # df %>% View()
# df.rNA %>% summary
# df.rNA %>% summary


#########################################
# Descriptive Plots


# Test code to makes stack bar plots
df.Plot <- df.rNA

# Making stacked bar plots via loop
df.stackPlots <- lapply(names(df.Plot)[-1], function(name){
  plot <- ggplot(df.Plot, aes(x = !!sym(name), fill = Class)) + 
    geom_bar() + 
    theme(legend.position = "none",
          text = element_text(size=6),
          axis.title = element_text(size=9),
          axis.text = element_text(size=9))
  return(plot)
})
# df.stackPlots[[3]]

df.Plot %>% filter(Odor == "FOUL" )
summary(df.Plot)
# Arranging all stacked bar plots in a grid
p <- df.stackPlots
grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]],
             p[[6]], p[[7]], p[[8]], p[[9]], p[[10]],
             p[[11]], p[[12]], p[[13]], p[[14]], p[[15]],
             p[[16]], p[[17]], p[[18]], p[[19]], p[[20]],
             p[[21]], p[[22]],
             nrow = 5, ncol = 5)

grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]],
             p[[6]], p[[7]], p[[8]], p[[9]], p[[10]],
             p[[11]],
             nrow = 3, ncol = 4,
             top = "Counts of Edible (Red) and Poisonous (Blue) Mushrooms")

grid.arrange(p[[12]], p[[13]], p[[14]], p[[15]],
             p[[16]], p[[17]], p[[18]], p[[19]], p[[20]],
             p[[21]], p[[22]],
             nrow = 3, ncol = 4,
             top = "Counts of Edible (Red) and Poisonous (Blue) Mushrooms")

grid.arrange(p[[5]], p[[6]], p[[2]],
             nrow = 1, ncol = 3,
             top = "Counts of Edible (Red) and Poisonous (Blue) Mushrooms")



############################################
# Removing Covariates
# Veil Type can be removed. All observations have the same factor
# df.rNA$VeilType %>% summary
df.rNa.rVeil <- subset(df.rNA, select = -c(VeilType))

# We also remove odor since it has such good separation
df.rNa.rVeil.rOdor <- subset(df.rNA, select = -c(VeilType, VeilColor, Odor))

df.rNa.rVeilT.rVeilC <- subset(df.rNA, select = -c(VeilType, VeilColor))

df.raw.final <- df.rNa.rVeilT.rVeilC %>% as.data.frame()
names(df.raw.final)


#########################################
# Splitting for training and test set
# Try classfication techniques. 
# Random Forest with Classification

# Subset the data base. For testing purposes - using the entire data set may be too slow. Use entire data set after all methods worked out.
df.subset.speed.index <- sample(
  1:length(df.raw.final$Class), 0.1*length(df.raw.final$Class))
df.raw.subset.speed <- df.raw.final[df.subset.speed.index,]
df.raw.forA <- df.raw.subset.speed %>% as.data.frame()
df.raw.forA.num <- data.matrix(df.raw.forA)

# df.raw.forA %>% typeof()

# Training set
df.train.index <- sample(
  1:length(df.raw.forA$Class), 0.75*length(df.raw.forA$Class))

#########################################
# Association Rules

predArules <- function(df.Final,
                       fSup = 0.055, 
                       fCon = 0.8,
                       fmin = 3,
                       fmax = 12) {
  arules.params <- list(support=fSup,confidence=fCon,minlen=fmin,maxlen=fmax)
  params <- arules.params # Here becasue std_lift requires an environmentlabel called params
  
  # Looking at only Poisonous
  arules.app.poi <- list(rhs=c("Class=POISONOUS"))
  arules.fit.poi <- apriori(df.Final, 
                            parameter = arules.params,
                            appearance = arules.app.poi)
  quality(arules.fit.poi) <- std_lift(arules.fit.poi, 
                                      df.Final)
  arules.fit.poi.rul <- inspect(sort(arules.fit.poi, by = "slift") %>% head(n=10))
  # arules.fit.poi.rul
  # Looking at only Edible
  arules.app.edi <- list(rhs=c("Class=EDIBLE"))
  arules.fit.edi <- apriori(df.Final, 
                            parameter = arules.params,
                            appearance = arules.app.edi)
  quality(arules.fit.edi) <- std_lift(arules.fit.edi, 
                                      df.Final)
  arules.fit.edi.rul <- inspect(sort(arules.fit.edi, by = "slift") %>% head(n=10))
  
  return(setNames(list(arules.fit.poi.rul, arules.fit.edi.rul),
                  c("aRulPoi", "aRulEdi")))
}
# out.aRul <- predArules(df.raw.forA)
# out.aRul$aRulPoi
# out.aRul$aRulEdi
# ?apriori
##########################################
# Functions
#########################################
# Random Forest Prediction

# Trains a Random forest with the given input parameters.
# Tests the trained random forest
# Returns the confusion table, ARI from the test set and the RandomForest Object
predForest <- function(formula = as.formula("Class ~ ."),
                       df.Final,
                       df.Index.Train,
                       fType = "class",
                       fImportance = T,
                       fMtry = 6,
                       fClass = "Class",
                       fnTree = 500
                     ) {
  
  datForest <- randomForest(
    formula,
    data = df.Final,
    subset = df.Index.Train,
    mtry = fMtry,
    importance = fImportance,
    type = fType,
    ntree = fnTree
  )
  # ?randomForest
  datForest.pred <- predict(datForest,
                            df.Final[-df.Index.Train, ],
                            type = fType)
  
  dat.Forest.tab <- table(df.Final[-df.Index.Train, fClass],
                          datForest.pred)
  dat.Forest.ari <- classAgreement(dat.Forest.tab)$crand
  
  return(setNames(list(datForest, dat.Forest.tab, dat.Forest.ari), c("ForestMod","ForestTab", "ForestAri")))
}
out.Forest <- predForest(df.Final =  df.raw.forA, df.Index.Train =  df.train.index, fMtry = 3, fnTree = 1000)
# importance(out.Forest$ForestMod)
# varImpPlot(out.Forest$ForestMod)
out.Forest$ForestAri
# out.Forest$ForestTab


#########################################
# Bagging Prediction

# Bagging is just a Random foest that considers all variables instead of a subsection upon splitting
predBag <- function(formula = as.formula("Class ~ ."),
                    df.Final,
                    df.Index.Train,
                    fType = "class",
                    fImportance = T) {
  
  out.ForestToBag <- predForest(df.Final = df.Final,
                       df.Index.Train = df.Index.Train,
                       fMtry = length(names(df.Final)) - 1) #-1 because df.Final includes the true label
  
  names(out.ForestToBag) <- c("BagMod","BagTab", "BagAri")
  return(out.ForestToBag)
}
# out.Bag <- predBag(df.Final = df.raw.forA, df.Index.Train = df.train.index)
# importance(out.Bag$BagMod)
# varImpPlot(out.Bag$BagMod)
# out.Bag$BagAri
# out.Bag$BagTab




##############################
# Tree Prediction

predTree <- function(formula = as.formula("Class ~ ."),
                     df.Final,
                     df.Index.Train,
                     fType = "class",
                     fClass = "Class") {
  
  datTree <- tree(formula, 
                  data = df.Final, 
                  subset=df.Index.Train,
                  method=fType)
  
  datTree.pred <- predict(datTree,
                          df.Final[-df.Index.Train,],
                          type=fType)
  
  datTree.tab <- table(df.Final[-df.Index.Train, fClass],
                     datTree.pred)
  
  datTree.ari <- classAgreement(datTree.tab)$crand
  return(setNames(list(datTree, datTree.tab, datTree.ari), c("TreeMod","TreeTab", "TreeAri")))
}

# out.Tree <- predTree(df.Final = df.raw.forA, df.Index.Train = df.train.index)
# summary(out.Tree$TreeMod)
# plot(out.Tree$TreeMod)
# text(out.Tree$TreeMod)
# out.Tree$TreeAri
# out.Tree$TreeTab



##############################
# Boosting

# Boosting Prediction
predBoost <- function(formula = as.formula("Class ~ ."),
                       df.Final,
                       df.Index.Train,
                       fType = "response",
                       fClass = "Class",
                       fNtrees = 5000,
                       fInteraction.Depth = 4,
                       fDistribution = "multinomial") {
  
  datBoost <- gbm(
    formula,
    data = df.Final[df.Index.Train, ],
    distribution = fDistribution,
    n.trees = fNtrees,
    interaction.depth = fInteraction.Depth
  )
  
  datBoost.pred.soft <- predict(
    datBoost,
    newdata = df.Final[-df.Index.Train, ],
    n.trees = fNtrees,
    distribution = fDistribution,
    type = fType
  )
  
  # For each row, return the index of the column with the maximum number
  datBoost.pred.hard <- apply(datBoost.pred.soft, 1, function(xrow) {
    return(which.max(xrow))
  })
  
  datBoost.tab <- table(
    df.Final[-df.Index.Train,fClass],
    datBoost.pred.hard)
  
  datBoost.ari <- classAgreement(datBoost.tab)$crand
  return(setNames(list(datBoost, datBoost.tab, datBoost.ari), c("BoostMod","BoostTab", "BoostAri")))
}

# out.Boost <- predBoost(df.Final = df.raw.forA, df.Index.Train = df.train.index)
# out.Boost$BoostAri
# out.Boost$BoostTab


#########################################
#  Mclust Prediction

# Mclust Prediction
predMclust <- function(fIndexResponse = 1,
                       df.Final.Num,
                       df.Index.Train) {
  
  datMclustDA <- MclustDA(
    df.Final.Num[df.Index.Train, -fIndexResponse] , 
    df.Final.Num[df.Index.Train, fIndexResponse])
  
  # Regular Summary
  # datMclustDA.sum <- summary(datMclustDA, parameters = TRUE)

  # Summary with Prediction
  datMclustDA.predict <- summary(datMclustDA,
                                 newdata = df.Final.Num[-df.Index.Train,-fIndexResponse], 
                                 newclass = df.Final.Num[-df.Index.Train,fIndexResponse])
  
  datMclustDA.tab <- datMclustDA.predict$tab.newdata

  datMclustDA.ari <- classAgreement(datMclustDA.predict$tab.newdata)$crand
  return(setNames(list(datMclustDA, datMclustDA.tab, datMclustDA.ari), c("MclustMod","MclustTab", "MclustAri")))
  
}
# out.Mclust <- predMclust(df.Final.Num = df.raw.forA.num, df.Index.Train = df.train.index)
# summary(out.Mclust$MclustMod, parameters = T)
# summary(out.Mclust$MclustMod, params = T)
# out.Mclust$MclustTab
# out.Mclust$MclustAri




#########################################
# nnet

# nnet Prediction. 
predNnet <- function(formula = as.formula("Class ~ ."),
                     df.Final,
                     df.Index.Train,
                     fType = "class",
                     fClassName = "Class",
                     fSize = 5,
                     fMaxit = 500,
                     fDecay = 0) {
  
  nnet.sum <- nnet(formula,
                   data = df.Final[df.Index.Train,], 
                   size=5, maxit=500,decay=0, trace = F)
  
  nnet.pred <- predict(nnet.sum, 
                       df.Final[-df.Index.Train,], 
                       type=fType)
  nnet.tab <- table(df.Final[-df.Index.Train,1],
                    nnet.pred)
  
  nnet.ari <- classAgreement(nnet.tab)$crand
  
  return(setNames(list(nnet.pred, nnet.tab, nnet.ari), c("NnetMod","NnetTab", "NnetAri")))
  
}

# out.Nnet <- predNnet(df.Final = df.raw.forA, df.Index.Train = df.train.index)
# summary(out.Nnet$NnetMod)
# out.Nnet$NnetTab
# out.Nnet$NnetAri

#########################################

# Interface to use all methods
predAll <- function(df.Final, df.Index.Train, df.Final.Num, fMtry = ceiling((ncol(df.Final) - 1)/2),DoBagMethod = T) {
  
  o.Nnet <- predNnet(df.Final = df.Final, df.Index.Train = df.Index.Train)
  # o.Mclust <- predMclust(df.Final.Num = df.Final.Num, df.Index.Train = df.Index.Train)
  o.Mclust <- NA
  if (DoBagMethod) {
    o.Boost <- predBoost(df.Final = df.Final, df.Index.Train = df.Index.Train)
  } else {
    o.Boost <- NA
  }
  o.Tree <- NA
  o.Bag <- NA
  # o.Tree <- predTree(df.Final = df.Final, df.Index.Train = df.Index.Train)
  # o.Bag <- predBag(df.Final = df.Final, df.Index.Train = df.Index.Train)
  o.Forest <- predForest(df.Final = df.Final, df.Index.Train = df.Index.Train, fMtry = fMtry)  
  
  o.all <- c(o.Tree, o.Forest, o.Bag, o.Boost, o.Mclust, o.Nnet)
  return(o.all)
}

# df.Final <- df.raw.forA
# df.Index.Train <- df.train.index
# df.Final.Num <- df.raw.forA.num
# out.All <- predAll(df.raw.forA, df.train.index, df.raw.forA.num)

# summary(out.All$TreeMod)
# plot(out.All$TreeMod)
# text(out.All$TreeMod)
# out.All$TreeAri
# out.All$TreeTab
# 
# importance(out.All$ForestMod)
# varImpPlot(out.All$ForestMod)
# out.All$ForestAri
# out.All$ForestTab
# 
# 
# importance(out.All$BagMod)
# varImpPlot(out.All$BagMod)
# out.All$BagAri
# out.All$BagTab
# 
# out.All$BoostAri
# out.All$BoostTab
# 
# summary(out.All$MclustMod, parameters = T)
# summary(out.All$MclustMod, params = T)
# out.All$MclustTab
# out.All$MclustAri
# 
# summary(out.All$NnetMod)
# out.All$NnetTab
# out.All$NnetAri
######################################
# GRAPHS

# Stacked Bar graph of covariates
grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]],
             p[[6]], p[[7]], p[[8]], p[[9]], p[[10]],
             p[[11]], p[[12]], p[[13]], p[[14]], p[[15]],
             p[[16]], p[[17]], p[[18]], p[[19]], p[[20]],
             p[[21]], p[[22]],
             nrow = 5, ncol = 5)

# Stacked Bar Graph
ggplot(df.Plot, aes(x = StalkColorAboveRing, fill = Class)) + 
  geom_bar() + 
  theme(legend.position = "none",
        text = element_text(size=6),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7))

# Adding noise to a numerical version of the factor data
df.raw.forA.num.noise <- apply( df.raw.forA.num[,-1], 2, function(x) {
  return(x + rnorm(length(x),mean = 0, sd = 0.1))
})

df.raw.forA.num.noise <- cbind(df.raw.forA.num[,1], df.raw.forA.num.noise) %>% as.data.frame()
colnames(df.raw.forA.num.noise) <- names(df.raw.forA)
df.raw.forA.num.noise$Class <- df.raw.forA.num.noise$Class %>% as.factor()
names(df.raw.forA.num.noise)

ggparcoord(
  df.raw.forA.num.noise, columns=2:ncol(df.raw.forA.num.noise), 
  alphaLines=0.1, groupColumn="Class", 
  scale="globalminmax") + 
  xlab("Covariate") + 
  ylab("Level") 



############################################
# Removing Covariates, Finding minimal set
cores = detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer. Leaves for regular use
registerDoParallel(cl)

#Variables to Choose Initally
# Take covariates involved in the rule with the highest standardized lift.
# 1 covariates that indicates POisonous and 1covariate that indicates edible
df.raw.arul <- subset(df.rNA) %>% as.data.frame()
params <- list(support=0.055,confidence=0.8,minlen=2,maxlen=4)
out.aRul <- predArules(df.raw.arul, fmin = 2, fmax = 4)

# 2 Variables
df.raw.final <- subset(df.rNA, select = c(Class, Odor, StalkColorBelowRing)) %>% as.data.frame()
df.raw.arul <- subset(df.rNA, select = -c(Odor, StalkColorBelowRing)) %>% as.data.frame()

params <- list(support=0.055,confidence=0.8,minlen=2,maxlen=4)
out.aRul <- predArules(df.raw.arul, fmin = 2, fmax = 4)
# Getting ARI for all methods
set.seed(102032)

test.AllARI1 <-
  foreach(
    index = 1:30, # Loops over the choice of alleles
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("dplyr", "randomForest", "e1071", 
                  "mclust", "tree", "nnet", "gbm")
  ) %dopar% {
    df.raw.forA <- df.raw.final
    df.raw.forA.num <- data.matrix(df.raw.forA)
    
    # Training set
    df.train.index <- sample(
      1:length(df.raw.forA$Class), 0.50*length(df.raw.forA$Class))
    
    out.All <- predAll(df.raw.forA, df.train.index, df.raw.forA.num, DoBagMethod = F, fMtry = 4)
    out.A.Ari <- setNames(
      list(out.All$TreeAri, out.All$ForestAri, 
           out.All$BagAri, out.All$BoostAri,
           out.All$MclustAri, out.All$NnetAri),
      c("Tree", "Forest", "Bag", "Boost", "Mclust", "Nnet"))
    out.df.row <- out.A.Ari %>% t %>% data.frame()
    rownames(out.df.row) <- paste("ARI", index, sep = "")
    
    return(out.df.row)
  } %>% apply(1, function(x) unlist(x)) %>% t
test.AllARI1.withv <- cbind(test.AllARI1, rep(2, nrow(test.AllARI1)))
################
# BoxPlot ARI 4 Variables
## 4 Variables
df.raw.final <- subset(df.rNA, select = c(Class, 
                                          Odor, StalkColorBelowRing,
                                          Population, StalkColorAboveRing)) %>% as.data.frame()

df.raw.arul <- subset(df.rNA, select = -c(Odor, StalkColorBelowRing, 
                                          Population, StalkColorAboveRing)) %>% as.data.frame()
params <- list(support=0.055,confidence=0.8,minlen=2,maxlen=4)
out.aRul <- predArules(df.raw.arul, fmin = 2, fmax = 4)
# Getting ARI for all methods
set.seed(102032)

test.AllARI2 <-
  foreach(
    index = 1:30, # Loops over the choice of alleles
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("dplyr", "randomForest", "e1071", 
                  "mclust", "tree", "nnet", "gbm")
  ) %dopar% {
    df.raw.forA <- df.raw.final
    df.raw.forA.num <- data.matrix(df.raw.forA)
    
    # Training set
    df.train.index <- sample(
      1:length(df.raw.forA$Class), 0.50*length(df.raw.forA$Class))
    
    out.All <- predAll(df.raw.forA, df.train.index, df.raw.forA.num, DoBagMethod = F, fMtry = 4)
    out.A.Ari <- setNames(
      list(out.All$TreeAri, out.All$ForestAri, 
           out.All$BagAri, out.All$BoostAri,
           out.All$MclustAri, out.All$NnetAri),
      c("Tree", "Forest", "Bag", "Boost", "Mclust", "Nnet"))
    out.df.row <- out.A.Ari %>% t %>% data.frame()
    rownames(out.df.row) <- paste("ARI", index, sep = "")
    
    return(out.df.row)
  } %>% apply(1, function(x) unlist(x)) %>% t
test.AllARI2.withv <- cbind(test.AllARI2, rep(4, nrow(test.AllARI2)))

################
# BoxPlot ARI 7 Variables
## 7 Variables
df.raw.final <- subset(df.rNA, select = c(Class, 
                                          Odor, StalkColorBelowRing,
                                          Population, StalkColorAboveRing,
                                          StalkSurfaceBelowRing, GillSize, GillColor)) %>% as.data.frame()

df.raw.arul <- subset(df.rNA, select = -c(Odor, StalkColorBelowRing, 
                                          Population, StalkColorAboveRing,
                                          StalkSurfaceBelowRing, GillSize, GillColor)) %>% as.data.frame()
params <- list(support=0.055,confidence=0.8,minlen=2,maxlen=4)
out.aRul <- predArules(df.raw.arul, fmin = 2, fmax = 4)
# Getting ARI for all methods

set.seed(102032)

test.AllARI3 <-
  foreach(
    index = 1:30, # Loops over the choice of alleles
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("dplyr", "randomForest", "e1071", 
                  "mclust", "tree", "nnet", "gbm")
  ) %dopar% {
    df.raw.forA <- df.raw.final
    df.raw.forA.num <- data.matrix(df.raw.forA)
    
    # Training set
    df.train.index <- sample(
      1:length(df.raw.forA$Class), 0.50*length(df.raw.forA$Class))
    
    out.All <- predAll(df.raw.forA, df.train.index, df.raw.forA.num, DoBagMethod = F, fMtry = 4)
    out.A.Ari <- setNames(
      list(out.All$TreeAri, out.All$ForestAri, 
           out.All$BagAri, out.All$BoostAri,
           out.All$MclustAri, out.All$NnetAri),
      c("Tree", "Forest", "Bag", "Boost", "Mclust", "Nnet"))
    out.df.row <- out.A.Ari %>% t %>% data.frame()
    rownames(out.df.row) <- paste("ARI", index, sep = "")
    
    return(out.df.row)
  } %>% apply(1, function(x) unlist(x)) %>% t
test.AllARI3.withv <- cbind(test.AllARI3, rep(7, nrow(test.AllARI3)))

################
# BoxPlot ARI 10 Variables
## 10 Variables
df.raw.final <- subset(df.rNA, select = c(Class, 
                                          Odor, StalkColorBelowRing,
                                          Population, StalkColorAboveRing,
                                          StalkSurfaceBelowRing, GillSize, GillColor,
                                          RingColor, Bruises, StalkRoot)) %>% as.data.frame()

df.raw.arul <- subset(df.rNA, select = -c(Odor, StalkColorBelowRing, 
                                          Population, StalkColorAboveRing,
                                          StalkSurfaceBelowRing, GillSize, GillColor,
                                          RingColor, Bruises, StalkRoot)) %>% as.data.frame()
params <- list(support=0.055,confidence=0.8,minlen=2,maxlen=4)
out.aRul <- predArules(df.raw.arul, fmin = 2, fmax = 4)
# Getting ARI for all methods
set.seed(102032)

test.AllARI4 <-
  foreach(
    index = 1:30, # Loops over the choice of alleles
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("dplyr", "randomForest", "e1071", 
                  "mclust", "tree", "nnet", "gbm")
  ) %dopar% {
    df.raw.forA <- df.raw.final %>% as.data.frame()
    df.raw.forA.num <- data.matrix(df.raw.forA)
    
    # Training set
    df.train.index <- sample(
      1:length(df.raw.forA$Class), 0.50*length(df.raw.forA$Class))
    
    out.All <- predAll(df.raw.forA, df.train.index, df.raw.forA.num, DoBagMethod = F, fMtry = 4)
    out.A.Ari <- setNames(
      list(out.All$TreeAri, out.All$ForestAri, 
           out.All$BagAri, out.All$BoostAri,
           out.All$MclustAri, out.All$NnetAri),
      c("Tree", "Forest", "Bag", "Boost", "Mclust", "Nnet"))
    out.df.row <- out.A.Ari %>% t %>% data.frame()
    rownames(out.df.row) <- paste("ARI", index, sep = "")
    
    return(out.df.row)
  } %>% apply(1, function(x) unlist(x)) %>% t
test.AllARI4.withv <- cbind(test.AllARI4, rep(10, nrow(test.AllARI4)))


###################
# BoxPlot ARI for 13 Variables
## 13 Variables
df.raw.final <- subset(df.rNA, select = c(Class, 
                                          Odor, StalkColorBelowRing,
                                          Population, StalkColorAboveRing,
                                          StalkSurfaceBelowRing, GillSize, GillColor,
                                          RingColor, Bruises, StalkRoot,
                                          StalkSurfaceAboveRing, CapColor, StalkShape)) %>% as.data.frame()

df.raw.arul <- subset(df.rNA, select = -c(Odor, StalkColorBelowRing, 
                                          Population, StalkColorAboveRing,
                                          StalkSurfaceBelowRing, GillSize, GillColor,
                                          RingColor, Bruises, StalkRoot,
                                          StalkSurfaceAboveRing, CapColor, StalkShape)) %>% as.data.frame()
params <- list(support=0.055,confidence=0.8,minlen=2,maxlen=4)
out.aRul <- predArules(df.raw.arul, fmin = 2, fmax = 4)
# Getting ARI for all methods
set.seed(102032)

test.AllARI5 <-
  foreach(
    index = 1:30, # Loops over the choice of alleles
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("dplyr", "randomForest", "e1071", 
                  "mclust", "tree", "nnet", "gbm")
  ) %dopar% {
    df.raw.forA <- df.raw.final
    df.raw.forA.num <- data.matrix(df.raw.forA)
    
    # Training set
    df.train.index <- sample(
      1:length(df.raw.forA$Class), 0.50*length(df.raw.forA$Class))
    
    out.All <- predAll(df.raw.forA, df.train.index, df.raw.forA.num, DoBagMethod = F, fMtry = 4)
    out.A.Ari <- setNames(
      list(out.All$TreeAri, out.All$ForestAri, 
           out.All$BagAri, out.All$BoostAri,
           out.All$MclustAri, out.All$NnetAri),
      c("Tree", "Forest", "Bag", "Boost", "Mclust", "Nnet"))
    out.df.row <- out.A.Ari %>% t %>% data.frame()
    rownames(out.df.row) <- paste("ARI", index, sep = "")
    
    return(out.df.row)
  } %>% apply(1, function(x) unlist(x)) %>% t
test.AllARI5.withv <- cbind(test.AllARI5, rep(13, nrow(test.AllARI5)))

##################
# BoxPlot ARI for 16 Variables
## 16 Variables
df.raw.final <- subset(df.rNA, select = c(Class, 
                                          Odor, StalkColorBelowRing,
                                          Population, StalkColorAboveRing,
                                          StalkSurfaceBelowRing, GillSize, GillColor,
                                          RingColor, Bruises, StalkRoot,
                                          StalkSurfaceAboveRing, CapColor, StalkShape,
                                          SporePrintColor, GillSpacing, Habitat)) %>% as.data.frame()

df.raw.arul <- subset(df.rNA, select = -c(Odor, StalkColorBelowRing, 
                                          Population, StalkColorAboveRing,
                                          StalkSurfaceBelowRing, GillSize, GillColor,
                                          RingColor, Bruises, StalkRoot,
                                          StalkSurfaceAboveRing, CapColor, StalkShape,
                                          SporePrintColor, GillSpacing, Habitat)) %>% as.data.frame()
params <- list(support=0.055,confidence=0.8,minlen=2,maxlen=4)
out.aRul <- predArules(df.raw.arul, fmin = 2, fmax = 4)
# Getting ARI for all methods
set.seed(102032)



test.AllARI6 <-
  foreach(
    index = 1:30, # Loops over the choice of alleles
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("dplyr", "randomForest", "e1071", 
                  "mclust", "tree", "nnet", "gbm")
  ) %dopar% {
    df.raw.forA <- df.raw.final %>% as.data.frame()
    df.raw.forA.num <- data.matrix(df.raw.forA)
    
    # Training set
    df.train.index <- sample(
      1:length(df.raw.forA$Class), 0.50*length(df.raw.forA$Class))
    
    out.All <- predAll(df.raw.forA, df.train.index, df.raw.forA.num, DoBagMethod = F, fMtry = 4)
    out.A.Ari <- setNames(
      list(out.All$TreeAri, out.All$ForestAri, 
           out.All$BagAri, out.All$BoostAri,
           out.All$MclustAri, out.All$NnetAri),
      c("Tree", "Forest", "Bag", "Boost", "Mclust", "Nnet"))
    out.df.row <- out.A.Ari %>% t %>% data.frame()
    rownames(out.df.row) <- paste("ARI", index, sep = "")
    
    return(out.df.row)
  } %>% apply(1, function(x) unlist(x)) %>% t
test.AllARI6.withv <- cbind(test.AllARI6, rep(16, nrow(test.AllARI6)))

##############
# Combining all variable ARIs to form Box plot

AriVarTest <- (rbind(
  test.AllARI1.withv,
  test.AllARI2.withv, 
  test.AllARI3.withv,
  test.AllARI4.withv,
  test.AllARI5.withv,
  test.AllARI6.withv)) %>% as.data.frame()
colnames(AriVarTest) <- c("Forest", "Nnet", "VariablesUsed")
rownames(AriVarTest) <- c(1:nrow(AriVarTest))
AriVarTest 

dfCleanAri <- melt(AriVarTest, "VariablesUsed")
colnames(dfCleanAri) <- c("VariablesUsed", "Method", "ARI")
dfCleanAri$VariablesUsed <- dfCleanAri$VariablesUsed %>% as.factor()
dfPlot <- dfCleanAri

# 2 data points removed so box plot is more readable
dfPlot %>% filter(ARI <= 0.85)
dfPlotBox <- dfPlot %>% filter(ARI > 0.85)
ggplot(dfPlotBox, aes(x=VariablesUsed, y=ARI, fill = Method)) +
  geom_boxplot() + 
  theme_tufte() + theme(axis.line=element_line()) +
  ggtitle("Comparing ARI Performance of Neutral Nets and Random Forests")

###################
# Optimizing Parameters on 13 Variables 13 VARIABLES

## 13 Variables
df.raw.final <- subset(df.rNA, select = c(Class, 
                                          Odor, StalkColorBelowRing,
                                          Population, StalkColorAboveRing,
                                          StalkSurfaceBelowRing, GillSize, GillColor,
                                          RingColor, Bruises, StalkRoot,
                                          StalkSurfaceAboveRing, CapColor, StalkShape)) %>% as.data.frame()

mTryChoice <- seq(1, 13, by=1)
nTreeChoice <- seq(300, 600, by=50)
mTryTreeComboMat <- expand.grid(mTryChoice,nTreeChoice)
colnames(mTryTreeComboMat) <- c("mTry", "nTree")
mTryTreeComboMat


# index <- 1
cores = detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer. Leaves for regular use
registerDoParallel(cl)
set.seed(102032)

opt.ForestParam <-
  foreach(
    index = 1:nrow(mTryTreeComboMat), # Loops over the choice of parameters
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("dplyr", "randomForest", "e1071", 
                  "mclust", "tree", "nnet", "gbm")
  ) %dopar% {

    cmTry <- mTryTreeComboMat[index,1]
    cnTree <- mTryTreeComboMat[index,2]
    
    # Average of 25 ARI's
    ariVec <- sapply(c(1:25), function(x) {
      df.raw.forA <- df.raw.final %>% as.data.frame()
      
      # Training set
      df.train.index <- sample(
        1:length(df.raw.forA$Class), 0.50*length(df.raw.forA$Class))
      
      out.All <- predForest(df.Final = df.raw.forA, df.Index.Train = df.train.index, fMtry = cmTry, fnTree = cnTree)
      return(out.All$ForestAri)
    })

    out.df.row <- c(mean(ariVec), cmTry, cnTree)
    return(out.df.row)
  }
colnames(opt.ForestParam) <- c("AverageARI", "mTry", "nTree")
stopCluster(cl)
df.opt.For <- opt.ForestParam %>% as.data.frame()
df.opt.For.o <- df.opt.For[order(df.opt.For$AverageARI, decreasing = T),]
df.opt.For.o

df.opt.For.o.print <- rbind(df.opt.For.o[1:5,], df.opt.For.o[(nrow(opt.ForestParam) - 4):nrow(opt.ForestParam),])
df.opt.For.o.print

print(
  xtable(  df.opt.For.o.print, type = "latex", digits = c(1,4,0,0)), 
  file = "ttForParam.txt", 
  include.rownames=F
)

ggplot(df.opt.For.o, aes(x = mTry, y = nTree, fill = AverageARI)) +
  geom_point(size=2, shape=23) 


#############
# Neural Net Optimization for 13 VARIABLES
df.raw.final <- subset(df.rNA, select = c(Class, 
                                          Odor, StalkColorBelowRing,
                                          Population, StalkColorAboveRing,
                                          StalkSurfaceBelowRing, GillSize, GillColor,
                                          RingColor, Bruises, StalkRoot,
                                          StalkSurfaceAboveRing, CapColor, StalkShape)) %>% as.data.frame()
fSizeChoice <- seq(2, 7, by=1)
fMaxitChoice <- seq(100, 500, by=100)
fDecayChoice <- seq(0, 0.5, by=0.25)

mPComboNnetMat <- expand.grid(fSizeChoice,fMaxitChoice, fDecayChoice)
colnames(mPComboNnetMat) <- c("fSize", "fMaxit", "fDecay")
mPComboNnetMat

cores = detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer. Leaves for regular use
registerDoParallel(cl)

set.seed(102032)
opt.NnetParam <-
  foreach(
    index = 1:nrow(mPComboNnetMat), # Loops over the choice of parameters
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("dplyr", "randomForest", "e1071", 
                  "mclust", "tree", "nnet", "gbm")
  ) %dopar% {
    
    cmSize <- mPComboNnetMat[index,1]
    cnMax <- mPComboNnetMat[index,2]
    cnDecay <- mPComboNnetMat[index,3]
    
    # Average of 25 ARI's
    ariVec <- sapply(c(1:25), function(x) {
      df.raw.forA <- df.raw.final %>% as.data.frame()
      
      # Training set
      df.train.index <- sample(
        1:length(df.raw.forA$Class), 0.50*length(df.raw.forA$Class))
      
      out.All <- predNnet(df.Final = df.raw.forA, df.Index.Train = df.train.index,
                          fSize = cmSize, fMaxit = cnMax, fDecay = cnDecay )
      
      return(out.All$NnetAri)
    })

    out.df.row <- c(mean(ariVec), cmSize, cnMax, cnDecay)
    return(out.df.row)
  }
stopCluster(cl)
  colnames(opt.NnetParam) <- c("AverageARI", "Size", "Maxit", "Decay")
df.opt.Nnet <- opt.NnetParam %>% as.data.frame()
df.opt.Nnet.o <- df.opt.Nnet[order(df.opt.Nnet$AverageARI, decreasing = T),]
df.opt.Nnet.o
df.opt.Nnet.o.print <- rbind(df.opt.Nnet.o[1:5,], df.opt.Nnet.o[(nrow(opt.NnetParam) - 4):nrow(opt.NnetParam),])
df.opt.Nnet.o.print
print(
  xtable(  df.opt.Nnet.o.print, type = "latex", digits = c(1,4,0,0,2)), 
  file = "ttNnetParam.txt", 
  include.rownames=F
)
###############################
# OPTIMIZED RANDOM FOREST And NNET
# 13 VARIABLES

df.raw.final <- subset(df.rNA, select = c(Class, 
                                          Odor, StalkColorBelowRing,
                                          Population, StalkColorAboveRing,
                                          StalkSurfaceBelowRing, GillSize, GillColor,
                                          RingColor, Bruises, StalkRoot,
                                          StalkSurfaceAboveRing, CapColor, StalkShape)) %>% as.data.frame()
df.raw.forA <- df.raw.final 

set.seed(102034)
sum.opt <- apply(diag(25), 2, function(x) {
  df.train.index <- sample(
    1:length(df.raw.forA$Class), 0.50*length(df.raw.forA$Class))
  
  opt.pred.net <- predNnet(df.Final = df.raw.forA, df.Index.Train = df.train.index,
                           fSize = 7, fMaxit = 500, fDecay = 0 )
  
  opt.pred.for <- predForest(df.Final = df.raw.forA, 
                             df.Index.Train = df.train.index, fMtry = 4, fnTree = 350)
  
  # opt.pred.net$NnetAri
  # opt.pred.for$ForestAri
  A <- opt.pred.for$ForestTab
  B <- opt.pred.net$NnetTab
  # sum(A[row(A)!=col(A)])
  # 
  # opt.pred.for$ForestTab
  # opt.pred.net$NnetTab
  c(opt.pred.net$NnetAri, opt.pred.for$ForestAri, sum(A[row(A)!=col(A)]), sum(B[row(B)!=col(B)]))
}) %>% t
sum.opt
# importance(out.Forest$ForestMod)
# varImpPlot(out.Forest$ForestMod)
sum.opt %>% summary

sum.opt[,1] %>% median
sum.opt[,2] %>% median
sum.opt[,3] %>% median
sum.opt[,4] %>% median

sum.opt[,1] %>% var
sum.opt[,2] %>% var
sum.opt[,3] %>% var
sum.opt[,4] %>% var

# ARI Plot
colnames(sum.opt) <- c("Nnet", "Forest", "ForestM", "NnetM")
sun.opt.plot.ari <- melt(sum.opt[,c(1,2)])
colnames(sun.opt.plot.ari) <- c("rm", "Method", "ARI")
sun.opt.plot.ari

ggplot(sun.opt.plot.ari, aes(x=Method, y=ARI)) +
  geom_boxplot() + 
  theme_tufte() + theme(axis.line=element_line()) +
  ggtitle("ARI of Optimized Neutral Nets and Random Forests")+ 
  theme(legend.position="none")


# Misclassification Number Plot
colnames(sum.opt) <- c("NnetM", "ForestM", "Forest", "Nnet")
sun.opt.plot.mis <- melt(sum.opt[,-c(1,2)])
colnames(sun.opt.plot.mis) <- c("rm", "Method", "MisclassificationNumber")
sun.opt.plot.mis %>% summary

ggplot(sun.opt.plot.mis, aes(x=Method, y=MisclassificationNumber)) +
  geom_boxplot() + 
  theme_tufte() + theme(axis.line=element_line()) +
  ggtitle("Misclassification Number of Optimized Neutral Nets and Random Forests") + 
  theme(legend.position="none")


# Forest
out.Forest <- predForest(df.Final = df.raw.forA, 
                           df.Index.Train = df.train.index, fMtry = 4, fnTree = 350)
importance(out.Forest$ForestMod)
varImpPlot(out.Forest$ForestMod, main="Random Forest Variable Contribution")
?varImpPlot

#####################################

stopCluster(cl)

# Parallel Coordinate Plot of Reduced Covaraites
# Adding noise to a numerical version of the factor data
df.raw.final <- subset(df.rNA, select = c(Class, 
                                          Odor, StalkColorBelowRing,
                                          Population, StalkColorAboveRing,
                                          StalkSurfaceBelowRing, GillSize, GillColor,
                                          RingColor, Bruises, StalkRoot,
                                          StalkSurfaceAboveRing, CapColor, StalkShape,
                                          SporePrintColor, GillSpacing, Habitat)) %>% as.data.frame()
df.raw.forA <- df.raw.final
df.raw.forA.num <- data.matrix(df.raw.forA)

df.raw.forA.num.noise <- apply(df.raw.forA.num[,-1], 2, function(x) {
  return(x + rnorm(length(x),mean = 0, sd = 0.1))
})

df.raw.forA.num.noise <- cbind(df.raw.forA.num[,1], df.raw.forA.num.noise) %>% as.data.frame()
colnames(df.raw.forA.num.noise) <- names(df.raw.forA)
df.raw.forA.num.noise$Class <- df.raw.forA.num.noise$Class %>% as.factor()
names(df.raw.forA.num.noise)

ggparcoord(
  df.raw.forA.num.noise, columns=2:ncol(df.raw.forA.num.noise), 
  alphaLines=0.1, groupColumn="Class", 
  scale="globalminmax") + 
  xlab("Covariates") + 
  ylab("Level") +
  ggtitle("Separation of Edible and Poisonous on Chosen Covariates")



