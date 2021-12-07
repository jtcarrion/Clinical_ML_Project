###   311 FINAL PROJECT   ###
###   TEAM 1 ~ DATASET 1   ###


############################################################
 ###   SET DIRECTORY  ### 
###   LOAD LIBRARIES    ###

setwd("C:/Users/labuser/Desktop/")

install.packages("lattice")
install.packages("e1071")
install.packages("caret")
install.packages("pROC")
install.packages("plyr")
install.packages("rgl")
install.packages("rpart")

library(rgl)
library(plyr)
library(pROC)
library(class)
library(caret)
library(lattice)
library(e1071)
library(ggplot2)
library(rpart)
library(rpart.plot)
library(rpartScore)
############################################################
###   READ IN DATA   ###

data= read.csv("breastCancer.csv", header = T, sep = ",", stringsAsFactors = FALSE)
dim(data)
head(data)

features= data[,2:ncol(data)]
dim(features)
head(features)
############################################################
###   CLEAN DATA  ###

badSamples= which(features$bare_nucleoli=='?')  # FIND MISSING DATA

cleanFeatures= features[-badSamples,]   # REMOVE DATA
dim(cleanFeatures)

### SINCE '?' MADE COLUMN STRING VARS, RESTORE BACK TO INT VARS ###
bare_nuc= as.integer(cleanFeatures$bare_nucleoli)
newFeatures= cbind.data.frame(cleanFeatures[,-6],bare_nuc)  # REMOVE OLD COLUMN, ADD NEW 

dim(newFeatures)
summary(newFeatures) # newFeatures consists of no missing values and all int vars
                                        # DID NOT SCALE #
class= newFeatures[,9]
pureFeat= newFeatures[,-9]

############################################################
###   MEAN ANALYSIS OF 3 VARS   ###

agg = aggregate(newFeatures[,2:4], list(newFeatures$class), function(x) c(mean=mean(x), sd=sd(x))) 
agg

########################################################
###   HISTOGRAM   ###
malignant = subset(data , class == "2") #creates a subset of the data that only includes rows where class is 2 
dim(malignant)

benign = subset(data, class == "4") #creates a subset of the data that only hold rows where class is 4
dim(benign)


clumpthickness_mal = data.frame(Clump_Thickness= rnorm(458, 2.96, 1.67))
clumpthickness_ben = data.frame(Clump_Thickness= rnorm(241, 7.19, 2.43))


clumpthickness_mal$variables = 'Malignant_Clump_Thickness'
clumpthickness_ben$variables = 'Benign_Clump_Thickness'

thickness = rbind(clumpthickness_ben, clumpthickness_mal)

ggplot(thickness, aes(Clump_Thickness, fill= variables)) +
  geom_histogram(alpha = 0.5, position = 'identity', binwidth = 2) +
  labs(y = 'Count')
##################################################################

#################################################################
####   SPLIT SAMPLES 4:1   ###
set.seed(16)
dt = sort(sample(nrow(newFeatures), nrow(newFeatures)*.75))
#dt= createDataPartition(newFeatures, p= .75, list = FALSE)

train = newFeatures[dt,]
train = as.data.frame(train)

test<-newFeatures[-dt,]
test= as.data.frame(test)


dim(train)
head(train)
dim(test)
head(test)

test_class= test[,9]
test_data= test[,-9]

############################################################################
###   PCA ANALYSIS   ###

pca <- prcomp(train, scale = T, center= T, retx = T)
summary(pca)

###   1ST 3 Principle Components retinas 82% of the variance   ###
## 5 Principle Components are needed in order to retain 90% variance  ##

rgl.open()
bg3d("white")
plot3d(pca$x[which(train$class == 2), 1:3], col = "blue")
points3d(pca$x[which(train$class == 4), 1:3], col = "red")

###########################################################################
###   10-fold cross validation   ###
##############################################################################

test_class
train$class[train$class==2] <- "benign"
train$class[train$class==4] <- "malignant"

test$class[test$class==2] <- "benign"
test$class[test$class==4] <- "malignant"


specs <- trainControl(method = "cv", number = 10,
                      savePredictions = "all",
                      classProbs = TRUE) #cross validation of 10 folds, save all predictions and class probabilities

set.seed(16) #set seed for random assignment of each sample

#names(getModelInfo())

model5 <- train(class~., data = train,
  method = "rpart", trControl = specs, 
  control= rpart.control(minsplit = 5)) #class is regressed onto these attributes, data gotten from train, method is generalized linear model

model25 <- train(class~., data = train,
          method = "rpart", trControl = specs, 
          control= rpart.control(minsplit = 25))

model50 <- train(class~., data = train,
                method = "rpart", trControl = specs, 
                control= rpart.control(minsplit = 50))

model75 <- train(class~., data = train,
                method = "rpart", trControl = specs, 
                control= rpart.control(minsplit = 75))

model100 <- train(class~., data = train,
                method = "rpart", trControl = specs, 
                control= rpart.control(minsplit = 100))

model1D <- train(class ~. , data = train,
                method= 'rpart', trControl= specs,
                control = rpart.control(maxdepth = 1, minsplit = 25))

model3D <- train(class ~. , data = train,
                 method= 'rpart', trControl= specs,
                 control = rpart.control(maxdepth = 3, minsplit = 25))

model7D <- train(class ~. , data = train,
                  method= 'rpart', trControl= specs,
                  control = rpart.control(maxdepth = 7, minsplit = 25))


confusionMatrix(model5)
confusionMatrix(model25)
confusionMatrix(model50)
confusionMatrix(model75)
confusionMatrix(model100)
confusionMatrix(model1D)
confusionMatrix(model3D)
confusionMatrix(model7D)

##############################################################################
###    ROC CURVE   ###
par(mfrow=c(1,1))

model5$pred
model25$pred
model50$pred
model75$pred
model100$pred

plot.roc(model5$pred$obs,
         model5$pred$malignant, print.auc= T, main= "Min Split 5")

plot.roc(model7D$pred$obs,
         model7D$pred$malignant, print.auc= T, main= "Max Depth 7")

plot.roc(model25$pred$obs,
         model25$pred$malignant, print.auc= T, main= "Min Split 25")

plot.roc(model50$pred$obs,
         model50$pred$malignant, print.auc= T, main= "Min Split 50")

plot.roc(model75$pred$obs,
         model75$pred$malignant, print.auc= T, main= "Min Split 75")

plot.roc(model100$pred$obs,
         model100$pred$malignant, print.auc= T, main= "Min Split 100")


###   COMPARE   ###
plot.roc(model25$pred$obs,
         model25$pred$malignant, main= "Min Split vs. Max Depth", col= "red")

plot.roc(model7D$pred$obs,
         model7D$pred$malignant, add= T, col= "blue")

plot.roc(model1D$pred$obs,
         model1D$pred$malignant, add= T, col= "green")

legend(0.05,0.15, legend = c("Min Split 25", "Max Depth 7", "Max Depth 1"), 
       col = c('red', 'blue', 'green'), lty = 1)


###############################################################################
###   TEST CLASS   ###

pre5= predict(model5, test_data, type = "raw")
pre5

pre7D = predict(model7D, test_data, type = "raw")
pre7D

pre25= predict(model25, test_data, type = "raw")
pre25

pre50= predict(model50, test_data, type = "raw")

pre75= predict(model75, test_data, type = "raw")

pre100= predict(model100, test_data, type = "raw")

test$class[test$class==2] <- "benign"
test$class[test$class==4] <- "malignant"

maligVben = test$class
maligVben = as.factor(maligVben)
maligVben

str(pre5)
str(maligVben)
maligVben
pre5

confusionMatrix(pre5, maligVben)
confusionMatrix(pre25, maligVben)
confusionMatrix(pre50, maligVben)
confusionMatrix(pre75, maligVben)
confusionMatrix(pre100, maligVben)

confusionMatrix(pre7D, maligVben)

###   FINAL TREE   ###

goodTree= model7D$finalModel
prp(goodTree, extra= 4)

