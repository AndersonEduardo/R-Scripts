## pequeno tutorial de machine learning com pacote Caret (em R) com exemplo


## EXEMPLO 1 ##

library(caret)
library(kernlab)


##abrindo dados
data(spam)

##particionando dados
inTrain <- createDataPartition(y=spam$type, p=0.75, list=FALSE) #indices dos dados pra treino
training <- spam[inTrain,]
testing <- spam[-inTrain,]

##construindo modelos
set.seed(32343)

modelFit <- train(type ~., data=training, method="glm")
modelFit

modelFit$finalModel

##predicoes do modelo
predictions <-  predict(modelFit, newdata=testing)
predictions

##validando o modelo
confusionMatrix(predictions, testing$type)


## EXEMPLO 2 ##
## FONTE: https://statistik-dresden.de/archives/14967

library(MASS)

##abrindo dados
data(Boston)
str(Boston)

##particionando os dados
trainctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 10)

##construindo modelos
set.seed(2018)
rpart_tree <- train(medv ~ ., data = Boston, method = "rpart", trControl = trainctrl) #decision tree
rf_tree <- train(medv ~ ., data = Boston, method = "rf", trControl = trainctrl) #random forest

##comparando modelos
resamps <- resamples(list(rpart_tree = rpart_tree, randomForest = rf_tree))
summary(resamps)

dotplot(resamps, metric = "RMSE")

##fazendo GBM e refinando os parametros do modelo - so como um exemplo mesmo

##variando os parametros a serem testados
myGrid <- expand.grid(n.trees = c(150, 175, 200, 225),
                      interaction.depth = c(5, 6, 7, 8, 9),
                      shrinkage = c(0.075, 0.1, 0.125, 0.15, 0.2),
                      n.minobsinnode = c(7, 10, 12, 15))

##construindo modelos
gbm_tree_tune <- train(medv ~ ., data = Boston, method = "gbm", distribution = "gaussian",
                       trControl = trainctrl, verbose = FALSE,
                       tuneGrid = myGrid)
gbm_tree_tune

myGrid <- gbm_tree_tune$bestTune #modelo com melhor parametrizacao

gbm_tree <- train(medv ~ ., data = Boston, method = "gbm",
                  trControl = trainctrl,
                  tuneGrid = myGrid, verbose = FALSE) #calibrando GBM otimizado
gbm_tree

##comparando modelos
resamps <- resamples(list(rpart_tree = rpart_tree, randomForest = rf_tree, gbm_tree_auto = gbm_tree_auto,gbm_tree_user = gbm_tree))
summary(resamps)

dotplot(resamps, metric = "RMSE", main = "Modell-Vergleich")
