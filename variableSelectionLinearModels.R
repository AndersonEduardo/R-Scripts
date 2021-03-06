# Stepwise Regression
library(MASS)
#fit <- lm(y~x1+x2+x3,data=dataSet)
fullModel = glm(pres ~ bioclim_10 +  I(bioclim_10^2) + bioclim_11 +  I(bioclim_11^2) + bioclim_15 +  I(bioclim_15^2) + bioclim_16 + I(bioclim_16^2) + bioclim_17 +  I(bioclim_17^2), family=binomial(link=logit), data=dataSet)
#
manualBest = glm(pres ~ bioclim_10 + bioclim_11 + bioclim_15 + bioclim_16 + I(bioclim_16^2) + bioclim_17, family=binomial(link=logit), data=dataSet)
#
RbestModel = glm(pres ~ bioclim_11 + bioclim_17 + I(bioclim_17^2) + bioclim_10, family=binomial(link=logit), data=dataSet)
    
stepWise <- stepAIC(fullModel, direction="both")
stepWise$anova # display results

# All Subsets Regression
library(leaps)
attach(mydata)
leaps<-regsubsets(y~x1+x2+x3+x4,data=mydata,nbest=10)
leaps<-regsubsets(pres ~ bioclim_10 +  I(bioclim_10^2) + bioclim_11 +  I(bioclim_11^2) + bioclim_15 +  I(bioclim_15^2) + bioclim_16 + I(bioclim_16^2) + bioclim_17 +  I(bioclim_17^2),data=dataSet,nbest=5)
# view results
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps,scale="r2")
# plot statistic by subset size
library(car)
subsets(leaps, statistic="rsq") 
