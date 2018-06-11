## Script para realizacao de DTT (dissmilarity through time)
## 06/junho/2018

##install.packages('geiger')
library(geiger)


##abrindo dados
data(geospiza)
attach(geospiza)

##DTT
disp.data = dtt(geospiza.tree, geospiza.data) #Evaluates disparity-through-time for either a single data set or multiple data sets



