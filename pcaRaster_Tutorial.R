## script para implmentacao de PCA (analise de componentes principais) sobre raters
## Junho/2016

#libraries related to maps
library(raster)
library(RStoolbox)# rasterPCA function

##abrir variaveis preditoras usadas no SDM
load("/home/anderson/Projetos/Darlan - macroevolucao caviomorfa/Variaveis selecionadas/uncorrelatedPredictorsSouthAmerica.R")

plot(predictors) #inspacao grafica


##now we can do pca
pcamap = rasterPCA(predictors, spca=TRUE)

##check loadings and eigenvalues
knitr::kable(round(pcamap$model$loadings[,1:3],3)) # top 3 loadings

##visualaizando os PCs no espaco
plot(pcamap$map)

##salvando rasters dos PCs
setwd('/home/anderson/Projetos/Darlan - macroevolucao caviomorfa/PCA')
writeRaster(x=pcamap$map,filename=names(pcamap$map), bylayer=TRUE, format="raster")
