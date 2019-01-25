## Funcao (simples) para criacao de especie virtual
## a partir de um conjunto das variaveis preditoras
## bioclim 01 e bioclim 12 (construcao de nicho climatico)
## Anderson A. eduardo
## 30/fev/2018

require(makeSpeciesSuitabilityIterator)

makeSpeciesSuitability = function(predictorsData){
  
  ##definindo parametros e variaveis locais
  ##matriz##
  if (is.list(predictorsData)){
    datMat = lapply( seq(length(predictorsData)), function(x) cbind(as.data.frame(predictorsData[[x]], xy=TRUE, na.rm=TRUE), fSp=0) ) #transformando raster em data.frame cada item da lista
  }else{
    datMat = cbind(as.data.frame(predictorsData, xy=TRUE, na.rm=TRUE), fSp=0) #transformando raster em data.frame  
  }
  
  #datMat = data.frame(datMat, fSp=0)
  #names(datMat) = c('lon', 'lat', 'bio1', 'bio12', paste('sp',i,sep='')) #ajustando os nomes das colunas do data.frame
  
  ## x = 1:100
  ## a=1 #altura do pico
  ## b=10 #posicao do centro
  ## c=1 #largura
  ## #
  ## fx = exp(-((x-b)^2/(2*c^2)))
  ## plot(fx~x,ylim=c(0,1.2))
  
  if ( is.list(datMat) ){
    spsDist = lapply( seq(length(datMat)), function(x) iterateFunc(x = datMat[[x]]) )
  }else{
    spsDist = iterateFunc(datMat)
  }
  
  ##mapa da distribuicao (presenca/ausencia) com automato celular
  ##source("/home/anderson/R/R-Scripts/rangeByAC.R")
  ##source("J:/Pesquisadorxs/Anderson_Eduardo/high_quality_PA/rangeByAC.R")
  ##source('D:/Anderson_Eduardo/rangeByAC.R')
  
  ##rasterSpDistributionAC = rangeByACcontinuous(rasterSpDistribution)     
  
  ## output da funcao
  return(spsDist)
  
}
