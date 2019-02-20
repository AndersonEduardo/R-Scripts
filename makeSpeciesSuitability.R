## Funcao (simples) para criacao de especie virtual
## a partir de um conjunto das variaveis preditoras
## bioclim 01 e bioclim 12 (construcao de nicho climatico)
## Anderson A. eduardo
## 30/fev/2018


makeSpeciesSuitability = function(predictorsData){
  
  ##definindo parametros e variaveis locais
  if (is.list(predictorsData)){
    #matriz
    datMat = lapply( seq(length(predictorsData)), function(i) cbind(as.data.frame(predictorsData[[i]], xy=TRUE, na.rm=TRUE), fSp=0) ) #transformando raster em data.frame cada item da lista
    datMat = lapply(datMat, setNames, c('lon','lat','bio1','bio12','fSp'))
    #df para funcao de resposta
    datMatRespFunc = datMat[[length(datMat)]]
    #area de fundo de referencia
    bgArea = predictorsData[[1]]*0
  }else{
    #matriz
    datMat = list(cbind(as.data.frame(predictorsData, xy=TRUE, na.rm=TRUE), fSp=0)) #transformando raster em data.frame  
    datMat = lapply(datMat, setNames, c('lon','lat','bio1','bio12','fSp'))
    #df para funcao de resposta
    datMatRespFunc = datMat[[1]]
    #area de fundo de referencia
    bgArea = predictorsData[[1]]*0
  }
  
  betaBio1 = runif(n=1, min=0.001, max=1)*sample(x=c(-1,1), size=1) #parametro para cada equacao de cada especie
  betaBio12 = runif(n=1, min=0.001, max=1)*sample(x=c(-1,1), size=1) #parametro para cada equacao de cada especie
  ##
  alphaBio1 = runif(n=1, min=quantile(datMatRespFunc$bio1, probs=0.25, na.rm=TRUE), max=quantile(datMatRespFunc$bio1, probs=0.75, na.rm=TRUE)) #parametro para cada equacao de cada especie
  alphaBio12 = runif(n=1, min=quantile(datMatRespFunc$bio12, probs=0.25, na.rm=TRUE), max=quantile(datMatRespFunc$bio12, probs=0.75, na.rm=TRUE)) #parametro para cada equacao de cada especie
  
  
  #criando distribuicao
  spsDist = lapply(seq(length(datMat)), function(i){
    datMatCurrent = datMat[[i]] #dados
    varBio1 = datMatCurrent$bio1 #variavel ambiental bioclim01
    varBio12 = datMatCurrent$bio12 #variavel ambiental bioclim12
    fBio1Sp_i = 1/(1+exp(betaBio1*(varBio1-alphaBio1))) #solucao da equacao com output binario ("suitability")
    fBio12Sp_i = 1/(1+exp(-betaBio12*(varBio12-alphaBio12))) #solucao da equacao com output binario ("suitability")
    datMatCurrent[,'fSp'] = fBio1Sp_i*fBio12Sp_i ##*fElevSp_i
    SpDistribution = datMatCurrent[,c('lon','lat','fSp')] #extraindo lon/lat e suitability (ou pres-aus) de cada especie
    coordinates(SpDistribution) = ~lon+lat #definindo colunas das coordenadas
    gridded(SpDistribution) = TRUE #definindo gridded
    proj4string(SpDistribution) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0' #definindo CRS
    SpDistribution = raster(SpDistribution) #rasterizing
    SpDistribution = merge(SpDistribution, bgArea) #quality control of resolution and other stuff
    
    return(SpDistribution)
  })
  
  ## output da funcao
  if(length(spsDist)==1){spsDist=spsDist[[1]]}
  
  return(spsDist)
  
}
