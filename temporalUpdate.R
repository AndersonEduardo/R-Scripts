## funcao para atualizar o gridfile do range de uma especie usando um novo mapa de suitability
## Anderson A. Eduardo
## 25/jan/2019

temporalUpdate = function(currentSpsRange, newSuitabilityMap, iter = 100){
  
  ## atualizando distribuicoes atraves do tempo ##
  # currentSpsRange: current binary distribution (sps range)
  # newSuitabilityMap: new suitability gridfile for reprojection of sps range
  # iter: number of iterations to be performed
  
  #defining initial parameters
  #suitT1 = rasterSpDistribution
  #NewSuit = rasterSpDistribution
  spsRange = currentSpsRange
  NewSuit = newSuitabilityMap
  
  SpDistDF = raster::as.data.frame(spsRange, xy=TRUE, na.rm=FALSE)
  NewSuitDF = raster::as.data.frame(NewSuit, xy=TRUE, na.rm=FALSE)
  
  occProb = runif(nrow(SpDistDF))
  
  #reduzindo a distribuicao usando o suitability novo (obs.: aqui eh somente a parte de reducao da distrib.)
  SpDistDF$DistUpdated = sapply(seq(nrow(SpDistDF)), function(i) ifelse(NewSuitDF[i,3] > occProb[i], SpDistDF[i,3], 0))
  
  #rastering the reduced range
  spRange = SpDistDF[,c('x','y','DistUpdated')] #extraindo lon/lat e suitability (ou pres-aus) de cada especie
  coordinates(spRange) = ~x+y #definindo colunas das coordenadas
  gridded(spRange) = TRUE #definindo gridded
  proj4string(spRange) = '+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84' #definindo proj
  spRangeUpdated = raster(spRange) #criando objeto raster
  #plot(spRangeUpdated)
  
  #filling NAs with zeros (help algorithm for updating sps range)
  values(spRangeUpdated)[ is.na(values(spRangeUpdated)) ] = 0
  
  #stating object for iterations
  spRange_i = spRangeUpdated
  
  #iterations
  for (i in seq(iter)){
    
    occupiedCells = Which(spRange_i == 1, cells = TRUE) #identifying occupied cells
    
    rangeBorder = boundaries(spRange_i, asNA=FALSE, classes=TRUE, type='outer', directions=sample(c(4,8),1)) #get the borders of sps range
    
    cellsBorder =  Which(rangeBorder == 1, cells = TRUE) #identifying boundary cells
    
    probOcc = runif(n = length(cellsBorder)) #prob of occupancy
    
    borderOcc = sapply(seq(length(cellsBorder)), function(i) ifelse(NewSuit[cellsBorder[i]] > probOcc[i], 1, 0) ) #probabilistic occupancy of borders
    
    spRange_i[c(occupiedCells, cellsBorder)] = c(rep(1,length(occupiedCells)), borderOcc) #updating sps range - filling cells
    
  }
  
  #adjusting final gridfile
  spRangeUpdated = spRange_i
  maskHelper = NewSuit
  maskHelper[values(maskHelper)>0] = 1
  spRangeUpdated = spRangeUpdated*maskHelper
  
  #output - parabin - parabem - paraboom!
  return(spRangeUpdated)
  
}