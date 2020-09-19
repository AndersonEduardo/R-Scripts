##algoritmo para simulacao da distribuicao geografica de uma especie usando um automato celular simples
##Anderson A. Eduardo
##26/jan/2018

library(raster)
library(dismo)

# rangeByACcontinuous = function(envAreas){
#     
#     ##spRange = raster('/home/anderson/Documentos/Projetos/Improved pseudo-absences_TESTE/virtual species/sp2.asc')
#     ##spRange = rasterSpDist
#     spRange = spSuit = envAreas
#     
#     ##tranformando 0 em NA
#     values(spSuit)[values(spSuit) == 0] = NA
#     values(spRange)[values(spRange) == 0] = NA
#     
#     ##trandformando 1 em 0 (i.e. preenchendo com 'ausencias' as areas climaticamente adequadas para a Sp)
#     values(spRange)[values(spRange) > 0] = 0
#     
#     ##semente (i.e. ponto de origem para o crescimento do range da Sp)
#     seedCoords = randomPoints(mask=spSuit, n=1, prob=TRUE)
#     
#     ##crescendo
#     spRange_i = SpatialPoints(seedCoords)
#     crs(spRange_i) = CRS('+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84 ')
#     
#     for (i in 1:10){
#         
#         ##criando Buffer
#         rangeBuff = rgeos::gBuffer(spRange_i, width=1) #capStyle="SQUARE"
#         
#         ##pegando celulas = zero dentro do Buffer e executando ocupacao usando suitability como probabilidade de occ.
#         vals = extract(spSuit, rangeBuff, cellnumbers=TRUE)[[1]]
#         vals2 = vals[which(vals[,2] > 0),]
#         probOcc = runif(n=nrow(vals2))
#         vals3 = vals2[vals2[,2] > probOcc, 1]
#         spRange[vals3] = 1 #; plot(spRange)
#         
#         ##fazendo Buffer sobre Buffer
#         #spRange_i = rangeBuff
#         
#         ##pontos para o novo Buffer
#         ptsDistUpadatedRaw = extract(spRange, rangeBuff, cellnumbers=TRUE)[[1]]
#         ptsDistUpadated = ptsDistUpadatedRaw[which(ptsDistUpadatedRaw[,2]==1),]
#         ptsMigrationTotal = xyFromCell(spRange,ptsDistUpadated[,1])
#         ptsMigrationSelected = ptsMigrationTotal[sample(x=nrow(ptsMigrationTotal), size=ceiling(0.01*nrow(ptsMigrationTotal))),]
#         spRange_i = SpatialPoints(ptsMigrationSelected)
#         
#         #plot(spRange)
#         
#     }
#     
#     return(spRange)
#     
# }


rangeByAC = function(envAreas, movRes, iter=100){
  
  ## atualizando distribuicoes atraves do tempo ##
  # envAreas: suitability gridfile for sps range
  # movRes: raster layer with data for resistence to species moviment throughout the landscape (migration)
  # iter: number of iterations to be performed
  
  ## validation of input data
  if(!compareRaster(envAreas, movRes, stopiffalse=FALSE)){
    stop('Os rasters de entrada precisam ter mesma resolucao e origem.')
  }
  
  #defining working gridfiles
  spRange = spSuit = envAreas
  movRes = movRes
  
  ##tranformando 0 em NA
  values(spSuit)[values(spSuit) == 0] = NA
  
  ##transformando tudo o que nao for NA em zero
  values(spRange)[ !is.na(values(spRange)) ] = 0
  
  ##semente (i.e. ponto de origem para o crescimento do range da Sp)
  seedCoords = randomPoints(mask=spSuit, n=1, ext = NULL, prob=TRUE)
  # extSeed = c(seedCoords[,1]-2, 
  #             seedCoords[,1]+2, 
  #             seedCoords[,2]-2, 
  #             seedCoords[,2]+2)
  # seedCoords = randomPoints(mask=spSuit, n=10, ext = extSeed, prob=TRUE)
  
  ##growing up sps range
  spRange_i = SpatialPoints(seedCoords)
  crs(spRange_i) = CRS('+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84')
  spRange_i = rasterize(x = spRange_i, y = spRange)
  values(spRange_i)[ is.na(values(spRange_i)) ] = 0
  
  #iterations
  for (i in 1:iter){
    
    occupiedCells = Which(spRange_i == 1, cells = TRUE) #identifying occupied cells
    
    rangeBorder = boundaries(spRange_i, asNA=FALSE, classes=TRUE, type='outer', directions=sample(c(4,8),1)) #get the borders os sps range
    
    cellsBorder =  Which(rangeBorder == 1, cells = TRUE) #identifying boundary cells
  
    probOcc = runif(n = length(cellsBorder)) #prob of occupancy
    
    borderOcc = as.integer(spSuit[cellsBorder] > probOcc) #probabilistic occupancy of borders
    #borderOcc = sapply(seq(length(cellsBorder)), function(i) ifelse(spSuit[cellsBorder[i]] > probOcc[i], 1, 0) ) #probabilistic occupancy of borders
    
    probMigr = runif(n = length(cellsBorder)) #prob of migration to cell i
    
    borderMig = as.integer(movRes[cellsBorder] > probMigr) #probabilistic migration to the borders
    
    borderOcc = borderOcc*borderMig #occupancy of borders accounting for environmental niche AND environmental resistance for species moviment
    
    #updating sps range - filling up cells
    spRange_i[c(occupiedCells, cellsBorder)] = c(rep(1,length(occupiedCells)), borderOcc) 
    
  }
  
  #adjusting final gridfile
  spRange = spRange_i
  maskHelper = spSuit
  maskHelper[values(maskHelper)>=0.1] = 1
  maskHelper[values(maskHelper)<0.1] = NA
  spRange = spRange*maskHelper
  
  #output - parabin - parabem - paraboom!
  return(spRange)
  
}