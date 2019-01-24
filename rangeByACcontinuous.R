##algoritmo para simulacao de um automato celular bastante simples
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


rangeByACcontinuous = function(envAreas){
  
  ##spRange = raster('/home/anderson/Documentos/Projetos/Improved pseudo-absences_TESTE/virtual species/sp2.asc')
  ##spRange = rasterSpDist
  spRange = spSuit = envAreas
  
  ##tranformando 0 em NA
  values(spSuit)[values(spSuit) == 0] = NA
  values(spRange)[ !is.na(values(spRange)) ] = 0
  
  ###trandformando 1 em 0 (i.e. preenchendo com 'ausencias' as areas climaticamente adequadas para a Sp)
  #values(spRange)[values(spRange) > 0] = 0
  
  ##semente (i.e. ponto de origem para o crescimento do range da Sp)
  seedCoords = randomPoints(mask=spSuit, n=1, prob=TRUE)
  
  ##crescendo
  spRange_i = SpatialPoints(seedCoords)
  crs(spRange_i) = CRS('+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84 ')
  spRange_i = rasterize(x = spRange_i, y = spRange)
  values(spRange_i)[ is.na(values(spRange_i)) ] = 0
  
  
  for (i in 1:50){
    
    ##raster do seed
    rangeBorder = boundaries(spRange_i, asNA=FALSE, classes=TRUE, type='outer')
    
    cellsBorder =  Which(rangeBorder == 1, cells = TRUE) 
    
    probOcc = runif(n = length(cellsBorder))
    
    borderOcc = sapply(seq(length(cellsBorder)), function(i) ifelse(spSuit[cellsBorder[i]] > probOcc[i], 1, 0) )
    
    rangeBorder[cellsBorder] = borderOcc
    
    rangeUp = sum(stack(spRange_i, rangeBorder), na.rm = TRUE)
    
    spRange_i = rangeUp
    
  }
    
    
    
    
    
    
    ##criando Buffer
    rangeBuff = rgeos::gBuffer(spRange_i, width=1) #capStyle="SQUARE"
    
    ##pegando celulas = zero dentro do Buffer e executando ocupacao usando suitability como probabilidade de occ.
    vals = extract(spSuit, rangeBuff, cellnumbers=TRUE)[[1]]
    vals2 = vals[which(vals[,2] > 0),]
    probOcc = runif(n=nrow(vals2))
    vals3 = vals2[vals2[,2] > probOcc, 1]
    spRange[vals3] = 1 #; plot(spRange)
    
    ##fazendo Buffer sobre Buffer
    #spRange_i = rangeBuff
    
    ##pontos para o novo Buffer
    ptsDistUpadatedRaw = extract(spRange, rangeBuff, cellnumbers=TRUE)[[1]]
    ptsDistUpadated = ptsDistUpadatedRaw[which(ptsDistUpadatedRaw[,2]==1),]
    ptsMigrationTotal = xyFromCell(spRange,ptsDistUpadated[,1])
    ptsMigrationSelected = ptsMigrationTotal[sample(x=nrow(ptsMigrationTotal), size=ceiling(0.01*nrow(ptsMigrationTotal))),]
    spRange_i = SpatialPoints(ptsMigrationSelected)
    
    #plot(spRange)
    
  }
  
  return(spRange)
  
}

