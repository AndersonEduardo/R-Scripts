##algoritmo para simulacao de um automato celular bastante simples
##Anderson A. Eduardo
##26/jan/2018

library(raster)
library(dismo)

rangeByAC = function(envAreas){

    #spRange = raster('/home/anderson/Documentos/Projetos/Improved pseudo-absences_TESTE/virtual species/sp2.asc')
    spRange = rasterSpDist
    #spRange = envAreas

    ##tranformando 0 em NA
    values(spRange)[values(spRange) == 0] = NA

    ##trandformando 1 em 0 (i.e. preenchendo com 'ausencias' as areas climaticamente adequadas para a Sp)
    values(spRange)[values(spRange) == 1] = 0

    ##semente (i.e. ponto de origem para o crescimento do range da Sp)
    seedCoords = dismo::randomPoints(mask=spRange, n=1)
    
    ##crescendo
    spRange_i = SpatialPoints(seedCoords)
    
    for (i in 1:20){

        ##criando Buffer
        rangeBuff = rgeos::gBuffer(spRange_i, width=1) #capStyle="SQUARE"
        
        ##pegando celulas = zero dentro do Buffer e igualando a 1
        vals = extract(spRange, rangeBuff, cellnumbers=TRUE)[[1]]
        vals2 = vals[which(vals[,2]==0),]
        spRange[vals2[,1]] = 1 #; plot(spRange)

        ##fazendo Buffer sobre Buffer
        #spRange_i = rangeBuff
        
        ##fazendo Buffer ponto a ponto
        pts = xyFromCell(spRange,vals2[,1])
        pts = pts[sample(x=nrow(pts), size=20),]
        spRange_i = SpatialPoints(pts)

    }
    
    return(spRange)
}
