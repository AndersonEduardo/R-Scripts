library(raster)
library(dismo)

rangeByAC = function(envAreas){

    ##spRange = raster('/home/anderson/Documentos/Projetos/Improved pseudo-absences_TESTE/virtual species/sp2.asc')
    spRange = envAreas

    ##tranformando 0 em NA
    values(spRange)[values(spRange) == 0] = NA

    ##trandformando 1 em 0 (i.e. preenchendo com 'ausencias' as areas climaticamente adequadas para a Sp)
    values(spRange)[values(spRange) == 1] = 0

    ##semente (i.e. ponto de origem para o crescimento do range da Sp)
    seedCoords = dismo::randomPoints(mask=spRange, n=1)
    
    ##crescendo
    spRange_i = SpatialPoints(seedCoords)
    
    for (i in 1:30){
        rangeBuff = rgeos::gBuffer(spRange_i, capStyle="SQUARE")
        ##
        vals = extract(spRange,rangeBuff,cellnumbers=TRUE)[[1]]
        vals2 = vals[which(vals[,2]==0),]
        ##
        spRange[vals2[,1]] = 1 ; plot(spRange)
        ##
        pts = xyFromCell(spRange,vals2[,1])
        ##
        spRange_i = SpatialPolygons(list(
            Polygons(list(
                Polygon(pts)),ID=1
                ))
            )
    }

    return(spRange)

}
