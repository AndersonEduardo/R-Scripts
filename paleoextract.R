##Funcao para extrair valores das variaveis ambientais nas coordenadas de registros de ocorrencias com diferentes idades
##Anderson A. Eduardo
##23/Ago/2018

paleoextract = function(x, y) {

    if( ncol(x)!=3 | class(x) != "data.frame" | any(names(x)!=c('lon','lat','age')) ){
        stop("O conjunto de dados de entrada deve ser um data.frame contendo as respectivas \n \t colunas: lon (longitude, lat (latitude), age (idade)")
    }

    if( any( is.na( match(unique(x$min), as.integer(list.files(y))) ) ) ){
        stop("As idades no conjunto de dados devem corresponder à pastas com variaveis ambientais para cada uma das idades.")
    }

    currentDataSet = x
    agesVector = unique(currentDataSet$age)
    predictorsData = data.frame()
    agesNA = vector()

    for (age_i in agesVector){

        if (!age_i %in% as.numeric(list.files(y))){
            agesNA = append(x=agesNA, values=age_i)
            warning("\n ATENÇÃO: alguns dos pontos de ocorrência não possuem dados ambientais. NAs produzidos. \n")
            next
        }
            
        currentPredictors = stack( list.files(file.path(y,age_i),pattern='asc', full.names=TRUE) )

        ##problema com meus dados
        if("landmask" %in% names(currentPredictors)){
            currentPredictors = dropLayer(currentPredictors, grep(pattern='landmask',x=names(currentPredictors)))
        }
        ##
        
        crs(currentPredictors) = crs(raster())
        currentVals = extract(x=currentPredictors, y=currentDataSet[ currentDataSet$age == age_i, c('lon','lat')])
        currentVals = data.frame(age=age_i,currentVals)
        predictorsData = rbind(predictorsData, currentVals)

    }

    if(length(agesNA) >0){
        dataNA = data.frame(age=agesNA, matrix(data=NA,nrow=length(agesNA),ncol=ncol(predictorsData)-1))
        names(dataNA) = names(predictorsData)
        predictorsData = rbind(predictorsData, dataNA)
        predictorsData = predictorsData[order(predictorsData$age),]
    }

    return(predictorsData)
    
}
