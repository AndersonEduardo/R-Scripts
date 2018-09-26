##Funcao para extrair valores das variaveis ambientais nas coordenadas de registros de ocorrencias com diferentes idades
##Anderson A. Eduardo
##23/Ago/2018

paleoextract = function(x, cols=names(x), path) {

    if( ncol(x) < 3 | class(x) != "data.frame" | sum(cols %in% c('lon','lat','age')) < 3 ){
        stop("O conjunto de dados de entrada deve ser um data.frame contendo (no mínimo) as respectivas \n \t colunas: lon (longitude), lat (latitude), age (idade)")
    }

    if( any( is.na( match(unique(x$age), as.integer(list.files(path))) ) ) ){
        warning("ATENÇÃO: as idades no conjunto de dados devem corresponder à pastas com variaveis ambientais para cada uma das idades. NAs produzidos. \n")
    }

    if("ID" %in% names(x)){
        stop("Desculpe, mas infelizmente o nome de coluna 'ID' não é permitido para esta funcao. Por favor, renomeie ou exclua do dataset de entrada desta funcao.")
    }

    ##variaveis locais
    currentDataSet = x
    path = path
    predictorsData = data.frame()
    age = unique(currentDataSet$age)
    idNA = vector()

    ##helper para tratar as idades que sao NA
    currentDataSet$ID = seq(nrow(currentDataSet))
    idNA = sapply( seq(nrow(currentDataSet)), function(x) currentDataSet[x,]$age  %in% as.numeric(list.files(path)) )
    ageNA = sapply( seq(length(age)), function(x) age[x]  %in% as.numeric(list.files(path)) )

    ##loop para coletar os dados ambientais dos pontos    
    infoData = do.call("rbind", lapply( X=age[ageNA], FUN=function(age_i)
        currentDataSet[grep(paste(age_i, collapse ='|'), currentDataSet$age), ] ) )
    
    predictorsData = do.call("rbind", lapply( X=age[ageNA], FUN=function(age_i)
        extract(x = stack(list.files(file.path(path, age_i), pattern='.asc', full.names=TRUE)),
                y = currentDataSet[grep(age_i, currentDataSet$age), c('lon','lat')]) ) )
    
    predictorsData = cbind(infoData, predictorsData)
  
    ## for ( i in seq(length(age)) ){
        
    ##     if ( ageNA[i] == FALSE ){ ##se nao houver pasta com dados ambientais da idade i, entao pular iteracao
            
    ##         next
            
    ##     }else{
            
    ##         currentPredictors = stack( list.files(file.path(path,  age[i]), pattern='.asc', full.names=TRUE) ) #abrindo dados
    ##         ##problema com meus dados
    ##         if("landmask" %in% names(currentPredictors)){
    ##             currentPredictors = dropLayer(currentPredictors, grep(pattern='landmask',x=names(currentPredictors)))
    ##         }
    ##         ##
    ##         crs(currentPredictors) = crs(raster()) #ajuste de projecao
    ##         currentVals = extract(x=currentPredictors, y=currentDataSet[which(currentDataSet$age == age[i]), c('lon','lat')]) #extraindo dados ambientais
    ##         predictorsData = rbind(predictorsData,
    ##                                data.frame(currentDataSet[which(currentDataSet$age == age[i]),], currentVals)) #dataset com dados ambientais
    ##     }
    ## }
    
    if ( sum(idNA==FALSE)>0 ){ #caso tenha ocorrido dados de ocorrencia sem dados ambientais para a idade deles

        envirDataCols = grep(paste(c(cols,"ID"), collapse="|"), names(predictorsData), value=TRUE, invert=TRUE) #colunas com variaveis ambientais
        matrixNA = matrix(NA, nrow=nrow(currentDataSet[!idNA,]), ncol=length(envirDataCols)) #matriz deNAs
        naData = cbind(currentDataSet[!idNA,], as.data.frame(matrixNA)) #dados associados a matriz de NAs
        names(naData) = names(predictorsData) #ajuste de nomes das colunas
        predictorsData = rbind(predictorsData, naData) #juntando NAs e dados obtidos ("nao-NAs")
    }

    if( nrow(predictorsData) > 1 ){
        predictorsData = predictorsData[order(predictorsData$ID),] #recuperando a organizacao original
    }
    
    predictorsData$ID = NULL #apagando ID helper

    ##output da funcao
    return(predictorsData)
    
}
