##Funcao para extrair valores das variaveis ambientais nas coordenadas de registros de ocorrencias com diferentes idades
##Anderson A. Eduardo
##23/Ago/2018

paleoextract = function(x, cols=names(x), path) {
  
  ##validacao dos dados de entrada
  if( ncol(x) < 3 | class(x) != "data.frame" | sum(cols %in% c('lon','lat','age')) < 3 ){
    stop("O conjunto de dados de entrada deve ser um data.frame contendo (no mínimo) as respectivas \n \t colunas: lon (longitude), lat (latitude), age (idade)")
  }
  
  if( any( is.na( match(unique(x$age), as.integer(list.files(path))) ) ) ){
    warning("ATENÇÃO: as idades no conjunto de dados devem corresponder à pastas com variaveis ambientais para cada uma das idades. NAs produzidos. \n")
  }
  
  if("ID" %in% names(x)){
    stop("Desculpe, mas infelizmente o nome de coluna 'ID' não é permitido para esta funcao. Por favor, renomeie ou exclua do dataset de entrada desta funcao.")
  }
  
  ##variaveis locais da funcao
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
  predictorsData = data.frame()
  
  for (age_i in age[ageNA]){
    
    currentLines = currentDataSet[which(age_i==currentDataSet$age), c('lon','lat','id','age','ID')]
    currentPredictors = stack(list.files(file.path(path, age_i), pattern='.asc', full.names=TRUE))
    varNames = basename(list.files(file.path(path, age_i), pattern='.asc', full.names=TRUE))
    varNames = gsub('\\..*', '', varNames)
    ##problema com meus dados
    if("landmask" %in% names(currentPredictors)){
      currentPredictors = dropLayer(currentPredictors, grep(pattern='landmask',x=names(currentPredictors)))
    }
    ##
    currentCoords =  currentLines[,c("lon","lat")] #currentDataSet[which(age_i==currentDataSet$age), c('lon','lat')]
    varValues = extract(x = currentPredictors, y = currentCoords)
    varValues = data.frame(varValues)
    names(varValues) = varNames
    predictorsData = rbind(predictorsData, cbind(currentLines, varValues))
  }
  
  ##caso tenha ocorrido dados de ocorrencia sem dados ambientais para a idade deles
  if ( sum(idNA==FALSE) > 0 ){ 
    
    envirDataCols = grep(paste(c(cols,"ID"), collapse="|"), names(predictorsData), value=TRUE, invert=TRUE) #colunas com variaveis ambientais
    matrixNA = matrix(NA, nrow=nrow(currentDataSet[!idNA,]), ncol=length(envirDataCols)) #matriz deNAs
    naData = cbind(currentDataSet[!idNA,], as.data.frame(matrixNA)) #dados associados a matriz de NAs
    names(naData) = names(predictorsData) #ajuste de nomes das colunas
    predictorsData = rbind(predictorsData, naData) #juntando NAs e dados obtidos ("nao-NAs")
  }

  ##recuperando a organizacao original
  if( nrow(predictorsData) > 1 ){
    predictorsData = predictorsData[order(predictorsData$ID),]
  }
  
  ##apagando ID helper
  predictorsData$ID = NULL 
  
  ##output da funcao
  return(predictorsData)
  
}
