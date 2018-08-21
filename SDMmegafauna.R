
##abrindo pacores necessarios
library(raster)


##definindo parametros e variaveis globais
projectFolder = "/home/anderson/Projetos/SDM megafauna Sul-Americana"
envFolder = "/home/anderson/gridfiles/dados_projeto"
AmSulBorders = rgdal::readOGR('/home/anderson/shapefiles/Am_Sul/borders.shp')


##abrindo e tratando o banco de dados

##arquivo do banco de dados
dataSetRaw = read.csv(file='/home/anderson/Projetos/SDM megafauna Sul-Americana/dataset_clean.csv',header=TRUE,dec='.',sep=',')

##subset do banco de dados
pts = dataSetRaw[,c('Species','Longitude','Latitude','Cal..Mean','Min.','Max.')]

##ajustando coordenadas
ptsCoords = apply(pts[,c('Longitude','Latitude')],1, as.numeric) #transforma informacao de texto em NA (ex.: gruta azul -> NA)
ptsCoords = data.frame(t(ptsCoords)) #transformando em data.frame
pts = data.frame(sp=pts$Species, lon=ptsCoords[,1], lat=ptsCoords[,2], mean=pts$Cal..Mean, min=pts$Min, max=pts$Max) #consolidando os dados

##ajustando as idades
ptsAge = apply(pts[,c("mean","min","max")], 1, as.numeric) #transforma informacao de texto em NA (ex.: pleistocene -> NA)
ptsAge = data.frame(t(ptsAge)) #consolidando os dados
ptsAgeMeanNA = apply( ptsAge, 1, function(x) ifelse(is.na(x[1]), mean(x[2:3]), x[1]) ) #qdo media=NA, obtem a partir do intervalo (max e min)
ptsAgeMinNA = apply( ptsAge, 1, function(x) ifelse(is.na(x[2]) & is.na(x[3]), x[1], x[2]) ) #qdo min=NA, min=mean (i.e. data altamente precisa)
ptsAgeMaxNA = apply( ptsAge, 1, function(x) ifelse(is.na(x[2]) & is.na(x[3]), x[1], x[3]) ) #qdo max=NA, max=mean (i.e. data altamente precisa)
ptsAge = data.frame(cbind(ptsAgeMeanNA,ptsAgeMinNA,ptsAgeMaxNA)) #consolidando os dados de idade

##consolidando dos dados
dataSet = data.frame( sp=pts$sp, lon=pts$lon, lat=pts$lat, mean=ptsAge$ptsAgeMeanNA, min=ptsAge$ptsAgeMinNA, max=ptsAge$ptsAgeMaxNA )
dataSet[,c('lat','lon')] = round(dataSet[,c('lat','lon')], 2) #arredondando para duas casas decimais
dataSet = dataSet[complete.cases(dataSet),] #retirando dados incompletos
dataSet = unique(dataSet) #retirando dados repetidos

## plot(AmSulBorders)
## points(pts[,c('lon','lat')], pch=20, cex=1.5, col=as.factor(pts$sp))
## table(pts$sp)


##funcao para testar sensibilidade do hipervolume do nicho ao intervalo de erro nas idades dos registros fosseis


sp_i = dataSet$sp[2]

currentDataSet = dataSet[grep(pattern=sp_i, x=dataSet$sp), ]
currentDataSet$age = apply( currentDataSet[,c('min','max')], 1, function(x) runif(1, min=x[1], max=x[2]) )
currentDataSet$age = round(currentDataSet$age/1000)
currentDataSet = currentDataSet[!base::duplicated(currentDataSet[,c('lon','lat','age')]), ]

paleoextract = function(x, y) {

    if( ncol(x)!=3 | class(x) != "data.frame" | any(names(x)!=c('lon','lat','age')) ){
        stop("O conjunto de dados de entrada deve ser um data.frame contendo as respectivas \n \t colunas: lon (longitude, lat (latitude), age (idade)")
    }

    if( any( is.na( match(unique(x$min), as.integer(list.files(y))) ) ) ){
        stop("As idades no conjunto de dados devem corresponder à pastas com variaveis ambientais para cada uma das idades.")
    }

    currentDataSet = x

    agesVector = sort( unique(currentDataSet$age) )
    predictorsData = data.frame()

    for (age_i in agesVector){

        currentPredictors = stack( list.files(file.path(y,i),pattern='asc', full.names=TRUE) )
        crs(currentPredictors) = crs(raster())
        currentVals = extract(x=currentPredictors, y=currentDataSet[ currentDataSet$age == age_i, c('lon','lat')])
        currentVals = data.frame(currentVals)

        predictorsData = rbind(predictorsData, currentVals)

    
        
}



