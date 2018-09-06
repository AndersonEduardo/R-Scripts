
##abrindo pacores necessarios
library(raster)
library(hypervolume)
source('/home/anderson/R-Scripts/paleoextract.R')
source('/home/anderson/R-Scripts/strings2na.R')
source('/home/anderson/R-Scripts/dataInstance.R')


##definindo parametros e variaveis globais
projectFolder = "/home/anderson/Projetos/SDM megafauna Sul-Americana"
envFolder = "/home/anderson/gridfiles/dados_projeto"
AmSulBorders = rgdal::readOGR('/home/anderson/shapefiles/Am_Sul/borders.shp')


##abrindo e tratando o banco de dados

##arquivo do banco de dados
dataSetRaw = read.csv(file='/home/anderson/Projetos/SDM megafauna Sul-Americana/dataset_clean.csv',header=TRUE,dec='.',sep=',')

##subset do banco de dados
pts = dataSetRaw[,c('Species','Longitude','Latitude','Cal..Mean','Min.','Max.')]

sp_i = 'Catonyx cuvieri'

pts = pts[which(pts$Species == sp_i),]

##ajustando dados
pts = strings2na(pts, 'Species') #transformando strings ao longo do dateset em NA
pts = dataInstance(pts, c('Cal..Mean','Min.','Max.'), n=10) #criando instancias de dados (i.e. definindo idades a partir dos itervalos)

##consolidando dos dados
##pts = as.data.frame(pts[[1]]) #pegando somente a primeira instancia dos dados
pts = lapply(seq(length(pts)), function(x) data.frame(sps=pts[[x]]$Species, round(pts[[x]][,c('Longitude','Latitude')], 2), age=pts[[x]]$age )) #arredondando para duas casas decimais
pts = lapply( seq(length(pts)), function(x) pts[[x]][complete.cases(pts[[x]]),] ) #retirando dados incompletos
pts = lapply( seq(length(pts)), function(x) unique(pts[[x]]) ) #retirando dados incompletos

plot(AmSulBorders)
points(pts[[1]][,c('Longitude','Latitude')], pch=20, cex=1.5, col=as.factor(pts[[1]]$sps))


##funcao para testar sensibilidade do hipervolume do nicho ao intervalo de erro nas idades dos registros fosseis



dataSetList = list()

currentDataSet = pts
currentDataSet = lapply( seq(length(currentDataSet)), function(x) data.frame(sps=currentDataSet[[x]]$sps, lon=currentDataSet[[x]]$Longitude, lat=currentDataSet[[x]]$Latitude, age=currentDataSet[[x]]$age) )

dataSet_i =  paleoextract( currentDataSet[[1]][c('lon','lat','age')], envFolder )

dataSet_i = dataSet_i[complete.cases(dataSet_i), ]

dataSetList[[1]] = dataSet_i

##global hypervolume
globaldataSet = do.call("rbind", dataSetList)
globaldataSet = unique(globaldataSet)
globalPC = prcomp(globaldataSet[,-1], center=TRUE, scale.=TRUE)
idx = which(summary(globalPC)$importance['Cumulative Proportion',] >= 0.95)[1]
globalHypvol = hypervolume_gaussian(globalPC$x[,1:idx])

PC = prcomp(dataSet_i[,-1], center=TRUE, scale.=TRUE)
idx = which(summary(PC)$importance['Cumulative Proportion',] >= 0.95)[1]
hypvol = hypervolume_gaussian(PC$x[,1:idx])

marginality = hypervolume_distance(globalHypvol, hypvol)
propSize =  get_volume(hypvol) / get_volume(globalHypvol)

outputData = rbind( outputData,
                   data.frame(marginality = marginality, propSize=propSize, matrix(dataSet_i$age,1,length(dataSet_i$age))) )

names(outputData) = c('marginality','propSize', paste('point_',seq(length(dataSet_i$age)),sep=''))
rownames(outputData) = NULL


#### PROXIMO PASSO: 1o) voltar la no comeco e dar um jeito para que as informacoes de cada ponto (ID, coordenadas e idade) possam ser recuperadas no output final; 2o) implementar o teste de sensibilidade do nicho (ou hipervolume) a incerteza em cada um dos pontos. ####




