
##abrindo pacores necessarios
library(raster)
library(hypervolume)
source('/home/anderson/R-Scripts/paleoextract.R')
source('/home/anderson/R-Scripts/strings2na.R')
source('/home/anderson/R-Scripts/dataInstance.R')
source('/home/anderson/R-Scripts/uplot.R')
source('/home/anderson/R-Scripts/uniche.R')


##definindo parametros e variaveis globais
projectFolder = "/home/anderson/Projetos/SDM megafauna Sul-Americana"
envFolder = "/home/anderson/gridfiles/dados_projeto"
AmSulBorders = rgdal::readOGR('/home/anderson/shapefiles/Am_Sul/borders.shp')


##abrindo e tratando o banco de dados

##arquivo do banco de dados
dataSetRaw = read.csv(file='/home/anderson/Projetos/SDM megafauna Sul-Americana/dataset_clean.csv', header=TRUE, dec='.', sep=',')

##subset do banco de dados
pts = dataSetRaw[,c('Species','Longitude','Latitude','Cal..Mean','Min.','Max.')]

sp_i = 'Catonyx cuvieri'

pts = pts[which(pts$Species == sp_i),]

##ajustando dados
pts = strings2na(pts, 'Species') #transformando strings ao longo do dateset em NA
pts = dataInstance(pts, c('Cal..Mean','Min.','Max.'), n=10) #criando instancias de dados (i.e. definindo idades a partir dos itervalos)
pts = cleanData(pts, c('Longitude','Latitude','age'), 2) #arredondando coords, limpando NAs e dados duplicados

## plot(AmSulBorders)
## points(pts[[1]][,c('Longitude','Latitude')], pch=20, cex=1.5, col=as.factor(pts[[1]]$sps))


##analise de incerteza e sensibilidade##


xxdat = lapply( seq(length(pts)), function(x) pts[[x]][,2:4] ) #data.frames somente com colunas longitude, latitude e idade (necessario pra funcao)

xx = uniche(xxdat) #analise de incerteza e sensibilidade
uplot(xx, AmSulBorders) #inspecao visual
