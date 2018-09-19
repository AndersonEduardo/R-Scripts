
##abrindo pacotes e funcoes necessarias##
library(raster)
library(hypervolume)
source('/home/anderson/R-Scripts/paleoextract.R')
source('/home/anderson/R-Scripts/strings2na.R')
source('/home/anderson/R-Scripts/dataInstance.R')
source('/home/anderson/R-Scripts/uplot.R')
source('/home/anderson/R-Scripts/uniche.R')
source('/home/anderson/R-Scripts/cleanData.R')


##definindo parametros e variaveis globais##
projectFolder = "/home/anderson/Projetos/SDM megafauna Sul-Americana" #pasta de trabalho do projeto
envFolder = "/home/anderson/gridfiles/dados_projeto" #pasta das variaveis ambientais
AmSulBorders = rgdal::readOGR('/home/anderson/shapefiles/Am_Sul/borders.shp') #shapefile da Am. do Sul


##abrindo e tratando o banco de dados##

##arquivo do banco de dados
dataSetRaw = read.csv(file='/home/anderson/Projetos/SDM megafauna Sul-Americana/dataset_clean.csv', header=TRUE, dec='.', sep=',')

##subset do banco de dados
pts = dataSetRaw[,c('Species','Longitude','Latitude','Cal..Mean','Min.','Max.')]
sp_i = 'Catonyx cuvieri' #'Notiomastodon platensis'
pts = pts[which(pts$Species == sp_i),]

##ajustando dados
pts = strings2na(pts, 'Species') #transformando strings ao longo do dateset em NA
pts = dataInstance(pts, c('Cal..Mean','Min.','Max.'), n=2*nrow(pts)) #criando instancias de dados (i.e. definindo idades a partir dos itervalos)
pts = cleanData(pts, c('Longitude','Latitude','age'), p=2) #arredondando coords, limpando NAs e dados duplicados
pts = lapply( X=sample(x=seq(length(pts)), size=2*nrow(pts[[1]])), FUN=function(x) pts[[x]] ) #ajustando tamanho baseado no numero de pontos reais

##inspeção visual dos dados
plot(AmSulBorders)
points(pts[[1]][,c('lon','lat')], pch=20, cex=1.5, col='red')

##analise de incerteza e sensibilidade##
analise.incerteza = uniche(pts, envFolder=envFolder) #analise de incerteza e sensibilidade

uplot(analise.incerteza, AmSulBorders, legend=FALSE) #inspecao visual - opcao 1

plot(AmSulBorders) #inspecao visual - opcao 2
uplot(x=xx, shape=NULL, niche.metric='marginality')



#### rodando para todas as espécies do banco de dados - somente um teste ainda ####



##abrindo e tratando o banco de dados##

##arquivo do banco de dados
dataSetRaw = read.csv(file='/home/anderson/Projetos/SDM megafauna Sul-Americana/dataset_clean.csv', header=TRUE, dec='.', sep=',')

##subset do banco de dados
dataSetRaw = dataSetRaw[,c('Species','Longitude','Latitude','Cal..Mean','Min.','Max.')]

##lista dos nomes das especies
sps = as.character( unique(dataSetRaw$Species) )

for (i in seq(length(sps))){

    ##dataset da especie atual
    pts = dataSetRaw[which(dataSetRaw$Species == sps[i]),]

    if (nrow(pts) <= 5 | nrow(pts) > 50){
        cat("\n ATENÇÃO: A espécie", sps[i], "possui menos que 5 registros, por isso não foi analisada. \n")
        next
    }else{
        cat("\n Rodando para a espécie", sps[i], "...\n")
    }

    ##ajustando dados
    pts = strings2na(pts, 'Species') #transformando strings ao longo do dateset em NA
    pts = dataInstance(pts, c('Cal..Mean','Min.','Max.'), n=2*nrow(pts)) #criando instancias de dados (i.e. definindo idades a partir dos itervalos)
    pts = cleanData(pts, c('Longitude','Latitude','age'), p=2) #arredondando coords, limpando NAs e dados duplicados

    ## ##inspeção visual dos dados
    ## plot(AmSulBorders)
    ## points(pts[[1]][,c('lon','lat')], pch=20, cex=1.5, col='red')

    ##analise de incerteza e sensibilidade
    analise.incerteza = uniche(x=pts, envFolder=envFolder) #analise de incerteza e sensibilidade

    ##salvando output da analise
    save(analise.incerteza, file=paste(projectFolder,'/Analise de sensibilidade e incerteza/',sps[i],'.R',sep=''))

    ##salvando output grafico
    jpeg(paste(projectFolder,'/Analise de sensibilidade e incerteza/',sps[i],'.jpeg',sep=''), 600, 600)
    uplot(analise.incerteza, AmSulBorders, legend=FALSE)
    dev.off()

    gc()

}


## ###testes


## for (i in seq(length(currentDataSet))){
##     cat(" Rodando iteração: ", i, "\n")
##     teste = paleoextract(currentDataSet[[i]], cols=c('lon','lat','age'), path=envFolder)
## }





