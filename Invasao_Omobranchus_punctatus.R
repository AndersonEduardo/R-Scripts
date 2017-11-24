## Script para analise do nicho de uma especie invasora, comparando nicho do ambiente nativo e do ambiente invadido.
## Anderson A. Eduardo
## novembro/2017


##abrindo pacotes necessarios
library(ecospat)
library(raster)
library(maps)

##definindo variaveis e parametros
projectFolder = "/home/anderson/Documentos/Projetos/Invasao_Omobranchus_punctatus" #pasta do projeto
#
envVarFolder = "???variaveis de oceano???" #pasta com as variaveis ambientais
predictors = stack(list.files(path=envVarFolder, full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis
predictors = predictors[[c('bioclim_01','bioclim_12')]] #selecionando as variaveis usadas
projection(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
#
spData = read.csv(file.path(projectFolder,'spOcc.csv'),header=TRUE) #dados de ocorrencia ambiente nativo
names(spData) = c('lon','lat')


## Separando dados de ocorrencia na area nativa e area invadida


wrld = rgdal::readOGR('/home/anderson/PosDoc/shapefiles/continents/continent.shp') #mapa mundi
plot(wrld) #criando imagem
points(spData,col='red') #plotando os pontos de ocorencia da sp
drawExtent() #delimitando area no mapa
natAreaExtent = extent(19.53604, 184.1086, -58.39874, 53.9756) #area distr. nativa

areaNat = crop(wrld,natAreaExtent) #recortando area selecionada  do mapa mundi
plot(wrld) #plotando mapa mundi
plot(areaNat,ad=TRUE,col='red') #verificando

sp::coordinates(spData) = ~lon+lat
crs(spData) = crs(areaNat)
#spData <- spTransform(spData, CRS(proj4string(areaNat))) # transform CRS

spNat = spData[areaNat, ]

plot(wrld) #plotando mapa mundi
plot(areaNat,ad=TRUE,col='red') #verificando
points(spNat, cex=2,col='blue')
points(spData, pch=20, cex=0.5, col='yellow')

## Analise de nicho

##The PCA is calibrated on all the sites of the study area
pca.env <- dudi.pca(rbind(spNatData,spInvData)[,c('bioclim_01','bioclim_12')],scannf=F,nf=2)
#ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) #grafico

##PCA scores for the whole study area
scores.globclim <- pca.env$li

##PCA scores for the species native distribution
scores.sp.nat <- suprow(pca.env,spNatData[which(spNatData[,'pres']==1),c('bioclim_01','bioclim_12')])$li

##PCA scores for the species invasive distribution
scores.sp.inv <- suprow(pca.env,spInvData[which(spInvData[,'pres']==1),c('bioclim_01','bioclim_12')])$li

##PCA scores for the whole native study area
scores.clim.nat <-suprow(pca.env,spNatData[,c('bioclim_01','bioclim_12')])$li

##PCA scores for the whole invaded study area
scores.clim.inv <- suprow(pca.env,spInvData[,c('bioclim_01','bioclim_12')])$li

##gridding the native niche
grid.clim.nat <-ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.nat, sp=scores.sp.nat, R=100, th.sp=0)

##gridding the invasive niche
grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.inv, sp=scores.sp.inv, R=100, th.sp=0)

##equivalencia de nicho
##OBS: Niche equivalency test H1: Is the overlap between the native and invaded niche higher than two random niches
eq.test <- ecospat.niche.equivalency.test(grid.clim.nat, grid.clim.inv, rep=10, alternative="greater")

Dobs = eq.test$obs$D #indice D observado
Iobs = eq.test$obs$I #indice I observado
DpValue = eq.test$p.D #p-valor indice D
IpValue = eq.test$p.I #p-valor indice I


