## Script para analise do nicho de uma especie invasora, comparando nicho do ambiente nativo e do ambiente invadido.
## Anderson A. Eduardo
## novembro/2017




## PARTE 1: carregando pacotes, parâmetros e variaveis




## limpando a memória (qdo necessario)
rm(list=ls())

## abrindo pacotes necessarios
library(ecospat)
library(rgdal)
library(raster)
library(ENMeval)
library(usdm)
library(iSDM)
library(pROC)
library(rJava)
Sys.setenv(JAVA_HOME = "/usr/lib/jvm/java-7-openjdk-amd64")
#Windows#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_91') # for 64-bit version
options(java.parameters = "Xmx7g")


## definindo variaveis e parametros
projectFolder = "/home/anderson/Projetos/Invasao_Omobranchus_punctatus" #pasta do projeto
envVarFolder = "/home/anderson/gridfiles/Bio-ORACLE" #pasta com as variaveis ambientais
predictors = stack(list.files(path=envVarFolder, full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis
spData = read.csv(file.path(projectFolder,'spOcc.csv'),header=TRUE) #dados de ocorrencia ambiente nativo
names(spData) = c('lon','lat')
wrld = readOGR('/home/anderson/shapefiles/ne_50m_ocean/ne_50m_ocean.shp') #mapa mundi




## PARTE 2: definindo dados da area nativa e area invadida




## separando area nativa e area invadida
plot(wrld) #criando imagem
points(spData,col='red') #plotando os pontos de ocorencia da sp
drawExtent() #delimitando area no mapa
natAreaExtent = extent(19.53604, 184.1086, -58.39874, 53.9756) #area distr. nativa
invAreaExtent = extent(-92.2889 , 18.83274, -88.59934, 24.47734) #area distr. invadida

## verificando area invadida e nativa
areaNat = crop(wrld,natAreaExtent) #recortando area nativa
areaInv = crop(wrld,invAreaExtent) #recortando area invadida
plot(wrld) #plotando mapa mundi
plot(areaNat,ad=TRUE,col='blue') #verificando
plot(areaInv,ad=TRUE,col='red') #verificando

## recortando pontos para area nativa
spNat = spData[which(spData$lon > natAreaExtent[1] &
                     spData$lon < natAreaExtent[2] &
                     spData$lat > natAreaExtent[3] &
                     spData$lat < natAreaExtent[4]),]

## recortando pontos para area invadida
spInv = spData[which(spData$lon > invAreaExtent[1] &
                     spData$lon < invAreaExtent[2] &
                     spData$lat > invAreaExtent[3] &
                     spData$lat < invAreaExtent[4]),]

## inspecao visual
plot(wrld) #plotando mapa mundi
plot(areaNat,ad=TRUE,col='blue') #verificando
plot(areaInv,ad=TRUE,col='red') #verificando
points(spNat, pch=20, col='red') #pontos da area nativa
points(spInv, pch=20, col='blue') #pontos da area invadida


##area acessivel para a especie na area nativa##

##minimo poligono convexo##
occMat = as.matrix(spNat[,c('lon','lat')]) #sps occ matrix

##centroide da distribuicao espacial

SpsCentroid = kmeans(occMat, 1)
SpsCentroid = SpsCentroid$centers

##espacializacao dos pontos
occMatSpatial = data.frame(occMat)
coordinates(occMatSpatial) = ~lon+lat
projection(occMatSpatial) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

##convex hull
cHull = dismo::convHull(occMatSpatial)

##distancia entre pontos e centroide
euc.dist <- function(x1) sqrt(sum((x1 - SpsCentroid) ^ 2))
meanDistCentroid = mean( apply(occMat, 1, euc.dist) )

##buffer - hipotese para area de alcance da especie
SpsBufferNativ = raster::buffer(polygons(cHull), width=meanDistCentroid)

##salvando no HD
save(SpsBufferNativ,file=paste(projectFolder,'/SpsBufferNativ.R',sep=''))


##area acessivel para a especie na area invadida##

##minimo poligono convexo##
occMat = as.matrix(spInv[,c('lon','lat')]) #sps occ matrix

##centroide da distribuicao espacial
SpsCentroid = kmeans(occMat, 1)
SpsCentroid = SpsCentroid$centers

##espacializacao dos pontos
occMatSpatial = data.frame(occMat)
coordinates(occMatSpatial) = ~lon+lat
projection(occMatSpatial) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

##convex hull
cHull = dismo::convHull(occMatSpatial)

##distancia entre pontos e centroide
euc.dist <- function(x1) sqrt(sum((x1 - SpsCentroid) ^ 2))
meanDistCentroid = mean( apply(occMat, 1, euc.dist) )

##buffer - hipotese para area de alcance da especie
SpsBufferInv = raster::buffer(polygons(cHull), width=meanDistCentroid)

##salvando no HD
save(SpsBufferInv,file=paste(projectFolder,'/SpsBufferInv.R',sep=''))

##inspecao visual
plot(wrld) #plotando mapa mundi
points(spNat, pch=19)
plot(SpsBufferNativ, add=T)
points(spInv, pch=19)
plot(SpsBufferInv, add=T)


##recortando as variaveis preditoras para a area acessivel para as especies##

predAreaNat = crop(x=stack(predictors), y=extent(SpsBufferNativ)) #recortando area selecionada
predAreaInv = crop(x=stack(predictors), y=extent(SpsBufferInv)) #recortando area selecionada


##analisando correlacao das variaveis na area nativa##

predictorsForVif = predAreaNat

vif(predictorsForVif)
predictorsVif1 = vifcor(predictorsForVif, th=0.7)
predictorsVif1

predictorsVif2 <- vifstep(predictorsForVif, th=10) # identify collinear variables that should be excluded
predictorsVif2

##comparando
predictorsVif1@results$Variables
predictorsVif2@results$Variables

##definindo variaveis preditoras a serem usadas nos modelos
predAreaNat = predAreaNat[[ predictorsVif1@results$Variables ]]
predAreaInv = predAreaInv[[ predictorsVif1@results$Variables ]]


##salvando no HD
writeRaster(x=predAreaNat, filename=paste(projectFolder,'/variaveis_ambientais/predAreaNat.asc',sep=''), overwrite=TRUE, bylayer=TRUE, suffix=names(predAreaNat))
writeRaster(x=predAreaInv, filename=paste(projectFolder,'/variaveis_ambientais/predAreaInv.asc',sep=''), overwrite=TRUE, bylayer=TRUE, suffix=names(predAreaInv))


##KDE para vies amostral a partir de ocorrencias no GBIF para o genero Omobranchus##

##abrindo arquivo salvo
omobranchusDataset = read.csv(file='/home/anderson/Projetos/Invasao_Omobranchus_punctatus/GBIF - genero Omobranchus - 13-set-2018/0000608-180412121330197.csv', header=TRUE, sep='\t', dec='.', stringsAsFactors=FALSE, na.strings="")
omobranchusOcc = omobranchusDataset[,c('decimallongitude','decimallatitude')]
names(omobranchusOcc) = c('lon','lat')
omobranchusOcc = omobranchusOcc[complete.cases(omobranchusOcc),]

## ##convertendo de data.frame para staialPointsDataFrame
## coordinates(omobranchusOcc) <- ~lon+lat
## proj4string(omobranchusOcc) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')

## ##abrindo o shapefile que 'cortara' os pontos
## wrld = readOGR("/home/anderson/shapefiles/ne_50m_land/ne_50m_land.shp")

## ##ajustando o sistema de coordenadas geograficas entre pontos e shapefile
## omobranchusOcc <- spTransform(omobranchusOcc, CRS(proj4string(wrld)))
## omobranchusOcc <- omobranchusOcc[wrld, ]

## ##sjuatando para toda a area da america do sul e grid file final
## #SOAextent = extent(-81.57551,-34.03384,-57.13385,12.99115)

## SAbg = predictors[[1]]*0 ##America do Sul como 'pano de fundo'
## crs(SAbg) = crs(raster())

## ##inspecionando pontos
## plot(omobranchusOcc)
## plot(wrld, add=TRUE)

## ##convertendo spatialPoints em data.frame para o KDE2D
## reduviidaeOcc = as.data.frame(reduviidaeOcc)


##kernel density estimation com o pacote MASS
dens = MASS::kde2d(x=omobranchusOcc$lon, y=omobranchusOcc$lat, n=100)
densRas = raster(dens)
densRas = mask(x=densRas, mask=wrld)
densRas = crop(x=densRas, y=extent(predAreaNat))


##ajustando projecao e fundindo com pano de fundo
predAreaNatBG = predAreaNat[[1]]*0
densRas = projectRaster(densRas, crs=proj4string(predAreaNatBG), res=res(predAreaNatBG), method="bilinear")
densRas = merge(densRas, predAreaNatBG, tolerance=0.5)


##salvado raster
writeRaster(x=densRas, file='/home/anderson/Projetos/Invasao_Omobranchus_punctatus/omobranchusBiasLayer.grd', overwrite=TRUE)


##pseudo-ausencia na area nativa com o mesmo vies dos dados de ocorrencia
biasLayer = raster(paste(projectFolder,'/omobranchusBiasLayer.grd',sep=''))
values(biasLayer)[values(biasLayer)<=0] = 0
bgPoints = dismo::randomPoints(mask=biasLayer, n=10000, p=spNat, prob=TRUE)
bgPoints = data.frame(lon=bgPoints[,1], lat=bgPoints[,2])


##consolindando dados de presenca e background para AREA NATIVA
dataSet = data.frame(lon = c(spNat$lon, spInv$lon, bgPoints$lon),
                     lat = c(spNat$lat, spInv$lat, bgPoints$lat),
                     occ = c(rep(1,nrow(rbind(spNat,spInv))),rep(0,nrow(bgPoints))),
                     area = c(rep('nat',nrow(spNat)), rep('inv',nrow(spInv)),rep('bg',nrow(bgPoints))))
#arredondando para garantir
dataSet[,c('lon','lat')] = round(dataSet[,c('lon','lat')], 2)
##variaveis ambientais
dataSetVars = extract(x=predictors[[predictorsVif1@results$Variables]], y=dataSet[,c('lon','lat')], method='bilinear', na.rm=TRUE) #extraindo variaeis ambientais
dataSet = data.frame(dataSet, dataSetVars) #juntando ao dataset
dataSet = dataSet[complete.cases(dataSet),] #retirando dados errados
dataSet = dataSet[!duplicated(dataSet[,c('lon','lat')]),] #retirando pontos numa mesma celula
##salvando dataset na memoria do HD
write.csv(dataSet,paste(projectFolder,'/dataSet.csv',sep=''),row.names=FALSE)


##mapa dos pontos de ocorrencia de back ground
jpeg(paste(projectFolder,'/mapaOccDataAreaNat.jpg',sep=''), width=600, height=600)
par(oma=c(0,0,0.5,7), cex=0.5)
plot(predAreaNat[[1]]*0, col='lightblue', legend=FALSE,  bty="n", box=FALSE)
points(dataSet[dataSet$area=='bg',c('lon','lat')], pch=20, cex=0.5, col=rgb(0.5,0.5,0.5,0.3))
points(dataSet[dataSet$area=='nat',c('lon','lat')], pch=20, cex=1, col=rgb(1,0,0,0.5))
legend('bottomleft', legend=c('Occorrence points','Background points'), pch=20, col=c('red','gray'), bty='n', inset=c(1,0.7),xpd=NA)
dev.off()



## PARTE 3: Analise de similaridade de equivalencia de nicho




##abrindo variaveis ambientais
predAreaNat = stack(list.files(path=paste(projectFolder,'/variaveis_ambientais',sep=''),pattern='predAreaNat',full.names=TRUE))
predAreaInv = stack(list.files(path=paste(projectFolder,'/variaveis_ambientais',sep=''),pattern='predAreaInv',full.names=TRUE))
#predAreaInv = predAreaInv[[gsub(pattern='predictors_',replacement='',x=names(predAreaNat))]]

## criando pontos de background
bgAreaNat = randomPoints(mask=predAreaNat[[1]], n=1000, p=spNat, prob=FALSE)
bgAreaNat = data.frame(bgAreaNat)
names(bgAreaNat) = names(spNat)
bgAreaInv = dismo::randomPoints(mask=predAreaInv[[1]], n=1000, p=spInv, prob=FALSE)
bgAreaInv = data.frame(bgAreaInv)
names(bgAreaInv) = names(spInv)

## definindo dados de pres a ausencia
spNatData = data.frame(rbind(spNat,bgAreaNat),occ=c(rep(1,nrow(spNat)),rep(0,nrow(bgAreaNat))))
spInvData = data.frame(rbind(spInv,bgAreaInv),occ=c(rep(1,nrow(spInv)),rep(0,nrow(bgAreaInv))))

## extarindo dados das variaveis ambientais
spNatDataEnv = extract(predAreaNat, spNatData[,c('lon','lat')],method='bilinear',na.romove=TRUE)
spInvDataEnv = extract(predAreaInv, spInvData[,c('lon','lat')],method='bilinear',na.romove=TRUE)

## juntando dados de ocorrencia e dados das variaveis ambientais
spNatData = data.frame(spNatData,spNatDataEnv)
spInvData = data.frame(spInvData,spInvDataEnv)

## garantindo dados completos
spNatData = spNatData[complete.cases(spNatData),]
spInvData = spInvData[complete.cases(spInvData),]

## deixando os nomes iguais (necessario para o metodo)
names(spNatData) = gsub(pattern='predAreaNat_', replacement='', x=names(spNatData))
names(spInvData) = gsub(pattern='predAreaInv_', replacement='', x=names(spInvData))

## The PCA is calibrated on all the sites of the study area
pca.env <- dudi.pca(rbind(spNatData,spInvData)[,grep(pattern='2o5', x=names(spNatData), value=TRUE)],scannf=F,nf=2)
## ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) #grafico

## PCA scores for the whole study area
scores.globclim <- pca.env$li

## PCA scores for the species native distribution
scores.sp.nat <- suprow(pca.env,spNatData[which(spNatData[,'occ']==1),grep(pattern='2o5', x=names(spNatData), value=TRUE)])$li

## PCA scores for the species invasive distribution
scores.sp.inv <- suprow(pca.env,spInvData[which(spInvData[,'occ']==1),grep(pattern='2o5', x=names(spInvData), value=TRUE)])$li

## PCA scores for the whole native study area
scores.clim.nat <-suprow(pca.env,spNatData[,grep(pattern='2o5', x=names(spNatData), value=TRUE)])$li

## PCA scores for the whole invaded study area
scores.clim.inv <- suprow(pca.env,spInvData[,grep(pattern='2o5', x=names(spInvData), value=TRUE)])$li

## gridding the native niche
grid.clim.nat <-ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.nat, sp=scores.sp.nat, R=100, th.sp=0)

## gridding the invasive niche
grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.inv, sp=scores.sp.inv, R=100, th.sp=0)

## equivalencia de nicho
equi.test <- ecospat.niche.equivalency.test(grid.clim.nat, grid.clim.inv, rep=100, alternative="greater")

## similaridade de nicho
simi.test <- ecospat.niche.similarity.test(grid.clim.nat, grid.clim.inv, rep=100, alternative="greater")

## tabela de resultados
nicheOverlapTab = data.frame(teste=c('similaridade','equivalencia'),
                             Schoener=c(simi.test$obs$D,equi.test$obs$D),
                             Schoener_p_value=c(simi.test$p.D,equi.test$p.D),
                             Hellinger=c(simi.test$obs$I,equi.test$obs$I),
                             Hellinger_p_value=c(simi.test$p.I,equi.test$p.I))

write.csv(nicheOverlapTab, paste(projectFolder,'/nicheOverlapTab.csv',sep=''), row.names=FALSE)




## PARTE 4: modelagem da distribuicao da especie




## ajustando o diretorio de trabalho (pois o biomod roda e salva tudo simultaneamente)
if(!file.exists(file.path(projectFolder,'maxent',sep=''))){
    dir.create(file.path(projectFolder,'maxent',sep=''),recursive=TRUE)
}
setwd(file.path(projectFolder,'maxent'))


##abrindo as variaveis ambientais
predAreaNat = stack(list.files(path=paste(projectFolder,'/variaveis_ambientais',sep=''),pattern='predAreaNat',full.names=TRUE))
predAreaInv = stack(list.files(path=paste(projectFolder,'/variaveis_ambientais',sep=''),pattern='predAreaInv',full.names=TRUE))


##abrindo dataSet (qdo ja estiver pronto)
dataSet = read.csv(paste(projectFolder,'/dataSet.csv',sep=''), header=TRUE)


##abrindo poligonos das distribuicoes nativa e invadida
load(paste(projectFolder,"/SpsBufferInv.R",sep=''))
load(paste(projectFolder,"/SpsBufferNativ.R",sep=''))



##ENMeval##



ENMblock = get.block(occ=dataSet[dataSet$occ==1 & dataSet$area=='nat',c('lon','lat')], bg.coords=dataSet[dataSet$occ==0 & dataSet$area=='bg',c('lon','lat')]) #dividindo em blocos

ENMblockUser = get.user(occ.grp=ENMblock$occ.grp, bg.grp=ENMblock$bg.grp) #usuario, por causa dos dados proprios de bg-points gerados com vies

rm(list=c('predAreaInv','predictors','spData','spInv','spNat','wrld','bgPoints','invAreaExtent','natAreaExtent','dataSetVars','areaInv','areaNat'))
gc()

SDMeval <- ENMevaluate(occ = dataSet[dataSet$occ==1 & dataSet$area=='nat',c('lon','lat')],
                       env = predAreaNat,
                       bg.coords = dataSet[dataSet$occ==0,c('lon','lat')],
                       occ.grp = ENMblockUser$occ.grp,
                       bg.grp = ENMblockUser$bg.grp,
                       method = 'user',
                       RMvalues = c(1,2,4),
                       fc = c("L", "LQ", "LQH", "LQHPT"),
                       clamp = FALSE,
                       parallel = FALSE)
                       
                       ## parallel = TRUE,
                       ## numCores = 2)


##melhor modelo
bestModel = SDMeval@results[SDMeval@results$AICc==min(SDMeval@results$AICc, na.rm=TRUE),]
bestModel = bestModel[complete.cases(bestModel),]
bestModel = bestModel[1,]


##importancia das variaveis
modelPars = SDMeval@models[[bestModel$settings]]
write.csv(var.importance(modelPars), paste(projectFolder,'/maxent/omobranchus_variablesImportance.csv',sep=''), row.names=FALSE)


##excluindo variaveis sem imporntacia para o modelo (ver importancia das variaveis do SDMeval)
dataSet = dataSet[,c('lon','lat','occ','area','bathy_2o5m','biogeo01_2o5m','biogeo02_2o5m','biogeo05_2o5m','biogeo06_2o5m','biogeo08_2o5m')]
predAreaNat = predAreaNat[[c('predAreaNat_bathy_2o5m','predAreaNat_biogeo01_2o5m','predAreaNat_biogeo02_2o5m','predAreaNat_biogeo05_2o5m','predAreaNat_biogeo06_2o5m','predAreaNat_biogeo08_2o5m')]]
predAreaInv = predAreaInv[[c('predAreaInv_bathy_2o5m','predAreaInv_biogeo01_2o5m','predAreaInv_biogeo02_2o5m','predAreaInv_biogeo05_2o5m','predAreaInv_biogeo06_2o5m','predAreaInv_biogeo08_2o5m')]]



##MAXENT##
SDMmaxent = maxent(x = dataSet[,grep(pattern='2o5m',x=names(dataSet),value=TRUE)],
                   p = dataSet[,'occ'],
                   args = c(unlist(make.args(RMvalues=bestModel$rm, fc=bestModel$features, labels=FALSE)),'threads=2'))


##predicao AREA NATIVA
SDMpredNativ = predict(predAreaNat, SDMmaxent) #projecao espacial
crs(SDMpredNativ) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'


##predicao AREA INVADIDA
SDMpredInv = predict(predAreaInv, SDMmaxent) #projecao espacial
crs(SDMpredInv) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'


##salvando gridfiles
writeRaster(SDMpredNativ, paste(projectFolder,'/maxent/SDMpredNativ_Suitability.asc',sep=''), overwrite=TRUE) #area nativa
writeRaster(SDMpredInv,paste(projectFolder,'/maxent/SDMpredInv_Suitability.asc',sep=''), overwrite=TRUE) #area nativa


##calculo do threshold##
##valores preditos pelo SDM em cada ponto do dataset na AREA NATIVA
pred =  extract(x=SDMpredNativ,
                y=dataSet[dataSet$area!='inv',c('lon','lat')],
                method='bilinear',
                na.rm=TRUE)
threshDF = data.frame(occ=dataSet[dataSet$area!='inv',]$occ, pred=pred) #juntando predicoes e observacoes
threshDF = threshDF[complete.cases(threshDF),] #retirando possiveis NAs

##OBS.: omissao -> falso negativo -> false negative rate (FNR) = 1 - TPR (True Positive Rate)
##TPR = sensitividade
myRoc = roc(predictor=threshDF$pred, response=threshDF$occ, positive=1) #curva ROC
rocDF = data.frame( FNR = 1-myRoc$sensitivities, thresholds = myRoc$thresholds ) #data.frame com omissao e thresholds
thre = max(rocDF[which(rocDF$FNR == max(rocDF[rocDF$FNR <= 0.05 ,]$FNR)),]$thresholds) #thresold de 5% omissao


##calculo do TSS para com o threshold -- OBS: TSS = sensitivity + specificity - 1
TSSvalue = myRoc$sensitivities[which(myRoc$thresholds == thre)] + myRoc$specificities[which(myRoc$thresholds == thre)] - 1


##salvando output do ENMeval
SDMevalOutput = as.data.frame(SDMeval@results)
write.csv(SDMevalOutput, paste(projectFolder,'/maxent/omobranchus_SDMeval.csv',sep=''), row.names=FALSE)


##salvando graficos do ENMeval
jpeg(paste(projectFolder,'/maxent/omobranchus_SDMeval.jpg',sep=''), width=800)
par(mfrow=c(1,2), mar=c(5,5,10,5))
eval.plot(SDMeval@results)
eval.plot(SDMeval@results, 'Mean.AUC', var='Var.AUC')
dev.off()


##salvando imagens jpeg no HD

mapInv = raster(paste(projectFolder,'/maxent/SDMpredInv_Suitability.asc',sep='')) #abrindo arquivo salvo

jpeg(paste(projectFolder,'/mapInv.jpg',sep=''), width=900, height=900)
par(oma=c(0,0,0,8), cex=1)
plot(mapInv > thre, legend=FALSE, bty="n", box=FALSE)
points(dataSet[dataSet$area=='inv',c('lon','lat')], pch=20, cex=1, col='red')
legend('bottomleft', legend=c('Suitable habitat','Non-suitable habitat', 'Occurrence points'), col=c('green','lightgrey', 'red'), pch=c(15,15,20), bty='n', inset=c(1,0.9),xpd=NA, cex=1.3)
dev.off()

mapNat = raster(paste(projectFolder,'/maxent/SDMpredNativ_Suitability.asc',sep='')) #abrindo arquivo salvo
jpeg(paste(projectFolder,'/mapNat.jpg',sep=''), width=900, height=900)
par(oma=c(0,0,0,8), cex=1)
plot(mapNat > thre, legend=FALSE, bty="n", box=FALSE)
points(dataSet[dataSet$area=='nat',c('lon','lat')], pch=20, cex=1, col='red')
legend('bottomleft', legend=c('Suitable habitat','Non-suitable habitat', 'Occurrence points'), col=c('green','lightgrey', 'red'), pch=c(15,15,20), bty='n', inset=c(1,0.9),xpd=NA, cex=1.3)
dev.off()




## PARTE 4: Analise de pDLA - iSDM  ##NAO FUNCIONOU...  :(




##abrindo as variaveis ambientais
predAreaNat = stack(list.files(path=paste(projectFolder,'/variaveis_ambientais',sep=''),pattern='predAreaNat',full.names=TRUE))
predAreaInv = stack(list.files(path=paste(projectFolder,'/variaveis_ambientais',sep=''),pattern='predAreaInv',full.names=TRUE))


##abrindo dataSet (qdo ja estiver pronto)
dataSet = read.csv(paste(projectFolder,'/dataSet.csv',sep=''), header=TRUE)


##variaveis ambientais
envDataInv = predAreaInv
envDataNativ = predAreaNat
names(envDataNativ) = gsub(pattern='predAreaNat_', replacement='', x=names(envDataNativ))
names(envDataInv) = names(envDataNativ) #necessario nome siguais para a funcao pDLA
proj4string(envDataNativ) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
proj4string(envDataInv) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

##occorencias
occInv = dataSet[which(dataSet$area == 'inv'),c('lon','lat','occ')] #rendera um spatial data.frame
occInvBG = randomPoints(mask=envDataInv, n=1000)
occInvBG = round(occInvBG, digits=2)
occInvBG = data.frame(lon=occInvBG[,1], lat=occInvBG[,2], occ=0)
occInv = rbind(occInv, occInvBG)
occNat = dataSet[which(dataSet$area == 'nat'),c('lon','lat')] #rendera um spatial points 


##espacializando pontos de ocorrencia
coordinates(occInv) = ~lon+lat
coordinates(occNat) = ~lon+lat
proj4string(occInv) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
proj4string(occNat) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"


##pDLA

pobability = pDLA(occData=occInv,
                  envData=envDataInv,
                  longlat=TRUE,
                  occNative=occNat,
                  envNative=envDataNativ)

##Display the results
par(mfrow=c(1,2),mar=c(2,2.5,2,2.5))
plot(realized.dist$occupied.area,main="Realized distribution")
plot(occData,col=ifelse(occData$SP==1,2,1),add=TRUE,pch=19,cex=0.8)

plot(SDMpredInv, main="Potential distribution")

scatterCol<-function(x){
    x <- (x-min(x, na.rm=TRUE))/(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))
    colorFunction <- colorRamp(colorRamps::matlab.like(100))
    zMatrix <- colorFunction(x)
    zColors <- rgb(zMatrix[,1], zMatrix[,2], zMatrix[,3], maxColorValue=255)
    return(zColors)
}


points(probability,pch=21, col=1,bg=scatterCol(probability@data[,"PDLA"]),cex=1)


