## Script para analise do nicho de uma especie invasora, comparando nicho do ambiente nativo e do ambiente invadido.
## Anderson A. Eduardo
## novembro/2017




## PARTE 1: carregando pacotes, parâmetros e variaveis




## limpando a memória (qdo necessario)
rm(list=ls(all=TRUE))
gc()

## abrindo pacotes necessarios
#Sys.setenv(JAVA_HOME = "/usr/lib/jvm/java-7-openjdk-amd64")
#Windows#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_91') # for 64-bit version
options(java.parameters = "Xmx7g")
library(ecospat)
library(rgdal)
library(raster)
library(ENMeval)
library(usdm)
library(iSDM)
library(pROC)
library(rJava)



# ## definindo variaveis e parametros - Notebook
# projectFolder = "/home/anderson/Projetos/Invasao_Omobranchus_punctatus" #pasta do projeto
# envVarFolder = "/home/anderson/Projetos/Invasao_Omobranchus_punctatus/variaveis_ambientais" #pasta com as variaveis ambientais
# predictors = stack(list.files(path=envVarFolder, full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis
# predictors = predictors[[grep(pattern=paste(c('Chlorophyll','Phytoplankton','Silicate','*.Mean$','*.Max$','*.Min$'),collapse='|'), names(predictors), value=FALSE,invert=TRUE)]]
# spData = read.csv(file.path(projectFolder,'spOcc.csv'),header=TRUE) #dados de ocorrencia ambiente nativo
# names(spData) = c('lon','lat')
# wrld = readOGR('/home/anderson/shapefiles/ne_50m_ocean/ne_50m_ocean.shp') #mapa mundi


## definindo variaveis e parametros - Iavana
projectFolder = "D:/Anderson_Eduardo/Invasao_Omobranchus_punctatus" #pasta do projeto
envVarFolder = "D:/Anderson_Eduardo/gridfiles/Bio-ORACLE - variaveis bentonicas" #pasta com as variaveis ambientais
predictors = stack(list.files(path=envVarFolder, full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis
predictors = predictors[[grep(pattern=paste(c('Chlorophyll','Phytoplankton','Silicate','*.Mean$','*.Max$','*.Min$'),collapse='|'), names(predictors), value=FALSE,invert=TRUE)]]
spData = read.csv(file.path(projectFolder,'spOcc.csv'),header=TRUE) #dados de ocorrencia ambiente nativo
names(spData) = c('lon','lat')
wrld = readOGR('D:/Anderson_Eduardo/shapefiles/ne_50m_ocean/ne_50m_ocean.shp') #mapa mundi




## PARTE 2: definindo dados da area nativa e area invadida




## separando area nativa e area invadida
#plot(wrld) #criando imagem
#points(spData,col='red') #plotando os pontos de ocorencia da sp
#drawExtent() #delimitando area no mapa
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

predAreaNat = crop(x=stack(predictors), y=natAreaExtent) #recortando area selecionada
predAreaInv = crop(x=stack(predictors), y=invAreaExtent) #recortando area selecionada


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
predAreaNat = predAreaNat[[ grep(pattern=paste(predictorsVif2@results$Variables,collapse='|'), x=names(predAreaNat), value=TRUE) ]]
predAreaInv = predAreaInv[[ grep(pattern=paste(predictorsVif2@results$Variables,collapse='|'), x=names(predAreaInv), value=TRUE) ]]

##salvando no HD
writeRaster(x=predAreaNat, filename=paste(projectFolder,'/variaveis_ambientais/predAreaNat.asc',sep=''), overwrite=TRUE, bylayer=TRUE, suffix=names(predAreaNat))
writeRaster(x=predAreaInv, filename=paste(projectFolder,'/variaveis_ambientais/predAreaInv.asc',sep=''), overwrite=TRUE, bylayer=TRUE, suffix=names(predAreaInv))


##KDE para vies amostral a partir de ocorrencias no GBIF para o genero Omobranchus##

##abrindo arquivo salvo
##omobranchusDataset = read.csv(file='D:/Anderson_Eduardo/Invasao_Omobranchus_punctatus/GBIF - genero Omobranchus - 13-set-2018/0000608-180412121330197.csv', header=TRUE, sep='\t', dec='.', stringsAsFactors=FALSE, na.strings="")
omobranchusDataset = read.csv(file='/home/anderson/Projetos/Invasao_Omobranchus_punctatus/GBIF - genero Omobranchus - 13-set-2018/0000608-180412121330197.csv', header=TRUE, sep='\t', dec='.', stringsAsFactors=FALSE, na.strings="")
omobranchusOcc = omobranchusDataset[,c('decimallongitude','decimallatitude')]
names(omobranchusOcc) = c('lon','lat')
omobranchusOcc = omobranchusOcc[complete.cases(omobranchusOcc),]

##convertendo de data.frame para staialPointsDataFrame
coordinates(omobranchusOcc) <- ~lon+lat
proj4string(omobranchusOcc) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')

## ##abrindo o shapefile que 'cortara' os pontos
## wrld = readOGR("/home/anderson/shapefiles/ne_50m_land/ne_50m_land.shp")

##ajustando o sistema de coordenadas geograficas entre pontos e shapefile
omobranchusOcc <- spTransform(omobranchusOcc, CRS(proj4string(wrld)))
omobranchusOcc <- omobranchusOcc[wrld, ]

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
values(densRas)[values(densRas)<=0] = 0
densRas = mask(x=densRas, mask=wrld)
densRas = crop(x=densRas, y=natAreaExtent)


##ajustando projecao e fundindo com pano de fundo
predAreaNatBG = raster('D:/Anderson_Eduardo/Invasao_Omobranchus_punctatus/variaveis_ambientais/predAreaNat_Present.Benthic.Mean.Depth.Primary.productivity.Lt.min.asc')
crs(predAreaNatBG) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
#predAreaNatBG = crop(x=predAreaNatBG, y=natAreaExtent)
predAreaNatBG = predAreaNatBG*0
densRas = projectRaster(densRas, crs=proj4string(predAreaNatBG), res=res(predAreaNatBG), method="bilinear")
densRas = merge(densRas, predAreaNatBG, tolerance=0.3)

extent(densRas) = extent(predAreaNatBG)
densRas = crop(x=densRas, y=predAreaNatBG)
res(densRas) =  res(predAreaNatBG)
densRas = extend(x=densRas, y=extent(predAreaNatBG), value=NA)

values(densRas)[values(densRas)<=0] = 0
densRas = mask(x=densRas, mask=wrld)

##salvado raster
writeRaster(x=densRas, file=paste(projectFolder,'/omobranchusBiasLayer.asc',sep=''), overwrite=TRUE)


##pseudo-ausencia na area nativa com o mesmo vies dos dados de ocorrencia
biasLayer = raster(paste(projectFolder,'/omobranchusBiasLayer.asc',sep=''))
bgPoints = randomPoints(mask=biasLayer, n=10000, p=spNat, prob=TRUE)
bgPoints = data.frame(lon=bgPoints[,1], lat=bgPoints[,2])


##consolindando dados de presenca e background para AREA NATIVA
dataSet = data.frame(lon = c(spNat$lon, spInv$lon, bgPoints$lon),
                     lat = c(spNat$lat, spInv$lat, bgPoints$lat),
                     occ = c(rep(1,nrow(rbind(spNat,spInv))),rep(0,nrow(bgPoints))),
                     area = c(rep('nat',nrow(spNat)), rep('inv',nrow(spInv)),rep('bg',nrow(bgPoints))))
#arredondando para garantir
dataSet[,c('lon','lat')] = round(dataSet[,c('lon','lat')], 2)
##variaveis ambientais
dataSetVars = extract(x=predictors[[ grep(pattern=paste(predictorsVif2@results$Variables,collapse='|'), x=names(predictors), value=TRUE) ]],
                      y=dataSet[,c('lon','lat')], 
                      method='bilinear', 
                      na.rm=TRUE) #extraindo variaeis ambientais
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
pca.env <- dudi.pca(rbind(spNatData,spInvData)[,grep(pattern='Present.Benthic.Mean.Depth.', x=names(spNatData), value=TRUE)],scannf=F,nf=2)
## ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) #grafico

## PCA scores for the whole study area
scores.globclim <- pca.env$li

## PCA scores for the species native distribution
scores.sp.nat <- suprow(pca.env,spNatData[which(spNatData[,'occ']==1),grep(pattern='Present.Benthic.Mean.Depth.', x=names(spNatData), value=TRUE)])$li

## PCA scores for the species invasive distribution
scores.sp.inv <- suprow(pca.env,spInvData[which(spInvData[,'occ']==1),grep(pattern='Present.Benthic.Mean.Depth.', x=names(spInvData), value=TRUE)])$li

## PCA scores for the whole native study area
scores.clim.nat <-suprow(pca.env,spNatData[,grep(pattern='Present.Benthic.Mean.Depth.', x=names(spNatData), value=TRUE)])$li

## PCA scores for the whole invaded study area
scores.clim.inv <- suprow(pca.env,spInvData[,grep(pattern='Present.Benthic.Mean.Depth.', x=names(spInvData), value=TRUE)])$li

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

gc()

SDMeval <- ENMevaluate(occ = dataSet[dataSet$occ==1 & dataSet$area=='nat',c('lon','lat')],
                       env = predAreaNat,
                       bg.coords = dataSet[dataSet$occ==0,c('lon','lat')],
                       occ.grp = ENMblockUser$occ.grp,
                       bg.grp = ENMblockUser$bg.grp,
                       method = 'user',
                       RMvalues = c(0.5,1,2,4),
                       fc = c("L", "LQ", "LQH", "LQHPT"),
                       clamp = FALSE,                       
                       parallel = TRUE,
                       numCores = 2)


##melhor modelo
bestModel = SDMeval@results[SDMeval@results$AICc==min(SDMeval@results$AICc, na.rm=TRUE),]
bestModel = bestModel[complete.cases(bestModel),]
bestModel = bestModel[1,]


##importancia das variaveis
modelPars = SDMeval@models[[bestModel$settings]]
##write.csv(var.importance(modelPars), paste(projectFolder,'/maxent/omobranchus_variablesImportance.csv',sep=''), row.names=FALSE)
varNames = as.vector(var.importance(modelPars)[which(var.importance(modelPars)$permutation.importance > 0),'variable'])
varNames = gsub(pattern='predAreaNat_', replacement='', x=varNames)


##excluindo variaveis sem imporntacia para o modelo (ver importancia das variaveis do SDMeval)
dataSet = dataSet[,c('lon','lat','occ','area',grep(pattern=paste(varNames, collapse='|'), x=names(dataSet), value=TRUE))]
write.csv(dataSet,paste(projectFolder,'/dataSet.csv',sep=''),row.names=FALSE) #atualiza o dataset salvo no hd


##MAXENT##
SDMmaxent = maxent(x = dataSet[which(dataSet$area!='inv'),grep(pattern=paste(as.character(varNames),collapse="|"),x=names(dataSet),value=TRUE)],
                   p = dataSet[which(dataSet$area!='inv'),'occ'],
                   args = c(unlist(make.args(RMvalues=bestModel$rm, fc=bestModel$features, labels=FALSE)),
                            'replicates=100',
                            'randomtestpoints=75',
                            'replicatetype=Subsample',
                            'responsecurves=TRUE',
                            'writemess=TRUE',
                            'pictures=TRUE',
                            'extrapolate=FALSE',
                            'doclamp=FALSE',
                            'fadebyclamping=FALSE',
                            'plots=TRUE',
                            'threads=2'))


##salvando modelo no HD
save(SDMmaxent, file=paste(projectFolder,'/maxent/SDMmaxent.R',sep=''))


##variableimportance final model
varImportOutput = rowMeans(SDMmaxent@results[grep("importance", rownames(SDMmaxent@results)), ])
varImportTable = data.frame(names(varImportOutput), permutation.importance=varImportOutput)
write.csv(varImportTable, paste(projectFolder,'/maxent/omobranchus_variablesImportance_finalModel.csv',sep=''), row.names=FALSE)


##predicao AREA NATIVA
names(predAreaNat) = gsub( pattern='^predAreaNat_', replacement='', x=names(predAreaNat) ) #retirando parte do nome para adequar ao modelo
predAreaNat = predAreaNat[[ grep(pattern=paste(varNames,collapse='|'), x=names(predAreaNat), value=TRUE) ]] #novo objeto, correto
SDMpredNativ = dismo::predict(SDMmaxent, predAreaNat, progress='text') #projecao espacial
crs(SDMpredNativ) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'


##predicao AREA INVADIDA
names(predAreaInv) = gsub( pattern='^predAreaInv_', replacement='', x=names(predAreaInv) ) #retirando parte do nome para adequar ao modelo
predAreaInv = predAreaInv[[ grep(pattern=paste(varNames,collapse='|'), x=names(predAreaNat), value=TRUE) ]] #novo objeto, correto
SDMpredInv = dismo::predict(SDMmaxent, predAreaInv, progress='text') #projecao espacial
crs(SDMpredInv) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'


##salvando gridfiles
writeRaster(calc(SDMpredNativ,mean), paste(projectFolder,'/maxent/SDMpredNativ_Suitability_Mean.asc',sep=''), overwrite=TRUE) #area nativa
writeRaster(calc(SDMpredNativ,sd), paste(projectFolder,'/maxent/SDMpredNativ_Suitability_SD.asc',sep=''), overwrite=TRUE) #area nativa
##
writeRaster(calc(SDMpredInv,mean),paste(projectFolder,'/maxent/SDMpredInv_Suitability_Mean.asc',sep=''), overwrite=TRUE) #area nativa
writeRaster(calc(SDMpredInv,sd),paste(projectFolder,'/maxent/SDMpredInv_Suitability_SD.asc',sep=''), overwrite=TRUE) #area nativa


##calculo do threshold##
##valores preditos pelo SDM em cada ponto do dataset na AREA NATIVA
pred =  extract(x=calc(SDMpredNativ,mean),
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

mapInv = raster(paste(projectFolder,'/maxent - Yavanna/SDMpredInv_Suitability_Mean.asc',sep='')) #abrindo arquivo salvo

##binario
jpeg(paste(projectFolder,'/mapInvBIN.jpg',sep=''), width=900, height=900)
par(oma=c(0,0,0,8), cex=1)
plot(mapInv > thre, legend=FALSE, bty="n", box=FALSE)
points(dataSet[dataSet$area=='inv',c('lon','lat')], pch=20, cex=1, col='red')
plot(wrld,lwd=0.1, add=TRUE)
legend('bottomleft', legend=c('Suitable habitat','Non-suitable habitat', 'Occurrence points'), col=c('green','lightgrey', 'red'), pch=c(15,15,20), bty='n', inset=c(1,0.9),xpd=NA, cex=1.3)
dev.off()

##suitability
jpeg(paste(projectFolder,'/mapInv.jpg',sep=''), width=900, height=900)
plot(mapInv, legend=TRUE, bty="n", box=FALSE)
points(dataSet[dataSet$area=='inv',c('lon','lat')], pch=20, cex=1, col='red')
plot(wrld,lwd=0.1, add=TRUE)
dev.off()


mapNat = raster(paste(projectFolder,'/maxent - Yavanna/SDMpredNativ_Suitability_Mean.asc',sep='')) #abrindo arquivo salvo

##binario
jpeg(paste(projectFolder,'/mapNatBIN.jpg',sep=''), width=900, height=900)
par(oma=c(0,0,0,8), cex=1)
plot(mapNat > thre, legend=FALSE, bty="n", box=FALSE)
points(dataSet[dataSet$area=='nat',c('lon','lat')], pch=20, cex=1, col='red')
plot(wrld,lwd=0.1, add=TRUE)
legend('bottomleft', legend=c('Suitable habitat','Non-suitable habitat', 'Occurrence points'), col=c('green','lightgrey', 'red'), pch=c(15,15,20), bty='n', inset=c(1,0.9),xpd=NA, cex=1.3)
dev.off()

##suitability
jpeg(paste(projectFolder,'/mapNat.jpg',sep=''), width=900, height=900)
plot(mapNat, legend=FALSE, bty="n", box=FALSE)
points(dataSet[dataSet$area=='nat',c('lon','lat')], pch=20, cex=1, col='red')
plot(wrld,lwd=0.1, add=TRUE)
dev.off()




## PARTE 4: Analise da correlacao do suitability com o Global Ocean Health Index (GOHI)




##mapa de suitability na area invadida
mapInv = raster(paste(projectFolder,'/maxent - Yavanna/SDMpredInv_Suitability_Mean.asc',sep=''), crs=CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'))
mapInv = crop(x=mapInv, y=areaInv) #so pra garantir

##indice (fazer apenas uma vez, salvar e usar o salvo - da mto trabalho ajustar esse negocio...)
gohi = raster('D:/Anderson_Eduardo/gridfiles/Global Ocean Health Index/cumulative_impact_one_2013_global_cumul_impact_2013_mol_20150714053146/global_cumul_impact_2013_all_layers.tif') ##global_cumul_impact_2013_all_layers_WGS.grd')
crs(mapInv) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
gohiInviAdjusted = projectRaster(gohi, mapInv)
#writeRaster(gohiInviAdjusted, filename=paste(projectFolder,'/global_cumul_impact_2013_all_layers_WGS_AREA_INVADIDA.asc',sep='' ))
gohiInviAdjusted = raster('D:/Anderson_Eduardo/gridfiles/Global Ocean Health Index/cumulative_impact_one_2013_global_cumul_impact_2013_mol_20150714053146/global_cumul_impact_2013_all_layers_WGS_AREA_INVADIDA.asc')

##extraindo valores de background para correlacao
##KDE area invadida (para background com mesmo vies amostral que as ocorrencias)
omobranchusDataset = read.csv(file=paste(projectFolder,'/GBIF - genero Omobranchus - 13-set-2018/0000608-180412121330197.csv',sep=''), header=TRUE, sep='\t', dec='.', stringsAsFactors=FALSE, na.strings="")
omobranchusOcc = omobranchusDataset[,c('decimallongitude','decimallatitude')]
names(omobranchusOcc) = c('lon','lat')
omobranchusOcc = omobranchusOcc[complete.cases(omobranchusOcc),]
coordinates(omobranchusOcc) <- ~lon+lat
proj4string(omobranchusOcc) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
##ajustando o sistema de coordenadas geograficas entre pontos e shapefile
omobranchusOcc <- spTransform(omobranchusOcc, CRS(proj4string(wrld)))
omobranchusOcc <- omobranchusOcc[wrld, ]
##kernel density estimation com o pacote MASS
dens = MASS::kde2d(x=omobranchusOcc$lon, y=omobranchusOcc$lat, n=100)
densRas = raster(dens)
values(densRas)[values(densRas)<=0] = 0
densRas = mask(x=densRas, mask=wrld)
densRasInv = crop(x=densRas, y=invAreaExtent)
##plot(densRasInv)
##plot(wrld,add=T)
##amostrando os pontos de bg
bgPts = randomPoints(mask=densRasInv, n=nrow(spInv)*2, p=spInv, prob=TRUE)
bgPts = data.frame(bgPts)
names(bgPts) = names(spInv)

##dataframe de pres/brackground
datasetGohi = data.frame(
  rbind(spInv,bgPts),
  pres = c(rep(1,nrow(spInv)),rep(0,nrow(bgPts))))
##alguns ajustes 
datasetGohi = datasetGohi[complete.cases(datasetGohi),]
datasetGohi = round(datasetGohi, digits=2)
datasetGohi = unique(datasetGohi)

##extraindo valores do indice GOHI e do suitability nos pontos (presenca/background)
gohiPts = extract(x=gohiInviAdjusted, y=datasetGohi[,c('lon','lat')])
suitabilityPts = extract(x=mapInv, y=datasetGohi[,c('lon','lat')])

##consolidando dataset
datasetGohi = cbind(datasetGohi,gohiPts,suitabilityPts)
datasetGohi = datasetGohi[complete.cases(datasetGohi),]


##modelos estatisticos##


##modelos estatisticos - GLM para presenca/background X GOHI##
modelGohi = glm(pres ~ gohiPts, data=datasetGohi, family=binomial) #modelo linear
summary(modelGohi)

modelGohiQuad = glm(pres ~ gohiPts + I(gohiPts^2), data=datasetGohi, family=binomial) #modelo quadratico
summary(modelGohiQuad)

aicPresXgohi =  AIC(modelGohi, modelGohiQuad) #comparando por AIC
anovaPresXgohi = anova(modelGohi, modelGohiQuad, test='Chisq') #comparando por Chisq

##visualizacao grafica
jpeg(filename=paste(projectFolder,'/presXgohi.jpg',sep=''))
plot(datasetGohi$pres ~ datasetGohi$gohiPts, xlim=range(datasetGohi$gohiPts), pch=20, ylab='Occurrrence', xlab='GOHI') 
points(x = c(-100:100/10),
       y = predict(modelGohi, newdata=data.frame(gohiPts=c(-100:100/10)), type='response'),
       col='red', type='l')
points(x = c(-100:100/10),
       y = predict(modelGohiQuad, newdata=data.frame(gohiPts=c(-100:100/10)), type='response'),
       col='blue', type='l')
legend(x=3.3, y=0.95, legend=c('Linear model','Quadratic model'), lty=1, col=c('red','blue'))
dev.off()


## ##correlacao (linear) usando pontos de presenca/background na area de invasao
## corMatrixPts = cor( datasetGohi[which(datasetGohi$pres==1),'gohiPts'], datasetGohi[which(datasetGohi$pres==1),'suitabilityPts'] )
## save(corMatrixPts,file=paste(projectFolder,'/corMatrixPts.R',sep=''))

## ##correlacao (linear) usando os rasters completos na area invadida
## stackedRasters = stack(gohiInviAdjusted, mapInv)
## names(stackedRasters) = c('gohi','suitability')
## valuesMatrix = data.frame(na.omit(values(stackedRasters)))
## ##
## corMatrixRasters = cor(valuesMatrix, use="complete.obs") #resultado: nao correlacionado (linearmente)
##save(corMatrixRasters,file=paste(projectFolder,'/corMatrixPts.R',sep=''))


# ##modelos estatisticos - GLM para suitability X GOHI
# modelSuitXGohi = glm(suitabilityPts ~ gohiPts,  family=binomial(link = "logit"), data=datasetGohi) #modelo linear
# summary(modelSuitXGohi)
# 
# modelSuitXGohiQuad = glm(suitabilityPts ~ gohiPts + I(gohiPts^2), family=binomial(link = "logit"), data=datasetGohi) #modelo quadratico
# summary(modelSuitXGohiQuad)
# 
# aicSuitXgohi =  AIC(modelSuitXGohi, modelSuitXGohiQuad) #comparando por AIC
# anovaSuitXgohi = anova(modelSuitXGohi, modelSuitXGohiQuad, test='Chisq') #comparando por Chisq
# 
# ##visualizacao grafica
# jpeg(filename=paste(projectFolder,'/suitabilityXgohi.jpg',sep=''))
# plot(datasetGohi$suitabilityPts ~ datasetGohi$gohiPts,  xlim=range(datasetGohi$gohiPts), pch=20, xlab='GOHI', ylab='Suitability') 
# points(x = c(-100:100/10),
#        y = predict(modelSuitXGohi, newdata=data.frame(gohiPts=c(-100:100/10)), type='response'),
#        col='red', type='l')
# points(x = c(-100:100/10),
#        y = predict(modelSuitXGohiQuad, newdata=data.frame(gohiPts=c(-100:100/10)), type='response'),
#        col='blue', type='l')
# legend(x=3.3,y=0.99, legend=c('Linear model','Quadratic model'), lty=1, col=c('red','blue'))
# dev.off()


##modelos estatisticos - GLM para presenca/background X GOHI##
modelSuit = glm(pres ~ suitabilityPts, data=datasetGohi, family=binomial) #modelo linear
summary(modelSuit)

modelSuitQuad = glm(pres ~ suitabilityPts + I(suitabilityPts^2), data=datasetGohi, family=binomial) #modelo quadratico
summary(modelSuitQuad)

aicPresXsuit =  AIC(modelSuit, modelSuitQuad) #comparando por AIC
anovaPresXsuit = anova(modelSuit, modelSuitQuad, test='Chisq') #comparando por Chisq


##modelos estatisticos - regressao multipla para occ ~ suitability + GOHI
modelOccXSuitXGohi = glm(pres ~ suitabilityPts + gohiPts, family=binomial(link = "logit"), data=datasetGohi) #modelo linear
summary(modelOccXSuitXGohi)


##visualizacao grafica
jpeg(filename=paste(projectFolder,'/OccXSuitXGohi.jpg',sep=''),  width=800, height=480)
par(mfrow=c(1,2), cex=1.5)
plot(datasetGohi$pres ~ datasetGohi$gohiPts,  xlim=range(datasetGohi$gohiPts), pch=20, xlab='GOHI', ylab='Occurrence') 
points(x = c(-100:100/10),
       y = predict(modelOccXSuitXGohi, newdata=data.frame(gohiPts=c(-100:100/10), suitabilityPts=c(-100:100/10)), type='response'),
       col='red', type='l', lwd=2)
##
plot(datasetGohi$pres ~ datasetGohi$suitabilityPts,  xlim=range(datasetGohi$suitabilityPts), pch=20, xlab='Suitability', ylab='Occurrence') 
points(x = c(-100:100/10),
       y = predict(modelSuitXGohiQuad, newdata=data.frame(gohiPts=c(-100:100/10), suitabilityPts=c(-100:100/10)), type='response'),
       col='red', type='l', lwd=2)
#legend(x=3.3,y=0.99, legend=c('Linear model','Quadratic model'), lty=1, col=c('red','blue'))
dev.off()



##tabela dos resultados dos modelos##

AIC(modelGohi, modelSuit, modelOccXSuitXGohi)
anova(modelSuit, modelGohi, test='Chisq') 
anova(modelSuit, modelOccXSuitXGohi, test='Chisq') 
anova(modelGohi, modelOccXSuitXGohi, test='Chisq') 

##pseudo R-squared para GLM
modelSuitR2 = DescTools::PseudoR2(modelSuit)
modelGohiR2 = DescTools::PseudoR2(modelGohi)
modelOccXSuitXGohiR2 = DescTools::PseudoR2(modelOccXSuitXGohi)



statModelsOutput = data.frame( responseVariable = c('occurrence','occurrence','occurrence'),
                         predictorVariable = c('GOHI','suitability','GOHI+suitability'),
                         degree_freedom = c(modelGohi$df.residual, modelSuit$df.residual, modelOccXSuitXGohi$df.residual),
                         deviance = c(modelGohi$deviance, modelSuit$deviance, modelOccXSuitXGohi$deviance),
                         PseudoR2 = c(modelGohiR2, modelSuitR2, modelOccXSuitXGohiR2),
                         aic = c(modelGohi$aic, modelSuit$aic, modelOccXSuitXGohi$aic) )

write.csv(statModelsOutput, paste(projectFolder,'/statModelsOutput.csv',sep='' ), row.names=FALSE)




## PARTE 5: Analise de pDLA - iSDM  ##NAO FUNCIONOU...  :(




##pacote
library(iSDM)

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

probability = pDLA(occData=occInv,
                   envData=envDataInv,
                   longlat=TRUE,
                   occNative=occNat,
                   envNative=envDataNativ)

##Display the results
par(mfrow=c(1,2),mar=c(2,2.5,2,2.5))
plot(realized.dist$occupied.area,main="Realized distribution")
plot(occData,col=ifelse(occData$SP==1,2,1),add=TRUE,pch=19,cex=0.8)

plot(SDMpredInv[[1]], main="Potential distribution")

scatterCol<-function(x){
  x <- (x-min(x, na.rm=TRUE))/(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))
  colorFunction <- colorRamp(colorRamps::matlab.like(100))
  zMatrix <- colorFunction(x)
  zColors <- rgb(zMatrix[,1], zMatrix[,2], zMatrix[,3], maxColorValue=255)
  return(zColors)
}


points(probability,pch=21, col=1,bg=scatterCol(probability@data[,"PDLA"]),cex=1)