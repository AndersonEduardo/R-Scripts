## Script para avaliar sobreposicao de nicho, diretamente no espaco ambiental
## junho/2018


##pacotes
library(raster)
library(ecospat)


## definindo variaveis e parametros - Notebook
projectFolder = "/home/anderson/Projetos/Invasao_Omobranchus_punctatus" #pasta do projeto
envVarFolder = "/home/anderson/Projetos/Invasao_Omobranchus_punctatus/variaveis_ambientais" #pasta com as variaveis ambientais
predictors = stack(list.files(path=envVarFolder, full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis
predictors = predictors[[grep(pattern=paste(c('Chlorophyll','Phytoplankton','Silicate','*.Mean$','*.Max$','*.Min$'),collapse='|'), names(predictors), value=FALSE,invert=TRUE)]]
spData = read.csv(file.path(projectFolder,'spOcc.csv'),header=TRUE) #dados de ocorrencia ambiente nativo
names(spData) = c('lon','lat')
wrld = readOGR('/home/anderson/shapefiles/ne_50m_ocean/ne_50m_ocean.shp') #mapa mundi


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

