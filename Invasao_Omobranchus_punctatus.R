## Script para analise do nicho de uma especie invasora, comparando nicho do ambiente nativo e do ambiente invadido.
## Anderson A. Eduardo
## novembro/2017

## PARTE 1: carregando pacotes, parâmetros e variaveis

## abrindo pacotes necessarios
library(ecospat)
library(raster)
library(maps)
library(biomod2)

## definindo variaveis e parametros
projectFolder = "/home/anderson/Documentos/Projetos/Invasao_Omobranchus_punctatus" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/Marspec" #pasta com as variaveis ambientais
predictors = stack(list.files(path=envVarFolder, full.names=TRUE, pattern='.tif')) #predictors com todas as variaveis
spData = read.csv(file.path(projectFolder,'spOcc.csv'),header=TRUE) #dados de ocorrencia ambiente nativo
names(spData) = c('lon','lat')


## PARTE 2: definindo dados da area nativa e area invadida

## separando area nativa e area invadida
wrld = rgdal::readOGR('/home/anderson/PosDoc/shapefiles/continents/continent.shp') #mapa mundi
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

## recortando variaveis ambientais para area de distribuicao nativa
predAreaNat = crop(predictors,natAreaExtent) #recortando area selecionada

## ##avaliando a correlacao nas variaveis ambientais
## predAreaNatMatrix = getValues(predAreaNat)
## corTest = as.data.frame(cor(predAreaNatMatrix,use="complete.obs"))
## write.csv(corTest,paste(projectFolder,'/corTestAreaNat.csv',sep=''))

## definindo variaveis preditoras nao correlacionadas 
predAreaNat = predAreaNat[[c('bathy_2o5m',
                             'biogeo01_2o5m',
                             'biogeo02_2o5m',
                             'biogeo05_2o5m',
                             'biogeo06_2o5m',
                             'biogeo07_2o5m',
                             'biogeo08_2o5m',
                             'biogeo11_2o5m',
                             'biogeo13_2o5m',
                             'biogeo16_2o5m')]] #selecionando as variaveis usadas

## recortando variaveis ambientais para area de distribuicao invadida
predAreaInv = crop(predictors,invAreaExtent) #recortando area selecionada

## ##avaliando a correlacao nas variaveis ambientais
## predAreaInvMatrix = getValues(predAreaInv)
## corTest = as.data.frame(cor(predAreaInvMatrix,use="complete.obs"))
## write.csv(corTest,paste(projectFolder,'/corTestAreaInv.csv',sep=''))

## definindo variaveis preditoras nao correlacionadas 
predAreaInv = predAreaInv[[c('bathy_2o5m',
                             'biogeo01_2o5m',
                             'biogeo02_2o5m',
                             'biogeo05_2o5m',
                             'biogeo06_2o5m',
                             'biogeo07_2o5m',
                             'biogeo08_2o5m',
                             'biogeo11_2o5m',
                             'biogeo13_2o5m',
                             'biogeo16_2o5m')]] #selecionando as variaveis usadas


## PARTE 3: Analise de similaridade de equivalencia de nicho

## criando pontos de background
bgAreaNat = dismo::randomPoints(mask=predAreaNat[[1]], n=1000, p=spNat, prob=FALSE)
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

## The PCA is calibrated on all the sites of the study area
pca.env <- dudi.pca(rbind(spNatData,spInvData)[,c('bathy_2o5m',
                                                  'biogeo01_2o5m',
                                                  'biogeo02_2o5m',
                                                  'biogeo05_2o5m',
                                                  'biogeo06_2o5m',
                                                  'biogeo07_2o5m',
                                                  'biogeo08_2o5m',
                                                  'biogeo11_2o5m',
                                                  'biogeo13_2o5m',
                                                  'biogeo16_2o5m')],scannf=F,nf=2)
## ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) #grafico

## PCA scores for the whole study area
scores.globclim <- pca.env$li

## PCA scores for the species native distribution
scores.sp.nat <- suprow(pca.env,spNatData[which(spNatData[,'occ']==1),c('bathy_2o5m',
                                                                         'biogeo01_2o5m',
                                                                         'biogeo02_2o5m',
                                                                         'biogeo05_2o5m',
                                                                         'biogeo06_2o5m',
                                                                         'biogeo07_2o5m',
                                                                         'biogeo08_2o5m',
                                                                         'biogeo11_2o5m',
                                                                         'biogeo13_2o5m',
                                                                         'biogeo16_2o5m')])$li

## PCA scores for the species invasive distribution
scores.sp.inv <- suprow(pca.env,spInvData[which(spInvData[,'occ']==1),c('bathy_2o5m',
                                                                         'biogeo01_2o5m',
                                                                         'biogeo02_2o5m',
                                                                         'biogeo05_2o5m',
                                                                         'biogeo06_2o5m',
                                                                         'biogeo07_2o5m',
                                                                         'biogeo08_2o5m',
                                                                         'biogeo11_2o5m',
                                                                         'biogeo13_2o5m',
                                                                         'biogeo16_2o5m')])$li

## PCA scores for the whole native study area
scores.clim.nat <-suprow(pca.env,spNatData[,c('bathy_2o5m',
                                              'biogeo01_2o5m',
                                              'biogeo02_2o5m',
                                              'biogeo05_2o5m',
                                              'biogeo06_2o5m',
                                              'biogeo07_2o5m',
                                              'biogeo08_2o5m',
                                              'biogeo11_2o5m',
                                              'biogeo13_2o5m',
                                              'biogeo16_2o5m')])$li

## PCA scores for the whole invaded study area
scores.clim.inv <- suprow(pca.env,spInvData[,c('bathy_2o5m',
                                               'biogeo01_2o5m',
                                              'biogeo02_2o5m',
                                              'biogeo05_2o5m',
                                              'biogeo06_2o5m',
                                              'biogeo07_2o5m',
                                              'biogeo08_2o5m',
                                              'biogeo11_2o5m',
                                              'biogeo13_2o5m',
                                              'biogeo16_2o5m')])$li

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

##variaveis e parametros locais especificos para o biomod2 OBS: os dados de 
myRespName = 'Omobranchus_punctatus' #nome
myExpl = predAreaNat  #variavel preditora (para biomod2)
myResp = spNat # pontos de presenca (ja foram criados na PARTE 2)
coordinates(spNat) = ~lon+lat
crs(myResp) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #transformando em spatialPoints

##ajuste de dados de entrada para biomod2
myBiomodData = BIOMOD_FormatingData(resp.var = spNat,
                                     expl.var = myExpl,
                                     resp.name = myRespName,
                                     PA.nb.rep = 1)

## ##inspecionando o objeto gerado pela funcao do biomod2
## myBiomodData
## plot(myBiomodData)

##parametrizando os modelos
myBiomodOption <- BIOMOD_ModelingOptions(
    MAXENT.Phillips=list(
        path_to_maxent.jar='/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java',
        maximumiterations=1000,
        linear=TRUE,
        quadratic=TRUE,
        product=FALSE,
        threshold=FALSE,
        hinge=FALSE,
        maximumiterations=1000,
        convergencethreshold=1.0E-5,
        threads=2))

##rodando o(s) algoritmo(s) (i.e. SDMs)
myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    models = c('MAXENT.Phillips'),
    models.options = myBiomodOption,
    NbRunEval = 10,
    DataSplit = 75,
    VarImport = 10,
    models.eval.meth = c('TSS','ROC'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = FALSE,
    modeling.id = myRespName)

##My output data
evaluationScores = get_evaluations(myBiomodModelOut)

##gravando estatistcas basicas do modelo
statResults = data.frame(sp = 'Omobranchus_punctatus',
                         AUC = mean(evaluationScores['ROC','Testing.data',,,]),
                         TSS = mean(evaluationScores['TSS','Testing.data',,,]))

write.csv(statResults,file=paste(projectFolder,'/maxent/statResults.csv',sep=''),row.names=FALSE)
                                            
##selecionando o melhor modelo para projecao
whichModel = names(evaluationScores['TSS','Testing.data',,,][which(evaluationScores['TSS','Testing.data',,,]== max(evaluationScores['TSS','Testing.data',,,]) )])
modelName = grep(pattern=whichModel, myBiomodModelOut@models.computed, value=TRUE)
                        
##rodando algoritmo de projecao para AREA NATIVA
myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = predAreaNat, #variaveis preditoras para area nativa (feito na PARTE 2)
    proj.name = 'AreaNativa',
    selected.models = modelName,
    binary.meth = 'TSS',
    compress = 'TRUE',
    build.clamping.mask = 'TRUE',
    output.format = '.tif')

##rodando algoritmo de projecao para AREA INVADIDA
myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = predAreaInv, #variaveis preditoras para area invadida (feito na PARTE 2)
    proj.name = 'AreaInvadida',
    selected.models = modelName,
    binary.meth = 'TSS',
    compress = 'TRUE',
    build.clamping.mask = 'TRUE',
    output.format = '.tif')

