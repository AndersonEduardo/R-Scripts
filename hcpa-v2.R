## Este script contem algoritmos para gerar 'n' sps artificiais aleatoriamente (Parte 1), modelar ausencias (Parte 2), modelar distribuicao com pseudoausencias melhoradas (Parte 3) e analise dos resultados atraves de graficos simples (Parte 4)
## Anderson A. Eduardo, outubro/2014


##abrindo pacotes necessarios
#Sys.setenv(JAVA_HOME = "/usr/lib/jvm/java-7-openjdk-amd64")
#options(java.parameters = "Xmx7g")
#library(rJava)
library(raster)
library(biomod2)
library(dismo)

##definindo prametros e variaveis globais (NOTEBOOK)
projectFolder = "/home/anderson/Projetos/HCPA" #pasta do projeto
envVarFolder = "/home/anderson/gridfiles/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("/home/anderson/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
maxentFolder = '/home/anderson/R/x86_64-pc-linux-gnu-library/3.4/dismo/java/maxent.jar' #pasta para resultados do maxent
## spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
## sdmTypes = c('normal','optimized')
sampleSizes = c(10,20,40,80,160)
NumRep = 10 #numero de replicas (de cada cenario amostral)
vies_levels = 5
##variaveis preditoras
## elevation = raster('/home/anderson/PosDoc/dados_ambientais/DEM/DEM.tif')
predictors = stack(list.files(path=envVarPaths[1],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
predictors = predictors[[c('bioclim_01','bioclim_12')]]
predictors = stack(mask(x=predictors, mask=AmSulShape))
crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
SDMreplicates = 10 #numero de replicas para calibracao/validacao dos SDMs
statResultsSDMnormal = data.frame() #tabela de estatisticas basicas do modelo
statResultsSDMimproved = data.frame()

## ##definindo prametros e variaveis globais (LORIEN)
## projectFolder = "J:/Pesquisadorxs/Anderson_Eduardo/high_quality_PA" #pasta do projeto
## envVarFolder = "J:/Pesquisadorxs/Anderson_Eduardo/dados_projeto/000" #pasta com as variaveis ambientais
## envVarPaths = list.files(path=envVarFolder, pattern='.asc', full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
## AmSulShape = rgdal::readOGR("J:/Pesquisadorxs/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
## maxentFolder = 'C:/Users/WS/Documents/R/win-library/3.4/dismo/java/maxent.jar' #pasta para resultados do maxent
## ## spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
## ## sdmTypes = c('normal','optimized')
## sampleSizes = c(10,20,40,80,160)
## NumRep = 10 #numero de replicas (de cada cenario amostral)
## vies_levels = 5
## ##variaveis preditoras
## elevation = raster('J:/Pesquisadorxs/Anderson_Eduardo/DEM/DEM.tif')
## predictors = stack(envVarPaths,elevation) #predictors com todas as variaveis (presente)
## predictors = predictors[[c('bioclim_01','bioclim_12')]]
## predictors = stack(mask(x=predictors, mask=AmSulShape))
## crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
## Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
## SDMreplicates = 10 #numero de replicas para calibracao/validacao dos SDMs
## statResultsSDMnormal = data.frame() #tabela de estatisticas basicas do modelo
## statResultsSDMimproved = data.frame()

# ##definindo prametros e variaveis globais (Yavanna)
# projectFolder = "D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas" #pasta do projeto
# envVarFolder = "D:/Anderson_Eduardo/gridfiles/variaveis ambientais AmSul 120kyr/000" #pasta com as variaveis ambientais
# envVarPaths = list.files(path=envVarFolder, pattern='.asc', full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
# AmSulShape = rgdal::readOGR("D:/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
# maxentFolder = 'D:/Anderson_Eduardo/maxent/maxent.jar' #pasta para resultados do maxent
# sampleSizes = c(20,40,80,160)
# NumRep = 10 #numero de replicas (de cada cenario amostral)
# vies_levels = 5
# ##variaveis preditoras
# ##elevation = raster('J:/Pesquisadorxs/Anderson_Eduardo/DEM/DEM.tif')
# predictors = stack(envVarPaths) #predictors com todas as variaveis (presente)
# predictors = predictors[[c('bioclim_01','bioclim_12')]]
# predictors = stack(mask(x=predictors, mask=AmSulShape))
# crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
# Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
# SDMreplicates = 10 #numero de replicas para calibracao/validacao dos SDMs
# statResultsSDMnormal = data.frame() #tabela de estatisticas basicas do modelo
# statResultsSDMimproved = data.frame()


## funcoes autorais :-) #NOTEBOOK
source('/home/anderson/R-Scripts/makeSpeciesSuitability.R')
source('/home/anderson/R-Scripts/rangeByAC.R')
source('/home/anderson/R-Scripts/makeOutput.R')
source('/home/anderson/R-Scripts/bestModel.R')

## ## funcoes autorais :-) #LORIEN
## source('D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas/makeSpecies.R')
## source('D:/Anderson_Eduardo/rangeByAC.R')
## source('D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas/makeOutput.R')
## source('D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas/bestModel.R')

# ## funcoes autorais :-) #Yavanna
# source('D:/Anderson_Eduardo/R/makeSpecies.R')
# source('D:/Anderson_Eduardo/R/rangeByAC.R')
# #source('D:/Anderson_Eduardo/R/makeOutput.R')
# source('D:/Anderson_Eduardo/R/bestModel.R')




##PARTE 1: CRIANDO AS ESPECIES ARTIFICIAIS




for(i in 1:Nsp){
  
  ##diretorio para o biomod2 salvar resultados para SDMnormal
  setwd(file.path(projectFolder))
  
  ##verifica e cria diretorio para salvar resultados da especie atual
  if (file.exists('virtual species')){
    setwd(paste(projectFolder,'/virtual species',sep=''))
  } else {
    dir.create('virtual species')
    setwd(paste(projectFolder,'/virtual species',sep=''))
  }
  
  SpDistAC = makeSpeciesSuitability(predictors) #criando nicho e projetando sua distribuicao
  SpDistAC = rangeByAC(SpDistAC) #criando a area de distribuicao da sps (usando automato celular)
  
  ##criando imagem da distribuicao de cada especie
  jpeg(filename=paste(projectFolder,'/virtual species/sp',i,'.jpeg',sep=''))
  plot(SpDistAC)
  dev.off()
  ##
  writeRaster(x=SpDistAC, filename=paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''), overwrite=TRUE)
  rm(SpDistAC) ##teste do bug persistente
  
}




##PARTE 2: SIMULANDO OS DATASETS




##arquivo de log
write.table(x=NULL, file=paste(projectFolder,'/logfileSDMnormalDatasets.txt',sep='')) #criando arquivo
cat('Log file - Started at: ', as.character(Sys.time()), "\n \n", file=paste(projectFolder,'/logfileSDMnormalDatasets.txt',sep=''), append=TRUE) #gravando no arquivo

for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    for(current_vies_level in 1:vies_levels){
      tryCatch({
        
        ##logfile
        cat('STARTING SCENARIO: sp = ', i, '/ sampleSizes = ', sampleSizes[j], '/ current_vies_level = ', current_vies_level, ' / time: ', as.character(Sys.time()), "\n \n", file = paste(projectFolder,'/logfileSDMnormalDatasets.txt',sep=''), append = TRUE) #gravando no arquivo
        
        ##diretorio de trabalho do projeto
        setwd(file.path(projectFolder))
        
        ##verifica e cria diretorio para salvar datasets 
        if (file.exists(paste(projectFolder,'/datasets', sep=''))){
          setwd(paste(projectFolder,'/datasets', sep=''))
        } else {
          dir.create(paste(projectFolder,'/datasets',sep=''), recursive=TRUE)
          setwd(paste(projectFolder,'/datasets', sep=''))
        }
        
        ##verifica e cria diretorio para salvar bias layers
        if (!file.exists(paste(projectFolder,'/biasLayers', sep=''))){
          dir.create(paste(projectFolder,'/biasLayers', sep=''))
        }
        
        
        ##definindo variaveis e parametros locais
        ##occPoints = read.csv(paste(mainSampleFolder,sdmTypes[h],'/',spsTypes[i],'/occ',sampleSizes[j],'.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia
        SpDistAC = raster(paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''))
        values(SpDistAC)[values(SpDistAC)==0] = NA
        vies_i = sample(x=sequence(vies_levels), size=current_vies_level) #aqui: escolhendo quais dos centroids do Kmeans sera usado (obs.: 'level' de 'bias' eh a quantidade de centroides)
        
        
        ##bias layer da iteracao atual - sps range
        SpDistAC = raster(paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''))     
        crs(SpDistAC) =  CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
        biasLayerBGpts =  data.frame( dismo::randomPoints(mask=SpDistAC, n=100, prob=TRUE) )
        ptsCluster = kmeans(x=biasLayerBGpts, centers=5, nstart=5)
        biasLayerBGpts$kcluster = ptsCluster$cluster
        biasLayerBGpts = biasLayerBGpts[ sapply(biasLayerBGpts$kcluster, function(x) x  %in% vies_i), ]
        
        
        ##KDE
        bgArea = SpDistAC*0 #predictors[[1]]*0
        crs(bgArea) = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
        biasLayer = data.frame(biasLayerBGpts[,1:2])
        coordinates(biasLayer) = ~x+y
        proj4string(biasLayer) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
        biasLayerBGptsKDE = MASS::kde2d(x=biasLayer$x, y=biasLayer$y, n=100, lims=c(xl=extent(bgArea)@xmin, xu=extent(bgArea)@xmax, yl=extent(bgArea)@ymin, yu=extent(bgArea)@ymax ))
        ##raster
        biasLayerBGptsKDE = raster(biasLayerBGptsKDE)
        biasLayerBGptsKDE = mask(biasLayerBGptsKDE, AmSulShape)
        extent(biasLayerBGptsKDE) = extent(SpDistAC)
        ##ajustando projecao, etc (com cuidado para nao quebrar com o grau de tolerancia)
        biasLayerBGptsKDE = projectRaster(biasLayerBGptsKDE, crs=proj4string(bgArea), res=res(bgArea), method="bilinear")
        tolerance = 0.1
        biasLayerBGptsKDEtry = try( merge(biasLayerBGptsKDE, bgArea, tolerance=tolerance), silent=TRUE )
        while('try-error' %in% class(biasLayerBGptsKDEtry)){
          tolerance = tolerance+0.5
          biasLayerBGptsKDEtry = merge(biasLayerBGptsKDE, bgArea, tolerance=tolerance)
        }
        biasLayerBGptsKDE = biasLayerBGptsKDEtry
        biasLayerBGptsKDE = crop(x=biasLayerBGptsKDE, y=bgArea)
        extent(biasLayerBGptsKDE) = extent(bgArea)
        biasLayerBGptsKDE = biasLayerBGptsKDE/biasLayerBGptsKDE@data@max #ajustando entre 0 e 1
        values(biasLayerBGptsKDE)[which(values(biasLayerBGptsKDE)<=0)] = 0
        ##salvando
        writeRaster(biasLayerBGptsKDE, paste(projectFolder,'/biasLayers/','sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.asc', sep=''), overwrite=TRUE)
        
        ##dados de ocorrencia
        occPoints = dismo::randomPoints(mask=SpDistAC*biasLayerBGptsKDE, n=sampleSizes[j], prob=TRUE) #sorteando pontos da distribuicao modelada
        ##rm(SpDistAC) ##teste do bug persistente
        occPoints = data.frame(lon=occPoints[,1],lat=occPoints[,2])
        write.csv(occPoints, paste(projectFolder,'/datasets/','occ_sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.csv', sep=''), row.names = FALSE)
        
        ##backgrownd points
        bgPoints = dismo::randomPoints(mask=biasLayerBGptsKDE, p=occPoints, n=10000, prob=TRUE) #sorteando pontos da distribuicao modelada
        bgPoints = data.frame(lon=bgPoints[,1], lat=bgPoints[,2])
        write.csv(bgPoints, paste(projectFolder,'/datasets/','bg_sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.csv', sep=''), row.names = FALSE)
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n \n", file=paste(projectFolder,'/logfileSDMnormalDatasets.txt',sep=''), append=TRUE)})
      
    }
  }
}   




##PARTE 3: SDMs COM BIOMOD2




##arquivo de log
write.table(x=NULL, file=paste(projectFolder,'/logfileSDMnormalBiomod2.txt',sep='')) #criando arquivo
cat('Log file - Started at: ', as.character(Sys.time()), "\n \n", file=paste(projectFolder,'/logfileSDMnormalBiomod2.txt',sep=''), append=TRUE) #gravando no arquivo

for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    for(current_vies_level in 1:vies_levels){
      tryCatch({
        
        ##logfile
        cat('STARTING SCENARIO: sp = ', i, '/ sampleSizes = ', sampleSizes[j], '/ current_vies_level = ', current_vies_level, ' / time: ', as.character(Sys.time()), "\n \n", file = paste(projectFolder,'/logfileSDMnormalBiomod2.txt',sep=''), append = TRUE) #gravando no arquivo
        
        ##diretorio de trabalho do projeto
        setwd(file.path(projectFolder))
        
        ##verifica e cria diretorio para salvar datasets 
        if (file.exists(paste(projectFolder,'/SDMnormal', sep=''))){
          setwd(paste(projectFolder,'/SDMnormal', sep=''))
        } else {
          dir.create(paste(projectFolder,'/SDMnormal',sep=''), recursive=TRUE)
          setwd(paste(projectFolder,'/SDMnormal', sep=''))
        }
        
        ##verifica e cria diretorio para salvar bias layers
        if (!file.exists(paste(projectFolder,'/biasLayers', sep=''))){
          dir.create(paste(projectFolder,'/biasLayers', sep=''))
        }
        
        ##dataset
        occPoints = read.csv(paste(projectFolder,'/datasets/','occ_sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.csv', sep=''))
        occPoints$pres = 1
        bgPoints =  read.csv(paste(projectFolder,'/datasets/','bg_sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.csv', sep=''))
        bgPoints$pres = 0
        dataSet = rbind(occPoints, bgPoints)
        
        ##preditoras
        
        #####################################
        ## FAZER SELECAO DE VARIAVEIS AQUI ##
        #####################################
        
        # myResp <- data.frame(lon=occPoints[,1], lat=occPoints[,2])
        # coordinates(myResp) <- ~ lon + lat #transformando em spatialPoints
        # crs(myResp) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #transformando em spatialPoints
        
        ##variaveis e parametros locais especificos para o biomod2
        myRespName <- paste('sp',i,sep='')  # nome do cenario atual (para biomod2)
        myResp <- dataSet[,c('pres')] # variavel resposta (para biomod2)
        myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
        myExpl = extract(predictors, dataSet[,c('lon','lat')]) #variavel preditora (para biomod2)
        
        ##ajuste de dados de entrada para biomod2
        # myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
        #                                      expl.var = myExpl,
        #                                      resp.name = paste(myRespName,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal',sep=''),
        #                                      PA.nb.rep = 0,
        #                                      PA.nb.absences = 10000)
        
        myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                             expl.var = myExpl,
                                             resp.xy = myRespXY,
                                             resp.name = paste(myRespName,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal',sep=''))
        
        ##parametrizando os modelos
        myBiomodOption <- BIOMOD_ModelingOptions(
          MAXENT.Phillips = list(path_to_maxent.jar = maxentFolder,
                                 jackknife = FALSE,
                                 randomseed = FALSE,
                                 randomtestpoints = 0,
                                 betamultiplier = 1,
                                 writebackgroundpredictions = FALSE,
                                 linear = TRUE,
                                 quadratic = TRUE,
                                 product = FALSE,
                                 threshold = FALSE,
                                 hinge = FALSE,
                                 extrapolate = FALSE,
                                 doclamp = FALSE,
                                 fadebyclamping = FALSE,
                                 threads = 2),
          GLM = list( type = 'quadratic',
                      interaction.level = 0,
                      myFormula = NULL,
                      test = 'AIC',
                      family = binomial(link = 'logit'),
                      mustart = 0.5,
                      control = glm.control(epsilon = 1e-07, maxit = 100, trace = FALSE)),
          GAM = list( algo = 'GAM_mgcv',
                      type = 's_smoother',
                      k = -1,
                      interaction.level = 0,
                      myFormula = NULL,
                      family = binomial(link = 'logit'),
                      method = 'GCV.Cp',
                      optimizer = c('outer','newton'),
                      select = FALSE,
                      knots = NULL,
                      paraPen = NULL,
                      control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07,
                                     maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15,
                                     rank.tol = 1.49011611938477e-08,
                                     nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0),
                                     optim = list(factr=1e+07),
                                     newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0),
                                     outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE,
                                     scale.est = 'fletcher', edge.correct = FALSE)),
          MARS = list( type = 'simple',
                       interaction.level = 0,
                       myFormula = NULL,
                       nk = NULL,
                       penalty = 2,
                       thresh = 0.001,
                       nprune = NULL,
                       pmethod = 'backward'),
          CTA = list( method = 'class',
                      parms = 'default',
                      cost = NULL,
                      control = list(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 25)),
          GBM = list( distribution = 'bernoulli',
                      n.trees = 2500,
                      interaction.depth = 7,
                      n.minobsinnode = 5,
                      shrinkage = 0.001,
                      bag.fraction = 0.5,
                      train.fraction = 1,
                      cv.folds = 3,
                      keep.data = FALSE,
                      verbose = FALSE,
                      perf.method = 'cv'),
          RF = list( do.classif = TRUE,
                     ntree = 500,
                     mtry = 'default',
                     nodesize = 5,
                     maxnodes = NULL)
        )                        
        
        ##rodando o(s) algoritmo(s) (i.e. SDMs)
        myBiomodModelOut <- BIOMOD_Modeling(
          data = myBiomodData,
          models = c('MAXENT.Phillips','GLM','GAM','MARS','CTA','GBM','RF'),
          models.options = myBiomodOption,
          NbRunEval = SDMreplicates,
          DataSplit = 80,
          models.eval.meth = c('TSS','ROC'),
          do.full.models = FALSE,
          modeling.id = paste(myRespName,'_sample',sampleSizes[j],'_SDMnormal',sep=''))
        
        ##My output data
        evaluationScores = get_evaluations(myBiomodModelOut, as.data.frame=TRUE)
        
        ##outputs 
        statResultsSDMnormal = rbind(statResultsSDMnormal,
                                     data.frame(SDM='normal', sp=paste('sp',i,sep=''), sampleSize=sampleSizes[j], biaslevel=current_vies_level, evaluationScores))
        write.csv(statResultsSDMnormal, file=paste(projectFolder,'/SDMnormal/StatisticalResults_SDMnormal.csv',sep=''), row.names=FALSE)
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n \n", file=paste(projectFolder,'/logfileSDMnormalBiomod2.txt',sep=''), append=TRUE)})
    }
  }
}      




##PARTE 4: IMPLEMENTANDO SDMs COM OS MODELOS MAIS SIMPLES (MAHALANOBIS, BIOCLIM, DOMAIN)




##arquivo de log
write.table(x=NULL, file=paste(projectFolder,'/logfileSDMnormalOutrosModelos.txt',sep='')) #criando arquivo
cat('Log file - Started at: ', as.character(Sys.time()), "\n \n", file=paste(projectFolder,'/logfileSDMnormalOutrosModelos.txt',sep=''), append=TRUE) #gravando no arquivo

##data.frame para ajudar a computar os outputs
evaluationScores = data.frame()

for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    for(current_vies_level in 1:vies_levels){
      tryCatch({
        
        ##logfile
        cat('STARTING SCENARIO: sp = ', i, '/ sampleSizes = ', sampleSizes[j], '/ current_vies_level = ', current_vies_level, ' / time: ', as.character(Sys.time()), "\n \n", file=paste(projectFolder,'/logfileSDMnormalOutrosModelos.txt',sep=''), append = TRUE) #gravando no arquivo
        
        ##diretorio para o biomod2 salvar resultados para SDMnormal
        setwd(file.path(projectFolder,'SDMnormal'))
        
        ##dataset
        occPoints = read.csv(paste(projectFolder,'/datasets/','occ_sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.csv', sep=''))
        bgPoints =  read.csv(paste(projectFolder,'/datasets/','bg_sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.csv', sep=''))
        
        ##preditoras
        
        #################################################################
        ## ABRIR AS VARIAVEIS JA SELECIONAAS E SALVAS NA ULTIMA ETAPA  ##
        #################################################################
        
        ##iterando replicas (ATENCAO: tem que ser o mesmo numero que no biomod2)
        SDMlist = list()
        
        for (iter in 1:SDMreplicates){
          
          trainIdxOcc = kfold(x = occPoints,k = 3) != 1
          testIdxOcc = !trainIdxOcc
          trainIdxBg = kfold(x = bgPoints,k = 3) != 1
          testIdxBg = !trainIdxBg
          curent_occData_train = subset( x=occPoints, subset=trainIdxOcc) #occ para teste
          curent_bgData_train = subset( x=bgPoints, subset=trainIdxBg) #bg para train
          curent_occData_test = subset( x=occPoints, subset=testIdxOcc) #occ para teste #occ para teste
          curent_bgData_test = subset( x=bgPoints, subset=testIdxBg) #bg para teste
          
          ##rodando Mahalanobis, Bioclim e Domain
          mahalanobisSDM = mahal( x=predictors, p=curent_occData_train )
          bioclimSDM = bioclim( x=predictors, p=curent_occData_train )
          domainSDM = domain( x=predictors, p=curent_occData_train )
          
          ##lista de objetos com os modelos
          SDMlist = list('mahalanobis'=mahalanobisSDM,'bioclim'=bioclimSDM,'domain'=domainSDM)
          names(SDMlist) = c(paste('mahalanobis_RUN',iter,sep=''), 
                             paste('bioclim_RUN',iter,sep=''), 
                             paste('domain_RUN',iter,sep=''))
          
          ##salvando os modelos na memoria fisica do PC
          assign(paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal_AllData_RUN',iter,'_mahalanobis', sep=''), mahalanobisSDM)
          save(list = paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal_AllData_RUN',iter,'_mahalanobis', sep=''), 
               file=paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal/models/sp',i,'_sample',sampleSizes[j],'_SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal_AllData_RUN',iter,'_mahalanobis.RData',sep=''))
          
          assign(paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal_AllData_RUN',iter,'_bioclim', sep=''), bioclimSDM)
          save(list = paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal_AllData_RUN',iter,'_bioclim', sep=''), 
               file=paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal/models/sp',i,'_sample',sampleSizes[j],'_SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal_AllData_RUN',iter,'_bioclim.RData',sep=''))          
          
          assign(paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal_AllData_RUN',iter,'_domain', sep=''), domainSDM)
          save(list = paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal_AllData_RUN',iter,'_domain', sep=''), 
               file=paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal/models/sp',i,'_sample',sampleSizes[j],'_SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal_AllData_RUN',iter,'_domain.RData',sep=''))          
          
          ##computando metricas dos modelos
          evaluationScoresList = lapply( seq(length(SDMlist)), function(i){
            SDMeval = evaluate(model = SDMlist[[i]],
                               p = curent_occData_test,
                               a = curent_bgData_test,
                               x = predictors) #algoritmo de avaliacao do modelo
            
            ##TSS
            tssVector = SDMeval@TPR + SDMeval@TNR - 1
            tssVal = max( SDMeval@TPR + SDMeval@TNR - 1 )
            ## Obs1: formula do TSS = true positive ratio + true negative ratio - 1
            ## FONTE: ALOUCHE & KADMON (2006), Journal of Applied Ecology.
            ## Obs2: outra forma de fazer a mesma coisa seria:
            ## tssVal = tssVector[ which(SDMeval@t == threshold(SDMeval, 'spec_sens')) ]
            ##threshold (maximizando o TSS)
            current_thre_maximizingTSS = threshold(SDMeval, 'spec_sens')
            ## Obs: outra forma de fazer a mesma coisa seria:
            ## tssVal = SDMeval@t[which(tssVector==max(tssVector))]
            ##AUC
            aucVal = SDMeval@auc
            current_thre_maximizingROC = max( SDMeval@t[which(
              sqrt( (1-SDMeval@TPR)^2 + (1-SDMeval@TNR)^2 ) == min(sqrt( (1-SDMeval@TPR)^2 + (1-SDMeval@TNR)^2 ))
            )] )
            ## Obs1: esse eh o procedimento realizado pelo Biomod. FONTE: http://lists.r-forge.r-project.org/pipermail/biomod-commits/2010-August/000129.html
            ## Obs2: Formula pitagorica retirada de Jimenez-Valverde (2012) Global Ecology and Biogreography.
            
            ##tabela com as metricas calculadas
            evaluationScoresOut = data.frame(Model.name = rep(paste(names(SDMlist)[i],'_AllData',sep=''), 2),
                                             Eval.metric = c('TSS','ROC'),
                                             Testing.data = c(tssVal, aucVal),
                                             Evaluating.data = c(NA,NA),
                                             Cutoff = c(current_thre_maximizingTSS*1000, current_thre_maximizingROC*1000),
                                             Sensitivity = c(SDMeval@TPR[which(SDMeval@t==current_thre_maximizingTSS)]*100, SDMeval@TPR[which(SDMeval@t==current_thre_maximizingROC)]*100),
                                             Specificity = c(SDMeval@TNR[which(SDMeval@t==current_thre_maximizingTSS)]*100, SDMeval@TNR[which(SDMeval@t==current_thre_maximizingROC)]*100))
            return(evaluationScoresOut)
          })
          
          evaluationScoresList = do.call('rbind', evaluationScoresList)
          evaluationScores = rbind(evaluationScores, evaluationScoresList)
          
        }
        
        ##outputs 
        statResultsSDMnormal = read.csv(paste(projectFolder,'/SDMnormal/StatisticalResults_SDMnormal.csv',sep=''), header = TRUE)
        statResultsSDMnormal = rbind(statResultsSDMnormal,
                                     data.frame(SDM='normal', sp=paste('sp',i,sep=''), sampleSize=sampleSizes[j], biaslevel=current_vies_level, evaluationScores))
        write.csv(statResultsSDMnormal, file=paste(projectFolder,'/SDMnormal/StatisticalResults_SDMnormal.csv',sep=''), row.names=FALSE)
        
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n \n", file=paste(projectFolder,'/logfileSDMnormalOutrosModelos.txt',sep=''), append=TRUE)})
    }
  }
}




##PARTE 5: SELECIONANDO MELHOR MODELO E IMPLEMENTANDO PROJECOES DOS SDMs




##arquivo de log
write.table(x=NULL, file=paste(projectFolder,'/logfileSDMnormalProjecoes.txt',sep='')) #criando arquivo
cat('Log file - Started at: ', as.character(Sys.time()), "\n \n", file=paste(projectFolder,'/logfileSDMnormalProjecoes.txt',sep=''), append=TRUE) #gravando no arquivo

for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    for(current_vies_level in 1:vies_levels){
      tryCatch({
        
        ##logfile
        cat('STARTING SCENARIO: sp = ', i, '/ sampleSizes = ', sampleSizes[j], '/ current_vies_level = ', current_vies_level, ' / time: ', as.character(Sys.time()),  "\n \n", file = paste(projectFolder,'/logfileSDMnormalProjecoes.txt',sep=''), append = TRUE) #gravando no arquivo
        
        ##diretorio para o biomod2 salvar resultados para SDMnormal
        setwd(file.path(projectFolder,'SDMnormal'))
        
        ##abrindo dados dos modelos
        statResultsSDMnormal = read.csv(paste(projectFolder,'/SDMnormal/StatisticalResults_SDMnormal.csv',sep=''), header = TRUE)
        evaluationScores = subset(statResultsSDMnormal, (statResultsSDMnormal$SDM == 'normal') & 
                                    (statResultsSDMnormal$sp == paste('sp',i,sep='')) & 
                                    (statResultsSDMnormal$sampleSize == sampleSizes[j]) & 
                                    (statResultsSDMnormal$biaslevel == current_vies_level))
        
        ##selecao do modelo de maior sensibilidade            
        bestModel = evaluationScores[which(evaluationScores$Specificity == max(evaluationScores$Specificity, na.rm = TRUE)),]
        modelNames = as.character(bestModel$Model.name)
        modelNames = strsplit(x = modelNames, split = "_")
        modelNames = lapply(seq(length(modelNames)), function(i) rev(modelNames[[i]]) )
        modelNames = lapply(seq(length(modelNames)), function(i) paste(modelNames[[i]], collapse = '_', sep='') )
        modelNames = unlist(modelNames)
        
        ##abrindo os modelos
        modelFiles = list.files(paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal/models/sp',i,'_sample',sampleSizes[j],'_SDMnormal',sep=''), pattern = paste(modelNames,collapse = '|'), full.names = TRUE)
        lapply(modelFiles, load, .GlobalEnv)
        
        ##implementado projecoes
        modelFiles = basename(modelFiles)
        if (sum(grep(pattern = '.RData', x = modelFiles)) > 0){modelFiles = gsub(pattern = '.RData', replacement = '', x = modelFiles)}
        modelNames = lapply(modelFiles, get)
        modelProj = lapply(seq(length(modelNames)), function(i)  predict(modelNames[[i]], predictors))
        
        ##salvando
        if (length(modelProj)==1){
          modelProj = modelProj[[1]]
          names(modelProj) = modelFiles
          crs(modelProj) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
          writeRaster(modelProj, bylayer=FALSE, filename=paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal/models/', modelFiles, sep=''), format='ascii', overwrite=TRUE)
        }else{
          modelProj = stack(modelProj)
          names(modelProj) = modelFiles
          crs(modelProj) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')  
          writeRaster(modelProj, bylayer=TRUE, filename=paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal/models/', modelFiles, sep=''), format='ascii', overwrite=TRUE)
        }
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n \n", file=paste(projectFolder,'/logfileSDMnormalProjecoes.txt',sep=''), append=TRUE)})
    }
  }
}




##PARTE 6: implementando datasets com pseudoausencias melhoradas (HCPA)




##arquivo de log
write.table(x=NULL, file=paste(projectFolder,'/logfileSDMimprovedDatasets.txt',sep='')) #criando arquivo
cat('Log file - Started at: ', as.character(Sys.time()), "\n \n", file=paste(projectFolder,'/logfileSDMimprovedDatasets.txt',sep=''), append=TRUE) #gravando no arquivo


for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    for(current_vies_level in 1:vies_levels){
      tryCatch({
        
        ##logfile
        cat('STARTING SCENARIO: sp = ', i, '/ sampleSizes = ', sampleSizes[j], '/ current_vies_level = ', current_vies_level, ' / time: ', as.character(Sys.time()),  "\n \n", file = paste(projectFolder,'/logfileSDMimprovedDatasets.txt',sep=''), append = TRUE) #gravando no arquivo
                
        ##diretorio para o biomod2 salvar resultados para SDMnormal
        setwd(file.path(projectFolder))
        
        ##verifica e cria diretorio 
        if (file.exists('SDMimproved')){
          setwd(paste(projectFolder,'/SDMimproved',sep=''))
        } else {
          dir.create('SDMimproved')
          setwd(paste(projectFolder,'/SDMimproved',sep=''))
        }
        
        ##abrindo dados de ocorrecias
        occPoints = read.csv(paste(projectFolder,'/datasets/','occ_sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.csv', sep=''))
        
        ##abrindo dados dos modelos - SDM normal
        statResultsSDMnormal = read.csv(paste(projectFolder,'/SDMnormal/StatisticalResults_SDMnormal.csv',sep=''), header = TRUE)
        evaluationScores = subset(statResultsSDMnormal, (statResultsSDMnormal$SDM == 'normal') & 
                                    (statResultsSDMnormal$sp == paste('sp',i,sep='')) & 
                                    (statResultsSDMnormal$sampleSize == sampleSizes[j]) & 
                                    (statResultsSDMnormal$biaslevel == current_vies_level))
        
        ##selecao do modelo de maior sensibilidade            
        bestModel = evaluationScores[which(evaluationScores$Specificity == max(evaluationScores$Specificity, na.rm = TRUE)),]
        modelNames = as.character(bestModel$Model.name)
        modelNames = strsplit(x = modelNames, split = "_")
        modelNames = lapply(seq(length(modelNames)), function(i) rev(modelNames[[i]]) )
        modelNames = lapply(seq(length(modelNames)), function(i) paste(modelNames[[i]], collapse = '_', sep='') )
        modelNames = unlist(modelNames)
        
        ##abrindo as projecoes dos modelos - SDM normal
        modelFiles = list.files(paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal/models',sep=''), pattern = paste(modelNames,collapse = '|'), full.names = TRUE)
        SDMnormalProjs = stack(modelFiles)
        crs(SDMnormalProjs) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
        #SDMnormalProjs = calc(SDMnormalProjs, mean, na.rm=TRUE) ## obs: usando media aritmetica simples das projecoes dos modelos
        #pegando os nomes das camadas pra garantir consistencia com tabela de melhores modelos
        layerNames = sub(".*normal_", "", names(SDMnormalProjs))
        layerNames = strsplit(x = layerNames, split = "_")
        layerNames = lapply(seq(length(layerNames)), function(i) rev(layerNames[[i]]) )
        layerNames = lapply(seq(length(layerNames)), function(i) paste(layerNames[[i]], collapse = '_', sep='') )
        layerNames = unlist(layerNames)
        idxTable = match( layerNames,  bestModel$Model.name)
        SDMnormalProjs = SDMnormalProjs > (bestModel[idxTable,'Cutoff']/1000.0)
        if (length(modelFiles)>1){
          SDMnormalProjs = sum(SDMnormalProjs, na.rm = TRUE)  
        }
        
        #deixando apenas areas de ausencia no mapa
        values(SDMnormalProjs)[values(SDMnormalProjs) != 0] = NA
        
        #################################################
        ##obtendo as pseudoausencias melhoradas - HCPA!##
        hcpa = randomPoints(mask=SDMnormalProjs, n=10000, p=occPoints)
        hcpa = data.frame(lon=hcpa[,1], lat=hcpa[,2])
        #################################################
        
        ##salvando HCPA dataset
        write.csv(hcpa, paste(projectFolder,'/datasets/','hcpa_sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.csv', sep=''), row.names = FALSE)

      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n \n", file=paste(projectFolder,'/logfileSDMimprovedDatasets.txt',sep=''), append=TRUE)})
    }
  }
}




##PARTE 7: SDMs USANDO HCPA COM BIOMOD2 - SDMimproved




##arquivo de log
write.table(x=NULL, file=paste(projectFolder,'/logfileSDMimprovedBiomod2.txt',sep='')) #criando arquivo
cat('Log file - Started at: ', as.character(Sys.time()), "\n \n", file=paste(projectFolder,'/logfileSDMimprovedBiomod2.txt',sep=''), append=TRUE) #gravando no arquivo

for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    for(current_vies_level in 1:vies_levels){
      tryCatch({
        
        ##logfile
        cat('STARTING SCENARIO: sp = ', i, '/ sampleSizes = ', sampleSizes[j], '/ current_vies_level = ', current_vies_level, ' / time: ', as.character(Sys.time()),  "\n \n", file=paste(projectFolder,'/logfileSDMimprovedBiomod2.txt',sep=''), append = TRUE) #gravando no arquivo
        
        ##diretorio para o biomod2 salvar resultados para SDMnormal
        ##verifica e cria diretorio para salvar resultados da especie atual
        setwd(projectFolder)
        
        if (file.exists('SDMimproved')){
          setwd(paste(projectFolder,'/SDMimproved',sep=''))
        } else {
          dir.create('SDMimproved')
          setwd(paste(projectFolder,'/SDMimproved',sep=''))
        }
        
        ##dataset
        occPoints = read.csv(paste(projectFolder,'/datasets/','occ_sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.csv', sep=''))
        occPoints$pres = 1
        bgPoints = read.csv(paste(projectFolder,'/datasets/','hcpa_sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.csv', sep=''))
        bgPoints$pres = 0
        dataSet = rbind(occPoints, bgPoints)
        
        ##preditoras
        
        #####################################
        ## FAZER SELECAO DE VARIAVEIS AQUI ##
        #####################################
        
        # myResp <- data.frame(lon=occPoints[,1], lat=occPoints[,2])
        # coordinates(myResp) <- ~ lon + lat #transformando em spatialPoints
        # crs(myResp) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #transformando em spatialPoints
        
        ##variaveis e parametros locais especificos para o biomod2
        myRespName <- paste('sp',i,sep='')  # nome do cenario atual (para biomod2)
        myResp <- dataSet[,c('pres')] # variavel resposta (para biomod2)
        myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
        myExpl = extract(predictors, dataSet[,c('lon','lat')]) #variavel preditora (para biomod2)
        
        ##ajuste de dados de entrada para biomod2
        # myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
        #                                      expl.var = myExpl,
        #                                      resp.name = paste(myRespName,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal',sep=''),
        #                                      PA.nb.rep = 0,
        #                                      PA.nb.absences = 10000)
        
        myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                             expl.var = myExpl,
                                             resp.xy = myRespXY,
                                             resp.name = paste(myRespName,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMimproved',sep=''))
        
        ##parametrizando os modelos
        myBiomodOption <- BIOMOD_ModelingOptions(
          MAXENT.Phillips = list(path_to_maxent.jar = maxentFolder,
                                 jackknife = FALSE,
                                 randomseed = FALSE,
                                 randomtestpoints = 0,
                                 betamultiplier = 1,
                                 writebackgroundpredictions = FALSE,
                                 linear = TRUE,
                                 quadratic = TRUE,
                                 product = FALSE,
                                 threshold = FALSE,
                                 hinge = FALSE,
                                 extrapolate = FALSE,
                                 doclamp = FALSE,
                                 fadebyclamping = FALSE,
                                 threads = 2),
          GLM = list( type = 'quadratic',
                      interaction.level = 0,
                      myFormula = NULL,
                      test = 'AIC',
                      family = binomial(link = 'logit'),
                      mustart = 0.5,
                      control = glm.control(epsilon = 1e-07, maxit = 100, trace = FALSE)),
          GAM = list( algo = 'GAM_mgcv',
                      type = 's_smoother',
                      k = -1,
                      interaction.level = 0,
                      myFormula = NULL,
                      family = binomial(link = 'logit'),
                      method = 'GCV.Cp',
                      optimizer = c('outer','newton'),
                      select = FALSE,
                      knots = NULL,
                      paraPen = NULL,
                      control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07,
                                     maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15,
                                     rank.tol = 1.49011611938477e-08,
                                     nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0),
                                     optim = list(factr=1e+07),
                                     newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0),
                                     outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE,
                                     scale.est = 'fletcher', edge.correct = FALSE)),
          MARS = list( type = 'simple',
                       interaction.level = 0,
                       myFormula = NULL,
                       nk = NULL,
                       penalty = 2,
                       thresh = 0.001,
                       nprune = NULL,
                       pmethod = 'backward'),
          CTA = list( method = 'class',
                      parms = 'default',
                      cost = NULL,
                      control = list(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 25)),
          GBM = list( distribution = 'bernoulli',
                      n.trees = 2500,
                      interaction.depth = 7,
                      n.minobsinnode = 5,
                      shrinkage = 0.001,
                      bag.fraction = 0.5,
                      train.fraction = 1,
                      cv.folds = 3,
                      keep.data = FALSE,
                      verbose = FALSE,
                      perf.method = 'cv'),
          RF = list( do.classif = TRUE,
                     ntree = 500,
                     mtry = 'default',
                     nodesize = 5,
                     maxnodes = NULL)
        )                        
        
        ##rodando o(s) algoritmo(s) (i.e. SDMs)
        myBiomodModelOut <- BIOMOD_Modeling(
          data = myBiomodData,
          models = c('MAXENT.Phillips','GLM','GAM','MARS','CTA','GBM','RF'),
          models.options = myBiomodOption,
          NbRunEval = SDMreplicates,
          DataSplit = 80,
          models.eval.meth = c('TSS','ROC'),
          do.full.models = FALSE,
          modeling.id = paste(myRespName,'_sample',sampleSizes[j],'_SDMimproved',sep=''))
        
        ##My output data
        evaluationScores = get_evaluations(myBiomodModelOut, as.data.frame=TRUE)
        
        ##outputs 
        statResultsSDMimproved = rbind(statResultsSDMimproved,
                                     data.frame(SDM='improved', sp=paste('sp',i,sep=''), sampleSize=sampleSizes[j], biaslevel=current_vies_level, evaluationScores))
        write.csv(statResultsSDMimproved, file=paste(projectFolder,'/SDMimproved/StatisticalResults_SDMimproved.csv',sep=''), row.names=FALSE)
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n \n", file=paste(projectFolder,'/logfileSDMimprovedBiomod2.txt',sep=''), append=TRUE)})
    }
  }
}      




##PARTE 8: IMPLEMENTANDO SDMs COM OS MODELOS MAIS SIMPLES (MAHALANOBIS, BIOCLIM, DOMAIN)




##arquivo de log
write.table(x=NULL, file=paste(projectFolder,'/logfileSDMimprovedOutrosModelos.txt',sep='')) #criando arquivo
cat('Log file - Started at: ', as.character(Sys.time()), "\n \n", file=paste(projectFolder,'/logfileSDMimprovedOutrosModelos.txt',sep=''), append=TRUE) #gravando no arquivo

##data.frame para ajudar a computar os outputs
evaluationScores = data.frame()

for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    for(current_vies_level in 1:vies_levels){
      tryCatch({
        
        ##logfile
        cat('STARTING SCENARIO: sp = ', i, '/ sampleSizes = ', sampleSizes[j], '/ current_vies_level = ', current_vies_level, ' / time: ', as.character(Sys.time()), "\n \n", file=paste(projectFolder,'/logfileSDMimprovedOutrosModelos.txt',sep=''), append = TRUE) #gravando no arquivo
        
        ##diretorio
        setwd(paste(projectFolder,'/SDMimproved',sep=''))
        
        ##dataset
        occPoints = read.csv(paste(projectFolder,'/datasets/','occ_sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.csv', sep=''))
        bgPoints =  read.csv(paste(projectFolder,'/datasets/','bg_sp',i,'_sampleSizes',sampleSizes[j],'_biasLevel',current_vies_level,'.csv', sep=''))
        
        ##preditoras
        
        #################################################################
        ## ABRIR AS VARIAVEIS JA SELECIONAAS E SALVAS NA ULTIMA ETAPA  ##
        #################################################################
        
        ##iterando replicas (ATENCAO: tem que ser o mesmo numero que no biomod2)
        SDMlist = list()
        
        for (iter in 1:SDMreplicates){
          
          trainIdxOcc = kfold(x = occPoints,k = 3) != 1
          testIdxOcc = !trainIdxOcc
          trainIdxBg = kfold(x = bgPoints,k = 3) != 1
          testIdxBg = !trainIdxBg
          curent_occData_train = subset( x=occPoints, subset=trainIdxOcc) #occ para teste
          curent_bgData_train = subset( x=bgPoints, subset=trainIdxBg) #bg para train
          curent_occData_test = subset( x=occPoints, subset=testIdxOcc) #occ para teste #occ para teste
          curent_bgData_test = subset( x=bgPoints, subset=testIdxBg) #bg para teste
          
          ##rodando Mahalanobis, Bioclim e Domain
          mahalanobisSDM = mahal( x=predictors, p=curent_occData_train )
          bioclimSDM = bioclim( x=predictors, p=curent_occData_train )
          domainSDM = domain( x=predictors, p=curent_occData_train )
          
          ##lista de objetos com os modelos
          SDMlist = list('mahalanobis'=mahalanobisSDM,'bioclim'=bioclimSDM,'domain'=domainSDM)
          names(SDMlist) = c(paste('mahalanobis_RUN',iter,sep=''), 
                             paste('bioclim_RUN',iter,sep=''), 
                             paste('domain_RUN',iter,sep=''))
          
          ##salvando os modelos na memoria fisica do PC
          assign(paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved_AllData_RUN',iter,'_mahalanobis', sep=''), mahalanobisSDM)
          save(list = paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved_AllData_RUN',iter,'_mahalanobis', sep=''), 
               file = paste(projectFolder,'/SDMimproved/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved/models/sp',i,'_sample',sampleSizes[j],'_SDMimproved/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved_AllData_RUN',iter,'_mahalanobis.RData',sep=''))
          
          assign(paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved_AllData_RUN',iter,'_bioclim', sep=''), bioclimSDM)
          save(list = paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved_AllData_RUN',iter,'_bioclim', sep=''), 
               file=paste(projectFolder,'/SDMimproved/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved/models/sp',i,'_sample',sampleSizes[j],'_SDMimproved/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved_AllData_RUN',iter,'_bioclim.RData',sep=''))          
          
          assign(paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved_AllData_RUN',iter,'_domain', sep=''), domainSDM)
          save(list = paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved_AllData_RUN',iter,'_domain', sep=''), 
               file=paste(projectFolder,'/SDMimproved/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved/models/sp',i,'_sample',sampleSizes[j],'_SDMimproved/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved_AllData_RUN',iter,'_domain.RData',sep=''))          
          
          
          
          
          ##computando metricas dos modelos
          evaluationScoresList = lapply( seq(length(SDMlist)), function(i){
            SDMeval = evaluate(model = SDMlist[[i]],
                               p = curent_occData_test,
                               a = curent_bgData_test,
                               x = predictors) #algoritmo de avaliacao do modelo
            
            ##TSS
            tssVector = SDMeval@TPR + SDMeval@TNR - 1
            tssVal = max( SDMeval@TPR + SDMeval@TNR - 1 )
            ## Obs1: formula do TSS = true positive ratio + true negative ratio - 1
            ## FONTE: ALOUCHE & KADMON (2006), Journal of Applied Ecology.
            ## Obs2: outra forma de fazer a mesma coisa seria:
            ## tssVal = tssVector[ which(SDMeval@t == threshold(SDMeval, 'spec_sens')) ]
            ##threshold (maximizando o TSS)
            current_thre_maximizingTSS = threshold(SDMeval, 'spec_sens')
            ## Obs: outra forma de fazer a mesma coisa seria:
            ## tssVal = SDMeval@t[which(tssVector==max(tssVector))]
            ##AUC
            aucVal = SDMeval@auc
            current_thre_maximizingROC = max( SDMeval@t[which(
              sqrt( (1-SDMeval@TPR)^2 + (1-SDMeval@TNR)^2 ) == min(sqrt( (1-SDMeval@TPR)^2 + (1-SDMeval@TNR)^2 ))
            )] )
            ## Obs1: esse eh o procedimento realizado pelo Biomod. FONTE: http://lists.r-forge.r-project.org/pipermail/biomod-commits/2010-August/000129.html
            ## Obs2: Formula pitagorica retirada de Jimenez-Valverde (2012) Global Ecology and Biogreography.
            
            ##tabela com as metricas calculadas
            evaluationScoresOut = data.frame(Model.name = rep(paste(names(SDMlist)[i],'_AllData',sep=''), 2),
                                             Eval.metric = c('TSS','ROC'),
                                             Testing.data = c(tssVal, aucVal),
                                             Evaluating.data = c(NA,NA),
                                             Cutoff = c(current_thre_maximizingTSS*1000, current_thre_maximizingROC*1000),
                                             Sensitivity = c(SDMeval@TPR[which(SDMeval@t==current_thre_maximizingTSS)]*100, SDMeval@TPR[which(SDMeval@t==current_thre_maximizingROC)]*100),
                                             Specificity = c(SDMeval@TNR[which(SDMeval@t==current_thre_maximizingTSS)]*100, SDMeval@TNR[which(SDMeval@t==current_thre_maximizingROC)]*100))
            return(evaluationScoresOut)
          })
          
          evaluationScoresList = do.call('rbind', evaluationScoresList)
          evaluationScores = rbind(evaluationScores, evaluationScoresList)
          
        }
        
        ##outputs 
        statResultsSDMimproved = read.csv(paste(projectFolder,'/SDMimproved/StatisticalResults_SDMimproved.csv',sep=''), header = TRUE)
        statResultsSDMimproved = rbind(statResultsSDMimproved,
                                     data.frame(SDM='improved', sp=paste('sp',i,sep=''), sampleSize=sampleSizes[j], biaslevel=current_vies_level, evaluationScores))
        write.csv(statResultsSDMimproved, file=paste(projectFolder,'/SDMimproved/StatisticalResults_SDMimproved.csv',sep=''), row.names=FALSE)
        
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n \n", file="logfileSDMimprovedOutrosModelos.txt", append=TRUE)})
    }
  }
}




##PARTE 5: SELECIONANDO MELHOR MODELO-HCPA E IMPLEMENTANDO PROJECOES DOS SDMs-HCPA




##arquivo de log
write.table(x=NULL, file=paste(projectFolder,'/logfileSDMimprovedProjecoes.txt',sep='')) #criando arquivo
cat('Log file - Started at: ', as.character(Sys.time()), "\n \n", file=paste(projectFolder,'/logfileSDMimprovedProjecoes.txt',sep=''), append=TRUE) #gravando no arquivo

for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    for(current_vies_level in 1:vies_levels){
      tryCatch({
        
        ##logfile
        cat('STARTING SCENARIO: sp = ', i, '/ sampleSizes = ', sampleSizes[j], '/ current_vies_level = ', current_vies_level, ' / time: ', as.character(Sys.time()),  "\n \n", file = paste(projectFolder,'/logfileSDMimprovedProjecoes.txt',sep=''), append = TRUE) #gravando no arquivo
        
        ##diretorio para o biomod2 salvar resultados para SDMimproved
        setwd(file.path(projectFolder,'SDMimproved'))
        
        ##abrindo dados dos modelos
        statResultsSDMimproved = read.csv(paste(projectFolder,'/SDMimproved/StatisticalResults_SDMimproved.csv',sep=''), header = TRUE)
        evaluationScores = subset(statResultsSDMimproved, (statResultsSDMimproved$SDM == 'improved') & 
                                    (statResultsSDMimproved$sp == paste('sp',i,sep='')) & 
                                    (statResultsSDMimproved$sampleSize == sampleSizes[j]) & 
                                    (statResultsSDMimproved$biaslevel == current_vies_level))
        
        ##selecao do modelo de maior sensibilidade            
        bestModel = evaluationScores[which(evaluationScores$Specificity == max(evaluationScores$Specificity, na.rm = TRUE)),]
        modelNames = as.character(bestModel$Model.name)
        modelNames = strsplit(x = modelNames, split = "_")
        modelNames = lapply(seq(length(modelNames)), function(i) rev(modelNames[[i]]) )
        modelNames = lapply(seq(length(modelNames)), function(i) paste(modelNames[[i]], collapse = '_', sep='') )
        modelNames = unlist(modelNames)
        
        ##abrindo os modelos
        modelFiles = list.files(paste(projectFolder,'/SDMimproved/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved/models/sp',i,'_sample',sampleSizes[j],'_SDMimproved',sep=''), pattern = paste(modelNames,collapse = '|'), full.names = TRUE)
        lapply(modelFiles, load, .GlobalEnv)
        
        ##implementado projecoes
        modelFiles = basename(modelFiles)
        if (sum(grep(pattern = '.RData', x = modelFiles)) > 0){modelFiles = gsub(pattern = '.RData', replacement = '', x = modelFiles)}
        modelNames = lapply(modelFiles, get)
        modelProj = lapply(seq(length(modelNames)), function(i)  predict(modelNames[[i]], predictors))
        
        ##salvando
        if (length(modelProj)==1){
          modelProj = modelProj[[1]]
          names(modelProj) = modelFiles
          crs(modelProj) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
          writeRaster(modelProj, bylayer=FALSE, filename=paste(projectFolder,'/SDMimproved/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved/models/', modelFiles, sep=''), format='ascii', overwrite=TRUE)
        }else{
          modelProj = stack(modelProj)
          names(modelProj) = modelFiles
          crs(modelProj) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')  
          writeRaster(modelProj, bylayer=TRUE, filename=paste(projectFolder,'/SDMimproved/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved/models/', modelFiles, sep=''), format='ascii', overwrite=TRUE)
        }
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n \n", file=paste(projectFolder,'/logfileSDMimprovedProjecoes.txt',sep=''), append=TRUE)})
    }
  }
}
