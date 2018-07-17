## Este script contem algoritmos para gerar 'n' sps artificiais aleatoriamente (Parte 1), modelar ausencias (Parte 2), modelar distribuicao com pseudoausencias melhoradas (Parte 3) e analise dos resultados atraves de graficos simples (Parte 4)
## Anderson A. Eduardo, outubro/2014


##abrindo pacotes necessarios
#Sys.setenv(JAVA_HOME = "/usr/lib/jvm/java-7-openjdk-amd64")
options(java.parameters = "Xmx7g")
library(rJava)
library(raster)
library(biomod2)
library(dismo)

# ##definindo prametros e variaveis globais (NOTEBOOK)
# projectFolder = "/home/anderson/Projetos/HCPA" #pasta do projeto
# envVarFolder = "/home/anderson/gridfiles/dados_projeto" #pasta com as variaveis ambientais
# envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
# AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
# maxentFolder = '/home/anderson/R/x86_64-pc-linux-gnu-library/3.4/dismo/java/maxent.jar' #pasta para resultados do maxent
# ## spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
# ## sdmTypes = c('normal','optimized')
# sampleSizes = c(10,20,40,80,160)
# NumRep = 10 #numero de replicas (de cada cenario amostral)
# vies_levels = 5
# ##variaveis preditoras
# ## elevation = raster('/home/anderson/PosDoc/dados_ambientais/DEM/DEM.tif')
# predictors = stack(list.files(path=envVarPaths[1],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
# predictors = predictors[[c('bioclim_01','bioclim_12')]]
# predictors = stack(mask(x=predictors, mask=AmSulShape))
# crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
# Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
# statResultsSDMnormal = data.frame() #tabela de estatisticas basicas do modelo
# statResultsSDMimproved = data.frame()

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
## statResultsSDMnormal = data.frame() #tabela de estatisticas basicas do modelo
## statResultsSDMimproved = data.frame()

##definindo prametros e variaveis globais (Yavanna)
projectFolder = "D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas" #pasta do projeto
envVarFolder = "D:/Anderson_Eduardo/gridfiles/variaveis ambientais AmSul 120kyr/000" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, pattern='.asc', full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("D:/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
maxentFolder = 'D:/Anderson_Eduardo/maxent/maxent.jar' #pasta para resultados do maxent
sampleSizes = c(20,40,80,160)
NumRep = 10 #numero de replicas (de cada cenario amostral)
vies_levels = 5
##variaveis preditoras
##elevation = raster('J:/Pesquisadorxs/Anderson_Eduardo/DEM/DEM.tif')
predictors = stack(envVarPaths) #predictors com todas as variaveis (presente)
predictors = predictors[[c('bioclim_01','bioclim_12')]]
predictors = stack(mask(x=predictors, mask=AmSulShape))
crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
SDMreplicates = 10 #numero de replicas para calibracao/validacao dos SDMs
statResultsSDMnormal = data.frame() #tabela de estatisticas basicas do modelo
statResultsSDMimproved = data.frame()


# ## funcoes autorais :-) #NOTEBOOK
# source('/home/anderson/R-Scripts/makeSpecies.R')
# source('/home/anderson/R-Scripts/rangeByAC.R')
# source('/home/anderson/R-Scripts/makeOutput.R')
# source('/home/anderson/R-Scripts/bestModel.R')

## ## funcoes autorais :-) #LORIEN
## source('D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas/makeSpecies.R')
## source('D:/Anderson_Eduardo/rangeByAC.R')
## source('D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas/makeOutput.R')
## source('D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas/bestModel.R')

## funcoes autorais :-) #Yavanna
source('D:/Anderson_Eduardo/R/makeSpecies.R')
source('D:/Anderson_Eduardo/R/rangeByAC.R')
#source('D:/Anderson_Eduardo/R/makeOutput.R')
source('D:/Anderson_Eduardo/R/bestModel.R')

##definindo diretorio de trabalho (importante porque o biomod2 salva tudo automaticamente)
setwd(projectFolder)




##PARTE 1: criando as especies artificiais


##registrando hora do incio
timeOne = Sys.time()


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
  
  
  SpDistAC = makeSpecies(predictors, i)
  
  ##criando imagem da distribuicao de cada especie
  jpeg(filename=paste(projectFolder,'/virtual species/sp',i,'.jpeg',sep=''))
  plot(SpDistAC)
  dev.off()
  ##
  writeRaster(x=SpDistAC, filename=paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''), overwrite=TRUE)
  rm(SpDistAC) ##teste do bug persistente
  
}


##calculando tempo gasto para processamento
Sys.time() - timeOne




##PARTE 2: modelando ausencias 




##registrando hora do incio
timeOne = Sys.time()

##arquivo de log
write.table(x=NULL, file=paste(projectFolder,'/logfileSDMnormal.txt',sep='')) #criando arquivo
cat('Log file - Started at: ', Sys.time(), "\n \n", file="logfileSDMnormal.txt", append=TRUE) #gravando no arquivo


for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    for(current_vies_level in 1:vies_levels){
      tryCatch({

        ##logfile
        cat('STARTING SCENARIO: sp = ', i, '/ sampleSizes = ', sampleSizes[j], '/ current_vies_level = ', current_vies_level,  "\n \n", file = "logfileSDMnormal.txt", append = TRUE) #gravando no arquivo
        
        ##diretorio para o biomod2 salvar resultados para SDMnormal
        setwd(file.path(projectFolder))
        
        ##verifica e cria diretorio para salvar resultados da especie atual
        if (file.exists('SDMnormal')){
          setwd(paste(projectFolder,'/SDMnormal',sep=''))
        } else {
          dir.create('SDMnormal')
          setwd(paste(projectFolder,'/SDMnormal',sep=''))
        }
      
        
        ##definindo variaveis e parametros locais
        ##occPoints = read.csv(paste(mainSampleFolder,sdmTypes[h],'/',spsTypes[i],'/occ',sampleSizes[j],'.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia
        SpDistAC = raster(paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''))
        values(SpDistAC)[values(SpDistAC)==0] = NA
        vies_i = sample(x=sequence(vies_levels), size=current_vies_level)
        
        
        ##pontos de ocorrencia com diferentes niveis de vies##
        
        
        ##bias layer da iteracao atual - sps range
        SpDistAC = raster(paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''))     
        crs(SpDistAC) =  CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
        biasLayerBGpts =  data.frame( dismo::randomPoints(mask=SpDistAC, n=100, prob=TRUE) )
        ptsCluster = kmeans(x=biasLayerBGpts, centers=5, nstart=5)
        biasLayerBGpts$kcluster = ptsCluster$cluster
        biasLayerBGpts = biasLayerBGpts[ sapply(biasLayerBGpts$kcluster, function(x) x  %in% vies_i), ]
        
        ## ##inspecao visual
        ## plot(AmSulShape)
        ## points(biasLayerBGpts, pch=19, col=ptsCluster$cluster)
        ## ptsNew = biasLayerBGpts[ptsCluster$cluster!=2,]
        ## points(ptsNew, pch='x', cex=2, col=ptsCluster$cluster)
        
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
        
        occPoints = dismo::randomPoints(mask=SpDistAC*biasLayerBGptsKDE, n=sampleSizes[j], prob=TRUE) #sorteando pontos da distribuicao modelada

        ##rm(SpDistAC) ##teste do bug persistente
        occPoints = data.frame(lon=occPoints[,1],lat=occPoints[,2])
    
        
        # ## amostrando pontos back ground
        # bgPoints = dismo::randomPoints(mask=biasLayerBGptsKDE+0.001, n=5000, p=occPoints, prob=TRUE) #sorteando pontos da distribuicao modelada
        # 
        # ##betterPseudoDF = extract(projStackBIN[[k]], be tterPseudoPoints) #distinguindo entre occ e ausencia
        # bgPointsDF = data.frame(lon=bgPoints[,1], lat=bgPoints[,2], occ=0) #distinguindo entre occ e ausencia
        # ##betterPseudo[[k]] =  data.frame(lon=betterPseudoPoints[,1], lat=betterPseudoPoints[,2], occ=betterPseudoDF) #data.frame
        # ##betterPseudo[[k]] = betterPseudo[[k]][which(betterPseudo[[k]]$occ==0),] #excluindo as presencas
        # bgPointsVar = extract(predictors, bgPointsDF[,c('lon','lat')]) #obtendo as variaveis preditoras nos pontos
        # bgPointsDF = data.frame(bgPointsDF, bgPointsVar) #motando dataset
        # 
        # ## plot(projStackBIN[[k]])
        # ## points(betterPseudo[[k]][,c('lon','lat')])
        # 
        # ##definindo variaveis e parametros locais para o biomod2 (que rodara a seguir)
        # occPointsDF = data.frame(occPoints, occ=1, extract(predictors,occPoints)) #assumindo que o stack predictors contem apenas as variaveis empregadas no projeto atual
        # 
        # ##agrupando ocorrencias e pseudo-ausencias melhoradas
        # dataSet = data.frame(rbind(occPointsDF, bgPointsDF)) #planilha de dados no formato SWD
        # ##variaveis e parametros locais especificos para o biomod2
        # myRespName <- paste('sp',i,sep='')  # nome do cenario atual (para biomod2)
        # myResp <- dataSet[,c('occ')] # variavel resposta (para biomod2)
        # myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
        # myExpl = dataSet[,c('bioclim_01','bioclim_12')]  #variavel preditora (para biomod2)
        # 
        # 
        # ##ajuste de dados de entrada para biomod2
        # myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
        #                                      expl.var = myExpl,
        #                                      resp.xy = myRespXY, 
        #                                      resp.name = paste(myRespName,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal',sep=''))
        # 
        
        myResp <- data.frame(lon=occPoints[,1], lat=occPoints[,2])
        coordinates(myResp) <- ~ lon + lat #transformando em spatialPoints
        crs(myResp) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #transformando em spatialPoints

        ##variaveis e parametros locais especificos para o biomod2
        myRespName <- paste('sp',i,sep='')  # nome do cenario atual (para biomod2)
        ##myResp <- dataSet[,c('pres')] # variavel resposta (para biomod2)
        ##myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
        myExpl = predictors  #variavel preditora (para biomod2)

        ##ajuste de dados de entrada para biomod2
        myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                             expl.var = myExpl,
                                             resp.name = paste(myRespName,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal',sep=''),
                                             PA.nb.rep = 1,
                                             PA.nb.absences = 5000)
        
        ## ##inspecionando o objeto gerado pela funcao do biomod2
        ## myBiomodData
        ## plot(myBiomodData)
        
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
        
        
        ##iterando replicas (ATENCAO: tem que ser o mesmo numero que no biomod2)
        
        SDMlist = list()
        
        for (iter in 1:SDMreplicates){
          
          curent_occData_train = subset( x=myBiomodData@coord[which(!is.na(myBiomodData@data.species) & myBiomodData@data.species==1),], subset=kfold(myBiomodData@coord[!is.na(myBiomodData@data.species) & myBiomodData@data.species==1,], k=3)!=1 ) #occ para teste
          curent_bgData_train = myBiomodData@coord[which(is.na(myBiomodData@data.species) | myBiomodData@data.species!=1),] #bg para train
          curent_occData_test = subset( x=myBiomodData@coord[which(!is.na(myBiomodData@data.species) & myBiomodData@data.species==1),], subset=kfold(myBiomodData@coord[!is.na(myBiomodData@data.species) & myBiomodData@data.species==1,], k=3)==1 ) #occ para teste
          curent_bgData_test = myBiomodData@coord[which(is.na(myBiomodData@data.species) | myBiomodData@data.species!=1),] #bg para teste
          
          ##rodando Mahalanobis, Bioclim e Domain
          mahalanobisSDM = mahal( x=predictors, p=curent_occData_train )
          bioclimSDM = bioclim( x=predictors, p=curent_occData_train )
          domainSDM = domain( x=predictors, p=curent_occData_train )
          
          ##lista de objetos com os modelos
          SDMlistRaw = list('mahalanobisSDM'=mahalanobisSDM,'bioclimSDM'=bioclimSDM,'domainSDM'=domainSDM)
          names(SDMlistRaw) = c(paste('mahalanobisSDM_RUN',iter,sep=''), 
                                paste('bioclimSDM_RUN',iter,sep=''), 
                                paste('domainSDM_RUN',iter,sep=''))
          SDMlist = append(SDMlist, SDMlistRaw)
          
          for (model_i in names(SDMlistRaw)){
            
            currentSDM = SDMlist[[model_i]] #modelo da iteracao atual
            
            SDMeval = evaluate(model = currentSDM,
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
            
            evaluationScores = rbind(evaluationScores,
                                     data.frame(Model.name = rep(paste(model_i,'_PA1',sep=''),2),
                                                Eval.metric = c('TSS','ROC'),
                                                Testing.data = c(tssVal, aucVal),
                                                Evaluating.data = c(NA,NA),
                                                Cutoff = c(current_thre_maximizingTSS*1000, current_thre_maximizingROC*1000),
                                                Sensitivity = c(SDMeval@TPR[which(SDMeval@t==current_thre_maximizingTSS)]*100, SDMeval@TPR[which(SDMeval@t==current_thre_maximizingROC)]*100),
                                                Specificity = c(SDMeval@TNR[which(SDMeval@t==current_thre_maximizingTSS)]*100, SDMeval@TNR[which(SDMeval@t==current_thre_maximizingROC)]*100)))
          }
        }
        
        ##construindo tabela de outputs
        ## statResultsSDMnormal = makeOutput(evaluationScores, statResultsSDMnormal, i, j, 'normal', sampleSizes[j])
        statResultsSDMnormal = rbind(statResultsSDMnormal,
                                     data.frame(SDM='normal', sp=paste('sp',i,sep=''), sampleSize=sampleSizes[j], biaslevel=current_vies_level, evaluationScores))
        
        write.csv(statResultsSDMnormal, file=paste(projectFolder,'/StatisticalResults_SDMnormal.csv',sep=''), row.names=FALSE)
        
        ##selecao do modelo de maior sensibilidade            
        modelNames = bestModel(evaluationScores, myBiomodModelOut)
        #rm(evaluationScores)
        
        
        ##rodando algortmo de projecao (i.e. rodando a projecao)##
        
        ##predicao espacial para SDMs implemnetados pelo BIOMOD2
        if ( length(grep(pattern=paste(modelNames,collapse='|'), x=myBiomodModelOut@models.computed, value=FALSE)) > 0 ){
          
          myExpl = predictors
          
          myBiomodProj <- BIOMOD_Projection(
            modeling.output = myBiomodModelOut,
            new.env = myExpl,
            binary.meth = c('ROC','TSS'),
            proj.name = paste('sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal',sep=''),
            selected.models = grep(pattern=paste(modelNames,collapse='|'),x=myBiomodModelOut@models.computed, value=TRUE),
            compress = 'TRUE',
            build.clamping.mask = 'FALSE')
          
          # ##pegando os gridfiles
          # myCurrentProj <- get_predictions(myBiomodProj) 
          # 
          # ##media e desvio padrao das projecoes
          # for(myCurrentProj_i in names(myCurrentProj)){
          #   ##media e desvio padrao
          #   rasterLayer_i_mean = calc( x=myCurrentProj[[myCurrentProj_i]], fun=mean )
          #   names(rasterLayer_i_mean) = myCurrentProj_i
          #   rasterLayer_i_SD = calc( x=myCurrentProj[[myCurrentProj_i]], fun=sd )
          #   names(rasterLayer_i_SD) = myCurrentProj_i
          #   ##salvando no HD
          #   writeRaster(rasterLayer_i_mean, file=paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal/','sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal_Mean.asc',sep=''), overwrite=TRUE)
          #   writeRaster(rasterLayer_i_SD, file=paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal/','sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal_SD.asc',sep=''), overwrite=TRUE)
            
          # }
        }
        
        
        ##predicao espacial para SDMs NAO implemnetados pelo BIOMOD2 (bioclim, mahalanobis, domain)
        ##verifica e cria diretorio para salvar resultados da especie atual
        currentProjFolder = paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal', sep='')
        if (!file.exists(currentProjFolder, recursive=TRUE)){
           dir.create(currentProjFolder, recursive=TRUE)
         }
        
        ##nomes dos modelos
        bimodNames = gsub( paste('sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal_PA1_',sep=''), '',myBiomodModelOut@models.computed )
        nonBiomod = grep(pattern=paste(bimodNames,collapse='|'), x=modelNames, value=TRUE, invert=TRUE)
        
        if ( length(nonBiomod) > 0 ){
          currentSDMpred = list()
          
          modelNamesAdjusted = vector()
          for(modName_i in 1:length(nonBiomod)){
            modelNamesRawSplit = strsplit(x=nonBiomod[[modName_i]], split='_')
            modelNamesAdjusted = append(x=modelNamesAdjusted, values=paste(modelNamesRawSplit[[1]][2],modelNamesRawSplit[[1]][1],sep='_'))
          }
          
          for (model_i in modelNamesAdjusted){ #rodando predict para cada um dos modelos indentificados
            for (repl_i in 1:SDMreplicates){ #replicas para mapear areas de incerteza
              ##dataset
              curent_occData_test = subset( x=myBiomodData@coord[!is.na(myBiomodData@data.species),], subset=kfold(myBiomodData@coord[!is.na(myBiomodData@data.species),], k=3)==1 ) #occ para teste
              ##predicao espacial
              currentSDM = SDMlist[[model_i]] #modelo da iteracao atual
              currentSDMpred = append(currentSDMpred, dismo::predict(predictors, currentSDM))
            }
            
            rasterLayer_i_mean = calc( x=stack(currentSDMpred), fun=mean ) #media das projecoes
            rasterLayer_i_SD = calc( x=stack(currentSDMpred), fun=sd ) #desvio padrao das projecoes
            ##salvando no HD - media e variancia do mapa de suitability
            writeRaster(rasterLayer_i_mean, 
                        file=paste(projectFolder,
                                   '/SDMnormal/sp',
                                   i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal',
                                   '/proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal/', 
                                   'proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal_sp',i,'.sample',sampleSizes[j],'.',model_i,'_biaslevel',current_vies_level,'_Mean_SDMnormal.asc',sep=''), overwrite=TRUE)
            
            writeRaster(rasterLayer_i_SD, 
                        file=paste(projectFolder,
                                   '/SDMnormal/sp',
                                   i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal',
                                   '/proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal/', 
                                   'proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal_sp',i,'.sample',sampleSizes[j],'.',model_i,'_biaslevel',current_vies_level,'_SD_SDMnormal.asc',sep=''), overwrite=TRUE)
            
            
            ##mapas binarios
            current_thres_TSS = statResultsSDMnormal[ statResultsSDMnormal$Model.name %in% grep(pattern=paste(unlist(strsplit(x=model_i, split='_')), collapse='*.*'), x=statResultsSDMnormal[,'Model.name'], value=TRUE) & statResultsSDMnormal$Eval.metric =='TSS' ,'Cutoff']/1000
            current_thres_TSS = mean(current_thres_TSS)
            current_thres_ROC = statResultsSDMnormal[ statResultsSDMnormal$Model.name %in% grep(pattern=paste(unlist(strsplit(x=model_i, split='_')), collapse='*.*'), x=statResultsSDMnormal[,'Model.name'], value=TRUE) & statResultsSDMnormal$Eval.metric =='ROC' ,'Cutoff']/1000
            current_thres_ROC = mean(current_thres_ROC)
            ##
            writeRaster(rasterLayer_i_mean > current_thres_TSS, 
                        file=paste(projectFolder,
                                   '/SDMnormal/sp',
                                   i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal',
                                   '/proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal/',
                                   'proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal_sp',i,'.sample',sampleSizes[j],'.',model_i,'_biaslevel',current_vies_level,'_Mean_SDMnormal_TSSbin.asc',sep=''), 
                        overwrite=TRUE)
            writeRaster(rasterLayer_i_mean > current_thres_ROC, 
                        file=paste(projectFolder,
                                   '/SDMnormal/sp',
                                   i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal',
                                   '/proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal/',
                                   'proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal_sp',i,'.sample',sampleSizes[j],'.',model_i,'_biaslevel',current_vies_level,'_Mean_SDMnormal_ROCbin.asc',sep=''), 
                        overwrite=TRUE) 
            
          }
        }
        
        
        ## ##rodando o algoritmo de consenso dos modelos (i.e. ensemble model)
        ## myBiomodEM = BIOMOD_EnsembleModeling(
        ##     modeling.output = myBiomodModelOut,
        ##     chosen.models = modelNames)
        
        ## ##forecasting com o consenso dos algoritmos (i.e. ensemble projection)
        ## myBiomodEF = BIOMOD_EnsembleForecasting(
        ##     EM.output = myBiomodEM,
        ##     binary.meth = c('TSS','ROC'),
        ##     projection.output = myBiomodProj)
        
        ##writeRaster(projStackBIN,file=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[j],'.replica',k,'/proj_',l,'kyr/proj_',i,'kyr','.sample',sampleSizes[j],'.replica',k,'_BIN.asc',sep=''),row.names=FALSE)
        
        ##projStackBIN = projStack>0.5  #BinaryTransformation(projStack,"10")
        
        
        ##projecoes distribuicao (i.e., AREAS DE PRESENCA) pelo SDM normal (as projecoes anteriores eram para AREAS DE AUSENCIA)##
        # nsps_distribution_folder = paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal/sps_distribution', sep='')
        # if (file.exists(nsps_distribution_folder, recursive=TRUE)){
        #   setwd(nsps_distribution_folder)
        # } else {
        #   dir.create(nsps_distribution_folder, recursive=TRUE)
        #   setwd(nsps_distribution_folder)
        # }
        
        
        #BIOMOD
        myBiomodProj <- BIOMOD_Projection(
          modeling.output = myBiomodModelOut,
          new.env = myExpl,
          binary.meth = c('ROC','TSS'),
          proj.name = paste('sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal_spsDistribution',sep=''),
          selected.models = 'all',
          compress = 'TRUE',
          build.clamping.mask = 'FALSE')
        
        
        #NONBIOMOD
        ##predicao espacial para SDMs NAO implemnetados pelo BIOMOD2 (bioclim, mahalanobis, domain)
        ##verifica e cria diretorio para salvar resultados da especie atual
        nsps_distribution_folder = paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal/proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal_spsDistribution',sep='')
        if (file.exists(nsps_distribution_folder, recursive=TRUE)){
          setwd(nsps_distribution_folder)
        } else {
          dir.create(nsps_distribution_folder, recursive=TRUE)
          setwd(nsps_distribution_folder)
        }
        
        currentSDMpred = list()
        nonBiomodNames = unique(gsub('_RUN.*', '', names(SDMlist)))
        
        for (model_i in nonBiomodNames){ #rodando predict para cada um dos modelos indentificados
          ids  = grep(pattern = model_i, x=names(SDMlist), value=FALSE)
          
          for (model_ids in ids){
            currentSDM = SDMlist[[model_ids]] #modelo da iteracao atual
            currentSDMpred = append(currentSDMpred, dismo::predict(predictors, currentSDM))
          }
          
          rasterLayer_i_mean = calc( x=stack(currentSDMpred), fun=mean ) #media das projecoes
          rasterLayer_i_SD = calc( x=stack(currentSDMpred), fun=sd ) #desvio padrao das projecoes
          ##salvando no HD - media e variancia do mapa de suitability
          writeRaster(rasterLayer_i_mean, 
                      file=paste('proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal_',model_i,'_spsDistribution.asc',sep=''), 
                      overwrite=TRUE)
          
          writeRaster(rasterLayer_i_SD, 
                      file=paste('proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal_',model_i,'_spsDistribution.asc',sep=''), 
                      overwrite=TRUE)
          
          
          ##mapas binarios
          current_thres_TSS = statResultsSDMnormal[ which( statResultsSDMnormal$Model.name %in% grep(pattern=model_i, x=statResultsSDMnormal[,'Model.name'], value=TRUE) & statResultsSDMnormal$Eval.metric =='TSS' ) , 'Cutoff']/1000
          current_thres_TSS = mean(current_thres_TSS)
          current_thres_ROC = statResultsSDMnormal[ statResultsSDMnormal$Model.name %in% grep(pattern=paste(unlist(strsplit(x=model_i, split='_')), collapse='*.*'), x=statResultsSDMnormal[,'Model.name'], value=TRUE) & statResultsSDMnormal$Eval.metric =='ROC' ,'Cutoff']/1000
          current_thres_ROC = mean(current_thres_ROC)
          ##
          writeRaster(rasterLayer_i_mean > current_thres_TSS, 
                      file=paste('proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal_',model_i,'_spsDistribution_TSSbin.asc',sep=''), 
                      overwrite=TRUE)
          writeRaster(rasterLayer_i_mean > current_thres_ROC, 
                      file=paste('proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal_',model_i,'_spsDistribution_ROCbin.asc',sep=''), 
                      overwrite=TRUE) 
        }
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n \n", file="logfileSDMnormal.txt", append=TRUE)})

    }
  }
}



##calculando tempo gasto para processamento
Sys.time() - timeOne




##PARTE 3: SDM com pseudoausencias melhoradas




##registrando hora do incio
timeOne = Sys.time()

##arquivo de log
write.table(x=NULL, file=paste(projectFolder,'/logfileSDMimproved.txt',sep='')) #criando arquivo
cat('Log file - Started at: ', timeOne, "\n \n", file="logfileSDMimproved.txt", append=TRUE) #gravando no arquivo


for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    for(current_vies_level in 1:vies_levels){
      tryCatch({
        
        ##logfile
        cat('STARTING SCENARIO: sp = ', i, '/ sampleSizes = ', sampleSizes[j], '/ current_vies_level = ', current_vies_level,  "\n \n", file = "logfileSDMimproved.txt", append = TRUE) #gravando no arquivo
        
        
        ##diretorio para o biomod2 salvar resultados para SDMnormal
        setwd(file.path(projectFolder))
        
        ##verifica e cria diretorio para salvar resultados da especie atual
        if (file.exists('SDMimproved')){
          setwd(paste(projectFolder,'/SDMimproved',sep=''))
        } else {
          dir.create('SDMimproved')
          setwd(paste(projectFolder,'/SDMimproved',sep=''))
        }
        
        ##definindo variaveis e parametros locais
        betterPseudo = list()
        betterPseudoVar = list()
        vies_i = sample(x=sequence(vies_levels), size=current_vies_level)
        
        
        ##pontos de ocorrencia com diferentes niveis de vies##
        
        
        ## ##bias layer da iteracao atual - America do Sul
        ## biasLayerBGpts =  data.frame( dismo::randomPoints(mask=predictors[[1]], n=10000) )
        ## ptsCluster = kmeans(x=biasLayerBGpts, centers=10, nstart=10)
        ## biasLayerBGpts$kcluster = ptsCluster$cluster
        ## biasLayerBGpts = biasLayerBGpts[ sapply(biasLayerBGpts$kcluster, function(x) x  %in% vies_i), ]
        
        ##bias layer da iteracao atual - sps range
        SpDistAC = raster(paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''))
        crs(SpDistAC) = CRS( '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0' )
        biasLayerBGpts =  data.frame( dismo::randomPoints(mask=SpDistAC, n=100, prob=TRUE) )
        ptsCluster = kmeans(x=biasLayerBGpts, centers=5, nstart=5)
        biasLayerBGpts$kcluster = ptsCluster$cluster
        biasLayerBGpts = biasLayerBGpts[ sapply(biasLayerBGpts$kcluster, function(x) x  %in% vies_i), ]
        
        ## ##inspecao visual
        ## plot(AmSulShape)
        ## points(biasLayerBGpts, pch=19, col=ptsCluster$cluster)
        ## ptsNew = biasLayerBGpts[ptsCluster$cluster!=2,]
        ## points(ptsNew, pch='x', cex=2, col=ptsCluster$cluster)
        
        ##KDE
        bgArea = SpDistAC*0 #predictors[[1]]*0
        crs(bgArea) = CRS( '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0' )
        biasLayer = data.frame(biasLayerBGpts[,1:2])
        coordinates(biasLayer) = ~x+y
        crs(biasLayer) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
        biasLayerBGptsKDE = MASS::kde2d(x=biasLayer$x, y=biasLayer$y, n=100, lims=c(xl=extent(bgArea)@xmin, xu=extent(bgArea)@xmax, yl=extent(bgArea)@ymin, yl=extent(bgArea)@ymax ))
        ##raster
        biasLayerBGptsKDE = raster(biasLayerBGptsKDE)
        biasLayerBGptsKDE = mask(biasLayerBGptsKDE, AmSulShape)
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
        
        
        ##potos de ocorrencia          
        ##SpDistAC = raster(paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''))
        ##biasLayerBGptsKDE = crop(biasLayerBGptsKDE,SpDistAC)
        ##extent(biasLayerBGptsKDE) =  extent(SpDistAC) 
        ## values(SpDistAC)[values(SpDistAC)==0] = NA
        occPoints = dismo::randomPoints( mask=SpDistAC*biasLayerBGptsKDE, n=sampleSizes[j], prob=TRUE ) #sorteando pontos da distribuicao modelada
        rm(SpDistAC) ##teste do bug persistente
        occPoints = data.frame(lon=occPoints[,1],lat=occPoints[,2])       
        
        ## binTSS = raster(paste(projectFolder,'/SDMnormal/','sp',i,'.sample',sampleSizes[j],'.SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal_sp',i,'.sample',sampleSizes[j],'.SDMnormal_ensemble_TSSbin.grd' ,sep=''))
        ## binAUC = raster(paste(projectFolder,'/SDMnormal/','sp',i,'.sample',sampleSizes[j],'.SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal_sp',i,'.sample',sampleSizes[j],'.SDMnormal_ensemble_ROCbin.grd' ,sep=''))
        
        ## binTSS = stack(list.files(path=paste(projectFolder,'/SDMnormal/','sp',i,'.sample',sampleSizes[j],'.SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal',sep=''), pattern='TSSbin', full.names=TRUE))
        ## binTSS = binTSS[[grep(pattern='_SD_', x=names(binTSS), invert=TRUE)]]
        ## binTSS = calc(x=binTSS, fun=sum)
        
        bin_files = list.files(path=paste(projectFolder,'/SDMnormal/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMnormal', sep=''), 
                               pattern='bin.grd|bin.asc', full.names=TRUE)
        bin_files = grep(pattern='*._SD_', x=bin_files, invert=TRUE, value=TRUE)
        
        projAbs = stack(bin_files)
        crs(projAbs) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
        
        ##excluindo projecoes que combrem o continente todo (i.e. erradas)
        
        verifyOverproj = function(x) length( values(x)[which(values(x)==0)] ) == 0 #funcao para verificar
        projAbsList = unstack(projAbs) #tranfromando o raster stack em lista (para a funcao)
        overProj = lapply(projAbsList, verifyOverproj) #lista dos raster com sobreprojecao
        projIdx = which(overProj != TRUE) #indices dos rasters com sobreprojecao
        
        if( length(projIdx) > 0 ){
          projAbs = projAbs[[projIdx]]
        }else{
          print(' ATENCAO! Essa iteracao falhou porque todos os rasters da etapa 1 estavam com sobreprojecao!')
          
          ##regitrando tabela de output
          statResultsSDMimproved = rbind(statResultsSDMimproved,
                                         data.frame(SDM = 'improved', 
                                                    sp=paste('sp',i,sep=''),
                                                    sampleSize = sampleSizes[j], 
                                                    biaslevel = current_vies_level,
                                                    Model.name = 'Erro de sobreprojecao na etapa 1!',
                                                    Eval.metric = 'Erro de sobreprojecao na etapa 1!',
                                                    Testing.data = NA,
                                                    Evaluating.data = NA,
                                                    Cutoff = NA,
                                                    Sensitivity = NA,
                                                    Specificity = NA))
          
          write.csv(statResultsSDMimproved, file=paste(projectFolder,'/StatisticalResults_SDMimproved','.csv',sep=''), row.names=FALSE)
          break
        }
        
        ## projStackBIN = stack(binTSS,binAUC) #empilhando mapas binarios (feitos com threshold a partir do AUC e TSS)
        ##projAbs = sum(projStackBIN) #somando (para depois pegar areas ausencia que ambos os thresolds concordam)
        projAbs = calc(projAbs, sum, na.rm=FALSE)
        
        ## amostrando pontos diretamente das areas de ausencia (abaixo do threshold) obtidas na Etapa 1 ##
        values(projAbs)[values(projAbs) != 0] = NA  #tranformando areas diferentes de zero em NA (retando somente os dados de ausencia)          
        betterPseudoPoints = dismo::randomPoints(mask=projAbs, n=5000) #sorteando pontos da distribuicao modelada
        
        ##betterPseudoDF = extract(projStackBIN[[k]], be tterPseudoPoints) #distinguindo entre occ e ausencia
        betterPseudoDF = data.frame(lon=betterPseudoPoints[,1], lat=betterPseudoPoints[,2], occ=0) #distinguindo entre occ e ausencia
        ##betterPseudo[[k]] =  data.frame(lon=betterPseudoPoints[,1], lat=betterPseudoPoints[,2], occ=betterPseudoDF) #data.frame
        ##betterPseudo[[k]] = betterPseudo[[k]][which(betterPseudo[[k]]$occ==0),] #excluindo as presencas
        betterPseudoVar = extract(predictors, betterPseudoDF[,c('lon','lat')]) #obtendo as variaveis preditoras nos pontos
        betterPseudoDF = data.frame(betterPseudoDF, betterPseudoVar) #motando dataset
        
        ## plot(projStackBIN[[k]])
        ## points(betterPseudo[[k]][,c('lon','lat')])
        
        ##definindo variaveis e parametros locais para o biomod2 (que rodara a seguir)
        occPointsDF = data.frame(occPoints, occ=1, extract(predictors,occPoints)) #assumindo que o stack predictors contem apenas as variaveis empregadas no projeto atual
        
        ##agrupando ocorrencias e pseudo-ausencias melhoradas
        dataSet = data.frame(rbind(occPointsDF,betterPseudoDF)) #planilha de dados no formato SWD
        ##variaveis e parametros locais especificos para o biomod2
        myRespName <- paste('sp',i,sep='')  # nome do cenario atual (para biomod2)
        myResp <- dataSet[,c('occ')] # variavel resposta (para biomod2)
        myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
        myExpl = dataSet[,c('bioclim_01','bioclim_12')]  #variavel preditora (para biomod2)
        
        
        ##ajuste de dados de entrada para biomod2
        myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                             expl.var = myExpl,
                                             resp.xy = myRespXY, 
                                             resp.name = paste(myRespName,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMimproved',sep=''))
        
        
        ## ##inspecionando o objeto gerado pela funcao do biomod2
        ## myBiomodData
        ## plot(myBiomodData)
        
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
        try(rm(myBiomodModelOut), silent=TRUE) #removendo informacoes da etapa anterior
        myBiomodModelOut <- BIOMOD_Modeling(
          data = myBiomodData,
          models = c('MAXENT.Phillips', 'GLM', 'GAM', 'MARS', 'CTA', 'GBM', 'RF'),
          models.options = myBiomodOption,
          NbRunEval = SDMreplicates,
          DataSplit = 75,
          models.eval.meth = c('TSS','ROC'),
          do.full.models = FALSE,
          modeling.id = myRespName) #paste(myRespName,'.sample',sampleSizes[j],'.SDMimproved',sep=''))
        
        
        ##My output data
        evaluationScores = get_evaluations(myBiomodModelOut, as.data.frame=TRUE)
        
        
        ##iterando 100 replicas (ATENCAO: tem que ser o mesmo numero que no biomod2)
        
        SDMlist = list()
        
        for (iter in 1:SDMreplicates){
          
          curent_occData_train = subset( x=myBiomodData@coord[which(!is.na(myBiomodData@data.species) & myBiomodData@data.species==1),], subset=kfold(myBiomodData@coord[!is.na(myBiomodData@data.species) & myBiomodData@data.species==1,], k=3)!=1 ) #occ para teste
          curent_bgData_train = myBiomodData@coord[which(is.na(myBiomodData@data.species) | myBiomodData@data.species!=1),] #bg para train
          curent_occData_test = subset( x=myBiomodData@coord[which(!is.na(myBiomodData@data.species) & myBiomodData@data.species==1),], subset=kfold(myBiomodData@coord[!is.na(myBiomodData@data.species) & myBiomodData@data.species==1,], k=3)==1 ) #occ para teste
          curent_bgData_test = myBiomodData@coord[which(is.na(myBiomodData@data.species) | myBiomodData@data.species!=1),] #bg para teste
          
          ##rodando Mahalanobis, Bioclim e Domain
          mahalanobisSDM = mahal( x=predictors, p=curent_occData_train )
          bioclimSDM = bioclim( x=predictors, p=curent_occData_train )
          domainSDM = domain( x=predictors, p=curent_occData_train )
          
          ##lista de objetos com os modelos
          SDMlistRaw = list('mahalanobisSDM'=mahalanobisSDM,'bioclimSDM'=bioclimSDM,'domainSDM'=domainSDM)
          names(SDMlistRaw) = c(paste('mahalanobisSDM_RUN',iter,sep=''), 
                                paste('bioclimSDM_RUN',iter,sep=''), 
                                paste('domainSDM_RUN',iter,sep=''))
          SDMlist = append(SDMlist, SDMlistRaw)
          
          for (model_i in names(SDMlistRaw)){
            
            currentSDM = SDMlist[[model_i]] #modelo da iteracao atual
            
            SDMeval = evaluate(model = currentSDM,
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
            
            evaluationScores = rbind(evaluationScores,
                                     data.frame(Model.name = rep(paste(model_i,'_PA1',sep=''),2),
                                                Eval.metric = c('TSS','ROC'),
                                                Testing.data = c(tssVal, aucVal),
                                                Evaluating.data = c(NA,NA),
                                                Cutoff = c(current_thre_maximizingTSS*1000, current_thre_maximizingROC*1000),
                                                Sensitivity = c(SDMeval@TPR[which(SDMeval@t==current_thre_maximizingTSS)]*100, SDMeval@TPR[which(SDMeval@t==current_thre_maximizingROC)]*100),
                                                Specificity = c(SDMeval@TNR[which(SDMeval@t==current_thre_maximizingTSS)]*100, SDMeval@TNR[which(SDMeval@t==current_thre_maximizingROC)]*100)))
          }
        }
        
        ##construindo tabela de outputs
        
        ## statResultsSDMimproved = makeOutput(evaluationScores, statResultsSDMimproved, i, j, 'improved', sampleSizes[j])
        ## rm(evaluationScores)
        
        statResultsSDMimproved = rbind(statResultsSDMimproved,
                                       data.frame(SDM='improved', 
                                                  sp=paste('sp',i,sep=''), 
                                                  sampleSize=sampleSizes[j], 
                                                  biaslevel=current_vies_level, 
                                                  evaluationScores))
        
        write.csv(statResultsSDMimproved, file=paste(projectFolder,'/StatisticalResults_SDMimproved','.csv',sep=''), row.names=FALSE)
        
        
        ##implementando projecoes dos modelos
        
        
        ##rodando algortmo de projecao (i.e. rodando a projecao)
        try(rm(myBiomodProj), silent=TRUE) #removendo informacao da etapa anterior
        myExpl = predictors
        
        myBiomodProj <- BIOMOD_Projection(
          modeling.output = myBiomodModelOut,
          new.env = myExpl,
          proj.name = paste(myRespName,'_sample',sampleSizes[j],'_SDMimproved',sep=''),
          selected.models = 'all',
          binary.meth = 'TSS',
          compress = 'TRUE',
          build.clamping.mask = 'TRUE',
          output.format = '.grd')
        
        
        #NON-BIOMOD
        ##predicao espacial para SDMs NAO implemnetados pelo BIOMOD2 (bioclim, mahalanobis, domain)
        ##verifica e cria diretorio para salvar resultados da especie atual
        nsps_distribution_folder = paste(projectFolder,'/SDMimproved/sp',i,'.sample',sampleSizes[j],'.biaslevel',current_vies_level,'.SDMimproved/proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMimproved_spsDistribution',sep='')
        if (file.exists(nsps_distribution_folder, recursive=TRUE)){
          setwd(nsps_distribution_folder)
        } else {
          dir.create(nsps_distribution_folder, recursive=TRUE)
          setwd(nsps_distribution_folder)
        }
        
        currentSDMpred = list()
        nonBiomodNames = unique(gsub('_RUN.*', '', names(SDMlist)))
        
        for (model_i in nonBiomodNames){ #rodando predict para cada um dos modelos indentificados
          ids  = grep(pattern = model_i, x=names(SDMlist), value=FALSE)
          
          for (model_ids in ids){
            currentSDM = SDMlist[[model_ids]] #modelo da iteracao atual
            currentSDMpred = append(currentSDMpred, dismo::predict(predictors, currentSDM))
          }
          
          rasterLayer_i_mean = calc( x=stack(currentSDMpred), fun=mean ) #media das projecoes
          rasterLayer_i_SD = calc( x=stack(currentSDMpred), fun=sd ) #desvio padrao das projecoes
          ##salvando no HD - media e variancia do mapa de suitability
          writeRaster(rasterLayer_i_mean, 
                      file=paste('proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMimproved_',model_i,'_spsDistribution.asc',sep=''), 
                      overwrite=TRUE)
          
          writeRaster(rasterLayer_i_SD, 
                      file=paste('proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMimproved_',model_i,'_spsDistribution.asc',sep=''), 
                      overwrite=TRUE)
          
          
          ##mapas binarios
          current_thres_TSS = statResultsSDMimproved[ which( statResultsSDMimproved$Model.name %in% grep(pattern=model_i, x=statResultsSDMimproved[,'Model.name'], value=TRUE) & statResultsSDMimproved$Eval.metric =='TSS' ) , 'Cutoff']/1000
          current_thres_TSS = mean(current_thres_TSS)
          current_thres_ROC = statResultsSDMimproved[ statResultsSDMimproved$Model.name %in% grep(pattern=paste(unlist(strsplit(x=model_i, split='_')), collapse='*.*'), x=statResultsSDMimproved[,'Model.name'], value=TRUE) & statResultsSDMimproved$Eval.metric =='ROC' ,'Cutoff']/1000
          current_thres_ROC = mean(current_thres_ROC)
          ##
          writeRaster(rasterLayer_i_mean > current_thres_TSS, 
                      file=paste('proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMimproved_',model_i,'_spsDistribution_TSSbin.asc',sep=''), 
                      overwrite=TRUE)
          writeRaster(rasterLayer_i_mean > current_thres_ROC, 
                      file=paste('proj_sp',i,'_sample',sampleSizes[j],'_biaslevel',current_vies_level,'_SDMimproved_',model_i,'_spsDistribution_ROCbin.asc',sep=''), 
                      overwrite=TRUE) 
        }
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n \n", file="logfileSDMimproved.txt", append=TRUE)})

    }
  }
}

##calculando tempo gasto para processamento
Sys.time() - timeOne




##PARTE 4: graficos e analise dos resultados




##ambiente para salvar os graficos
setwd(projectFolder)

##planilhas de dados
statResultsSDMimproved = read.csv(paste(projectFolder,'/StatisticalResults_SDMimproved','.csv',sep=''), header=TRUE)
statResultsSDMnormal = read.csv(paste(projectFolder,'/StatisticalResults_SDMnormal','.csv',sep=''), header=TRUE)

# statResultsSDMimproved = read.csv(paste(projectFolder,'/improved','.csv',sep=''), header=TRUE)
# statResultsSDMnormal = read.csv(paste(projectFolder,'/normal','.csv',sep=''), header=TRUE)

## deixando so os cenarios que rodaram para os dois SDMs
dim(statResultsSDMnormal)
dim(statResultsSDMimproved)

## Verificando atraves do grafico de barras simples ##
##igualando os nomes na coluna modelos (estam diferentes)
statResultsSDMimproved$Model.name = gsub(pattern='_AllData|_PA1' , replacement='' , x=as.character(statResultsSDMimproved$Model.name))
statResultsSDMnormal$Model.name = gsub(pattern='_PA1' , replacement='' , x=as.character(statResultsSDMnormal$Model.name))

##contando (pela coluna que optar para ser discriminada)
xxImp = table(statResultsSDMimproved$sampleSize)
xxNor = table(statResultsSDMnormal$sampleSize)

pooledData = cbind(xxImp,xxNor) ##juntando numa matrix (necessario para o barplot)

barplot( t(pooledData), beside=T, las=2, mar=c(12,5,5,5), legend.text=c('HCPA','Normal')) #grafico
####


##obs: usar o merge() para recortar a planilha que for maior
##exemplo:
##newTab=merge(tabB,tabA,by=c('infoA','infoB'))[,-3] ##a terceira coluna e da planilha menor

statResultsSDMnormal = merge(statResultsSDMimproved,
                             statResultsSDMnormal, 
                             all.y=TRUE,
                             by=c('sampleSize','biaslevel','Model.name','Eval.metric'))[,c('SDM.y',
                                                                 'sampleSize',
                                                                 'biaslevel',
                                                                 'Model.name', 
                                                                 'Eval.metric', 
                                                                 'Testing.data.y',
                                                                 "Evaluating.data.y",
                                                                 "Cutoff.y",
                                                                 "Sensitivity.y",
                                                                 "Specificity.y")]
 

##ajustando so nomes
names(statResultsSDMnormal) = names(statResultsSDMimproved)

modelNames = unique( gsub(pattern='_.*' , replacement='' , x=as.character(statResultsSDMnormal$Model.name)) )


##chacagem rapida
for(i in modelNames){
  print(paste(i, 
              ' ==>> SDMnormal = ', length(statResultsSDMnormal$Model.name[ grep(modelNames[1], statResultsSDMnormal$Model.name) ]),
              ' // ',
              'SDMimproved = ',length(statResultsSDMimproved$Model.name[ grep(modelNames[1], statResultsSDMimproved$Model.name) ]),
              sep=''))
}


##boxplot AUC
jpeg(filename='boxplotAUC.jpeg')
boxplot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data, statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data, ylim=c(0,1), names=c('SDM normal','SDM improved'), ylab='AUC')
dev.off()

idxNorm = which(is.finite(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data))
idxImp = which(is.finite(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data))
jpeg(filename='densidadeAUC.jpeg')
plot(density( statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data[idxNorm] ), ylim=c(0,40),col='blue', lwd=3, xlab='AUC values', main="")
lines(density( statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data[idxImp] ),col='red', lwd=3)
legend(x='topleft', legend=c('SDM normal', 'SDM improved'), lty=1, lwd=3, col=c('blue', 'red'))
dev.off()

##boxplot TSS
jpeg(filename='boxplotTSS.jpeg')
boxplot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data, statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data, ylim=c(0,1), names=c('SDM normal','SDM improved'), ylab='TSS')
dev.off()

idxNorm = which(is.finite(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data))
idxImp = which(is.finite(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data))
jpeg(filename='densidadeTSS.jpeg')
plot(density( statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data[idxNorm] ), ylim=c(0,40),col='blue', lwd=3, xlab='TSS values', main="")
lines(density( statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data[idxImp] ),col='red', lwd=3)
legend(x='topleft', legend=c('SDM normal', 'SDM improved'), lty=1, lwd=3, col=c('blue', 'red'))
dev.off()

##AUC x sample size
jpeg(filename='AUC_&_sampleSize.jpeg')
plot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$sampleSize, ylim=c(0,1), xlim=c(0,170), cex=2, pch=19, col=rgb(0,0,0,0.5), xlab='Sample size', ylab='AUC')
tendenciaSDMnormalAUC = lm(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$sampleSize)
abline(tendenciaSDMnormalAUC)
points(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$sampleSize, ylim=c(0,1), xlim=c(0,170), cex=1.5, pch=20, col=rgb(1,0,0,0.5))
tendenciaSDMimprovedAUC = lm(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$sampleSize)
abline(tendenciaSDMimprovedAUC, col='red')
legend('bottomleft',legend=c('SDM normal','SDM improved'), pch=c(19,20), col=c(rgb(0,0,0,0.5),rgb(1,0,0,0.5)))
dev.off()

##TSS x sample size
jpeg(filename='TSS_&_sampleSize.jpeg')
plot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$sampleSize, ylim=c(0,1), xlim=c(0,170), cex=2, pch=19, col=rgb(0,0,0,0.5), xlab='Sample size', ylab='TSS')
tendenciaSDMnormalTSS = lm(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$sampleSize)
abline(tendenciaSDMnormalTSS)
points(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$sampleSize, ylim=c(0,1), xlim=c(0,170), cex=1.5, pch=20, col=rgb(1,0,0,0.5))
tendenciaSDMimprovedTSS = lm(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$sampleSize)
abline(tendenciaSDMimprovedTSS, col='red')
legend('bottomright',legend=c('SDM normal','SDM improved'), pch=c(19,20), col=c(rgb(0,0,0,0.5),rgb(1,0,0,0.5)))
dev.off()


##teste de significancia

wilcox.test(statResultsSDMimproved$AUCvalue_bestModel,statResultsSDMnormal$AUCvalue_bestModel) #resultado: p<0.05

wilcox.test(statResultsSDMimproved$TSSvalue_bestModel,statResultsSDMnormal$TSSvalue_bestModel)  #resultado:p<<0.05


## TSS e AUC por algoritmo

jpeg(filename='TSS_por_algoritmo.jpeg', width=800)
par(mfrow=c(1,2), las=2, mar=c(8,5,5,1))
boxplot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Model.name, ylim=c(0,1), ylab='TSS', main='SDM normal')
boxplot(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Model.name, ylim=c(0,1), ylab='TSS', main='SDM improved')
dev.off()

jpeg(filename='AUC_por_algoritmo.jpeg', width=800)
par(mfrow=c(1,2), las=2, mar=c(12,5,5,1))
boxplot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Model.name, ylim=c(0,1), ylab='AUC', main='SDM normal')
boxplot(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Model.name, ylim=c(0,1), ylab='AUC', main='SDM improved')
dev.off()


##teste TSS
kruskal.test(TSSvalue_bestModel ~ model, data=statResultsSDMnormal) #resultado: p>0.05
kruskal.test(TSSvalue_bestModel ~ model, data=statResultsSDMimproved) #resultado: p<<0.05

##teste AUC
kruskal.test(AUCvalue_bestModel ~ model, data=statResultsSDMnormal) #resultado: p<0.05
kruskal.test(AUCvalue_bestModel ~ model, data=statResultsSDMimproved) #resultado: p<<0.05


## especificidade (escolha do melhor modelo)

jpeg(filename='especificidade_por_algoritmo.jpeg', width=800)
par(mfrow=c(1,2), las=2, mar=c(8,5,5,1))
boxplot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Specificity ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Model.name, ylab='Specificity', main='Maximization of AUC')
boxplot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Specificity ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Model.name, ylab='Specificity', main='Maximization of TSS')
dev.off()

#teste especificade
kruskal.test(maxTSSspecificity ~ model, data=statResultsSDMnormal) #rsultado: p>0.05
kruskal.test(maxAUCspecificity ~ model, data=statResultsSDMnormal) #resultado:p>0.05


## niveis de vies

jpeg(filename='aucXbiasleval.jpeg', width=800) #AUC
plot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$biaslevel, ylim=c(0,1), xlab='Bias level', ylab='AUC', pch=19, col=rgb(0,0,0,0.5))
tendenciaNorm = lm(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='ROC',]$biaslevel)
abline(tendenciaNorm)
points(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$biaslevel, ylim=c(0,1), pch=19, col=rgb(1,0,0,0.5))
tendenciaImprov = lm(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='ROC',]$biaslevel)
abline(tendenciaImprov, col='red')
dev.off()

jpeg(filename='tssXbiasleval.jpeg', width=800) #TSS
plot(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$biaslevel, ylim=c(0,1), xlab='Bias level', ylab='TSS', pch=19, col=rgb(0,0,0,0.5))
tendenciaNorm = lm(statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMnormal[statResultsSDMnormal$Eval.metric=='TSS',]$biaslevel)
abline(tendenciaNorm)
points(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$biaslevel, ylim=c(0,1), pch=19, col=rgb(1,0,0,0.5))
tendenciaImprov = lm(statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$Testing.data ~ statResultsSDMimproved[statResultsSDMimproved$Eval.metric=='TSS',]$biaslevel)
abline(tendenciaImprov, col='red')
dev.off()

