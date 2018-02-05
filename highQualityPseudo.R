## Este script contem algoritmos para gerar 'n' sps artificiais aleatoriamente (Parte 1), modelar ausencias (Parte 2), modelar distribuicao com pseudoausencias melhoradas (Parte 3) e analise dos resultados atraves de graficos simples (Parte 4)
## Anderson A. Eduardo, outubro/2014

##registrando hora do incio
timeOne = Sys.time()

##abrindo pacotes necessarios
library(raster)
library(rJava)
library(biomod2)
Sys.setenv(JAVA_HOME = "/usr/lib/jvm/java-7-openjdk-amd64")
options(java.parameters = "Xmx7g")

##definindo prametros e variaveis globais (NOTEBOOK)
projectFolder = "/home/anderson/Documentos/Projetos/Improved pseudo-absences" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
maxentFolder = '/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java/maxent.jar' #pasta para resultados do maxent
## spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
## sdmTypes = c('normal','optimized')
sampleSizes = 100  #c(10,20,40,80,160)
NumRep = 2 #10 #numero de replicas (de cada cenario amostral)
##variaveis preditoras
## elevation = raster('/home/anderson/PosDoc/dados_ambientais/DEM/DEM.tif')
predictors = stack(list.files(path=envVarPaths[1],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
predictors = predictors[[c('bioclim_01','bioclim_12')]]
predictors = stack(mask(x=predictors, mask=AmSulShape))
crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
statResultsSDMnormal = data.frame() #tabela de estatisticas basicas do modelo
statResultsSDMimproved = data.frame()

# ##definindo prametros e variaveis globais (LORIEN)
# projectFolder = "J:/Pesquisadorxs/Anderson_Eduardo/high_quality_PA" #pasta do projeto
# envVarFolder = "J:/Pesquisadorxs/Anderson_Eduardo/dados_projeto/000" #pasta com as variaveis ambientais
# envVarPaths = list.files(path=envVarFolder, pattern='.asc', full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
# AmSulShape = rgdal::readOGR("J:/Pesquisadorxs/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
# maxentFolder = 'C:/Users/WS/Documents/R/win-library/3.4/dismo/java/maxent.jar' #pasta para resultados do maxent
# ## spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
# ## sdmTypes = c('normal','optimized')
# sampleSizes = c(10,20,40,80,160)
# NumRep = 10 #numero de replicas (de cada cenario amostral)
# ##variaveis preditoras
# elevation = raster('J:/Pesquisadorxs/Anderson_Eduardo/DEM/DEM.tif')
# predictors = stack(envVarPaths,elevation) #predictors com todas as variaveis (presente)
# predictors = predictors[[c('bioclim_01','bioclim_12')]]
# predictors = stack(mask(x=predictors, mask=AmSulShape))
# crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
# Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
# statResultsSDMnormal = data.frame() #tabela de estatisticas basicas do modelo
# statResultsSDMimproved = data.frame()

## ##definindo prametros e variaveis globais (PC nupeg)
## projectFolder = "D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas" #pasta do projeto
## envVarFolder = "D:/Anderson_Eduardo/variaveis ambientais AmSul 120kyr/000" #pasta com as variaveis ambientais
## envVarPaths = list.files(path=envVarFolder, pattern='.asc', full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
## AmSulShape = rgdal::readOGR("D:/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
## maxentFolder = 'D:/Anderson_Eduardo/maxent/maxent.jar' #pasta para resultados do maxent
## sampleSizes = c(20,40,80,160)
## NumRep = 10 #numero de replicas (de cada cenario amostral)
## ##variaveis preditoras
## #elevation = raster('J:/Pesquisadorxs/Anderson_Eduardo/DEM/DEM.tif')
## predictors = stack(envVarPaths) #predictors com todas as variaveis (presente)
## predictors = predictors[[c('bioclim_01','bioclim_12')]]
## predictors = stack(mask(x=predictors, mask=AmSulShape))
## crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
## Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
## statResultsSDMnormal = data.frame() #tabela de estatisticas basicas do modelo
## statResultsSDMimproved = data.frame()

## funcoes autorais :-)
source('D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas/makeSpecies.R')
source('D:/Anderson_Eduardo/rangeByAC.R')
source('D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas/makeOutput.R')
source('D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas/bestModel.R')


##definindo diretorio de trabalho (importante porque o biomod2 salva tudo automaticamente)
setwd(projectFolder)


##PARTE 1: criando as especies artificiais


for(i in 1:Nsp){
  
  SpDistAC = makeSpecies(predictors, i)
  
  ##criando imagem da distribuicao de cada especie
  jpeg(filename=paste(projectFolder,'/virtual species/sp',i,'.jpeg',sep=''))
  plot(SpDistAC)
  dev.off()
  ##
  writeRaster(x=SpDistAC, filename=paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''), overwrite=TRUE)
  rm(SpDistAC) ##teste do bug persistente
  
}


##PARTE 2: modelando ausencias 

for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    tryCatch({
      
      ##diretorio para o biomod2 salvar resultados para SDMnormal
      setwd(file.path(projectFolder,'SDMnormal'))
      
      ##definindo variaveis e parametros locais
      ##occPoints = read.csv(paste(mainSampleFolder,sdmTypes[h],'/',spsTypes[i],'/occ',sampleSizes[j],'.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia
      SpDistAC = raster(paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''))
      values(SpDistAC)[values(SpDistAC)==0] = NA
      
      ##amostra de pontos
      occPoints = dismo::randomPoints(mask=SpDistAC, n=sampleSizes[j]) #sorteando pontos da distribuicao modelada
      rm(SpDistAC) ##teste do bug persistente
      occPoints = data.frame(lon=occPoints[,1],lat=occPoints[,2])
      
      myResp <- data.frame(lon=occPoints[,1], lat=occPoints[,2])
      coordinates(myResp) <- ~ lon + lat #transformando em spatialPoints
      crs(myResp) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #transformando em spatialPoints
      
      ##variaveis e parametros locais especificos para o biomod2
      myRespName <- paste('sp',i,sep='') # nome do cenario atual (para biomod2)
      ##myResp <- dataSet[,c('pres')] # variavel resposta (para biomod2)
      ##myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
      myExpl = predictors  #variavel preditora (para biomod2)
      
      ##ajuste de dados de entrada para biomod2
      myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                           expl.var = myExpl,
                                           resp.name = paste(myRespName,'_sample',sampleSizes[j],'_SDMnormal',sep=''),
                                           PA.nb.rep = 1)
      
      ## ##inspecionando o objeto gerado pela funcao do biomod2
      ## myBiomodData
      ## plot(myBiomodData)
      
      ##parametrizando os modelos
      myBiomodOption <- BIOMOD_ModelingOptions(
        MAXENT.Phillips = list(path_to_maxent.jar= maxentFolder,
                               linear=TRUE,
                               quadratic=TRUE,
                               product=FALSE,
                               threshold=FALSE,
                               hinge=FALSE,
                               maximumiterations=1000,
                               convergencethreshold=1.0E-5,
                               threads=2),
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
        models = c('MAXENT.Phillips', 'GLM', 'GAM', 'MARS', 'CTA', 'GBM', 'RF'),
        models.options = myBiomodOption,
        NbRunEval = 10,
        DataSplit = 75,
        models.eval.meth = c('TSS','ROC'),
        do.full.models = FALSE,
        modeling.id = paste(myRespName,'_sample',sampleSizes[j],'_SDMnormal',sep=''))
      
      ##My output data
      evaluationScores = get_evaluations(myBiomodModelOut)
      
      # ##especificidade
      # outputRawTSSspec = evaluationScores['TSS','Specificity',,,]
      # outputRawTSSspec = outputRawTSSspec[complete.cases(outputRawTSSspec),]
      # outputRawAUCspec = evaluationScores['ROC','Specificity',,,]
      # outputRawAUCspec = outputRawAUCspec[complete.cases(outputRawAUCspec),]
      # ##auc e tss
      # outputRawTSSvalue = evaluationScores['TSS','Testing.data',,,]
      # outputRawTSSvalue = outputRawTSSvalue[complete.cases(outputRawTSSvalue),]
      # outputRawAUCvalue = evaluationScores['ROC','Testing.data',,,]
      # outputRawAUCvalue = outputRawAUCvalue[complete.cases(outputRawAUCvalue),]
      # 
      # ##maior valore de especificidade de cada algoritmo implementado (tanto para TSS qto para AUC)
      # TSSspec = as.numeric(apply(outputRawTSSspec, 1, max, na.rm=TRUE))
      # AUCspec = as.numeric(apply(outputRawAUCspec, 1, max, na.rm=TRUE))
      # 
      # ##tabela auxiliar para obtencao das informacoes do melhor modelo
      # tabBestScoresTSS = data.frame(outputRawTSSspec, bestvalue=TSSspec)
      # tabBestScoresAUC = data.frame(outputRawAUCspec, bestvalue=AUCspec)
      # 
      # ##vetores vazios (para os nomes dos melhores modelos)
      # bestModelRunTSS = vector()
      # bestModelRunAUC = vector()
      # TSSvalues = vector()
      # AUCvalues = vector()
      # 
      # for (l in 1:nrow(tabBestScoresTSS)){
      #   bestRunNameTSS = names(tabBestScoresTSS)[match(tabBestScoresTSS$bestvalue[l], tabBestScoresTSS[l, 1:ncol(tabBestScoresTSS)-1])] #pegando o nome do melhor modelo (e.g. RUN2)
      #   TSSvalueBestModel = outputRawTSSvalue[l,bestRunNameTSS] #pegando o TSS do melhor modelo
      #   bestModelRunTSS = append(bestModelRunTSS, bestRunNameTSS) #empilhando os nomes dos melhores modelos
      #   TSSvalues = append(TSSvalues, TSSvalueBestModel)
      #   ##
      #   bestRunNameAUC = names(tabBestScoresAUC)[match(tabBestScoresAUC$bestvalue[l], tabBestScoresAUC[l, 1:ncol(tabBestScoresAUC)-1])] #pegando o nome do melhor modelo (e.g. RUN2)
      #   AUCvalueBestModel = outputRawAUCvalue[l,bestRunNameAUC] #pegando o AUC do melhor modelo
      #   bestModelRunAUC = append(bestModelRunAUC, bestRunNameAUC)
      #   AUCvalues = append(AUCvalues, AUCvalueBestModel)
      # }
      # 
      # ##medias
      # meanTSSValue = as.numeric(apply(outputRawTSSvalue, 1, mean, na.rm=TRUE))
      # meanTSSspecificity = as.numeric(apply(outputRawTSSspec, 1, mean, na.rm=TRUE))
      # meanAUCValue = as.numeric(apply(outputRawAUCvalue, 1, mean, na.rm=TRUE))
      # meanAUCspecificity = as.numeric(apply(outputRawAUCspec, 1, mean, na.rm=TRUE))
      # 
      # ##gravando estatisticas basicas do melhor modelo de cada algoritmo
      # statResultsSDMnormal = rbind(statResultsSDMnormal,
      #                                data.frame(SDM = 'normal',
      #                                           sampleSize = sampleSizes[j],
      #                                           sp = i,
      #                                           ##tss
      #                                           model = rownames(tabBestScoresTSS),
      #                                           meanTSSValue = meanTSSValue,
      #                                           meanTSSspecificity = meanTSSspecificity,
      #                                           maxTSSspecificity = TSSspec,
      #                                           bestModelTSS = bestModelRunTSS,
      #                                           TSSvalue_bestModel= TSSvalues,
      #                                           ##auc
      #                                           meanAUCValue = meanAUCValue,
      #                                           meanAUCspecificity = meanAUCspecificity,
      #                                           maxAUCspecificity = AUCspec,
      #                                           bestModelAUC = bestModelRunAUC,
      #                                           AUCvalue_bestModel= AUCvalues,
      #                                           stringsAsFactors = FALSE)
      # )
      
      statResultsSDMnormal = makeOutput(evaluationScores, statResultsSDMnormal, i, j, 'normal', sampleSizes[j])
      
      write.csv(statResultsSDMnormal, file=paste(projectFolder,'/StatisticalResults_SDMnormal.csv',sep=''), row.names=FALSE)
      
      ##selecao do modelo de maior sensibilidade
      
      ## tssMax = max(as.numeric(statResultsSDMnormal[which(statResultsSDMnormal$sp == i),]$TSSspec))
      ## aucMax = max(as.numeric(statResultsSDMnormal[which(statResultsSDMnormal$sp == i),]$AUCspec))
      ## bestAlgorithmTSS = statResultsSDMnormal[which(statResultsSDMnormal$TSSspec==tssMax),]$model
      ## bestAlgorithmAUC = statResultsSDMnormal[which(statResultsSDMnormal$TSSspec==tssMax),]$model
      ## ind = match(bestAlgorithmTSS, bestAlgorithmAUC)
      ## modelToProj = bestAlgorithmAUC[ind]
      ## modelNames = grep(pattern=paste(modelToProj,collapse='|'), x=myBiomodModelOut@models.computed, value=TRUE)
      
      # ##maior TSS e AUC para especie da iteracao atual
      # tssMax = max(TSSspec)
      # aucMax = max(AUCspec)
      # 
      # ##formacao de vetor com nome e RUN do melhor modelo (TSS)
      # bestAlgorithmTSS = statResultsSDMnormal[which(statResultsSDMnormal$maxTSSspecificity==tssMax),]$model
      # bestRunTSS = statResultsSDMnormal[which(statResultsSDMnormal$maxTSSspecificity==tssMax),]$bestModelTSS
      # patternsTSS = paste(bestRunTSS,bestAlgorithmTSS,sep='_')
      # 
      # ##formacao de vetor com nome e RUN do melhor modelo (AUC)
      # bestAlgorithmAUC = statResultsSDMnormal[which(statResultsSDMnormal$maxAUCspecificity==aucMax),]$model
      # bestRunAUC = statResultsSDMnormal[which(statResultsSDMnormal$maxAUCspecificity==aucMax),]$bestModelAUC
      # patternsAUC = paste(bestRunAUC,bestAlgorithmAUC,sep='_')
      # 
      # ##nomes dos melhores modelos
      # modelNames = grep(pattern=paste(c(patternsTSS,patternsAUC),collapse='|'), x=myBiomodModelOut@models.computed, value=TRUE)
      
      modelNames = bestModel(evaluationScores, statResultsSDMnormal, myBiomodModelOut)
      rm(evaluationScores)
      
      ## if (bestAlgorithmTSS == bestAlgorithmAUC){
      ##     modelToProj = bestAlgorithmTSS
      ## }else{
      ##     modelToProj = c(bestAlgorithmAUC, bestAlgorithmTSS)
      ## }
      ## ind = grep(pattern=modelToProj, x=myBiomodModelOut@models.computed) ##pegando os indices
      ## modelNames = myBiomodModelOut@models.computed[ind] ##pegando os nomes dos modelos para projecao (aqueles com maior especificidade)
      
      ##rodando algortmo de projecao (i.e. rodando a projecao)
      myBiomodProj <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = myExpl,
        proj.name = paste('sp',i,'_sample',sampleSizes[j],'_SDMnormal',sep=''),
        selected.models = modelNames,
        compress = 'TRUE',
        build.clamping.mask = 'TRUE')
      
      ##rodando o algoritmo de consenso dos modelos (i.e. ensemble model)
      myBiomodEM = BIOMOD_EnsembleModeling(
        modeling.output = myBiomodModelOut,
        chosen.models = modelNames)
      
      ##forecasting com o consenso dos algoritmos (i.e. ensemble projection)
      myBiomodEF = BIOMOD_EnsembleForecasting(
        EM.output = myBiomodEM,
        binary.meth = c('TSS','ROC'),
        projection.output = myBiomodProj)
      
      ##writeRaster(projStackBIN,file=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[j],'.replica',k,'/proj_',l,'kyr/proj_',i,'kyr','.sample',sampleSizes[j],'.replica',k,'_BIN.asc',sep=''),row.names=FALSE)
      
      ##projStackBIN = projStack>0.5  #BinaryTransformation(projStack,"10")
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "/n")})
  }
}



##PARTE 3: SDM com pseudoausencias melhoradas

for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    tryCatch({     
      
      ##diretorio para biomod salvar os resultados
      setwd(file.path(projectFolder,'SDMimproved'))
      
      ##definindo variaveis e parametros locais
      betterPseudo = list()
      betterPseudoVar = list()
      
      ##projecoes de ausencias do SDM (rodado na etapa 2, acima)
      binTSS = raster(paste(projectFolder,'/SDMnormal/','sp',i,'.sample',sampleSizes[j],'.SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal_sp',i,'.sample',sampleSizes[j],'.SDMnormal_ensemble_TSSbin.grd' ,sep=''))
      binAUC = raster(paste(projectFolder,'/SDMnormal/','sp',i,'.sample',sampleSizes[j],'.SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal_sp',i,'.sample',sampleSizes[j],'.SDMnormal_ensemble_ROCbin.grd' ,sep=''))
      projStackBIN = stack(binTSS,binAUC) #empilhando mapas binarios (feitos com threshold a partir do AUC e TSS)
      projAbs = sum(projStackBIN) #somando (para depois pegar areas ausencia que ambos os thresolds concordam)
      
      ## amostrando pontos diretamente das areas de ausencia (abaixo do threshold) obtidas na etapa 1 ##
      values(projAbs)[values(projAbs) != 0] = NA  #tranformando areas diferentes de zero em NA (retando somente os dados de ausencia)          
      betterPseudoPoints = dismo::randomPoints(mask=projAbs, n=1000) #sorteando pontos da distribuicao modelada
      
      ##betterPseudoDF = extract(projStackBIN[[k]], be tterPseudoPoints) #distinguindo entre occ e ausencia
      betterPseudoDF = data.frame(lon=betterPseudoPoints[,1], lat=betterPseudoPoints[,2], occ=0) #distinguindo entre occ e ausencia
      ##betterPseudo[[k]] =  data.frame(lon=betterPseudoPoints[,1], lat=betterPseudoPoints[,2], occ=betterPseudoDF) #data.frame
      ##betterPseudo[[k]] = betterPseudo[[k]][which(betterPseudo[[k]]$occ==0),] #excluindo as presencas
      betterPseudoVar = extract(predictors,betterPseudoDF[,c('lon','lat')]) #obtendo as variaveis preditoras nos pontos
      betterPseudoDF = data.frame(betterPseudoDF,betterPseudoVar) #motando dataset
      
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
                                           resp.name = paste(myRespName,'_sample',sampleSizes[j],'_SDMimproved',sep=''))
      
      ## ##inspecionando o objeto gerado pela funcao do biomod2
      ## myBiomodData
      ## plot(myBiomodData)
      
      ##parametrizando os modelos
      myBiomodOption <- BIOMOD_ModelingOptions(
        MAXENT.Phillips = list(path_to_maxent.jar=maxentFolder,
                               maximumiterations=1000,
                               linear=TRUE,
                               quadratic=TRUE,
                               product=FALSE,
                               threshold=FALSE,
                               hinge=FALSE,
                               maximumiterations=1000,
                               convergencethreshold=1.0E-5,
                               threads=2),
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
      rm(myBiomodModelOut) #removendo informacoes da etapa anterior
      myBiomodModelOut <- BIOMOD_Modeling(
        data = myBiomodData,
        models = c('MAXENT.Phillips','GLM', 'GAM', 'MARS', 'CTA', 'GBM', 'RF'),
        models.options = myBiomodOption,
        NbRunEval = 10,
        DataSplit = 75,
        models.eval.meth = c('TSS','ROC'),
        do.full.models = FALSE,
        modeling.id = paste(myRespName,'_sample',sampleSizes[j],'_SDMimproved',sep=''))
      
      
      ##My output data
      evaluationScores = get_evaluations(myBiomodModelOut)
      
      # ##especificidade
      # outputRawTSSspec = evaluationScores['TSS','Specificity',,,]
      # outputRawTSSspec = outputRawTSSspec[complete.cases(outputRawTSSspec),]
      # outputRawAUCspec = evaluationScores['ROC','Specificity',,,]
      # outputRawAUCspec = outputRawAUCspec[complete.cases(outputRawAUCspec),]
      # ##auc e tss
      # outputRawTSSvalue = evaluationScores['TSS','Testing.data',,,]
      # outputRawTSSvalue = outputRawTSSvalue[complete.cases(outputRawTSSvalue),]
      # outputRawAUCvalue = evaluationScores['ROC','Testing.data',,,]
      # outputRawAUCvalue = outputRawAUCvalue[complete.cases(outputRawAUCvalue),]
      # 
      # ##maior valore de especificidade de cada algoritmo implementado (tanto para TSS qto para AUC)
      # TSSspec = as.numeric(apply(outputRawTSSspec, 1, max, na.rm=TRUE))
      # AUCspec = as.numeric(apply(outputRawAUCspec, 1, max, na.rm=TRUE))
      # 
      # ##tabela auxiliar para obtencao das informacoes do melhor modelo
      # tabBestScoresTSS = data.frame(outputRawTSSspec, bestvalue=TSSspec)
      # tabBestScoresAUC = data.frame(outputRawAUCspec, bestvalue=AUCspec)
      # 
      # ##vetores vazios (para os nomes dos melhores modelos)
      # bestModelRunTSS = vector()
      # bestModelRunAUC = vector()
      # TSSvalues = vector()
      # AUCvalues = vector()
      # 
      # for (l in 1:nrow(tabBestScoresTSS)){
      #   bestRunNameTSS = names(tabBestScoresTSS)[match(tabBestScoresTSS$bestvalue[l], tabBestScoresTSS[l, 1:ncol(tabBestScoresTSS)-1])] #pegando o nome do melhor modelo (e.g. RUN2)
      #   TSSvalueBestModel = outputRawTSSvalue[l,bestRunNameTSS] #pegando o TSS do melhor modelo
      #   bestModelRunTSS = append(bestModelRunTSS, bestRunNameTSS) #empilhando os nomes dos melhores modelos
      #   TSSvalues = append(TSSvalues, TSSvalueBestModel)
      #   ##
      #   bestRunNameAUC = names(tabBestScoresAUC)[match(tabBestScoresAUC$bestvalue[l], tabBestScoresAUC[l, 1:ncol(tabBestScoresAUC)-1])] #pegando o nome do melhor modelo (e.g. RUN2)
      #   AUCvalueBestModel = outputRawAUCvalue[l,bestRunNameAUC] #pegando o AUC do melhor modelo
      #   bestModelRunAUC = append(bestModelRunAUC, bestRunNameAUC)
      #   AUCvalues = append(AUCvalues, AUCvalueBestModel)
      # }
      # 
      # ##medias
      # meanTSSValue = as.numeric(apply(outputRawTSSvalue, 1, mean, na.rm=TRUE))
      # meanTSSspecificity = as.numeric(apply(outputRawTSSspec, 1, mean, na.rm=TRUE))
      # meanAUCValue = as.numeric(apply(outputRawAUCvalue, 1, mean, na.rm=TRUE))
      # meanAUCspecificity = as.numeric(apply(outputRawAUCspec, 1, mean, na.rm=TRUE))
      # 
      # ##gravando estatisticas basicas do melhor modelo de cada algoritmo
      # statResultsSDMimproved = rbind(statResultsSDMimproved,
      #                                data.frame(SDM = 'improved',
      #                                           sampleSize = sampleSizes[j],
      #                                           sp = i,
      #                                           ##tss
      #                                           model = rownames(tabBestScoresTSS),
      #                                           meanTSSValue = meanTSSValue,
      #                                           meanTSSspecificity = meanTSSspecificity,
      #                                           maxTSSspecificity = TSSspec,
      #                                           bestModelTSS = bestModelRunTSS,
      #                                           TSSvalue_bestModel= TSSvalues,
      #                                           ##auc
      #                                           meanAUCValue = meanAUCValue,
      #                                           meanAUCspecificity = meanAUCspecificity,
      #                                           maxAUCspecificity = AUCspec,
      #                                           bestModelAUC = bestModelRunAUC,
      #                                           AUCvalue_bestModel= AUCvalues,
      #                                           stringsAsFactors = FALSE)
      # )
      
      statResultsSDMimproved = makeOutput(evaluationScores, statResultsSDMimproved, i, j, 'improved', sampleSizes[j])
      rm(evaluationScores)
      
      write.csv(statResultsSDMimproved, file=paste(projectFolder,'/StatisticalResults_SDMimproved','.csv',sep=''), row.names=FALSE)
      
      ##implementando projecoes do modelo
      
      ##rodando algortmo de projecao (i.e. rodando a projecao)
      rm(myBiomodProj) #removendo informacao da etapa anterior
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
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "/n")})
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

##obs: usar o merge() para recortar a planilha que for maior
##exemplo:
##newTab=merge(tabB,tabA,by=c('infoA','infoB'))[,-3] ##a terceira coluna e da planilha menor

statResultsSDMnormal = merge(statResultsSDMimproved,
           statResultsSDMnormal, 
           by=c('sampleSize','sp','model'))[,c('SDM.y',
                                               'sampleSize',
                                               'sp',
                                               'model', 
                                               'meanTSSValue.y', 
                                               'meanTSSspecificity.y', 
                                               'maxTSSspecificity.y', 
                                               'bestModelTSS.y', 
                                               'TSSvalue_bestModel.y', 
                                               'meanAUCValue.y', 
                                               'meanAUCspecificity.y', 
                                               'maxAUCspecificity.y', 
                                               'bestModelAUC.y', 
                                               'AUCvalue_bestModel.y')]

##ajustando so nomes
names(statResultsSDMnormal) = names(statResultsSDMimproved)

##chacagem rapida
for(i in c('MAXENT.Phillips','CTA',  'GAM',  'GBM',  'GLM',  'MARS', 'RF')){
  print(paste(i, 
    ' ==>> SDMnormal = ', dim(statResultsSDMnormal[statResultsSDMnormal$model==i,])[1], ' linhas e ', dim(statResultsSDMnormal[statResultsSDMnormal$model==i,])[2], ' colunas',
    ' // ',
    'SDMimproved = ', dim(statResultsSDMimproved[statResultsSDMimproved$model==i,])[1], ' linhas e ',dim(statResultsSDMimproved[statResultsSDMimproved$model==i,])[2], ' colunas',
    sep=''))
}


##boxplot AUC
jpeg(filename='boxplotAUC.jpeg')
boxplot(statResultsSDMnormal$AUCvalue_bestModel, statResultsSDMimproved$AUCvalue_bestModel, ylim=c(0,1), names=c('SDM normal','SDM improved'), ylab='AUC')
dev.off()

plot(density(statResultsSDMnormal$AUCvalue_bestModel, ylim=c(0,1), names=c('SDM normal','SDM improved'), ylab='TSS'), ylim=c(0,20),col='blue')
lines(density(statResultsSDMimproved$AUCvalue_bestModel, ylim=c(0,1), names=c('SDM normal','SDM improved'), ylab='TSS'),col='red')

##boxplot TSS
jpeg(filename='boxplotTSS.jpeg')
boxplot(statResultsSDMnormal$TSSvalue_bestModel, statResultsSDMimproved$TSSvalue_bestModel, ylim=c(0,1), names=c('SDM normal','SDM improved'), ylab='TSS')
dev.off()

plot(density(statResultsSDMnormal$TSSvalue_bestModel, ylim=c(0,1), names=c('SDM normal','SDM improved'), ylab='TSS'), ylim=c(0,11),col='blue')
lines(density(statResultsSDMimproved$TSSvalue_bestModel, ylim=c(0,1), names=c('SDM normal','SDM improved'), ylab='TSS'),col='red')

##AUC x sample size
jpeg(filename='AUC_&_sampleSize.jpeg')
plot(statResultsSDMnormal$AUCvalue_bestModel~statResultsSDMnormal$sampleSize, ylim=c(0,1), xlim=c(0,170), cex=2, pch=19, col=rgb(0,0,0,0.5), xlab='Sample size', ylab='AUC')
tendenciaSDMnormalAUC = lm(statResultsSDMnormal$AUCvalue_bestModel~statResultsSDMnormal$sampleSize)
abline(tendenciaSDMnormalAUC)
points(statResultsSDMimproved$AUCvalue_bestModel~statResultsSDMimproved$sampleSize, ylim=c(0,1), xlim=c(0,170), cex=1.5, pch=20, col=rgb(0,0,1,0.5))
tendenciaSDMimprovedAUC = lm(statResultsSDMimproved$AUCvalue_bestModel~statResultsSDMimproved$sampleSize)
abline(tendenciaSDMimprovedAUC, col='blue')
legend('bottomleft',legend=c('SDM normal','SDM improved'), pch=c(19,20), col=c(rgb(0,0,0,0.5),rgb(0,0,1,0.5)))
dev.off()

##TSS x sample size
jpeg(filename='TSS_&_sampleSize.jpeg')
plot(statResultsSDMnormal$TSSvalue_bestModel~statResultsSDMnormal$sampleSize, ylim=c(0,1), xlim=c(0,170), cex=2, pch=19, col=rgb(0,0,0,0.5), xlab='Sample size', ylab='TSS')
tendenciaSDMnormalTSS = lm(statResultsSDMnormal$TSSvalue_bestModel~statResultsSDMnormal$sampleSize)
abline(tendenciaSDMnormalTSS)
points(statResultsSDMimproved$TSSvalue_bestModel~statResultsSDMimproved$sampleSize, ylim=c(0,1), xlim=c(0,170), cex=1.5, pch=20, col=rgb(0,0,1,0.5))
tendenciaSDMimprovedTSS = lm(statResultsSDMimproved$TSSvalue_bestModel~statResultsSDMimproved$sampleSize)
abline(tendenciaSDMimprovedTSS, col='blue')
legend('bottomleft',legend=c('SDM normal','SDM improved'), pch=c(19,20), col=c(rgb(0,0,0,0.5),rgb(0,0,1,0.5)),bg='white')
dev.off()


##teste de significancia

wilcox.test(statResultsSDMimproved$AUCvalue_bestModel,statResultsSDMnormal$AUCvalue_bestModel) #resultado: p<0.05

wilcox.test(statResultsSDMimproved$TSSvalue_bestModel,statResultsSDMnormal$TSSvalue_bestModel)  #resultado:p<<0.05


## TSS e AUC por algoritmo

jpeg(filename='TSS_por_algoritmo.jpeg', width=800)
par(mfrow=c(1,2), las=2, mar=c(8,5,5,1))
boxplot(statResultsSDMnormal$TSSvalue_bestModel ~ statResultsSDMnormal$model, ylim=c(0,1), ylab='TSS', main='SDM normal')
boxplot(statResultsSDMimproved$TSSvalue_bestModel ~ statResultsSDMimproved$model, ylim=c(0,1), ylab='TSS', main='SDM improved')
dev.off()

jpeg(filename='AUC_por_algoritmo.jpeg', width=800)
par(mfrow=c(1,2), las=2, mar=c(8,5,5,1))
boxplot(statResultsSDMnormal$AUCvalue_bestModel ~ statResultsSDMnormal$model, ylim=c(0,1), ylab='AUC', main='SDM normal')
boxplot(statResultsSDMimproved$AUCvalue_bestModel ~ statResultsSDMimproved$model, ylim=c(0,1), ylab='AUC', main='SDM improved')
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
boxplot(statResultsSDMnormal$maxTSSspecificity ~ statResultsSDMnormal$model, ylab='Specificity', main='Maximization of TSS')
boxplot(statResultsSDMnormal$maxAUCspecificity ~ statResultsSDMnormal$model, ylab='Specificity', main='Maximization of AUC')
dev.off()

#teste especificade
kruskal.test(maxTSSspecificity ~ model, data=statResultsSDMnormal) #rsultado: p>0.05
kruskal.test(maxAUCspecificity ~ model, data=statResultsSDMnormal) #resultado:p>0.05


