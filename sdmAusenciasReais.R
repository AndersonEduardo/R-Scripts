##registrando hora do incio
timeOne = Sys.time()

##abrindo pacotes necessarios
library(raster)
library(rJava)
library(biomod2)
Sys.setenv(JAVA_HOME = "/usr/lib/jvm/java-7-openjdk-amd64")
options(java.parameters = "Xmx7g")

# ##definindo prametros e variaveis globais (NOTEBOOK)
# projectFolder = "/home/anderson/Documentos/Projetos/Improved pseudo-absences_TESTE" #pasta do projeto
# envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
# envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
# AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
# maxentFolder = '/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java/maxent.jar' #pasta para resultados do maxent
# ## spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
# ## sdmTypes = 'real'
# sampleSizes = 100  #c(10,20,40,80,160)
# NumRep = 2 #10 #numero de replicas (de cada cenario amostral)
# ##variaveis preditoras
# ## elevation = raster('/home/anderson/PosDoc/dados_ambientais/DEM/DEM.tif')
# predictors = stack(list.files(path=envVarPaths[1],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
# predictors = predictors[[c('bioclim_01','bioclim_12')]]
# predictors = stack(mask(x=predictors, mask=AmSulShape))
# crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
# Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
# statResultsSDMreal = data.frame() #tabela de estatisticas basicas do modelo

# ##definindo prametros e variaveis globais (LORIEN)
# projectFolder = "J:/Pesquisadorxs/Anderson_Eduardo/high_quality_PA" #pasta do projeto
# envVarFolder = "J:/Pesquisadorxs/Anderson_Eduardo/dados_projeto/000" #pasta com as variaveis ambientais
# envVarPaths = list.files(path=envVarFolder, pattern='.asc', full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
# AmSulShape = rgdal::readOGR("J:/Pesquisadorxs/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
# maxentFolder = 'C:/Users/WS/Documents/R/win-library/3.4/dismo/java/maxent.jar' #pasta para resultados do maxent
# ## spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
# ## sdmTypes = 'real'
# sampleSizes = c(10,20,40,80,160)
# NumRep = 10 #numero de replicas (de cada cenario amostral)
# ##variaveis preditoras
# elevation = raster('J:/Pesquisadorxs/Anderson_Eduardo/DEM/DEM.tif')
# predictors = stack(envVarPaths,elevation) #predictors com todas as variaveis (presente)
# predictors = predictors[[c('bioclim_01','bioclim_12')]]
# predictors = stack(mask(x=predictors, mask=AmSulShape))
# crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
# Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
# statResultsSDMreal = data.frame() #tabela de estatisticas basicas do modelo

##definindo prametros e variaveis globais (PC nupeg)
projectFolder = "D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas" #pasta do projeto
envVarFolder = "D:/Anderson_Eduardo/variaveis ambientais AmSul 120kyr/000" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, pattern='.asc', full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("D:/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
maxentFolder = 'D:/Anderson_Eduardo/maxent/maxent.jar' #pasta para resultados do maxent
sampleSizes = c(20,40,80,160)
NumRep = 10 #numero de replicas (de cada cenario amostral)
##variaveis preditoras
#elevation = raster('J:/Pesquisadorxs/Anderson_Eduardo/DEM/DEM.tif')
predictors = stack(envVarPaths) #predictors com todas as variaveis (presente)
predictors = predictors[[c('bioclim_01','bioclim_12')]]
predictors = stack(mask(x=predictors, mask=AmSulShape))
crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
statResultsSDMreal = data.frame() #tabela de estatisticas basicas do modelo

## funcoes autorais :-)
source('D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas/makeSpecies.R')
source('D:/Anderson_Eduardo/rangeByAC.R')
source('D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas/makeOutput.R')
source('D:/Anderson_Eduardo/SDM com pseudoausencias melhoradas/bestModel.R')


##definindo diretorio de trabalho (importante porque o biomod2 salva tudo automaticamente)
setwd(projectFolder)

##PARTE 2: modelando ausencias 

for(i in 1:Nsp){
  for(j in 1:length(sampleSizes)){
    tryCatch({
      
      ##diretorio para o biomod2 salvar resultados para SDMreal
      setwd(file.path(projectFolder,'ausencias_reais'))
      
      ##pontos de ocorrencia
      SpDistAC = raster(paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''))
      values(SpDistAC)[values(SpDistAC)==0] = NA
      occPoints = dismo::randomPoints(mask=SpDistAC, n=sampleSizes[j]) #sorteando pontos da distribuicao modelada
      occPoints = data.frame(lon=occPoints[,1],lat=occPoints[,2], occ=1)
      
      ## amostrando pontos das areas de ausencia real
      SpDistAC = raster(paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''))
      values(SpDistAC)[values(SpDistAC)==0] = NA
      occAreaPts = as.data.frame(SpDistAC, xy=TRUE, na.rm=TRUE)
      backgroundArea = raster(nrows=nrow(SpDistAC), ncols=ncol(SpDistAC),xmn=extent(SpDistAC)[1], xmx=extent(SpDistAC)[2], ymn = extent(SpDistAC)[3], ymx=extent(SpDistAC)[4] )
      values(backgroundArea) = 0
      backgroundArea = mask(x=backgroundArea, mask=AmSulShape)
      absencesPts = dismo::randomPoints(mask=backgroundArea, n=1000, p=occAreaPts[,1:2], excludep=TRUE)
      absencesPts = data.frame(absencesPts, occ=0)
      names(absencesPts) = names(occPoints)

      dataSet = data.frame(rbind(occPoints,absencesPts)) #planilha de dados no formato SWD
      dataSetVars = extract(predictors, dataSet[,c('lon','lat')], na.rm=TRUE)
      dataSet = data.frame(dataSet, dataSetVars)
      
      
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
        modeling.id = paste(myRespName,'_sample',sampleSizes[j],'_SDMreal',sep=''))
      
      ##My output data
      evaluationScores = get_evaluations(myBiomodModelOut)
      
      statResultsSDMreal = makeOutput(evaluationScores, statResultsSDMreal, i, j, 'real', sampleSizes[j])
      
      write.csv(statResultsSDMreal, file=paste(projectFolder,'/StatisticalResults_SDMreal.csv',sep=''), row.names=FALSE)
      
      modelNames = bestModel(evaluationScores, statResultsSDMreal, myBiomodModelOut)
      rm(evaluationScores)
      
      ##rodando algortmo de projecao (i.e. rodando a projecao)
      myExpl = predictors
      
      myBiomodProj <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = myExpl,
        proj.name = paste('sp',i,'_sample',sampleSizes[j],'_SDMreal',sep=''),
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

      xx = get_predictions(myBiomodEF)
      plot(xx)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "/n")})
  }
}


##PARTE 2: graficos e analise dos resultados


##ambiente para salvar os graficos
setwd("/home/anderson/Documentos/Projetos/Improved pseudo-absences/Resultados completos - versao final/graficos")
# statResultsSDMnormal = read.csv(paste(projectFolder,'/normal','.csv',sep=''), header=TRUE)

##planilhas de dados
statResultsSDMimproved = read.csv(paste(projectFolder,'/StatisticalResults_SDMimproved','.csv',sep=''), header=TRUE)
statResultsSDMnormal = read.csv(paste(projectFolder,'/StatisticalResults_SDMnormal','.csv',sep=''), header=TRUE)
statResultsSDMreal = read.csv(paste(projectFolder,'/StatisticalResults_SDMreal','.csv',sep=''), header=TRUE)

## deixando so os cenarios que rodaram para os dois SDMs
dim(statResultsSDMnormal)
dim(statResultsSDMimproved)
dim(statResultsSDMreal)

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

statResultsSDMreal = merge(statResultsSDMimproved,
                           statResultsSDMreal, 
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
names(statResultsSDMreal) = names(statResultsSDMimproved)

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
boxplot(statResultsSDMreal$AUCvalue_bestModel,statResultsSDMnormal$AUCvalue_bestModel, statResultsSDMimproved$AUCvalue_bestModel, ylim=c(0,1), names=c('SDM real','SDM normal','SDM improved'), ylab='AUC')
dev.off()

jpeg(filename='densidadeAUC.jpeg')
plot(density(statResultsSDMnormal$AUCvalue_bestModel), ylim=c(0,40),col='blue', lwd=3, xlab='AUC values', main="")
lines(density(statResultsSDMimproved$AUCvalue_bestModel),col='red', lwd=3)
lines(density(statResultsSDMreal$AUCvalue_bestModel),col='black', lwd=3)
legend(x='topleft', legend=c('SDM normal', 'SDM improved', 'SDM real'), lty=1, lwd=3, col=c('blue', 'red', 'black'))
dev.off()

##boxplot TSS
jpeg(filename='boxplotTSS.jpeg')
boxplot(statResultsSDMreal$TSSvalue_bestModel,statResultsSDMnormal$TSSvalue_bestModel, statResultsSDMimproved$TSSvalue_bestModel, ylim=c(0,1), names=c('SDM real','SDM normal','SDM improved'), ylab='TSS')
dev.off()

jpeg(filename='densidadeTSS.jpeg')
plot(density(statResultsSDMnormal$TSSvalue_bestModel), ylim=c(0,12), col='blue', lwd=3, xlab='TSS values', main="")
lines(density(statResultsSDMimproved$TSSvalue_bestModel),col='red', lwd=3)
lines(density(statResultsSDMreal$TSSvalue_bestModel),col='black', lwd=3)
legend(x='topleft', legend=c('SDM normal', 'SDM improved', 'SDM real'), lty=1, lwd=3, col=c('blue', 'red', 'black'))
dev.off()

##teste de significancia
dados = rbind(statResultsSDMreal, statResultsSDMnormal, statResultsSDMimproved)

kruskal.test(AUCvalue_bestModel~SDM,data=dados) #p << 0.05
kruskal.test(TSSvalue_bestModel~SDM,data=dados) #p << 0.05

##analise post-hoc - teste de Dunn (par a par)
FSA::dunnTest(AUCvalue_bestModel~SDM,data=dados[,]) ##diferencas signif entre todos os pares
FSA::dunnTest(TSSvalue_bestModel~SDM,data=dados[,]) ##diferencas signif entre todos os pares


##AUC x sample size
dados = rbind(statResultsSDMreal, statResultsSDMnormal, statResultsSDMimproved)

jpeg(filename='AUC_&_sampleSize.jpeg')
boxplot(AUCvalue_bestModel~SDM*sampleSize,col=c(rgb(0,0,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),data=dados,frame.plot=TRUE,axes=FALSE,ylim=c(0.2,1), xlab='Sample sizes', ylab='AUC')
axis(1,at=c(2,5,8,11),labels=c(20,40,60,160),cex.axis=0.7)
axis(2,at=c(0.2,0.4,0.6,0.8,1.0),labels=c(0.2,0.4,0.6,0.8,1.0),cex.axis=0.7)
abline(v=3.5); abline(v=6.5); abline(v=9.5)
legend('bottomright',legend=c('SDM real','SDM normal','SDM improved'),pch=15,col=c('black','blue','red'), cex=0.7)
dev.off()

jpeg(filename='AUC_&_sampleSize_densityPlot.jpeg')
par(mfrow=c(2,2), cex=1.1, mar=c(5,5,1,1))
#samplesize 20
plot( density(dados[dados$SDM=='real' & dados$sampleSize==20,]$AUCvalue_bestModel), col='black', xlab='AUC', main='Sample size 20')
lines( density(dados[dados$SDM=='improved' & dados$sampleSize==20,]$AUCvalue_bestModel), col='red')
lines( density(dados[dados$SDM=='normal' & dados$sampleSize==20,]$AUCvalue_bestModel) , col='blue')
legend('topleft',legend=c('SDM real','SDM normal','SDM improved'),pch=15,col=c('black','blue','red'), cex=0.5)
#samplesize 40
plot( density(dados[dados$SDM=='real' & dados$sampleSize==40,]$AUCvalue_bestModel), ylim=c(0,20), col='black', xlab='AUC', main='Sample size 40')
lines( density(dados[dados$SDM=='improved' & dados$sampleSize==40,]$AUCvalue_bestModel), col='red')
lines( density(dados[dados$SDM=='normal' & dados$sampleSize==40,]$AUCvalue_bestModel) , col='blue')
legend('topleft',legend=c('SDM real','SDM normal','SDM improved'),pch=15,col=c('black','blue','red'), cex=0.5)
#samplesize 80
plot( density(dados[dados$SDM=='real' & dados$sampleSize==80,]$AUCvalue_bestModel), ylim=c(0,40), col='black', xlab='AUC', main='Sample size 80')
lines( density(dados[dados$SDM=='improved' & dados$sampleSize==80,]$AUCvalue_bestModel), col='red')
lines( density(dados[dados$SDM=='normal' & dados$sampleSize==80,]$AUCvalue_bestModel) , col='blue')
legend('topleft',legend=c('SDM real','SDM normal','SDM improved'),pch=15,col=c('black','blue','red'), cex=0.5)
#samplesize 160
plot( density(dados[dados$SDM=='real' & dados$sampleSize==160,]$AUCvalue_bestModel), ylim=c(0,40), col='black', xlab='AUC', main='Sample size 160')
lines( density(dados[dados$SDM=='improved' & dados$sampleSize==160,]$AUCvalue_bestModel), col='red')
lines( density(dados[dados$SDM=='normal' & dados$sampleSize==160,]$AUCvalue_bestModel) , col='blue')
legend('topleft',legend=c('SDM real','SDM normal','SDM improved'),pch=15,col=c('black','blue','red'), cex=0.5)
dev.off()


##teste AUC x sample size
##sample size 20
kruskal.test(AUCvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==20),]) #resultado: p<0.05
FSA::dunnTest(AUCvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==20),]) #diferenca sig entre todos os pares
##sample size 40
kruskal.test(AUCvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==40),]) #resultado: p<0.05
FSA::dunnTest(AUCvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==40),]) #nao signf. entre normal e real
##sample size 80
kruskal.test(AUCvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==80),]) #resultado: p<0.05
FSA::dunnTest(AUCvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==80),]) #diferenca sig entre todos os pares
##sample size 20
kruskal.test(AUCvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==160),]) #resultado: p<0.05
FSA::dunnTest(AUCvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==160),]) #diferenca sig entre todos os pares

##SDMnormal
kruskal.test(AUCvalue_bestModel ~ sampleSize, data=dados[which(dados$SDM=='normal'),]) #resultado: p>0.05
FSA::dunnTest(AUCvalue_bestModel ~ sampleSize, data=dados[which(dados$SDM=='normal'),]) #nenhuma diferenca sig entre os pares
##SDMimproved
kruskal.test(AUCvalue_bestModel ~ sampleSize, data=dados[which(dados$SDM=='improved'),]) #resultado: p<0.05
FSA::dunnTest(AUCvalue_bestModel ~ sampleSize, data=dados[which(dados$SDM=='improved'),]) 
##SDMreal
kruskal.test(AUCvalue_bestModel ~ sampleSize, data=dados[which(dados$SDM=='real'),]) #resultado: p>0.05
FSA::dunnTest(AUCvalue_bestModel ~ sampleSize, data=dados[which(dados$SDM=='real'),]) 

##TSS x sample size
dados = rbind(statResultsSDMreal, statResultsSDMnormal, statResultsSDMimproved)

jpeg(filename='TSS_&_sampleSize.jpeg')
boxplot(TSSvalue_bestModel~SDM*sampleSize,col=c(rgb(0,0,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),data=dados,frame.plot=TRUE,axes=FALSE,ylim=c(0.2,1), xlab='Sample sizes', ylab='TSS')
axis(1,at=c(2,5,8,11),labels=c(20,40,60,160),cex.axis=0.7)
axis(2,at=c(0.2,0.4,0.6,0.8,1.0),labels=c(0.2,0.4,0.6,0.8,1.0),cex.axis=0.7)
abline(v=3.5); abline(v=6.5); abline(v=9.5)
legend('bottomright',legend=c('SDM real','SDM normal','SDM improved'),pch=15,col=c('black','blue','red'), cex=0.7)
dev.off()

##teste TSS x sample size
##sample size 20
kruskal.test(TSSvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==20),]) #resultado: p<0.05
FSA::dunnTest(TSSvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==20),]) #diferenca sig entre todos os pares
##sample size 40
kruskal.test(TSSvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==40),]) #resultado: p<0.05
FSA::dunnTest(AUCvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==40),]) #diferenca sig entre o SDM improved e os outros
##sample size 80
kruskal.test(TSSvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==80),]) #resultado: p>0.05
FSA::dunnTest(AUCvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==80),]) #diferenca sig entre todos os pares
##sample size 20
kruskal.test(TSSvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==160),]) #resultado: p>0.05
FSA::dunnTest(TSSvalue_bestModel ~ SDM, data=dados[which(dados$sampleSize==160),]) #diferenca sig entre todos os pares


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
FSA::dunnTest(TSSvalue_bestModel ~ model, data=statResultsSDMnormal, method="bh") 

kruskal.test(TSSvalue_bestModel ~ model, data=statResultsSDMimproved) #resultado: p<<0.05
##analise post-hoc - teste de Dunn (par a par)
FSA::dunnTest(TSSvalue_bestModel ~ model, data=statResultsSDMimproved, method="bh") 

##teste AUC
kruskal.test(AUCvalue_bestModel ~ model, data=statResultsSDMnormal) #resultado: p<0.05
##analise post-hoc - teste de Dunn (par a par)
FSA::dunnTest(AUCvalue_bestModel ~ model, data=statResultsSDMnormal, method="bh") 

kruskal.test(AUCvalue_bestModel ~ model, data=statResultsSDMimproved) #resultado: p<<0.05
##analise post-hoc - teste de Dunn (par a par)
FSA::dunnTest(AUCvalue_bestModel ~ model, data=statResultsSDMimproved, method="bh") 


## especificidade (escolha do melhor modelo)

jpeg(filename='especificidade_por_algoritmo.jpeg', width=800)
par(mfrow=c(1,2), las=2, mar=c(8,5,5,1))
boxplot(statResultsSDMnormal$maxTSSspecificity ~ statResultsSDMnormal$model, ylab='Specificity', main='Maximization of TSS')
boxplot(statResultsSDMnormal$maxAUCspecificity ~ statResultsSDMnormal$model, ylab='Specificity', main='Maximization of AUC')
dev.off()

#teste especificade
kruskal.test(maxTSSspecificity ~ model, data=statResultsSDMnormal) #rsultado: p>0.05
kruskal.test(maxAUCspecificity ~ model, data=statResultsSDMnormal) #resultado:p>0.05
