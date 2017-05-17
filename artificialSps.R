#########################################################################################
####SCRIPT PARA FAZER A DISTRIBUICAO REAL DAS ESPECIES EM VARIOS MOMENTOS DO TEMPO####
#########################################################################################
library(virtualspecies)
library(maptools)
library(dismo)
library(raster)
library(phyloclim) #para funcao niche.overlap()
source("/home/anderson/R/R-Scripts/TSSmaxent.R")get
source("/home/anderson/R/R-Scripts/AUCrand.R")

###PRIMEIRA PARTE: criando sps virtuais###

###Parametros necessarios###
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
caminhosCamadasTemp = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
############################

for (i in 1:length(caminhosCamadasTemp)){
    predictors = stack(paste(caminhosCamadasTemp[i],'/bioclim_01.asc',sep=''),paste(caminhosCamadasTemp[i],'/bioclim_12.asc',sep='')) #carregando as variaveis ambientais
    predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
    nameScenario = basename(caminhosCamadasTemp[i])
    
    #Especie de clima quente e umido
    parametersHW <- formatFunctions(bioclim_01=c(fun='betaFun',p1=200,p2=295,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=2000,p2=3500,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais

    #Especie de clima quente e seco
    parametersHD <- formatFunctions(bioclim_01=c(fun='betaFun',p1=200,p2=260,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=50,p2=1800,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais
 
    #Especie de clima frio e seco
    parametersCD <- formatFunctions(bioclim_01=c(fun='betaFun',p1=50,p2=220,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=50,p2=1800,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais

    spHW <- generateSpFromFun(predictors, parametersHW) #criando a especie artifical (clima quente e umido)
    spHD <- generateSpFromFun(predictors, parametersHD) #criando a especie artifical (clima quente e seco)
    spCD <- generateSpFromFun(predictors, parametersCD) #criando a especie artifical (clima frio e seco)

    auxVector=stack(c(spHW$suitab.raster,spHD$suitab.raster,spCD$suitab.raster))
    names(auxVector) = c('spHW', 'spHD', 'spCD')
    
    for(j in 1:dim(auxVector)[3]){
        projection(auxVector[[j]]) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
        writeRaster(auxVector[[j]], filename=paste(projectFolder,'NichoReal/',names(auxVector[[j]]),'/',nameScenario,'.asc',sep=""), overwrite=TRUE,prj=TRUE) #salvando o raster do mapa da sp
    }
}

###SEGUNDA PARTE: amostragem de pontos de ocorrencia em diferentes camadas de tempo###

envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto/" #pasta com as variaveis ambientais
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha com os pontos amostrados sera salva
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
biomodFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/biomod/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
sampleSizes = c(5,15,25,35,45,55,65,75,85,95) #tamanhos das amostras
NumRep = 10 #numero de replicas (de cada cenario amostral)
Tmax = 22 #idade maxima (no passado)
bgPoints = 1000 #numero de pontos de background
sampleData = data.frame()


#ocorrencias
for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
    for (sSize in sampleSizes){ #numero de pontos (registros, dados) na amostra
        sampledAges = vector()
        sampledAges = round(runif(sSize,0,Tmax)) #selecionando 'n' camadas de tempo aleatoriamente
        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da
        for (j in 1:NumRep){ #replicas do cenario amostral
            for (sAge in sampledAges){ #amostrando em cada camada de tempo que consta na amostra
                sampleData_i = sampleOccurrences(x=raster(nicheRealPath[sAge+1]),n=1,plot=FALSE)$sample.points[,1:2] #amostra d ponto
                scenarioName = basename(nicheRealPath[1:24][sAge+1]) #tempo vinculado ao cenario para variaveis ambientais
                scenarioName = gsub('.asc','',scenarioName) #retirando do nome o '.asc'
                layers_i = extract(
                    x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                    y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
                sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
            }
            names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
            write.csv(sampleData,paste(projectFolder,'Amostras/',spsTypes[i],'/occ',sSize,'pts', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
            sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
        }
    }
}

#background points
for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
    for (sSize in sampleSizes){ #numero de pontos de ocoorrencia (que parearah com a amostra de background dessa  iteracao)
        sampledAges = vector()
        sampledAges = round(runif(bgPoints,0,Tmax)) #selecionando 'n' camadas de tempo aleatoriamente ('n' = bgPoints)
#        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        for (j in 1:NumRep){ #replicas do cenario amostral
            for (sAge in sampledAges){ #amostrando em cada camada de tempo que consta na amostra

                envVarPath = list.files(path=paste(envVarFolder,list.files(path=paste(envVarFolder))[sAge+1],sep=''), full.names=TRUE, pattern='.asc') #lista com os enderecos das variaveis ambientais no tempo corresposndente a interacao 
                
                sampleData_i = randomPoints(mask=raster(envVarPath[1],crs=CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),n=1) #amostra do ponto
                scenarioName = list.files(path=paste(envVarFolder))[sAge+1] #nome do cenario
                layers_i = extract(
                    x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                    y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
                sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
            }
            names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
            write.csv(sampleData,paste(projectFolder,'Amostras/',spsTypes[i],'/bg',sSize,'pts', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
            sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
        }
    }
}

########
########
#backgroung points 2
for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
  for (j in 0:(length(envVarPaths[1:24])-1)){ #amostrando em cada camada de tempo que consta na amostra
      envVarPath = list.files(path=paste(envVarFolder,list.files(path=paste(envVarFolder))[j+1],sep=''), full.names=TRUE, pattern='.asc') #lista com os enderecos das variaveis ambientais no tempo corresposndente a interacao 
      sampleData_j = randomPoints(mask=raster(envVarPath[1],crs=CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),n=500) #amostra do ponto
      scenarioName = list.files(path=paste(envVarFolder))[j+1] #nome do cenario
      layers_j = extract(
        x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='asc', full.names=TRUE)),
        y=sampleData_j) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
        sampleData = rbind(sampleData, cbind(sampleData_j,layers_j,j)) #juntando com os dados das outras camadas de tempo amostradas
  }
  names(sampleData) = c('lon','lat',names(as.data.frame(layers_j)),'kyrBP') #ajustando os nomes
  write.csv(sampleData,paste(projectFolder,'Amostras/',spsTypes[i],'/bg.csv',sep=''),row.names=FALSE) #salvando
  sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
}
########
########

#ausencias reais
for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
  for (sSize in sampleSizes){ #numero de pontos de ocoorrencia (que parearah com a amostra de background dessa  iteracao)
    sampledAges = vector()
    sampledAges = round(runif(bgPoints,0,Tmax)) #selecionando 'n' camadas de tempo aleatoriamente ('n' = bgPoints)
    #        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
    for (j in 1:NumRep){ #replicas do cenario amostral
      for (sAge in sampledAges){ #amostrando em cada camada de tempo que consta na amostra
        envVarPath = list.files(path=paste(envVarFolder,list.files(path=paste(envVarFolder))[sAge+1],sep=''), full.names=TRUE, pattern='.asc') #lista com os enderecos das variaveis ambientais no tempo corresposndente a interacao 
        sampleData_i = sampleOccurrences(x=raster(nicheRealPath[sAge+1]),n=1000,type='presence-absence',sample.prevalence=0,plot=FALSE)$sample.points[,1:2] #amostra d ponto
        scenarioName = list.files(path=paste(envVarFolder))[sAge+1] #nome do cenario
        layers_i = extract(
          x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='asc', full.names=TRUE)),
          y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
        sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
      }
      names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
      write.csv(sampleData,paste(projectFolder,'Amostras/',spsTypes[i],'/ausences',sSize,'pts', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
      sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
    }
  }
}


###TERCEIRA PARTE: SDM usando de pontos de ocorrencia em diferentes camadas de tempo (do atual ate 120 kyr BP)###


#######################################################
####################### MAXENT ########################
#######################################################

options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha 
maxentFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
source("/home/anderson/R/R-Scripts/TSSmaxent.R")
evaluation = list()
TSSvector = data.frame()
sampleSizes = c(5,15,25,35,45,55,65,75,85,95)
NumRep = 10 #numero de replicas (de cada cenario amostral)

for (i in 1:length(spsTypes)){
    for (j in sampleSizes){
        for (k in 1:NumRep){ #loop sobre o numero de replicas 
        
            occPoints = read.csv(paste(mainSampleFolder,spsTypes[i],'/occ',j,'pts',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia
            ##backgroundPoints = read.csv(paste(mainSampleFolder,spsTypes[i],'/bg',j,'pts',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de background
            backgroundPoints = read.csv(paste(mainSampleFolder,spsTypes[i],'/bg.csv',sep=''),header=TRUE) #abrindo pontos de background
            backgroundPoints = backgroundPoints[sample(nrow(backgroundPoints),5000),]
            ##realAusences = read.csv(paste(mainSampleFolder,spsTypes[i],'/ausences',j,'pts',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de background
            names(backgroundPoints) = names(occPoints) #certificando que os nomes das colunas estão iguais (cuidado aqui...)
            dataSet = data.frame(cbind(rbind(occPoints,backgroundPoints),pres=c(rep(1,nrow(occPoints)),rep(0,nrow(backgroundPoints))))) #planilha de dados no formato SWD
            ##dataSetRealAusences = data.frame(cbind(rbind(occPoints,realAusences),pres=c(rep(1,nrow(occPoints)),rep(0,nrow(realAusences))))) #planilha de dados no formato SWD
            
            me = maxent(
                ##x=dataSet[,c("bioclim_01","bioclim_04","bioclim_10","bioclim_11","bioclim_12","bioclim_15","bioclim_16","bioclim_17")],
                x=dataSet[,c("bioclim_01","bioclim_12")],
                p=dataSet$pres,
                path=paste(maxentFolder,spsTypes[i],'sampleSize',j,sep=''),
                removeDuplicates = TRUE,
                args=c('responsecurves=TRUE',
                       'jackknife=TRUE',
                       'randomseed=TRUE',
                       'randomtestpoints=25',
                       'maximumbackground=5000',
                       'replicates=10',
                       'replicatetype=subsample',
                       'writebackgroundpredictions=TRUE',
                       'linear=TRUE',
                       'quadratic=TRUE',
                       'product=FALSE',
                       'threshold=FALSE',
                       'hinge=FALSE',
                       'maximumiterations=1000',
                       'convergencethreshold=1.0E-5',
                       'threads=2'
                       ))

            ##rodando a avaliacao do modelo
            TSSvector = TSSmaxent(paste(maxentFolder,spsTypes[i],'/sampleSize',j,'/',sep=''))
            
            thres = as.data.frame(read.csv(paste(maxentFolder,spsTypes[i],'/sampleSize',j,'/maxentResults.csv',sep=''),header=TRUE))$X10.percentile.training.presence.Cloglog.threshold[11] #threshold 10 percentile training occ. OBS: '11' e o numero da linha no arquivo maxntResults.csv em q esta a media
                  
            AUCrandVector = AUCrand(x=dataSet[,c("bioclim_01","bioclim_12")], p=dataSet$pres, path=paste(maxentFolder,spsTypes[i],sep=''), args=c('randomseed=TRUE','randomtestpoints=25','maximumbackground=5000','replicates=3','replicatetype=subsample','writebackgroundpredictions=TRUE','linear=TRUE','quadratic=TRUE','product=FALSE','threshold=FALSE','hinge=FALSE','maximumiterations=1000','convergencethreshold=1.0E-5','threads=6'))
            
            pValue = sum(as.numeric(AUCrandVector) >=  TSSvector$tesAUC) / length(AUCrandVector)
            
            statResults = rbind(statResults,cbind(sp=spsTypes[i],sampleSize=j,replicate=k,AUC=TSSvector$tesAUC,p_value=pValue,TSS=TSSvector$TSS,Threshold=thres,numbOfTimeLayers=length(unique(occPoints$kyrBP)),medianKyr=median(occPoints$kyrBP),minAge=min(occPoints$kyrBP),maxAge=max(occPoints$kyrBP)))
            
            write.csv(statResults,file=paste(projectFolder,'Maxent/',spsTypes[i],'/StatisticsResults-',spsTypes[i],'.csv',sep=''),row.names=FALSE)
            
            for (l in 1:length(envVarPaths[1:24])){
                predictors = stack(list.files(path=envVarPaths[l],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
                predictors = predictors[[c('bioclim_01','bioclim_12')]]
                predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
                crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
                proj = predict(me,predictors,crs=crs) #realizando projetacoes (para cada replica)
                writeRaster(mean(proj),paste(maxentFolder,spsTypes[i],'/projections/projection-Time',l-1,'kyrBP','-Replica',k,'-Sample',j,'.asc',sep=''),overwrite=TRUE) #salvando a projecao media
                
            }            
        }

        ##criando um mapa binario
        ## thresholdValues = NULL
        ## aucValues = NULL
        ## for (m in 1:length(evaluation)){
        ##     thresholdValues <- append(thresholdValues, threshold(evaluation[[m]],'spec_sens'))
        ##     aucValues = append(aucValues, evaluation[[m]]@auc)
        ## }
        ## aucMean = mean(aucValues)
        ## thresholdMean = mean(thresholdValues)
        ## TSSmean = mean(TSSvector$TSS)
        ## statResults = data.frame(sp=spsTypes[i],AUCmean=aucMean,TSSmean=TSSmean,ThresholdMean=thresholdMean)
        ## write.csv(statResults,file=paste(projectFolder,'Maxent/',spsTypes[i],'/StatisticsResults-',spsTypes[i],'.csv',sep=''),row.names=FALSE)
        
        ##esvaziando o vetor para a proxima especie
        evaluation = list()
        TSSvector = data.frame()
    }
}
#######################################################
#######################################################
#######################################################

###QUARTA PARTE: comparando projecao do SDM e a distribuicao espacial real do nicho da sp###

###
## clim1 = na.exclude(as.data.frame(predictors1,xy=TRUE)) #em que predictors e um stack com as variaveis bioclim01 e bioclim12 para 10 kyrBP
## clim2 = na.exclude(as.data.frame(predictors2,xy=TRUE)) #em que predictors e um stack com as variaveis bioclim01 e bioclim12 para 10 kyrBP
## clim12 = rbind(clim1,clim2) #dados ambientais para todo o espaco estudado
## spOcc1 = occPoints #pontos de ocorrencia com dados para as variaveis ambientais
## spOcc2 = occPoints #pontos de ocorrencia com dados para as variaveis ambientais

## scores.clim12.MAXENT <- data.frame(dismo::predict(object=me@models[[1]], x=clim12[,-c(1,2)]))
## scores.clim1.MAXENT <- data.frame(dismo::predict(object=me@models[[1]], x=clim1[,-c(1,2)]))
## scores.clim2.MAXENT <- data.frame(dismo::predict(object=me@models[[1]], x=clim2[,-c(1,2)]))
## scores.sp1.MAXENT <- data.frame(dismo::predict(object=me@models[[1]], x=spOcc1[,c("bioclim_01","bioclim_12")]))
## scores.sp2.MAXENT <- data.frame(dismo::predict(object=me@models[[1]], x=spOcc2[,c("bioclim_01","bioclim_12")]))

## R=100

## z1<- grid.clim(scores.clim12.MAXENT,scores.clim1.MAXENT,scores.sp1.MAXENT,R)
## z2<- grid.clim(scores.clim12.MAXENT,scores.clim2.MAXENT,scores.sp2.MAXENT,R)
## a<-niche.equivalency.test(z1,z2,rep=100)# test of niche equivalency and similarity according to Warren et al. 2008
## b<-niche.similarity.test(z1,z2,rep=100)
## b2<-niche.similarity.test(z2,z1,rep=100)

###

projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha com os pontos amostrados sera salva
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
sampleSizes = c(5,15,25,35,45,55,65,75,85,95) #aqui, deve ser igual ao usasado nas partes anteriores do script
NumRep = 9
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul

for (i in 1:length(spsTypes)){
  
  nicheRealFolder = paste(projectFolder,'/NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
  nicheRealPath = list.files(path=nicheRealFolder,pattern='asc',full.names=TRUE) #lista com os enderecos dos mapas de distribuicao da sp
  projectionsFolder = paste(projectFolder,'/maxent/',spsTypes[i],'/projections',sep='') #pasta com as projecoes do cenario
  projectionsPath = list.files(path=projectionsFolder, pattern='asc',full.names=T) #caminhos para os .asc na paste do cenario
  outputData = data.frame()

  
  for (l in 1:length(nicheRealPath[1:24])){ #loop sobre as cadamdas de tempo
    
    realNiche = nicheRealPath[l] #nicho real
    
    predictors = stack(list.files(path=envVarPaths[l],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
    predictors = predictors[[c('bioclim_01','bioclim_12')]]
    predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
    projection(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
    
    for (m in sampleSizes){ #loop sobre os tamanhos amostrais
      
      timeSampleData = list.files(path=projectionsFolder, pattern=glob2rx(paste('*Time',l-1,'*Sample',m,'.asc',sep='')),full.names=TRUE)
      #output_i = data.frame() #dataframe vazio para o loop abaixo
      
      for(n in 1:NumRep){ #loop sobre replicas de cada combinacao de tempo e tamanho amostral
        
        ##vetores vazios
        Ddistribution = numeric()
        Idistribution = numeric()
        
        sdmNiche = timeSampleData[[n]] #mapa de suitability gerado por SDM
        nicheOverlapObs = niche.overlap(c(sdmNiche,realNiche))
        Dobs = nicheOverlapObs[1,2]
        Iobs = nicheOverlapObs[2,1]
        
        ##aleatorizando ocorrencias para teste de significancia
        occPoints = read.csv(paste(mainSampleFolder,'/',spsTypes[i],'/occ',m,'pts',n,'rep.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia
        occPoints[occPoints==0] = NA
        occPoints = occPoints[complete.cases(occPoints),]
        occPoints = round(occPoints, digits=2)
        occPoints = occPoints[!duplicated(occPoints),]                 
        
        backgroundPoints = read.csv(paste(mainSampleFolder,'/',spsTypes[i],'/bg.csv',sep=''),header=TRUE) #abrindo pontos de background
        backgroundPoints = backgroundPoints[sample(nrow(backgroundPoints),5000),]
        backgroundPoints[backgroundPoints==0] = NA
        backgroundPoints = backgroundPoints[complete.cases(backgroundPoints),]
        backgroundPoints = round(backgroundPoints, digits=2)
        backgroundPoints = backgroundPoints[!duplicated(backgroundPoints),]                 
        
        names(backgroundPoints) = names(occPoints) #certificando que os nomes das colunas estão iguais (cuidado aqui...)
        dataSet = data.frame(cbind(rbind(occPoints,backgroundPoints),pres=c(rep(1,nrow(occPoints)),rep(0,nrow(backgroundPoints))))) #planilha de dados no formato SWD
        
        # for(o in 1:100){ #replicas da distribuicao nula de D e I
        #   
        #   presMix = dataSet[sample(nrow(dataSet)),c('pres')]
        #   sampleMix = dataSet
        #   sampleMix$pres = presMix
        #   
        #   me = maxent(
        #     x=sampleMix[,c("bioclim_01","bioclim_12")],
        #     p=sampleMix$pres,
        #     args=c('responsecurves=FALSE',
        #            'jackknife=FALSE',
        #            'randomseed=FALSE',
        #            'randomtestpoints=0',
        #            'maximumbackground=5000',
        #            'replicates=1',
        #            'writebackgroundpredictions=FALSE',
        #            'linear=TRUE',
        #            'quadratic=TRUE',
        #            'product=FALSE',
        #            'threshold=FALSE',
        #            'hinge=FALSE',
        #            'maximumiterations=1000',
        #            'convergencethreshold=1.0E-5',
        #            'threads=2'
        #     ))
        #   
        #   proj = predict(me,predictors,crs=crsInfo) #realizando projetacoes (para cada replica)                    
        #   
        #   sdmNicheMix = as(proj,'SpatialGridDataFrame')
        #   
        #   nicheOverlap_i= niche.overlap(c(sdmNicheMix,realNiche))
        #   Ddistribution = append(Ddistribution,nicheOverlap_i[1,2])
        #   Idistribution = append(Idistribution,nicheOverlap_i[2,1])
        #   
        # }
        
        DpValue = NA #sum(Ddistribution >= Dobs) / length(Ddistribution)
        IpValue = NA #sum(Idistribution >= Iobs) / length(Idistribution)
        
        outputData = rbind(outputData,cbind(kyrBP=l-1,sampleSize=m,replicate=n,numbOfTimeLayers=length(unique(occPoints$kyrBP)),medianKyr=median(unique(occPoints$kyrBP)),minAge=min(unique(occPoints$kyrBP)),maxAge=max(unique(occPoints$kyrBP)),Schoeners_D=Dobs,p_value=DpValue,Hellinger_I=Iobs,p_value=IpValue))
        
        write.csv(outputData, file=paste(projectFolder,'/maxent/',spsTypes[i],"/projections/NO.csv",sep=""),row.names=FALSE) #salvando os dados do cenario
        
      }
    }
  }
}


### QUINTA PARTE: construindo graficos dos resultados ###

spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
###abrindo as planilhas de dados
outputData = list()
vetor.nomes = vector()
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto

for (i in 1:length(spsTypes)){
    outputData[[i]] = read.csv(file=paste(projectFolder,'Maxent/',spsTypes[i],'/projections/NO.csv',sep=''),header=TRUE)
    vetor.nomes = append(vetor.nomes,paste(spsTypes[i],sep=''))
}
names(outputData) = vetor.nomes

###graficos

##Boxplots 

##HW
HWdataD = data.frame(kyrBP=outputData$spHW[outputData$spHW$sampleSize>=15,]$kyrBP,indexD=outputData$spHW[outputData$spHW$sampleSize>=15,]$Schoeners_D)
HWdataH = data.frame(kyrBP=outputData$spHW[outputData$spHW$sampleSize>=15,]$kyrBP,indexH=outputData$spHW[outputData$spHW$sampleSize>=15,]$Hellinger_I)

##HD
HDdataD = data.frame(kyrBP=outputData$spHD[outputData$spHD$sampleSize>=15,]$kyrBP,indexD=outputData$spHD[outputData$spHD$sampleSize>=15,]$Schoeners_D)
HDdataH = data.frame(kyrBP=outputData$spHD[outputData$spHD$sampleSize>=15,]$kyrBP,indexH=outputData$spHD[outputData$spHD$sampleSize>=15,]$Hellinger_I)

##CD
CDdataD = data.frame(kyrBP=outputData$spCD[outputData$spCD$sampleSize>=15,]$kyrBP,indexD=outputData$spCD[outputData$spCD$sampleSize>=15,]$Schoeners_D)
CDdataH = data.frame(kyrBP=outputData$spCD[outputData$spCD$sampleSize>=15,]$kyrBP,indexH=outputData$spCD[outputData$spCD$sampleSize>=15,]$Hellinger_I)

jpeg(file='/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/graficos/boxplots.jpg',width=900,height=600)
par(mfrow=c(1,2),mar=c(3,5,3,3))
##D
Ddata = rbind(data.frame(kyrBP=HWdataD$kyrBP,indexD=HWdataD$indexD,spsType='Hot and Wet'),data.frame(kyrBP=HDdataD$kyrBP,indexD=HDdataD$indexD,spsType='Hot and Dry'),data.frame(kyrBP=CDdataD$kyrBP,indexD=CDdataD$indexD,spsType='Cold and Dry'))
stripchart(indexD~spsType,data=Ddata,vertical=TRUE,method="jitter",pch=20,cex=1.5,col=rgb(0.5,0.5,0.5,0.05),ylim=c(0.2,1),ylab="Schoeners' D",cex.main=2,cex.lab=1.5,cex.axis=1.5)
boxplot(indexD~spsType,data=Ddata,boxwex=c(0.5,0.5,0.5),col=rgb(0,0,0,alpha=0.05),lwd=3,ylim=c(0.2,1),names=c('','',''),main="",cex.main=2,cex.lab=1.5,cex.axis=1.5,add=TRUE)
#H
Hdata = rbind(data.frame(kyrBP=HWdataH$kyrBP,indexH=HWdataH$indexH,spsType='Hot and Wet'),data.frame(kyrBP=HDdataH$kyrBP,indexH=HDdataH$indexH,spsType='Hot and Dry'),data.frame(kyrBP=CDdataH$kyrBP,indexH=CDdataH$indexH,spsType='Cold and Dry'))
stripchart(indexH~spsType,data=Hdata,vertical=TRUE,method="jitter",pch=20,cex=1.5,col=rgb(0.5,0.5,0.5,0.05),ylim=c(0.2,1),ylab="Hellinger I",cex.main=2,cex.lab=1.5,cex.axis=1.5,) 
boxplot(indexH~spsType,data=Hdata,boxwex=c(0.5,0.5,0.5),col=rgb(0,0,0,alpha=0.05),lwd=3,ylim=c(0.2,1),names=c('','',''),main="",cex.main=2,cex.lab=1.5,cex.axis=1.5,add=TRUE)
dev.off()

##teste de sidnnficancia (comparando as especies)

##kruskal-Wallis
##indice D
kruskal.test(indexD~spsType,data=Ddata)
##indice I
kruskal.test(indexH~spsType,data=Hdata)

##comparacoes par a par
##indice D
pairwise.wilcox.test(Ddata$indexD,Ddata$spsType,p.adjust.method='bonferroni')
##indice H
pairwise.wilcox.test(Hdata$indexH,Hdata$spsType,p.adjust.method='bonferroni')


##Schoeners' D ao longo do tempo##

jpeg(file='/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/graficos/metricsXtime.jpg',width=1400,height=750)
par(mar=c(5,6,5,5),mfrow=c(1,2))
#D
plot(indexD~as.factor(kyrBP),data=HWdataD,type='p',ylab="Schoeners' D",xlab="Time (kyr BP)",ylim=c(0,1),pch=1,col=rgb(0,0,0,alpha=0.5),cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
plot(indexD~as.factor(kyrBP),data=HDdataD,type='p',ylab="Schoeners' D",xlab="Time (kyr BP)",ylim=c(0,1),pch=2,col=rgb(0,0,1,alpha=0.5),cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2,add=TRUE)
plot(indexD~as.factor(kyrBP),data=CDdataD,type='p',ylab="Schoeners' D",xlab="Time (kyr BP)",ylim=c(0,1),pch=3,col=rgb(1,0,0,alpha=0.5),cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2,add=TRUE)
legend(x=20,y=0.2,legend=c('H&W','H&D','C&D'),pch=20,col=c('black','blue','red'),cex=1.8)
#H
plot(indexH~as.factor(kyrBP),data=HWdataH,type='p',ylab="Hellinger I",xlab="Time (kyr BP)",ylim=c(0,1),pch=1,col=rgb(0,0,0,alpha=0.5),cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
plot(indexH~as.factor(kyrBP),data=HDdataH,type='p',ylab="Hellinger I",xlab="Time (kyr BP)",ylim=c(0,1),pch=2,col=rgb(0,0,1,alpha=0.5),cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2,add=TRUE)
plot(indexH~as.factor(kyrBP),data=CDdataH,type='p',ylab="Hellinger I",xlab="Time (kyr BP)",ylim=c(0,1),pch=3,col=rgb(1,0,0,alpha=0.5),cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2,add=TRUE)
legend(x=20,y=0.2,legend=c('H&W','H&D','C&D'),pch=20,col=c('black','blue','red'),cex=1.8)
dev.off()

## ##testando significancia
## ##indice D
cor.test(HWdataD$kyrBP,HWdataD$indexD,method='spearman')
cor.test(HDdataD$kyrBP,HDdataD$indexD,method='spearman')
cor.test(CDdataD$kyrBP,CDdataD$indexD,method='spearman')
## ##indice I
cor.test(HWdataH$kyrBP,HWdataH$indexH,method='spearman')
cor.test(HDdataH$kyrBP,HDdataH$indexH,method='spearman')
cor.test(CDdataH$kyrBP,CDdataH$indexH,method='spearman')


##metricas X tamanho amostral

#maxent
jpeg(file='/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/graficos/metricsXsampleSize.jpg',width=1100,height=750)
par(mfrow=c(1,2),mar=c(5,6,4,4))
#D
plot(outputData$spHW$Schoeners_D~as.factor(outputData$spHW$sampleSize),type='p',ylim=c(0,1),xlab='Sample size',ylab="Shoeners' D",pch=1,col=rgb(0,0,0,alpha=0.5),cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
plot(outputData$spHD$Schoeners_D~as.factor(outputData$spHD$sampleSize),data=HDdataD,type='p',ylim=c(0,1),xlab=NULL,ylab=NULL,pch=2,col=rgb(0,0,1,alpha=0.5),cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2,add=TRUE)
plot(outputData$spCD$Schoeners_D~as.factor(outputData$spCD$sampleSize),data=CDdataD,type='p',ylim=c(0,1),xlab=NULL,ylab=NULL,pch=3,col=rgb(1,0,0,alpha=0.5),cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2,add=TRUE)
legend(x=7.5,y=0.2,legend=c('H&W','H&D','C&D'),pch=20,col=c('black','blue','red'),cex=1.7)
#H
plot(outputData$spHW$Hellinger_I~as.factor(outputData$spHW$sampleSize),type='p',ylim=c(0,1),xlab='Sample size',ylab="Hellinger I",pch=1,col=rgb(0,0,0,alpha=0.5),cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
plot(outputData$spHD$Hellinger_I~as.factor(outputData$spHD$sampleSize),data=HDdataD,type='p',ylim=c(0,1),xlab=NULL,ylab=NULL,pch=2,col=rgb(0,0,1,alpha=0.5),cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2,add=TRUE)
plot(outputData$spCD$Hellinger_I~as.factor(outputData$spCD$sampleSize),data=CDdataD,type='p',ylim=c(0,1),xlab=NULL,ylab=NULL,pch=3,col=rgb(1,0,0,alpha=0.5),cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2,add=TRUE)
legend(x=7.5,y=0.2,legend=c('H&W','H&D','C&D'),pch=20,col=c('black','blue','red'),cex=1.7)
dev.off()

##teste de significancia
##indice D
cor.test(outputData$spHW$sampleSize,outputData$spHW$Schoeners_D,method='spearman')
cor.test(outputData$spHD$sampleSize,outputData$spHD$Schoeners_D,method='spearman')
cor.test(outputData$spCD$sampleSize,outputData$spCD$Schoeners_D,method='spearman')
##indice I
cor.test(outputData$spHW$sampleSize,outputData$spHW$Hellinger_distances,method='spearman')
cor.test(outputData$spHD$sampleSize,outputData$spHD$Hellinger_distances,method='spearman')
cor.test(outputData$spCD$sampleSize,outputData$spCD$Hellinger_distances,method='spearman')


##distribuicao presente, inter e maximo glacial
library(raster)

##threshold para modelo
threHW5 = read.csv(paste(projectFolder,'Maxent/spHW/StatisticsResults-spHW.csv',sep=''),header=TRUE); threHW5 = threHW5[threHW5$sampleSize==5,]; threHW5 = mean(threHW5$Threshold)
threHW45 = read.csv(paste(projectFolder,'Maxent/spHW/StatisticsResults-spHW.csv',sep=''),header=TRUE); threHW45 = threHW45[threHW45$sampleSize==45,]; threHW45 = mean(threHW45$Threshold)
threHW95 = read.csv(paste(projectFolder,'Maxent/spHW/StatisticsResults-spHW.csv',sep=''),header=TRUE); threHW95 = threHW95[threHW95$sampleSize==95,]; threHW95 = mean(threHW95$Threshold)
##threshold paraa distribuicao real
HWcurrentReal = raster(paste(projectFolder,'NichoReal/spHW/000.asc',sep='')) > 0.1

HW22Real = raster(paste(projectFolder,'NichoReal/spHW/022.asc',sep='')) > 0.1
HWModel_0kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threHW5
HWModel_0kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threHW45
HWModel_0kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threHW95
HWModel_22kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threHW5
HWModel_22kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threHW45
HWModel_22kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threHW95

##threshold para modelo
threHD5 = read.csv(paste(projectFolder,'Maxent/spHD/StatisticsResults-spHD.csv',sep=''),header=TRUE); threHD5 = threHD5[threHD5$sampleSize==5,]; threHD5 = mean(threHD5$Threshold)
threHD45 = read.csv(paste(projectFolder,'Maxent/spHD/StatisticsResults-spHD.csv',sep=''),header=TRUE); threHD45 = threHD45[threHD45$sampleSize==45,]; threHD45 = mean(threHD45$Threshold)
threHD95 = read.csv(paste(projectFolder,'Maxent/spHD/StatisticsResults-spHD.csv',sep=''),header=TRUE); threHD95 = threHD95[threHD95$sampleSize==95,]; threHD95 = mean(threHD95$Threshold)
##threshold paraa distribuicao real
HDcurrentReal = raster(paste(projectFolder,'NichoReal/spHD/000.asc',sep='')) > 0.1

HD22Real = raster(paste(projectFolder,'NichoReal/spHD/022.asc',sep='')) > 0.1
HDModel_0kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threHD5
HDModel_0kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threHD45
HDModel_0kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threHD95
HDModel_22kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threHD5
HDModel_22kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threHD45
HDModel_22kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threHD95

##threshold para modelo
threCD5 = read.csv(paste(projectFolder,'Maxent/spCD/StatisticsResults-spCD.csv',sep=''),header=TRUE); threCD5 = threCD5[threCD5$sampleSize==5,]; threCD5 = mean(threCD5$Threshold)
threCD45 = read.csv(paste(projectFolder,'Maxent/spCD/StatisticsResults-spCD.csv',sep=''),header=TRUE); threCD45 = threCD45[threCD45$sampleSize==45,]; threCD45 = mean(threCD45$Threshold)
threCD95 = read.csv(paste(projectFolder,'Maxent/spCD/StatisticsResults-spCD.csv',sep=''),header=TRUE); threCD95 = threCD95[threCD95$sampleSize==95,]; threCD95 = mean(threCD95$Threshold)
##threshold paraa distribuicao real
CDcurrentReal = raster(paste(projectFolder,'NichoReal/spCD/000.asc',sep='')) > 0.1

CD22Real = raster(paste(projectFolder,'NichoReal/spCD/022.asc',sep='')) > 0.1
CDModel_0kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threCD5
CDModel_0kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threCD45
CDModel_0kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threCD95
CDModel_22kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threCD5
CDModel_22kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threCD45
CDModel_22kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threCD95

library(maptools)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp")

##sobreposicoes spHW
jpeg(filename=paste(projectFolder,'Maxent/graficos/sobreposicoesHW.jpg',sep=''), width = 1100 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,0))
plot(HWcurrentReal*1+HWModel_0kyrSample5*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample45*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample95*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample5*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample45*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample95*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
grid()
mtext('Hot & Wet sp.',outer=TRUE,cex=4)
dev.off()

##sobreposicoes spHD
jpeg(filename=paste(projectFolder,'Maxent/graficos/sobreposicoesHD.jpg',sep=''), width = 1100 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,0))
plot(HDcurrentReal*1+HDModel_0kyrSample5*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(HDcurrentReal*1+HDModel_0kyrSample45*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(HDcurrentReal*1+HDModel_0kyrSample95*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
grid()
plot(HD22Real*1+HDModel_22kyrSample5*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(HD22Real*1+HDModel_22kyrSample45*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(HD22Real*1+HDModel_22kyrSample95*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
grid()
mtext('Hot & Dry sp.',outer=TRUE,cex=4)
dev.off()

##sobreposicoes spCD
jpeg(filename=paste(projectFolder,'Maxent/graficos/sobreposicoesCD.jpg',sep=''), width = 1100 , height = 1100) 
par(mfrow=c(2,3)) 
plot(CDcurrentReal*1+CDModel_0kyrSample5*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample45*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample95*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample5*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample45*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample95*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
legend("bottomright", inset=c(-1,0),legend=c('Virtual species','Maxent projection','Overlap'),pch=20,col=c('blue','green','dark green'),cex=2)
grid()
mtext('Cold & Dry sp.',outer=TRUE,cex=4)
dev.off()
