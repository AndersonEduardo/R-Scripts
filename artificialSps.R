#########################################################################################
####SCRIPT PARA FAZER A DISTRIBUICAO REAL DAS ESPECIES EM VARIOS MOMENTOS DO TEMPO####
#########################################################################################
library(virtualspecies)
library(maptools)
library(dismo)
library(raster)
library(phyloclim) #para funcao niche.overlap()
source("/home/anderson/R/R-Scripts/TSSmaxent.R")

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
        writeRaster(auxVector[[j]], filename=paste(projectFolder,'NichoReal/',names(auxVector[[j]]),'/',nameScenario,'.asc',sep=""), overwrite=T,prj=T) #salvando o raster do mapa da sp
    }
}

###SEGUNDA PARTE: amostragem de pontos de ocorrencia em diferentes camadas de tempo###

envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto/" #pasta com as variaveis ambientais
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha com os pontos amostrados sera salva
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
biomodFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/biomod/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
sampleSizes = c(5,15,25,35,45,55,65,75,85,95)
Tmax = 22
amostra = data.frame()

#ocorrencias
for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
    for (sSize in sampleSizes){ #numero de pontos (registros, dados) na amostra
        sampledAges = vector()
        sampledAges = round(runif(sSize,0,Tmax)) #selecionando 'n' camadas de tempo aleatoriamente
        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da
        for (j in 1:5){ #replicas do cenario amostral
            for (sAge in sampledAges){ #amostrando em cada camada de tempo que consta na amostra
                amostra_i = sampleOccurrences(x=raster(nicheRealPath[sAge+1]),n=1,plot=FALSE)$sample.points[,1:2] #amostra d ponto
                scenarioName = basename(nicheRealPath[1:24][sAge+1]) #tempo vinculado ao cenario para variaveis ambientais
                scenarioName = gsub('.asc','',scenarioName) #retirando do nome o '.asc'
                layers_i = extract(
                    x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                    y=amostra_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
                amostra = rbind(amostra, cbind(amostra_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
            }
            names(amostra) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
            write.csv(amostra,paste(projectFolder,'Amostras/',spsTypes[i],'/occ',sSize,'pts', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
            amostra = data.frame() #devolvendo data.frame vazio para proxima rodada
        }
    }
}

#background
for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
    for (sSize in sampleSizes){ #numero de pontos (registros, dados) na amostra
        sampledAges = vector()
        sampledAges = round(runif(sSize,0,Tmax)) #selecionando 'n' camadas de tempo aleatoriamente
#        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        for (j in 1:5){ #replicas do cenario amostral
            for (sAge in sampledAges){ #amostrando em cada camada de tempo que consta na amostra

                envVarPath = list.files(path=paste(envVarFolder,list.files(path=paste(envVarFolder))[sAge+1],sep=''), full.names=TRUE, pattern='.asc') #lista com os enderecos das variaveis ambientais no tempo corresposndente a interacao 
                
                amostra_i = randomPoints(mask=raster(envVarPath[1]),n=1) #amostra d ponto
                scenarioName = list.files(path=paste(envVarFolder))[sAge+1] #nome do cenario
                layers_i = extract(
                    x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                    y=amostra_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
                amostra = rbind(amostra, cbind(amostra_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
            }
            names(amostra) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
            write.csv(amostra,paste(projectFolder,'Amostras/',spsTypes[i],'/bg',sSize,'pts', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
            amostra = data.frame() #devolvendo data.frame vazio para proxima rodada
        }
    }
}

## ###SEGUNDA PARTE: amostragem de pontos de ocorrencia em diferentes camadas de tempo

## caminhosCamadasTemp = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
## projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
## mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha com os pontos amostrados sera salva
## AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
## biomodFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/biomod/' #pasta para resultados do maxent
## spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies

## ############################
## for(g in 1:10){ #loop para replicar a 'amostragem de pontos fosseis'
##     for (h in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
##         nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[h],sep='') #pasta com os mapas de nicho real da sp
##         nicheRealPath = list.files(path=nicheRealFolder, full.names=T, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da sp
##         for (j in 1:length(Npres)){ #loop sobre os cenarios de amostragem para o presente
##             for (k in 1:length(Npass)){ #loop sobre os cenarios de amostragem para o passado
##                 amostra = data.frame()
##                 for (i in 1:length(nicheRealPath[1:24])){ #loop sobre o numero de camadas de tempo a serem amostradas
##                     realNicheMap = raster(nicheRealPath[i]) #abrindo o mapa de occ da sp
##                     ##nameScenario = names(realNicheMap)
##                     if (names(realNicheMap) != 'X000'){ #verificando se e o mapa da ocorrencia no presente
##                         amostra_i = sampleOccurrences(realNicheMap,Npass[k],plot=FALSE) #realizando a amostragem
##                         layers_i = stack(list.files(path=paste(envVarFolder,'/000',sep=''), pattern='asc', full.names=T)) ### stack all rasters in
##                         envVar_i = extract(layers_i, amostra_i[[1]][1:2], method='bilinear', buffer=NULL, fun=NULL, df=TRUE)
##                         amostra_i = cbind(amostra_i[[1]][1:2],envVar_i[,2:ncol(envVar_i)])
##                     }
##                     else{
##                         amostra_i = sampleOccurrences(realNicheMap,Npres[j],plot=FALSE) #realizando a amostragem
##                         scenarioName = basename(nicheRealPath[1:24][i])
##                         scenarioName = gsub('.asc','',scenarioName)
##                         layers_i = stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='asc', full.names=T)) ### stack all rasters in
##                         envVar_i = extract(layers_i, amostra_i[[1]][1:2], method='bilinear', buffer=NULL, fun=NULL, df=TRUE)
##                         amostra_i = cbind(amostra_i[[1]][1:2],envVar_i[,2:ncol(envVar_i)])
##                     }
##                     amostra = rbind(amostra,amostra_i)
##                 }

##                 setwd(paste(projectFolder,'Amostras/',spsTypes[h],'/',sep='')) #mudando a pasta de trabalho para os outputs
##                 names(amostra) = c('lon','lat',"bioclim_01","bioclim_04","bioclim_10","bioclim_11","bioclim_12","bioclim_15","bioclim_16","bioclim_17")
##                 write.csv(amostra,file=paste('occurrences_',g,'.csv',sep=''),row.names=FALSE)#salvando a planilha com os dados da amostra
##             }
##         }
##     }
## }

## ##Background points
## for (g in 1:10){ #loop para replicar a 'amostragem de pontos fosseis'
##     for (i in 1:length(spsTypes)){  #loop sobre os 'tipos de especies'
##         backgroundPoints = data.frame()
##         nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
##         nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da sp
##         for (j in 1:length(nicheRealPath[1:24])){ #loop sobre o numero de camadas de tempo a serem amostradas
##             predictors = stack(list.files(path=paste(caminhosCamadasTemp[j],sep=''),pattern='asc',full.names=TRUE)) #carregando as variaveis ambientais
##             #predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais #ARRUMAR: RETIRAR, P FICAR MAIS RAPIDO
            
##             #pooledOccPoints = read.csv(paste(mainSampleFolder,spsTypes[i],'/occurrences_',g,'.csv',sep=''),header=TRUE) #abrindo .csv de occ #ARRUMAR: RETIRAR ESTA LINHA (E MELHOR GERAR BACKGROUND CONSIDERANDO TUDO, E DEPOIS ELIMINAR)
##             #backgroundPoints_i<- randomPoints(mask=predictors[[1]],n=200,p=pooledOccPoints[,c("lon","lat")], excludep=TRUE) #sorteando coordenadas #ARRUMAR: RETIRAR O ARGUMENTO 'P' (VIDE COMENTARIO NA LINHA ANTERIOR)
##             backgroundPoints_i <- randomPoints(mask=predictors[[1]],n=100, excludep=TRUE) # 100 pontos de fundo por camada de tempo
##             colnames(backgroundPoints_i) <- c("lon", "lat")
            
##             ##extraindo dados da variavel climatica nos pontos de background
##             ausencesVars <- extract(predictors,backgroundPoints_i,method='bilinear',buffer=NULL,fun=NULL) #extraindo variaveis ambientais nas coordenadas para cada time slice 'j'
##             backgroundPoints = data.frame(rbind(backgroundPoints,data.frame(backgroundPoints_i,ausencesVars))) #dados completos dos background points
##         }
##         ##'limpando' os background points
##         backgroundPoints1 = round(backgroundPoints, digits=3) ##ARRUMAR: RETIRAR COLUNAS DAS COORDENADAS PARA 'LIMPAR'
## ##        backgroundPoints1 = round(backgroundPoints[,c("bioclim_01","bioclim_04","bioclim_10","bioclim_11","bioclim_12","bioclim_15","bioclim_16","bioclim_17")], digits=4)
##         backgroundPoints2 <- backgroundPoints1[!duplicated(backgroundPoints1[,1:2]),]
##         backgroundPoints3 <- backgroundPoints2[complete.cases(backgroundPoints2),]
##         backgroundPoints <- backgroundPoints3
##         setwd(paste(projectFolder,'Amostras/',spsTypes[i],'/',sep='')) 
##         write.csv(backgroundPoints,file=paste('background_',g,'.csv',sep=''),row.names=FALSE)
##     }
##}                                       

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

for (i in 1:length(spsTypes)){
    for (j in sampleSizes){
        for (k in 1:5){ #loop sobre o numero de replicas #10
        
            occPoints = read.csv(paste(mainSampleFolder,spsTypes[i],'/occ',j,'pts',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia
            backgroundPoints = read.csv(paste(mainSampleFolder,spsTypes[i],'/bg',j,'pts',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de background
            names(backgroundPoints) = names(occPoints) #certificando que os nomes das colunas estão iguais (cuidado aqui...)
            dataSet = data.frame(cbind(rbind(occPoints,backgroundPoints),pres=c(rep(1,nrow(occPoints)),rep(0,nrow(backgroundPoints))))) #planilha de dados no formato SWD
            
            me = maxent(
                ##x=dataSet[,c("bioclim_01","bioclim_04","bioclim_10","bioclim_11","bioclim_12","bioclim_15","bioclim_16","bioclim_17")],
                x=dataSet[,c("bioclim_01","bioclim_12")],
                p=dataSet$pres,
                path=paste(maxentFolder,spsTypes[i],sep=''),
                args=c('responsecurves=TRUE',
                       'jackknife=TRUE',
                       'randomseed=TRUE',
                       'randomtestpoints=25',
                       'maximumbackground=5000',
                       'replicates=1',
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
            TSSvector = rbind(TSSvector, TSSmaxent(paste(maxentFolder,spsTypes[i],'/',sep='')))
            evaluation = append(evaluation, evaluate(p=occPoints,a=backgroundPoints,model=me))
        
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
        thresholdValues = NULL
        aucValues = NULL
        for (m in 1:length(evaluation)){
            thresholdValues <- append(thresholdValues, threshold(evaluation[[m]],'spec_sens'))
            aucValues = append(aucValues, evaluation[[m]]@auc)
        }
        aucMean = mean(aucValues)
        thresholdMean = mean(thresholdValues)
        TSSmean = mean(TSSvector$TSS)
        statResults = data.frame(sp=spsTypes[i],AUCmean=aucMean,TSSmean=TSSmean,ThresholdMean=thresholdMean)
        write.csv(statResults,file=paste(projectFolder,'Maxent/',spsTypes[i],'/StatisticsResults-',spsTypes[i],'.csv',sep=''),row.names=FALSE)
        
        ##esvaziando o vetor para a proxima especie
        evaluation = list()
        TSSvector = data.frame()
    }
}
#######################################################
#######################################################
#######################################################

#######################################################
####################### GLM ###########################
#######################################################
    
options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
GLMfolder = '/home/anderson/Documentos/Projetos/Sps artificiais/GLM/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
                                        #
model1 = pres ~ bioclim_01 + bioclim_04 + bioclim_10 + bioclim_11 + bioclim_12 + bioclim_15 + bioclim_16 + bioclim_17 
model2 = pres ~ bioclim_01 + I(bioclim_01^2) + bioclim_04 + I(bioclim_04^2) + bioclim_10 + I(bioclim_10^2) + bioclim_11 + I(bioclim_11^2) + bioclim_12 + I(bioclim_12^2) + bioclim_15 + I(bioclim_15^2) + bioclim_16 + I(bioclim_16^2) + bioclim_17 + I(bioclim_17^2)
model3 = pres ~ bioclim_01 + bioclim_12
model4 = pres ~ bioclim_01 + I(bioclim_01^2) + bioclim_12 + I(bioclim_12^2)
model = c(model1,model2,model3,model4)
scenarioModel = c('8varLinearModel','8varQuadModel','2varLinearModel','2varQuadModel')
                                        #
varNames = list(c('bioclim_01','bioclim_04','bioclim_10','bioclim_11','bioclim_12','bioclim_15','bioclim_16','bioclim_17'),c('bioclim_01','bioclim_04','bioclim_10','bioclim_11','bioclim_12','bioclim_15','bioclim_16','bioclim_17'),c('bioclim_01','bioclim_12'),c('bioclim_01','bioclim_12'))


for (i in 1:length(spsTypes)){
    modelInfo = data.frame()
    for (j in 1:length(model)){
        
        ##preparando o conjunto de dados
        predictors = stack(list.files(path=envVarPaths[1],full.names=T, pattern='.asc')) #predictors com todas as variaveis (presente)
        predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
        sampleFolder = paste(mainSampleFolder,spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        samplePaths = list.files(path=sampleFolder, full.names=T, pattern='.csv') #lista com os enderecos dos mapas de distribuicao da sp
        sp.data = read.csv(paste(samplePaths[2],sep=''),header=TRUE) #abrindo a planilha de pontos de occ amostrados
        sp.data = data.frame(cbind(sp.data,pres=1))
        names(sp.data) = c('lon','lat',"bioclim_01","bioclim_04","bioclim_10","bioclim_11","bioclim_12","bioclim_15","bioclim_16","bioclim_17",'pres')
        pseudoausencia1.occ <- randomPoints(mask=predictors[[1]], n=1000, p=sp.data[1:2], excludep=TRUE) 
        pseudoausencia2.occ <- round(pseudoausencia1.occ[,1:2], digits=4)
        pseudoausencia3.occ <- pseudoausencia2.occ[!duplicated(pseudoausencia2.occ),]
        pseudoausencia4.occ <- pseudoausencia3.occ[complete.cases(pseudoausencia3.occ),]
        pseudoausencia.occ <- data.frame(pseudoausencia4.occ)
        colnames(pseudoausencia.occ) <- c("lon", "lat")
        pseudoausencia.clim <- extract(predictors, pseudoausencia.occ, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)
        pseudoausencia.data <- data.frame(cbind(pseudoausencia.occ,pseudoausencia.clim[2:ncol(pseudoausencia.clim)],pres=0))
        dataset = data.frame(rbind(sp.data,pseudoausencia.data))
        
        ##avaliando o modelo
        V <- numeric()#abrir un vector vazio 
        
        for (k in 1:10){
            ##reparando uma porcao dos dados de presenca e ausencia (background) para calibrar (treinar) o modelo
            rand = round(0.75*runif(nrow(sp.data)))
            presenciasTrain = sp.data[rand==0,]
            ausenciasTrain = pseudoausencia.data[rand==0,]
            
            ##juntando presencas e ausencias da calibracao
            presausTrainRaw <- rbind(presenciasTrain, ausenciasTrain)
            presausTrainRaw = data.frame(presausTrainRaw)
            presausTrainRaw = presausTrainRaw[!duplicated(presausTrainRaw[,1:2]),] #selecionar colunas de longitude e latitude
            presausTrainRaw<-presausTrainRaw[complete.cases(presausTrainRaw),]
            presausTrain = presausTrainRaw
            
            ##AJUSTANDO O MODELO##                
            GLM = glm(model[[j]], family=binomial(link=logit), data=dataset)
            
            ##TSSfunction(paste(maxentFolder,spsTypes[k],'/',names(feature)[i],'/',betamultiplier[j],sep='')) #TSS
            
            ##pegando a porcao dos dados separados para a avaliacao (validacao) do modelo
            presencias.evaluacion = sp.data[rand==1,]
            presencias.evaluacion <- cbind(presencias.evaluacion$lon,presencias.evaluacion$lat)
            pseudoausencias.evaluacion = pseudoausencia.data[rand==1,]
            pseudoausencias.evaluacion = cbind(pseudoausencias.evaluacion$lon,pseudoausencias.evaluacion$la)

            ##RODANDO A AVALIACAO DO MODELO##
            evaluacion=evaluate(presencias.evaluacion, pseudoausencias.evaluacion, GLM, predictors)
            
            ##registrando o valor de AUC em um objeto
            V[k]<-evaluacion@"auc" #sacamos el valor de auc (fíjate que es una @ en lugar de $ para mirar dentro de los slots)y guardamos en vector
        }

        auc<-mean(V, na.rm=T)#media de los vectores de las iternaciones de j
        
        modelInfo = rbind(modelInfo, data.frame(Scenario=scenarioModel[j],Deviance=GLM$deviance,Null_deviance=GLM$null.deviance,AIC=GLM$aic,AUC=auc))
        names(modelInfo) = c('Scenario','Deviance','Null_deviance','AIC','AUC')
        
        ##projecoes do modelo
        for (l in 1:length(envVarPaths[1:24])){
            predictors.proj = stack(list.files(path=envVarPaths[l],full.names=T, pattern='.asc')) #predictors com todas as variaveis (presente)
            predictors.proj = predictors.proj[[varNames[[j]]]] #carregando as variaveis ambientais
            nameScenario = basename(envVarPaths[l])
            crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
            proj <- predict(predictors.proj,GLM,type=c('response'),crs=crs) #rodando a projecao espacial do modelo
            names(proj) <- nameScenario #ajustando nome
            writeRaster(proj, filename=paste(GLMfolder,spsTypes[i],'/',scenarioModel[j],'/Projections/',nameScenario,'.asc',sep=''), overwrite=T) #salvando o mapa de suitability projetado
        }
    }
    write.csv(modelInfo,paste(GLMfolder,spsTypes[i],'/','modelInfo.csv',sep=''),row.names=F)
}

###############################################
###############################################
###############################################

###############################################
################### BIOMOD ####################
###############################################
library(biomod2)

options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha 
biomodFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/biomod/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies

for (g in 1:10){
    for (i in 1:length(spsTypes)){
        setwd(paste(biomodFolder,spsTypes[i],'/biomodOutputs',sep='')) #mudando a pasta de trabalho para os outputs do biomod
        ##
        if (!dir.exists(paste(g))){
            dir.create(paste(g,sep=''))
        }
        ##
        setwd(paste(g))
        ##
        occPoints = read.csv(paste(mainSampleFolder,spsTypes[i],'/occurrences_',g,'.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia
        backgroundPoints = read.csv(paste(mainSampleFolder,spsTypes[i],'/background_',g,'.csv',sep=''),header=TRUE) #abrindo pontos de background
        ##names(backgroundPoints) = names(occPoints) #certificando que os nomes das colunas estão iguais (cuidado aqui...)
        dataSet = data.frame(cbind(rbind(occPoints,backgroundPoints),pres=c(rep(1,nrow(occPoints)),rep(0,nrow(backgroundPoints))))) #planilha de dados no formato SWD
        ##
        ##DADOS DE ENTRADA PARA O BIOMOD2###
        myResp = dataSet$pres
        predictors = dataSet[,c("bioclim_01","bioclim_12")]
        myRespXY = dataSet[,c('lon','lat')] 
        myRespName = spsTypes[i]
        
        myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                             resp.xy = myRespXY,
                                             expl.var = predictors,
                                             resp.name = myRespName)
        
        myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips=list(path_to_maxent.jar="/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java",maximumiterations=5000,linear=TRUE,quadratic=TRUE,product=FALSE,threshold=FALSE,hinge=FALSE,threads=2))
        
        myBiomodModelOut <- BIOMOD_Modeling(
            myBiomodData,
            models = c('GLM','MAXENT.Phillips'),
            models.options = myBiomodOption,
            NbRunEval = 3,
            DataSplit = 75,
            VarImport = 3,
            models.eval.meth = c('TSS','ROC'),
            SaveObj = TRUE,
            rescal.all.models = TRUE,
            do.full.models = FALSE,
            modeling.id = paste(myRespName,"_Models",sep=""))
        
        ##PROJECOES##
        for (j in 1:length(envVarPaths[1:24])){
            setwd(paste(biomodFolder,spsTypes[i],'/biomodOutputs/',g,sep='')) #mudando a pasta de trabalho para os outputs do biomod
            predictors = stack(list.files(path=envVarPaths[j],full.names=T, pattern='.asc')) #predictors com todas as variaveis (presente)
            predictors = stack(c(predictors$bioclim_01,predictors$bioclim_12)) #recortando as variaveis ambientais
            crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
            
            myBiomodProj <- BIOMOD_Projection(
                modeling.output = myBiomodModelOut,
                new.env = stack(predictors),
                proj.name = paste(j-1,'kyrBP',sep=''),
                selected.models = paste(myBiomodModelOut@models.computed,sep=''),
                binary.meth = 'TSS',
                compress = FALSE,
                clamping.mask = TRUE,
                output.format = '.grd',
                on_0_1000 = FALSE)
            
            ##My output data
            projStack = get_predictions(myBiomodProj)
            varImportance = get_variables_importance(myBiomodModelOut)
            evaluationScores = get_evaluations(myBiomodModelOut)
            ##
            setwd(paste(biomodFolder,spsTypes[i],'/projections/',sep='')) #mudando a pasta de trabalho para os outputs do biomod
            if (!dir.exists(paste(g))){
                dir.create(paste(g,sep=''))
            }
            setwd(paste(g))
            writeRaster(projStack,filename=paste(biomodFolder,spsTypes[i],'/projections/',g,'/',names(projStack),'-',j-1,'kyrBP',sep=''),bylayer=TRUE,format='ascii',overwrite=TRUE)
            ##
            setwd(paste(biomodFolder,spsTypes[i],'/varImportance/',sep='')) #mudando a pasta de trabalho para os outputs do biomod
            if (!dir.exists(paste(g))){
                dir.create(paste(g,sep=''))
            }
            setwd(paste(g))
            write.csv(data.frame(varImportance),paste(biomodFolder,spsTypes[i],'/varImportance/',g,'/varImportance_',myRespName,'-',j-1,'kyrBP.csv',sep=''),row.names=TRUE)
            ##
            setwd(paste(biomodFolder,spsTypes[i],'/evaluationScores/',sep='')) #mudando a pasta de trabalho para os outputs do biomod
            if (!dir.exists(paste(g))){
                dir.create(paste(g,sep=''))
            }
            setwd(paste(g))
            write.csv(data.frame(evaluationScores),paste(biomodFolder,spsTypes[i],'/evaluationScores/',g,'/evaluationScores_',myRespName,'-',j-1,'kyrBP.csv',sep=''),row.names=TRUE)
        }
    }
}
###############################################
###############################################
###############################################


###QUARTA PARTE: comparando projecao do SDM e a distribuicao espacial real do nicho da sp###

projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
scenarioModel = c('8varLinearModel','8varQuadModel','2varLinearModel','2varQuadModel')
sampleSizes = c(5,15,25,35,45,55,65,75,85,95) #aqui, deve ser igual ao usasado nas partes anteriores do script

model1 = pres ~ bioclim_01 + bioclim_04 + bioclim_10 + bioclim_11 + bioclim_12 + bioclim_15 + bioclim_16 + bioclim_17 
model2 = pres ~ bioclim_01 + I(bioclim_01^2) + bioclim_04 + I(bioclim_04^2) + bioclim_10 + I(bioclim_10^2) + bioclim_11 + I(bioclim_11^2) + bioclim_12 + I(bioclim_12^2) + bioclim_15 + I(bioclim_15^2) + bioclim_16 + I(bioclim_16^2) + bioclim_17 + I(bioclim_17^2)
model3 = pres ~ bioclim_01 + bioclim_12
model4 = pres ~ bioclim_01 + I(bioclim_01^2) + bioclim_12 + I(bioclim_12^2)
model = c(model1,model2,model3,model4)

for (i in 1:length(spsTypes)){
    
    nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
    nicheRealPath = list.files(path=nicheRealFolder,pattern='asc',full.names=T) #lista com os enderecos dos mapas de distribuicao da sp
    projectionsFolder = paste(projectFolder,'Maxent/',spsTypes[i],'/projections',sep='') #pasta com as projecoes do cenario
    projectionsPath = list.files(path=projectionsFolder, pattern='asc',full.names=T) #caminhos para os .asc na paste do cenario
    outputData = data.frame(kyrBP=numeric(),sampleSize=numeric(),Schoeners_D=numeric(),Hellinger_distances=numeric())
    
    for (l in 1:length(nicheRealPath[1:24])){ #loop sobre as cadamdas de tempo

        realNiche = raster(nicheRealPath[l]) #nicho real
        realNiche.spgrid = as(realNiche,'SpatialGridDataFrame')
        
        for (m in sampleSizes){ #loop sobre os tamnhos amostrais

            timeSampleData =  stack(list.files(path=projectionsFolder, pattern=glob2rx(paste('*Time',l-1,'*Sample',m,'.asc',sep='')),full.names=TRUE))
            output_i = data.frame(kyrBP=numeric(),sapleSize=numeric(),Schoeners_D=numeric(),Hellinger_distances=numeric()) #dataframe vazio para o loop abaixo
            
            for(n in 1:5){
              
                sdmNiche = timeSampleData[[n]] #mapa de suitability gerado por SDM
                sdmNiche.spgrid = as(sdmNiche,'SpatialGridDataFrame')
                nicheOverlap = niche.overlap(c(realNiche.spgrid,sdmNiche.spgrid))
                output_i= rbind(output_i,cbind(kyrBP=l-1,sampleSize=m,Schoeners_D=nicheOverlap[1,2],Hellinger_distances=nicheOverlap[2,1]))
                
            }
        
            outputMeans = colMeans(output_i) #media das iteracoes (para a camada de tempo atual)
            outputData = data.frame(rbind(outputData,outputMeans))
            
        }
    }
    
    names(outputData) = c('kyrBP','sampleSize','Schoeners_D','Hellinger_distances')  
    write.csv(outputData, file=paste(projectionsFolder,"/NO.csv",sep=""),row.names=FALSE) #salvando os dados do cenario
    
}


### QUINTA PARTE: construindo graficos dos resultados ###

spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
###abrindo as planilhas de dados
outputData = list()
vetor.nomes = vector()
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto

####GLM
## for (i in 1:length(spsTypes)){
##     for (j in 1:length(model)){
##         k = j + length(model)*(i-1)
##         outputData[[k]] = read.csv(file=paste(projectFolder,'GLM/',spsTypes[i],'/',scenarioModel[j],'/Projections','/NO.csv',sep=""),header=TRUE)
##         vetor.nomes = append(vetor.nomes,paste(spsTypes[i],scenarioModel[j],sep=''))
##     }
## }
## names(outputData) = vetor.nomes

##maxent
for (i in 1:length(spsTypes)){
    outputData[[i]] = read.csv(file=paste(projectFolder,'Maxent/',spsTypes[i],'/projections/NO.csv',sep=''),header=TRUE)
    vetor.nomes = append(vetor.nomes,paste(spsTypes[i],sep=''))
}
names(outputData) = vetor.nomes

###graficos

##Boxplots 

##GLM##
##HW
HWdataD = data.frame(rbind(data.frame(indexD=outputData$spHW8varQuadModel$Schoeners_D,model='8var.&quadratic'),data.frame(indexD=outputData$spHW8varLinearModel$Schoeners_D,model='8var.&linear'),data.frame(indexD=outputData$spHW2varQuadModel$Schoeners_D,model='2var.&quadratic'),data.frame(indexD=outputData$spHW2varLinearModel$Schoeners_D,model='2var.&linear')))
HWdataH = data.frame(rbind(data.frame(indexH=outputData$spHW8varQuadModel$Hellinger_distances,model='8var.&quadratic'),data.frame(indexH=outputData$spHW8varLinearModel$Hellinger_distances,model='8var.&linear'),data.frame(indexH=outputData$spHW2varQuadModel$Hellinger_distances,model='2var.&quadratic'),data.frame(indexH=outputData$spHW2varLinearModel$Hellinger_distances,model='2var.&linear')))

##HD
HDdataD = data.frame(rbind(data.frame(indexD=outputData$spHD8varQuadModel$Schoeners_D,model='8var.&quadratic'),data.frame(indexD=outputData$spHD8varLinearModel$Schoeners_D,model='8var.&linear'),data.frame(indexD=outputData$spHD2varQuadModel$Schoeners_D,model='2var.&quadratic'),data.frame(indexD=outputData$spHD2varLinearModel$Schoeners_D,model='2var.&linear')))
HDdataH = data.frame(rbind(data.frame(indexH=outputData$spHD8varQuadModel$Hellinger_distances,model='8var.&quadratic'),data.frame(indexH=outputData$spHD8varLinearModel$Hellinger_distances,model='8var.&linear'),data.frame(indexH=outputData$spHD2varQuadModel$Hellinger_distances,model='2var.&quadratic'),data.frame(indexH=outputData$spHD2varLinearModel$Hellinger_distances,model='2var.&linear')))

##CD
CDdataD = data.frame(rbind(data.frame(indexD=outputData$spCD8varQuadModel$Schoeners_D,model='8var.&quadratic'),data.frame(indexD=outputData$spCD8varLinearModel$Schoeners_D,model='8var.&linear'),data.frame(indexD=outputData$spCD2varQuadModel$Schoeners_D,model='2var.&quadratic'),data.frame(indexD=outputData$spCD2varLinearModel$Schoeners_D,model='2var.&linear')))
CDdataH = data.frame(rbind(data.frame(indexH=outputData$spCD8varQuadModel$Hellinger_distances,model='8var.&quadratic'),data.frame(indexH=outputData$spCD8varLinearModel$Hellinger_distances,model='8var.&linear'),data.frame(indexH=outputData$spCD2varQuadModel$Hellinger_distances,model='2var.&quadratic'),data.frame(indexH=outputData$spCD2varLinearModel$Hellinger_distances,model='2var.&linear')))


##maxent##
##HW
HWdataD = data.frame(kyrBP=outputData$spHW$kyrBP,indexD=outputData$spHW$Schoeners_D)
HWdataH = data.frame(kyrBP=outputData$spHW$kyrBP,indexH=outputData$spHW$Hellinger_distances)

##HD
HDdataD = data.frame(kyrBP=outputData$spHW$kyrBP,indexD=outputData$spHD$Schoeners_D)
HDdataH = data.frame(kyrBP=outputData$spHW$kyrBP,indexH=outputData$spHD$Hellinger_distances)

##CD
CDdataD = data.frame(kyrBP=outputData$spHW$kyrBP,indexD=outputData$spCD$Schoeners_D)
CDdataH = data.frame(kyrBP=outputData$spHW$kyrBP,indexH=outputData$spCD$Hellinger_distances)

#glm
pdf(file='/home/anderson/Documentos/Projetos/Sps artificiais/GLM/boxplots.pdf')
par(mfrow=c(2,3))
##
boxplot(indexD~model,data=HWdataD,ylim=c(0.3,1.0),main='"Hot & Wet" species')
stripchart(indexD~model,data=HWdataD,vertical=TRUE,method="jitter",pch=20,cex=1,col=rgb(0.5,0.5,0.5,0.2),add=TRUE) 
##
boxplot(indexD~model,data=HDdataD,ylim=c(0.3,1.0),main='"Hot & Dry" species')
stripchart(indexD~model,data=HDdataD,vertical=TRUE,method="jitter",pch=20,cex=1,col=rgb(0.5,0.5,0.5,0.2),add=TRUE) 
##
boxplot(indexD~model,data=CDdataD,ylim=c(0.3,1.0),main='"Cold & Dry" species')
stripchart(indexD~model,data=CDdataD,vertical=TRUE,method="jitter",pch=20,cex=1,col=rgb(0.5,0.5,0.5,0.2),add=TRUE) 
##
boxplot(indexH~model,data=HWdataH,ylim=c(0.3,1.0),main='"Hot & Wet" species')
stripchart(indexH~model,data=HWdataH,vertical=TRUE,method="jitter",pch=20,cex=1,col=rgb(0.5,0.5,0.5,0.2),add=TRUE) 
##
boxplot(indexH~model,data=HDdataH,ylim=c(0.3,1.0),main='"Hot & Dry" species')
stripchart(indexH~model,data=HDdataH,vertical=TRUE,method="jitter",pch=20,cex=1,col=rgb(0.5,0.5,0.5,0.2),add=TRUE) 
##
boxplot(indexH~model,data=CDdataH,ylim=c(0.3,1.0),main='"Cold & Dry" species')
stripchart(indexH~model,data=CDdataH,vertical=TRUE,method="jitter",pch=20,cex=1,col=rgb(0.5,0.5,0.5,0.2),add=TRUE) 
dev.off()

#maxent
jpeg(file='/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/graficos/boxplots.jpg',width=900,height=600)
par(mfrow=c(1,2))
#D
Ddata = rbind(data.frame(kyrBP=HWdataD$kyrBP,indexD=HWdataD$indexD,spsType='Hot and Wet'),data.frame(kyrBP=HDdataD$kyrBP,indexD=HDdataD$indexD,spsType='Hot and Dry'),data.frame(kyrBP=CDdataD$kyrBP,indexD=CDdataD$indexD,spsType='Cold and Dry'))
boxplot(indexD~spsType,data=Ddata,ylim=c(0.2,1),names=c('H&W','H&D','C&D'),main="Schoeners' D",cex.main=2,cex.lab=1.5,cex.axis=1.5,lwd=3)
stripchart(indexD~spsType,data=Ddata,vertical=TRUE,method="jitter",pch=20,cex=1.5,col=rgb(0.5,0.5,0.5,0.2),add=TRUE) 
#H
Hdata = rbind(data.frame(kyrBP=HWdataH$kyrBP,indexH=HWdataH$indexH,spsType='Hot and Wet'),data.frame(kyrBP=HDdataH$kyrBP,indexH=HDdataH$indexH,spsType='Hot and Dry'),data.frame(kyrBP=CDdataH$kyrBP,indexH=CDdataH$indexH,spsType='Cold and Dry'))
boxplot(indexH~spsType,data=Hdata,ylim=c(0.2,1),names=c('H&W','H&D','C&D'),main='Hellinger',cex.main=2,cex.lab=1.5,cex.axis=1.5,lwd=3)
stripchart(indexH~spsType,data=Hdata,vertical=TRUE,method="jitter",pch=20,cex=1.5,col=rgb(0.5,0.5,0.5,0.2),add=TRUE) 
dev.off()

##Schoeners' D ao longo do tempo##

#GLM
pdf(file='/home/anderson/Documentos/Projetos/Sps artificiais/GLM/DxTime.pdf')
par(mfrow=c( 1,3))
##HW species
plot(outputData$spHW8varQuadModel$Schoeners_D[1:23]~outputData$spHW8varQuadModel$kyrBP[1:23],t='b',pch=20,cex=1.5,ylim=c(0,1),col='black',ylab="Schoeners' D",xlab='kyr BP',main='Hot & Wet species')
points(outputData$spHW8varLinearModel$Schoeners_D[1:23]~outputData$spHW8varLinearModel$kyrBP[1:23],t='b',pch=18,cex=1.5,ylim=c(0.5,1),col='blue')
points(outputData$spHW2varQuadModel$Schoeners_D[1:23]~outputData$spHW2varQuadModel$kyrBP[1:23],t='b',ylim=c(0.5,1),pch=17,cex=1.5,col='green')
points(outputData$spHW2varLinearModel$Schoeners_D[1:23]~outputData$spHW2varLinearModel$kyrBP[1:23],t='b',ylim=c(0.5,1),pch=16,cex=1.5,col='red')
legend('bottomright',legend=c('8 var. - quadratic','8 var. - linear','2 var. - quadratic','2 var. - linear'),pch=c(20,18,17,16),col=c('black','blue','green','red'))
##HD species
plot(outputData$spHD8varQuadModel$Schoeners_D[1:23]~outputData$spHD8varQuadModel$kyrBP[1:23],t='b',pch=20,cex=1.5,ylim=c(0,1),col='black',ylab="Schoeners' D",xlab='kyr BP',main='Hot & Dry species')
points(outputData$spHD8varLinearModel$Schoeners_D[1:23]~outputData$spHD8varLinearModel$kyrBP[1:23],t='b',pch=18,cex=1.5,ylim=c(0.5,1),col='blue')
points(outputData$spHD2varQuadModel$Schoeners_D[1:23]~outputData$spHD2varQuadModel$kyrBP[1:23],t='b',ylim=c(0.5,1),pch=17,cex=1.5,col='green')
points(outputData$spHD2varLinearModel$Schoeners_D[1:23]~outputData$spHD2varLinearModel$kyrBP[1:23],t='b',ylim=c(0.5,1),pch=16,cex=1.5,col='red')
legend('bottomright',legend=c('8 var. - quadratic','8 var. - linear','2 var. - quadratic','2 var. - linear'),pch=c(20,18,17,16), col=c('black','blue','green','red'))#
##CD species
plot(outputData$spCD8varQuadModel$Schoeners_D[1:23]~outputData$spCD8varQuadModel$kyrBP[1:23],t='b',pch=20,cex=1.5,ylim=c(0,1),col='black',ylab="Schoeners' D",xlab='kyr BP',main='Cold & Dry species')
points(outputData$spCD8varLinearModel$Schoeners_D[1:23]~outputData$spCD8varLinearModel$kyrBP[1:23],t='b',pch=18,cex=1.5,ylim=c(0.5,1),col='blue')
points(outputData$spCD2varQuadModel$Schoeners_D[1:23]~outputData$spCD2varQuadModel$kyrBP[1:23],t='b',ylim=c(0.5,1),pch=17,cex=1.5,col='green')
points(outputData$spCD2varLinearModel$Schoeners_D[1:23]~outputData$spCD2varLinearModel$kyrBP[1:23],t='b',ylim=c(0.5,1),pch=16,cex=1.5,col='red')
legend('bottomright',legend=c('8 var. - quadratic','8 var. - linear','2 var. - quadratic','2 var. - linear'),pch=c(20,18,17,16),col=c('black','blue','green','red'))#
dev.off()

##Distancia de Hellinger ao longo do tempo
pdf(file='/home/anderson/Documentos/Projetos/Sps artificiais/GLM/HellingerXtime.pdf')
par(mfrow=c( 1,3))
##HW species
plot(outputData$spHW8varQuadModel$Hellinger_distances[1:23]~outputData$spHW8varQuadModel$kyrBP[1:23],t='b',pch=20,cex=1.5,ylim=c(0,1),col='black',ylab="Hellinger distance",xlab='kyr BP',main='Hot & Wet species')
points(outputData$spHW8varLinearModel$Hellinger_distances[1:23]~outputData$spHW8varLinearModel$kyrBP[1:23],t='b',pch=18,cex=1.5,ylim=c(0.5,1),col='blue')
points(outputData$spHW2varQuadModel$Hellinger_distances[1:23]~outputData$spHW2varQuadModel$kyrBP[1:23],t='b',ylim=c(0.5,1),pch=17,cex=1.5,col='green')
points(outputData$spHW2varLinearModel$Hellinger_distances[1:23]~outputData$spHW2varLinearModel$kyrBP[1:23],t='b',ylim=c(0.5,1),pch=16,cex=1.5,col='red')
legend('bottomright',legend=c('8 var. - quadratic','8 var. - linear','2 var. - quadratic','2 var. - linear'),pch=c(20,18,17,16),col=c('black','blue','green','red'))
##HD species
plot(outputData$spHD8varQuadModel$Hellinger_distances[1:23]~outputData$spHD8varQuadModel$kyrBP[1:23],t='b',pch=20,cex=1.5,ylim=c(0,1),col='black',ylab="Hellinger distance",xlab='kyr BP',main='Hot & Dry species')
points(outputData$spHD8varLinearModel$Hellinger_distances[1:23]~outputData$spHD8varLinearModel$kyrBP[1:23],t='b',pch=18,cex=1.5,ylim=c(0.5,1),col='blue')
points(outputData$spHD2varQuadModel$Hellinger_distances[1:23]~outputData$spHD2varQuadModel$kyrBP[1:23],t='b',ylim=c(0.5,1),pch=17,cex=1.5,col='green')
points(outputData$spHD2varLinearModel$Hellinger_distances[1:23]~outputData$spHD2varLinearModel$kyrBP[1:23],t='b',ylim=c(0.5,1),pch=16,cex=1.5,col='red')
legend('bottomright',legend=c('8 var. - quadratic','8 var. - linear','2 var. - quadratic','2 var. - linear'),pch=c(20,18,17,16), col=c('black','blue','green','red'))#
##CD species
plot(outputData$spCD8varQuadModel$Hellinger_distances[1:23]~outputData$spCD8varQuadModel$kyrBP[1:23],t='b',pch=20,cex=1.5,ylim=c(0,1),col='black',ylab="Hellinger distance",xlab='kyr BP',main='Cold & Dry species')
points(outputData$spCD8varLinearModel$Hellinger_distances[1:23]~outputData$spCD8varLinearModel$kyrBP[1:23],t='b',pch=18,cex=1.5,ylim=c(0.5,1),col='blue')
points(outputData$spCD2varQuadModel$Hellinger_distances[1:23]~outputData$spCD2varQuadModel$kyrBP[1:23],t='b',ylim=c(0.5,1),pch=17,cex=1.5,col='green')
points(outputData$spCD2varLinearModel$Hellinger_distances[1:23]~outputData$spCD2varLinearModel$kyrBP[1:23],t='b',ylim=c(0.5,1),pch=16,cex=1.5,col='red')
legend('bottomright',legend=c('8 var. - quadratic','8 var. - linear','2 var. - quadratic','2 var. - linear'),pch=c(20,18,17,16),col=c('black','blue','green','red'))#
dev.off()

#maxent
jpeg(file='/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/graficos/metricsXtime.jpg',width=1100,height=750)
par(mfrow=c(1,2))
#D
plot(indexD~kyrBP,data=HWdataD,type='b',ylim=c(0,1),pch=1,col='black',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(indexD~kyrBP,data=HDdataD,type='b',ylim=c(0,1),pch=2,col='blue',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(indexD~kyrBP,data=CDdataD,type='b',ylim=c(0,1),pch=3,col='red',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
legend(x=19,y=1,legend=c('H&W','H&D','C&D'),pch=20,col=c('black','blue','red'),cex=1.5)
#H
plot(indexH~kyrBP,data=HWdataH,type='b',ylim=c(0,1),pch=1,col='black',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(indexH~kyrBP,data=HDdataH,type='b',ylim=c(0,1),pch=2,col='blue',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(indexH~kyrBP,data=CDdataH,type='b',ylim=c(0,1),pch=3,col='red',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
legend(x=19,y=1,legend=c('H&W','H&D','C&D'),pch=20,col=c('black','blue','red'),cex=1.5)
dev.off()

##metricas X tamanho amostral

#maxent
jpeg(file='/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/graficos/metricsXtime.jpg',width=1100,height=750)
par(mfrow=c(1,2))
#D
plot(outputData$spHW$Schoeners_D~outputData$spHW$sampleSize,type='p',ylim=c(0,1),pch=1,col='black',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(outputData$spHD$Schoeners_D~outputData$spHD$sampleSize,data=HDdataD,type='p',ylim=c(0,1),pch=2,col='blue',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(outputData$spCD$Schoeners_D~outputData$spCD$sampleSize,data=CDdataD,type='p',ylim=c(0,1),pch=3,col='red',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
legend(x=6,y=1,legend=c('H&W','H&D','C&D'),pch=20,col=c('black','blue','red'),cex=1.5)
#H
plot(outputData$spHW$Hellinger_distances~outputData$spHW$sampleSize,type='p',ylim=c(0,1),pch=1,col='black',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(outputData$spHD$Hellinger_distances~outputData$spHD$sampleSize,data=HDdataD,type='p',ylim=c(0,1),pch=2,col='blue',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(outputData$spCD$Hellinger_distances~outputData$spCD$sampleSize,data=CDdataD,type='p',ylim=c(0,1),pch=3,col='red',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
legend(x=6,y=1,legend=c('H&W','H&D','C&D'),pch=20,col=c('black','blue','red'),cex=1.5)
dev.off()


##distribuicao presente, inter e maximo glacial

##GLM
HWpresentReal = "/home/anderson/Documentos/Projetos/Sps artificiais/NichoReal/spHW/000.asc"
HWpresentSDM8VQ = "/home/anderson/Documentos/Projetos/Sps artificiais/GLM/spHW/8varQuadModel/Projections/000.asc"
HWpresentSDM8VL = "/home/anderson/Documentos/Projetos/Sps artificiais/GLM/spHW/8varLinearModel/Projections/000.asc"
HWpresentSDM2VQ = "/home/anderson/Documentos/Projetos/Sps artificiais/GLM/spHW/2varQuadModel/Projections/000.asc"
HWpresentSDM2VL = "/home/anderson/Documentos/Projetos/Sps artificiais/GLM/spHW/2varLinearModel/Projections/000.asc"
#
HDpresentReal = "/home/anderson/Documentos/Projetos/Sps artificiais/NichoReal/spHD/000.asc"
HDpresentSDM8VQ = "/home/anderson/Documentos/Projetos/Sps artificiais/GLM/spHD/8varQuadModel/Projections/000.asc"
HDpresentSDM8VL = "/home/anderson/Documentos/Projetos/Sps artificiais/GLM/spHD/8varLinearModel/Projections/000.asc"
HDpresentSDM2VQ = "/home/anderson/Documentos/Projetos/Sps artificiais/GLM/spHD/2varQuadModel/Projections/000.asc"
HDpresentSDM2VL = "/home/anderson/Documentos/Projetos/Sps artificiais/GLM/spHD/2varLinearModel/Projections/000.asc"
#
CDpresentReal = "/home/anderson/Documentos/Projetos/Sps artificiais/NichoReal/spCD/000.asc"
CDpresentSDM8VQ = "/home/anderson/Documentos/Projetos/Sps artificiais/GLM/spCD/8varQuadModel/Projections/000.asc"
CDpresentSDM8VL = "/home/anderson/Documentos/Projetos/Sps artificiais/GLM/spCD/8varLinearModel/Projections/000.asc"
CDpresentSDM2VQ = "/home/anderson/Documentos/Projetos/Sps artificiais/GLM/spCD/2varQuadModel/Projections/000.asc"
CDpresentSDM2VL = "/home/anderson/Documentos/Projetos/Sps artificiais/GLM/spCD/2varLinearModel/Projections/000.asc"

speciesLayers = stack(c(HWpresentReal,HDpresentReal,CDpresentReal,
                        HWpresentSDM8VQ,HDpresentSDM8VQ,CDpresentSDM8VQ,
                        HWpresentSDM8VL,HDpresentSDM8VL,CDpresentSDM8VL,
                        HWpresentSDM2VQ,HDpresentSDM2VQ,CDpresentSDM2VQ,
                        HWpresentSDM2VL,HDpresentSDM2VL,CDpresentSDM2VL))
speciesLayers = mask(speciesLayers,AmSulShape)

nomesSubgraficos = c('Hot&Wet sp','Hot&Dry sp','Cold&Dry sp','8 var./Quadratic','8 var./Quadratic','8 var./Quadratic','8 var./Linear','8 var./Linear','8 var./Linear','2 var./Quadratic','2 var./Quadratic','2 var./Quadratic','2 var./Linear','2 var./Linear','2 var./Linear')

pdf(file='/home/anderson/Documentos/Projetos/Sps artificiais/GLM/realXmodel.pdf')
rasterVis::levelplot(speciesLayers,scales=list(x=list(cex=0.6), y=list(cex=0.6)),between=list(x=1.8, y=0.25),par.strip.text=list(cex=0.6),layout=c(3,5), main='',names.attr=nomesSubgraficos,colorkey=list(space="right")) + layer(sp.polygons(AmSulShape))
dev.off()


##maxent
threHW = read.csv(paste(projectFolder,'Maxent/spHW/StatisticsResults-spHW.csv',sep=''),header=TRUE)$ThresholdMean
HWcurrentReal = raster(paste(projectFolder,'NichoReal/spHW/000.asc',sep='')) > 0.1
HWModel_0kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threHW
HWModel_0kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threHW
HWModel_0kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threHW
HWModel_22kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threHW
HWModel_22kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threHW
HWModel_22kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threHW

threHD = read.csv(paste(projectFolder,'Maxent/spHD/StatisticsResults-spHD.csv',sep=''),header=TRUE)$ThresholdMean
HDcurrentReal = raster(paste(projectFolder,'NichoReal/spHD/000.asc',sep='')) > 0.1
HDModel_0kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threHD
HDModel_0kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threHD
HDModel_0kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threHD
HDModel_22kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threHD
HDModel_22kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threHD
HDModel_22kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threHD

threCD = read.csv(paste(projectFolder,'Maxent/spCD/StatisticsResults-spCD.csv',sep=''),header=TRUE)$ThresholdMean
CDcurrentReal = raster(paste(projectFolder,'NichoReal/spCD/000.asc',sep='')) > 0.1
CDModel_0kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threCD
CDModel_0kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threCD
CDModel_0kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threCD
CDModel_22kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threCD
CDModel_22kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threCD
CDModel_22kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threCD

library(rasterVis)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp")

## levelplot(stack(c(HWcurrentReal*1+HWcurrentModel*2,HDcurrentReal*1+HDcurrentModel*2,CDcurrentReal*1+CDcurrentModel*2)),scales=list(x=list(cex=1), y=list(cex=1)),between=list(x=1.8, y=0.25),par.strip.text=list(cex=1),layout=c(3,1),col.regions=colorRampPalette(c("white","light green","blue","dark green")), main='',names.attr=c('Hot&Wet sp.','Hot&Dry sp.','Cold&Dry sp.'),colorkey=list(space="right",labels=list(cex=1.2))) + layer(sp.polygons(AmSulShape))

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
plot(HWcurrentReal*1+HWModel_22kyrSample5*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_22kyrSample45*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_22kyrSample95*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
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
grid()
plot(HDcurrentReal*1+HDModel_0kyrSample45*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
grid()
plot(HDcurrentReal*1+HDModel_0kyrSample95*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
grid()
plot(HDcurrentReal*1+HDModel_22kyrSample5*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
grid()
plot(HDcurrentReal*1+HDModel_22kyrSample45*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
grid()
plot(HDcurrentReal*1+HDModel_22kyrSample95*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
grid()
mtext('Hot & Dry sp.',outer=TRUE,cex=4)
dev.off()

##sobreposicoes spCD
jpeg(filename=paste(projectFolder,'Maxent/graficos/sobreposicoesCD.jpg',sep=''), width = 1100 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,0))
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
plot(CDcurrentReal*1+CDModel_22kyrSample5*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_22kyrSample45*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_22kyrSample95*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
grid()
mtext('Cold & Dry sp.',outer=TRUE,cex=4)
dev.off()

