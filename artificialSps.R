#########################################################################################
####ALGORITMO PARA FAZER A DISTRIBUICAO REAL DAS ESPECIES EM VARIOS MOMENTOS DO TEMPO####
#########################################################################################

library(virtualspecies)
library(maptools)
library(dismo)
library(raster)
source("/home/anderson/R/R-Scripts/TSSfunction.R")

###PRIMEIRA PARTE: criando sps virtuais###

###Parametros necessarios###
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
projectFolder = "/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/" #pasta do projeto
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
caminhosCamadasTemp = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
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

###SEGUNDA PARTE: amostragem de pontos de ocorrencia em diferentes camadas de tempo para fazer o pooled niche model###

###Parametros necessarios###
Npass = 1:10 #numero de pontos a serem amostrados para camadas do passado (pensando em pontos fosseis)
Npres = c(10,100,200,400,800) #numero de pontos a serem amostrados para camadas do presente
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
projectFolder = "/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/" #pasta do projeto
mainSampleFolder = '/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/Amostras/' #caminho para pasta onde a planilha com os pontos amostrados sera salva
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies

#abrindo e cortando camadas de variaveis ambientais para o presente
filesRaw <- stack(list.files(path=paste(envVarFolder,"dados_projeto/000",sep=''), pattern='asc', full.names=T)) ### stack all rasters in Bioclim folder
#files <- stack(list.files(path = "/home/anderson/R/PosDoc/dados_ambientais/bcmidbi_2-5m _asc/dados_ambientais_para_projeto", pattern='asc', full.names=T))
files = mask(filesRaw,AmSulShape) #cortando para Am. do Sul


############################

for (h in 1:length(spsTypes)){
    nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[h],sep='') #pasta com os mapas de nicho real da sp
    nicheRealPath = list.files(path=nicheRealFolder, full.names=T, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da sp
    for (j in 1:length(Npres)){
        for (k in 1:length(Npass)){
            amostra = data.frame()
            for (i in 1:length(nicheRealPath[1:24])){
                realNicheMap = raster(nicheRealPath[i]) #abrindo o mapa de occ da sp
                #nameScenario = names(realNicheMap)
                if (names(realNicheMap) != 'X000'){ #verificando se e o mapa da ocorrencia no presente
                    amostra_i = sampleOccurrences(realNicheMap,Npass[k],plot=FALSE) #realizando a amostragem
                    layers_i = stack(list.files(path=paste(envVarFolder,'/000',sep=''), pattern='asc', full.names=T)) ### stack all rasters in
                    envVar_i = extract(layers_i, amostra_i[[1]][1:2], method='bilinear', buffer=NULL, fun=NULL, df=TRUE)
                    amostra_i = cbind(amostra_i[[1]][1:2],envVar_i[,2:ncol(envVar_i)])
                }
                else{
                    amostra_i = sampleOccurrences(realNicheMap,Npres[j],plot=FALSE) #realizando a amostragem
                    scenarioName = basename(nicheRealPath[1:24][i])
                    scenarioName = gsub('.asc','',scenarioName)
                    layers_i = stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='asc', full.names=T)) ### stack all rasters in
                    envVar_i = extract(layers_i, amostra_i[[1]][1:2], method='bilinear', buffer=NULL, fun=NULL, df=TRUE)
                    amostra_i = cbind(amostra_i[[1]][1:2],envVar_i[,2:ncol(envVar_i)])
                }
                amostra = rbind(amostra,amostra_i)
            }
            write.csv(amostra,file=paste(mainSampleFolder,spsTypes[h],'/ptsPres',Npres[j],'ptsPass',Npass[k],'.csv',sep=''),row.names=FALSE)#salvando a planilha com os dados da amostra
        }
    }
}

08:33###TERCEIRA PARTE: SDM usando de pontos de ocorrencia em diferentes camadas de tempo (atual a 120 kyr BP)###

#######################################################
####################### GLM ###########################
#######################################################

options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder = "/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = '/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/Amostras/' #caminho para pasta onde a planilha spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
GLMfolder = '/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/GLM/' #pasta para resultados do maxent
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

            #pegando a porcao dos dados separados para a avaliacao (validacao) do modelo
            presencias.evaluacion = sp.data[rand==1,]
            presencias.evaluacion <- cbind(presencias.evaluacion$lon,presencias.evaluacion$lat)
            pseudoausencias.evaluacion = pseudoausencia.data[rand==1,]
            pseudoausencias.evaluacion = cbind(pseudoausencias.evaluacion$lon,pseudoausencias.evaluacion$la)

            ##RODANDO A AVALIACAO DO MODELO##
            evaluacion=evaluate(presencias.evaluacion, pseudoausencias.evaluacion, GLM, predictors)

            #registrando o valor de AUC em um objeto
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


###QUARTA PARTE: comparando projecao do SDM e a distribuicao espacial real do nicho da sp###


projectFolder = "/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/" #pasta do projeto
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
scenarioModel = c('8varLinearModel','8varQuadModel','2varLinearModel','2varQuadModel')

for (i in 1:length(spsTypes)){
    for (j in 1:length(model)){
        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        nicheRealPath = list.files(path=nicheRealFolder,pattern='asc',full.names=T) #lista com os enderecos dos mapas de distribuicao da sp
        projectionsFolder = paste(projectFolder,'GLM/',spsTypes[i],'/',scenarioModel[j],'/Projections',sep='') #pasta com as projecoes do cenario
        projectionsPath = list.files(path=projectionsFolder, pattern='asc',full.names=T) #caminhos para os .asc na paste do cenario
        outputData = data.frame()
        l=1
        for (l in 1:length(nicheRealPath[1:24])){
            realNiche = raster(nicheRealPath[l]) #nicho real
            realNiche.spgrid = as(realNiche,'SpatialGridDataFrame')
            sdmNiche = raster(projectionsPath[l]) #mapa de suitability gerado por SDM
            sdmNiche.spgrid = as(sdmNiche,'SpatialGridDataFrame')
            nicheOverlap = niche.overlap(c(realNiche.spgrid,sdmNiche.spgrid))
            output_i= c(paste(l,'kyrBP',sep=''),nicheOverlap[1,2],nicheOverlap[2,1])
            outputData = rbind(outputData,output_i)
        }
        names(outputData) = c('kyrBP','Schoeners_D','Hellinger_distances')
        

        
        
        
        
write.csv(cor.matrix, file=paste(projectionsFolder,"/NO.csv",sep=""),row.names=TRUE) #salvando os dados do cenario
    }
}



##########################################################################
#########################TESTANDO BIOMOD##################################

#myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
myRespXY = dataset[,c('lon','lat')]
#myResp = rep(1,nrow(myRespXY))
myResp = dataset[,'pres']
predictors = stack(predictors)

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = predictors,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    models = c('GLM','RF'),
    models.options = myBiomodOption,
    NbRunEval=3,
    DataSplit=50,
    Prevalence=0.5,
    VarImport=3,
    models.eval.meth = c('TSS','ROC'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = FALSE,
    modeling.id = paste(myRespName,"FirstModeling",sep=""))

myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = predictors,
    proj.name = 'current',
    selected.models = 'all',
    binary.meth = 'TSS',
    compress = 'xz',
    clamping.mask = F,
    output.format = '.grd'
)



##RODANDO SDMs com BIOMOD
projectFolder = "/home/anderson/Documentos/Minha produção bibliográfica/Sps artificiais/" #pasta do projeto
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies

for (i in 1:length(spsTypes)){

    sampleFolder = paste(mainSampleFolder,spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
    samplePaths = list.files(path=sampleFolder, full.names=T, pattern='.csv') #lista com os enderecos dos mapas de distribuicao da sp
    predictors = stack(list.files(path=envVarPaths[1],full.names=T, pattern='.asc')) #predictors com todas as variaveis
    predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
    sp.data = read.csv(paste(samplePaths[2],sep=''),header=TRUE) #abrindo a planilha de pontos de occ amostrados
    names(sp.data) = c('lon','lat',"bioclim_01","bioclim_04","bioclim_10","bioclim_11","bioclim_12","bioclim_15","bioclim_16","bioclim_17")
    pseudoausencia1.occ <- randomPoints(mask=predictors[[1]], n=nrow(sp.data), p=sp.data[,1:2], excludep=TRUE) #este sera usado no loop para gerar ausencias de teste, la embaixo
    pseudoausencia2.occ <- round(pseudoausencia1.occ[,1:2], digits=4)
    pseudoausencia3.occ <- pseudoausencia2.occ[!duplicated(pseudoausencia2.occ),]
    pseudoausencia4.occ <- pseudoausencia3.occ[complete.cases(pseudoausencia3.occ),]
    pseudoausencia.occ <- data.frame(pseudoausencia4.occ)
    colnames(pseudoausencia.occ) <- c("lon", "lat")
    pseudoausencia.clim <- extract(predictors, pseudoausencia.occ, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)
    pres = c(rep(1, nrow(sp.data)),rep(0, nrow(pseudoausencia.clim)))
    pseudoausencia.data <- cbind(pseudoausencia.occ,pseudoausencia.clim[2:ncol(pseudoausencia.clim)])
    dataset = data.frame(rbind(sp.data,pseudoausencia.data))
    dataset = cbind(dataset,pres)

    setwd(paste(projectFolder,'/BioMOD',sep=''))

    myRespName = paste(spsTypes[i])
    myRespXY = dataset[,c('lon','lat')]
    myResp = dataset[,'pres']
    predictors = stack(predictors)
    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = predictors,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

    myBiomodOption <- BIOMOD_ModelingOptions()

    myBiomodModelOut <- BIOMOD_Modeling(
        myBiomodData,
        models = c('GLM','RF'),
        models.options = myBiomodOption,
        NbRunEval=3,
        DataSplit=75,
        VarImport=3,
        models.eval.meth = c('TSS','ROC'),
        SaveObj = TRUE,
        rescal.all.models = TRUE,
        do.full.models = FALSE,
        modeling.id = paste(myRespName,"FirstModeling",sep=""))
    
    myBiomodProj <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = predictors,
        proj.name = 'current',
        selected.models = 'all',
        binary.meth = 'TSS',
        compress = 'xz',
        clamping.mask = F,
        output.format = '.grd'
    )
}


    
    
