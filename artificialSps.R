#########################################################################################
####SCRIPT PARA FAZER A DISTRIBUICAO REAL DAS ESPECIES EM VARIOS MOMENTOS DO TEMPO####
#########################################################################################

##pacotes necessarios
library(virtualspecies)
library(maptools)
library(dismo)
library(raster)
library(phyloclim) #para funcao niche.overlap()
source("/home/anderson/R/R-Scripts/TSSmaxent.R")
source("/home/anderson/R/R-Scripts/AUCrand.R")


###PRIMEIRA PARTE: criando sps virtuais###


###Parametros necessarios###
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
caminhosCamadasTemp = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
############################

for (i in 1:length(caminhosCamadasTemp)){

    ##variaveis preditoras
    predictors = stack(paste(caminhosCamadasTemp[i],'/bioclim_01.asc',sep=''),paste(caminhosCamadasTemp[i],'/bioclim_12.asc',sep='')) #carregando as variaveis ambientais
    predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
    nameScenario = basename(caminhosCamadasTemp[i])
    
    ##funcoes especie de clima quente e umido
    parametersHW <- formatFunctions(bioclim_01=c(fun='betaFun',p1=200,p2=295,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=2000,p2=3500,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais

    ##funcoes especie de clima quente e seco
    parametersHD <- formatFunctions(bioclim_01=c(fun='betaFun',p1=200,p2=260,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=50,p2=1800,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais
 
    ##funcoes especie de clima frio e seco
    parametersCD <- formatFunctions(bioclim_01=c(fun='betaFun',p1=50,p2=220,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=50,p2=1800,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais

    ##criando distribuicao geografica das sps
    spHW <- generateSpFromFun(predictors, parametersHW) #criando a especie artifical (clima quente e umido)
    spHD <- generateSpFromFun(predictors, parametersHD) #criando a especie artifical (clima quente e seco)
    spCD <- generateSpFromFun(predictors, parametersCD) #criando a especie artifical (clima frio e seco)

    ##empilhando distribuicoes geradas
    auxVector=stack(c(spHW$suitab.raster,spHD$suitab.raster,spCD$suitab.raster))
    names(auxVector) = c('spHW', 'spHD', 'spCD')

    ##ajustando o raster e salvando
    for(j in 1:dim(auxVector)[3]){
        projection(auxVector[[j]]) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
        writeRaster(auxVector[[j]], filename=paste(projectFolder,'NichoReal/',names(auxVector[[j]]),'/',nameScenario,'.asc',sep=""), overwrite=TRUE,prj=TRUE) #salvando o raster do mapa da sp
        
    }
}


###SEGUNDA PARTE: amostragem de pontos de ocorrencia em diferentes camadas de tempo###


##definindo variaveis e parametros
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto/" #pasta com as variaveis ambientais
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha com os pontos amostrados sera salva
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
crs(AmSulShape) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
biasLayer = raster('/home/anderson/Documentos/Projetos/Sps artificiais/biasLayer.grd')
biomodFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/biomod/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
sampleSizes = c(5,15,25,35,45,55,65,75,85,95) #tamanhos das amostras
NumRep = 10 #numero de replicas (de cada cenario amostral)
Tmax = 22 #idade maxima (no passado)
bgPoints = 1000 #numero de pontos de background


##PARA SDM MULTITEMPORAL E SEM VIES AMOSTRAL


##ocorrencias

sampleData = data.frame()

for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
    
    for (sSize in sampleSizes){ #numero de pontos (registros, dados) na amostra
        
        sampledAges = vector()
        sampledAges = round(runif(sSize,0,Tmax)) #selecionando 'n' camadas de tempo aleatoriamente
        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da
        
        for (j in 1:NumRep){ #replicas do cenario amostral

            for (sAge in sampledAges){ #amostrando em cada camada de tempo que consta na amostra

                sampleData_i = randomPoints(mask=raster(nicheRealPath[sAge+1]), n=1) #amostra d ponto
                scenarioName = basename(nicheRealPath[1:24][sAge+1]) #tempo vinculado ao cenario para variaveis ambientais
                scenarioName = gsub('.asc','',scenarioName) #retirando do nome o '.asc'
                layers_i = extract(
                    x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                    y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
                sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
                
            }
            
            names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
            write.csv(sampleData,paste(projectFolder,'Amostras/multitemporal/',spsTypes[i],'/occ_',sSize,'pts_multitemporal_', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
            sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
            
        }
    }
}

##background points

sampleData = data.frame()

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
            write.csv(sampleData,paste(projectFolder,'Amostras/multitemporal/',spsTypes[i],'/bg_',sSize, 'pts_multitemporal_', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
            sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
            
        }
    }
}


##PARA SDM MULTITEMPORAL E COM VIES AMOSTRAL


##ocorrencias

sampleData = data.frame()

for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
    
    for (sSize in sampleSizes){ #numero de pontos (registros, dados) na amostra
        
        sampledAges = vector()
        sampledAges = round(runif(sSize,0,Tmax)) #selecionando 'n' camadas de tempo aleatoriamente
        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da
        biasLayerAdjusted = projectRaster(biasLayer,raster(nicheRealPath[1])) #alinhando o biasLayer com os rasters do projeto
            
        for (j in 1:NumRep){ #replicas do cenario amostral

            for (sAge in sampledAges){ #amostrando em cada camada de tempo que consta na amostra

                sampleData_i = randomPoints(mask=raster(nicheRealPath[sAge+1])*biasLayerAdjusted, n=1, prob=TRUE) #amostra d ponto
                scenarioName = basename(nicheRealPath[1:24][sAge+1]) #tempo vinculado ao cenario para variaveis ambientais
                scenarioName = gsub('.asc','',scenarioName) #retirando do nome o '.asc'
                layers_i = extract(
                    x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                    y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
                sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
                
            }
            
            names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
            write.csv(sampleData,paste(projectFolder,'Amostras/multitemporal/',spsTypes[i],'/occ_',sSize,'pts_multitemporal_comVIES_', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
            sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada

        }
    }
}

##background points

sampleData = data.frame()

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
                    x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='.asc', full.names=TRUE)),
                    y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
                sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
                
            }
            
            names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
            write.csv(sampleData,paste(projectFolder,'Amostras/multitemporal/',spsTypes[i],'/bg_',sSize,'pts_multitemporal_comVIES_', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
            sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
            
        }
    }
}


##PARA SDM MONOTEMPORAL (SO PONTOS DO TEMPO PRESENTE) E SEM VIES AMOSTRAL


##ocorrencias
for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'

    ##criando uma pasta da especie, se nao exisitir
    if(!file.exists(file.path(projectFolder,'Amostras','monotemporal',spsTypes[i],sep='')))
        dir.create(file.path(projectFolder,'Amostras','monotemporal',spsTypes[i],sep=''))
    
    for (sSize in sampleSizes){ #numero de pontos (registros, dados) na amostra
        
        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da
        
        for (i in 1:NumRep){ #replicas do cenario amostral

            sampleData_i = randomPoints(mask=raster(nicheRealPath[1],crs=CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')), n=sSize) #amostra do ponto
            scenarioName = basename(nicheRealPath[1:24][1]) #tempo vinculado ao cenario para variaveis ambientais
            scenarioName = gsub('.asc','',scenarioName) #retirando do nome o '.asc'
            layers_i = extract(
                x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
            sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
            names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
            write.csv(sampleData,paste(projectFolder,'Amostras/monotemporal/',spsTypes[i],'/occ_',sSize,'pts_monotemporal_','rep',i,'.csv',sep=''),row.names=FALSE) #salvando
            sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
            
        }
    }
}

##background points
for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'

    ##criando uma pasta da especie, se nao exisitir
    if(!file.exists(file.path(projectFolder,'Amostras','monotemporal',spsTypes[i],sep='')))
        dir.create(file.path(projectFolder,'Amostras','monotemporal',spsTypes[i],sep=''))
    
    for (sSize in sampleSizes){ #numero de pontos (registros, dados) na amostra
        
        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da
        
        for (i in 1:NumRep){ #replicas do cenario amostral

            sampleData_i = randomPoints(mask=raster(envVarPath[1],crs=CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')), n=bgPoints) #amostra do ponto
            scenarioName = basename(nicheRealPath[1:24][1]) #tempo vinculado ao cenario para variaveis ambientais
            scenarioName = gsub('.asc','',scenarioName) #retirando do nome o '.asc'
            layers_i = extract(
                x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
            sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
            names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
            write.csv(sampleData,paste(projectFolder,'Amostras/monotemporal/',spsTypes[i],'/bg_',sSize,'pts_monotemporal_','rep',i,'.csv',sep=''),row.names=FALSE) #salvando
            sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
            
        }
    }
}


###TERCEIRA PARTE: SDM usando de pontos de ocorrencia em diferentes camadas de tempo (do atual ate 120 kyr BP)###


#######################################################
####################### MAXENT ########################
#######################################################

##pacotes
library(biomod2)

##definindo variaveis e parametros
options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha 
maxentFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
sdmTypes = c('multitemporal','monotemporal')
source("/home/anderson/R/R-Scripts/TSSmaxent.R")
sampleSize = c(5,15,25,35,45,55,65,75,85,95)
NumRep = 10 #numero de replicas (de cada cenario amostral)

##algoritmo da analise do projeto
for (h in 1:length(sdmTypes)){
    for (i in 1:length(spsTypes)){
        for (j in 1:length(sampleSizes)){
            for (k in 1:NumRep){ #loop sobre o numero de replicas 

                ##ajustando o diretorio de trabalho (pois o biomod roda e salva tudo simultaneamente)
                if(!file.exists(file.path(projectFolder,'maxent',sdmTypes[h], spsTypes[i],sep='')))
                    dir.create(file.path(projectFolder,'maxent',sdmTypes[h],spsTypes[i],sep=''))
                setwd(file.path(projectFolder,'maxent',sdmTypes[h],spsTypes[i]))
                
                ##definindo variaveis e parametros locais
                statResults = data.frame() #tabela de estatisticas basicas do modelo
                occPoints = read.csv(paste(mainSampleFolder,sdmTypes[h],'/',spsTypes[i],'/occ',sampleSizes[j],'pts',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia    
                backgroundPoints = read.csv(paste(mainSampleFolder,sdmTypes[h],'/',spsTypes[i],'/bg.csv',sep=''),header=TRUE) #abrindo pontos de background
                backgroundPoints = backgroundPoints[sample(nrow(backgroundPoints),2000),]
                ##agrupando ocorrencias e pseudo-ausencias
                names(backgroundPoints) = names(occPoints) #certificando que os nomes das colunas estão iguais (cuidado aqui...)
                dataSet = data.frame(cbind(rbind(occPoints,backgroundPoints),pres=c(rep(1,nrow(occPoints)),rep(0,nrow(backgroundPoints))))) #planilha de dados no formato SWD
                ##variaveis e parametros locais especificos para o biomod2
                myRespName <- paste(spsTypes[i],'_sample',sampleSizes[j],'_replica',k,sep='') # nome do cenario atual (para biomod2)
                myResp <- dataSet[,c('pres')] # variavel resposta (para biomod2)
                myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
                myExpl = dataSet[,c('bioclim_01','bioclim_12')]  #variavel preditora (para biomod2)

                ##ajuste de dados de entrada para biomod2
                myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                                     expl.var = myExpl,
                                                     resp.xy = myRespXY,
                                                     resp.name = myRespName)

                ## ##inspecionando o objeto gerado pela funcao do biomod2
                ## myBiomodData
                ## plot(myBiomodData)

                ##parametrizando os modelos
                myBiomodOption <- BIOMOD_ModelingOptions(
                    MAXENT.Phillips=list(
                        path_to_maxent.jar="/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java",
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
                    NbRunEval = 3,
                    DataSplit = 75,
                    VarImport = 5,
                    models.eval.meth = c('TSS','ROC'),
                    SaveObj = TRUE,
                    rescal.all.models = TRUE,
                    do.full.models = FALSE,
                    modeling.id = paste(myRespName))
                
                ##My output data
                evaluationScores = get_evaluations(myBiomodModelOut)

                ##gravando estatistcas basicas do modelo
                statResults = rbind(statResults,cbind(
                                                    sp = spsTypes[i],
                                                    sampleSize = sampleSizes[j],
                                                    replicate = k,
                                                    AUC = mean(evaluationScores['ROC','Testing.data',,,]),
                                                    TSS = mean(evaluationScores['TSS','Testing.data',,,]),
                                                    numbOfTimeLayers = length(unique(occPoints$kyrBP)),
                                                    medianKyr = median(occPoints$kyrBP),
                                                    minAge = min(occPoints$kyrBP),
                                                    maxAge = max(occPoints$kyrBP)))

                write.csv(statResults,file=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/StatisticalResults-',spsTypes[i],'.csv',sep=''),row.names=FALSE)

                ##implementando projecoes do modelo
                for (l in 1:length(envVarPaths[1:24])){

                    ##definindo variaveis e parametros internos
                    predictors = stack(list.files(path=envVarPaths[l],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
                    predictors = predictors[[c('bioclim_01','bioclim_12')]]
                                        #                predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
                    crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS

                    ##rodando algortmo de projecao (i.e. rodando a projecao)
                    myBiomodProj <- BIOMOD_Projection(
                        modeling.output = myBiomodModelOut,
                        new.env = predictors,
                        proj.name = paste(l,'kyr',sep=''),
                        selected.models = 'all',
                        compress = 'FALSE',
                        build.clamping.mask = 'FALSE',
                        output.format = '.grd')

                    ##gerando e salvando um mapa binario (threshold 10%)
                    projStack = get_predictions(myBiomodProj) #extrai as projecoes
                    projStackBIN = BinaryTransformation(stack(mean(projStack)),'10')
                    writeRaster(projStackBIN,file=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[j],'.replica',k,'/proj_',l,'kyr/proj_',i,'kyr','.sample',sampleSizes[j],'.replica',k,'_BIN.asc',sep=''),row.names=FALSE)

                    
                }
            }
        }
    }
}


###QUARTA PARTE: comparando projecao do SDM e a distribuicao espacial real do nicho da sp###


##abrindo pacotes necessarios
library(raster)
library(ecospat)

##definindo variaveis e parametros
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha com os pontos amostrados sera salva
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
sdmTypes = c('multitemporal','monotemporal')
sampleSizes = c(5,15,25,35,45,55,65,75,85,95) #aqui, deve ser igual ao usasado nas partes anteriores do script
NumRep = 10
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul

##algoritmo da analise do projeto
for (h in 1:length(sdmTypes))
    for (i in 1:length(spsTypes)){

        ##definindo variaveis e parametros locais
        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        nicheRealPath = list.files(path=nicheRealFolder,pattern='asc',full.names=TRUE) #lista com os enderecos dos mapas de distribuicao da sp
        outputData = data.frame() #tabela de dados de saida

        #loop sobre as cadamdas de tempo
        for (l in 1:length(nicheRealPath[1:24])){ 

            ##definindo variaveis e parametros locais
            realNiche = nicheRealPath[l] #nicho real

            ##amostrando pontos da distribuicao real para compracao dos SDMs
            binMap = raster(realNiche,crs=CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'))>0.1 #mapa binario do real
            realNicheDataOccCoord = dismo::randomPoints(binMap,1000) #amostrando 1000 pontos do binario real
            realNicheDataOccPres = extract(binMap,realNicheDataOccCoord,na.rm=TRUE) #definindo occ e ausencias para os pontos
            realNicheDataOcc = data.frame(cbind(realNicheDataOccCoord,realNicheDataOccPres)) #tabela lon, lat e pres
            names(realNicheDataOcc) = c('longitude','latitude','pres') #ajustando nomes das colunas
            predictors = stack(list.files(path=envVarPaths[l],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis
            predictors = predictors[[c('bioclim_01','bioclim_12')]] #selecionando as variaveis usadas
            predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
            projection(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
            realNicheDataPred = extract(x=predictors,y=realNicheDataOcc[,c('longitude','latitude')],na.rm=TRUE) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
            realNicheData = cbind(realNicheDataOcc, realNicheDataPred) #juntando com os dados das outras camadas de tempo amostradas
            
            #loop sobre os tamanhos amostrais
            for (m in sampleSizes){ 




##################################################################################################################################
##################################################################################################################################
C O N T I N U A R    D A Q U I ! ! !
##################################################################################################################################
##################################################################################################################################

                
#                timeSampleData = list.files(path=projectionsFolder, pattern=glob2rx(paste('*Time',l-1,'*Sample',m,'.asc',sep='')),full.names=TRUE)
                timeSampleData = list.files(path=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],sep=''), pattern=glob2rx(paste('*',l,'kyr','*sample',m,'*.asc',sep='')), recursive=TRUE, full.names=TRUE) #pasta com as projecoes do cenario

                ##loop sobre replicas de cada combinacao de tempo e tamanho amostral
                for(n in 1:NumRep){ 

                    ##definindo variaveis e parametros locais
                    projectionsFolder = paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',m,'.replica',n,sep='') #pasta com as projecoes do cenario
                    projectionsPath = list.files(path=projectionsFolder, pattern='.asc',recursive=TRUE,full.names=T) #caminhos para os .asc na paste do cenario


                    x = paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',m,'.replica',n,'/proj_',l,'kyr/','proj_',l,'kyr_',spsTypes[i],'.sample',m,'.replica',n,'.grd',sep='') #pasta com as projecoes do cenario
                    
                    
                    sdmNiche = x #timeSampleData[[n]] #mapa de suitability gerado por SDM
                    binMapSDM = raster(sdmNiche)>0.1
                    SDMDataOccCoord = dismo::randomPoints(binMapSDM,1000)
                    SDMDataOccPres = extract(binMapSDM,SDMDataOccCoord,na.rm=TRUE)
                    SDMDataOcc = data.frame(cbind(SDMDataOccCoord,SDMDataOccPres))
                    names(SDMDataOcc) = c('longitude','latitude','pres')      

                    SDMDataPred = extract(x=predictors,y=SDMDataOcc[,c('longitude','latitude')],na.rm=TRUE) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
                    SDMData = cbind(SDMDataOcc, SDMDataPred) #juntando com os dados das outras camadas de tempo amostradas

                    ##The PCA is calibrated on all the sites of the study area
                    pca.env <-dudi.pca(rbind(realNicheData,SDMData)[,c('bioclim_01','bioclim_12')],scannf=F,nf=2)
                    ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) #grafico

                    ##PCA scores for the whole study area
                    scores.globclim <- pca.env$li
                    ##PCA scores for the species native distribution
                    scores.sp.realNiche <- suprow(pca.env,realNicheData[which(realNicheData[,'pres']==1),c('bioclim_01','bioclim_12')])$li

                    ##PCA scores for the species invasive distribution
                    scores.sp.SDMniche <- suprow(pca.env,SDMData[which(SDMData[,'pres']==1),c('bioclim_01','bioclim_12')])$li

                    ##PCA scores for the whole native study area
                    scores.clim.realNiche <-suprow(pca.env,realNicheData[,c('bioclim_01','bioclim_12')])$li

                    ##PCA scores for the whole invaded study area
                    scores.clim.SDMniche <- suprow(pca.env,SDMData[,c('bioclim_01','bioclim_12')])$li

                    ##gridding the native niche
                    grid.clim.realNiche <-ecospat.grid.clim.dyn(glob=scores.globclim,glob1=scores.clim.realNiche,sp=scores.sp.realNiche, R=100,th.sp=0)

                    ##gridding the invasive niche
                    grid.clim.SDMniche <- ecospat.grid.clim.dyn(glob=scores.globclim,glob1=scores.clim.SDMniche,sp=scores.sp.SDMniche, R=100,th.sp=0)

                    ##equivalencia de nicho
                    ##OBS: Niche equivalency test H1: Is the overlap between the native and invaded niche higher than two random niches?
                    eq.test <- ecospat.niche.equivalency.test(grid.clim.realNiche, grid.clim.SDMniche,rep=10, alternative = "greater")

                    
                    ## nicheOverlapObs = niche.overlap(c(sdmNiche,realNiche))
                    ## Dobs = nicheOverlapObs[1,2]
                    ## Iobs = nicheOverlapObs[2,1]
                    
                    ## ##aleatorizando ocorrencias para teste de significancia
                    ## occPoints = read.csv(paste(mainSampleFolder,'/',spsTypes[i],'/occ',m,'pts',n,'rep.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia
                    ## occPoints[occPoints==0] = NA
                    ## occPoints = occPoints[complete.cases(occPoints),]
                    ## occPoints = round(occPoints, digits=2)
                    ## occPoints = occPoints[!duplicated(occPoints),]                 
                    
                    ## backgroundPoints = read.csv(paste(mainSampleFolder,'/',spsTypes[i],'/bg.csv',sep=''),header=TRUE) #abrindo pontos de background
                    ## backgroundPoints = backgroundPoints[sample(nrow(backgroundPoints),5000),]
                    ## backgroundPoints[backgroundPoints==0] = NA
                    ## backgroundPoints = backgroundPoints[complete.cases(backgroundPoints),]
                    ## backgroundPoints = round(backgroundPoints, digits=2)
                    ## backgroundPoints = backgroundPoints[!duplicated(backgroundPoints),]                 
                    
                    ## names(backgroundPoints) = names(occPoints) #certificando que os nomes das colunas estão iguais (cuidado aqui...)
                    ## dataSet = data.frame(cbind(rbind(occPoints,backgroundPoints),pres=c(rep(1,nrow(occPoints)),rep(0,nrow(backgroundPoints))))) #planilha de dados no formato SWD
                    
                    ## # for(o in 1:100){ #replicas da distribuicao nula de D e I
                    ## #   
                    ## #   presMix = dataSet[sample(nrow(dataSet)),c('pres')]
                    ## #   sampleMix = dataSet
                    ## #   sampleMix$pres = presMix
                    ## #   
                    ## #   me = maxent(
                    ## #     x=sampleMix[,c("bioclim_01","bioclim_12")],
                    ## #     p=sampleMix$pres,
                    ## #     args=c('responsecurves=FALSE',
                    ## #            'jackknife=FALSE',
                    ## #            'randomseed=FALSE',
                    ## #            'randomtestpoints=0',
                    ## #            'maximumbackground=5000',
                    ## #            'replicates=1',
                    ## #            'writebackgroundpredictions=FALSE',
                    ## #            'linear=TRUE',
                    ## #            'quadratic=TRUE',
                    ## #            'product=FALSE',
                    ## #            'threshold=FALSE',
                    ## #            'hinge=FALSE',
                    ## #            'maximumiterations=1000',
                    ## #            'convergencethreshold=1.0E-5',
                    ## #            'threads=2'
                    ## #     ))
                    ## #   
                    ## #   proj = predict(me,predictors,crs=crsInfo) #realizando projetacoes (para cada replica)                    
                    ## #   
                    ## #   sdmNicheMix = as(proj,'SpatialGridDataFrame')
                    ## #   
                    ## #   nicheOverlap_i= niche.overlap(c(sdmNicheMix,realNiche))
                    ## #   Ddistribution = append(Ddistribution,nicheOverlap_i[1,2])
                    ## #   Idistribution = append(Idistribution,nicheOverlap_i[2,1])
                    ## #   
                    ## # }

                    Dobs = eq.test$obs$D
                    Iobs = eq.test$obs$I
                    DpValue = eq.test$p.D  #sum(Ddistribution >= Dobs) / length(Ddistribution)
                    IpValue = eq.test$p.I #sum(Idistribution >= Iobs) / length(Idistribution)

                    ##abrindo planilha de pontos para extrair dados do cenario
                    occPoints = read.csv(paste(mainSampleFolder,'/',spsTypes[i],'/occ',m,'pts',n,'rep.csv',sep=''),header=TRUE) 
                    occPoints[occPoints==0] = NA
                    occPoints = occPoints[complete.cases(occPoints),]
                    occPoints = round(occPoints, digits=2)
                    occPoints = occPoints[!duplicated(occPoints),]                 
                    
                    outputData = rbind(outputData,cbind(kyrBP=l-1,
                                                        sampleSize=m,
                                                        replicate=n,
                                                        numbOfTimeLayers=length(unique(occPoints$kyrBP)),
                                                        medianKyr=median(unique(occPoints$kyrBP)),
                                                        minAge=min(unique(occPoints$kyrBP)),
                                                        maxAge=max(unique(occPoints$kyrBP)),
                                                        Schoeners_D=Dobs,
                                                        p_value=DpValue,
                                                        Hellinger_I=Iobs,
                                                        p_value=IpValue))
                    
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

###GRAFICOS

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


###CASO DA HIENA, Crocuta crocuta

#abrindo os pacotes necessarios

library(raster)
library(biomod2)

#fixando caminhos das pastas do projeto e carregando funcoes proprias

options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
maxentFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/' #pasta para resultados do maxent
source("/home/anderson/R/R-Scripts/TSSmaxent.R")

#obtendo pontos de ocorrencia (para o PRESENTE)

occDataRaw = spocc::occ(query='Crocuta crocuta',from='gbif',limit=10000) #baixando os dados
occDataRawDF = data.frame(occDataRaw$gbif$data$Crocuta_crocuta[,c('name','longitude','latitude','basisOfRecord')])
occDataRawDFcurrent = occDataRawDF[occDataRawDF$basisOfRecord!="FOSSIL_SPECIMEN",]
occDataRawDFcoord1 = occDataRawDFcurrent[complete.cases(occDataRawDFcurrent),] #retirando dados incompletos
occDataRawDFcoord2 = round(occDataRawDFcoord1[,c('longitude','latitude')], digits = 2) #arredondado para duas casas decimais
occDataRawDFcoord3 = occDataRawDFcoord2[!duplicated(occDataRawDFcoord2), ] #retirando dados redundantes (pontos repetidos)
occData = data.frame(cbind(species='Crocuta crocuta',occDataRawDFcoord3)) #tabela de dados de ocorrencia final

##inspecionado os pontos de ocorrencia

maps::map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0)) #plotando o mapa do mundo
occShape = rgdal::readOGR('/home/anderson/Documentos/Projetos/Sps artificiais/hiena/iucn_distribution/species_5674.shp') #shape IUCN
plot(occShape,add=TRUE,col=rgb(1,0,0,0.5)) #plotando shape da IUCN para a distribuição da especie
points(occData[,c('longitude','latitude')], col="blue", pch=16) #plotando os pontos de ocorrencia (meus dados)

##eliminando os pontos que estao fora do poligono da area de ocorrencia

pts = occData[,c('longitude','latitude')] #pegando apenas long e lat
coordinates(pts) <- ~ longitude + latitude #transformando em spatialPoints
proj4string(pts) = proj4string(occShape) #transformando em spatialPoints
ptsInside = pts[!is.na(sp::over(pts,as(occShape,"SpatialPolygons")))] #sobrepondo pontos e poligono

points(ptsInside,pch=20,cex=0.5) #plotando para checar (em cima do ultimo grafico!)

#obtendo pontos de ocorrencia (para o PASSADO)

fossilDataRaw = paleobioDB::pbdb_occurrences(limit="all", base_name="Crocuta crocuta",show=c("coords", "phylo", "ident")) #dados do PBDB

#inspecionando pontos fosseis

X11(width=13, height=7.8); paleobioDB::pbdb_map(fossilDataRaw,pch=19,col.point=c("pink","red"), col.ocean="light blue",main="Crocuta crocuta") #mapa dos pontos 
X11(width=13, height=7.8); paleobioDB::pbdb_map_occur(fossilDataRaw) #mapa do eforco amostral
paleobioDB::pbdb_temp_range(fossilDataRaw,rank="species") #range temporal dos registros

#'limpando' os dados de registros fosseis

occDataFossilLatLong = fossilDataRaw[,c("lng","lat","eag","lag")] #pegando apenas long e lat do PBDB data
occDataFossilRound = round(occDataFossilLatLong, digits=2) #arredondando lat e long para 2 casas decimais
occDataFossilClean1 = occDataFossilRound[complete.cases(occDataFossilRound),] #retirando dados incompletos
occDataFossilClean2 = occDataFossilClean1[!duplicated(occDataFossilClean1),] #retirando pontos em sobreposicao
occDataFossil = occDataFossilClean2 #gravando objeto com o conjundo de dados final
occDataFossil = cbind(occDataFossil[,c('lng','lat')],age=(occDataFossil$eag-occDataFossil$lag)/2)

##salvando dados de ocorrencia "tratados"

write.csv(ptsInside, file="/home/anderson/Documentos/Projetos/Sps artificiais/hiena/Crocuta_crocuta_Occ.csv",row.names=FALSE)
write.csv(occDataFossil, file="/home/anderson/Documentos/Projetos/Sps artificiais/hiena/Crocuta_crocuta_Fossil.csv",row.names=FALSE)


#############################################
#############################################
#############################################


##modelando a distribuição das espécies: MODELO 'MONOTEMPORAL'

##pacotes
library(biomod2)

##abrindo dados de ocorrencia "tratados"
occData = read.csv(file="/home/anderson/Documentos/Projetos/Sps artificiais/hiena/Crocuta_crocuta_Occ.csv",header=T,stringsAsFactors=FALSE)
occDataFossil = read.csv(file="/home/anderson/Documentos/Projetos/Sps artificiais/hiena/Crocuta_crocuta_Fossil.csv",header=T,stringsAsFactors=FALSE)    

# the name of studied species
myRespName <- 'Hiena'

# the XY coordinates of species data
myRespXY <- occData
coordinates(myRespXY) <- ~ longitude + latitude #transformando em spatialPoints
crs(myRespXY) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #transformando em spatialPoints


# load the environmental raster layers
myExpl = raster::stack(list.files(envVarPaths[1],pattern='bioclim',full.names=TRUE))
crs(myExpl) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #transformando em spatialPoints

myBiomodData <- BIOMOD_FormatingData(resp.var = myRespXY,
                                     expl.var = myExpl,
                                     resp.name = myRespName,
                                     PA.nb.rep = 1)

#inspecionando o objeto gerado pela funcao do biomod2
myBiomodData
plot(myBiomodData)



myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips=list(path_to_maxent.jar="/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java",maximumiterations=2000,memory_allocated=NULL))
    
myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    models = c('MAXENT.Phillips'),
    models.options = myBiomodOption,
    NbRunEval = 3,
    DataSplit = 75,
    VarImport = 3,
    models.eval.meth = c('TSS','ROC'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = FALSE,
    modeling.id = paste(myRespName,"FirstModeling",sep=""))
    
###PROJECAO PARA O PRESENTE###
    myBiomodProj <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = myExpl,
        proj.name = '000kyrBP',
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
