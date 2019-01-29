#########################################################################################
####SCRIPT PARA FAZER A DISTRIBUICAO REAL DAS ESPECIES EM VARIOS MOMENTOS DO TEMPO####
#########################################################################################

##pacotes necessarios
library(virtualspecies)
library(maptools)
library(dismo)
library(raster)
library(phyloclim) #para funcao niche.overlap()
#source("/home/anderson/R/R-Scripts/TSSmaxent.R")
#source("/home/anderson/R/R-Scripts/AUCrand.R")


###PRIMEIRA PARTE: criando sps virtuais###


###Parametros necessarios###
##Workstation
## envVarFolder = "J:/Anderson_Eduardo/dados_projeto" #pasta com as variaveis ambientais
## caminhosCamadasTemp = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
## projectFolder = "J:/Anderson_Eduardo/spsArtificiais/" #pasta do projeto
## AmSulShape = rgdal::readOGR("J:/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
##Meu notebook
envVarFolder = "/home/anderson/gridfiles/dados_projeto" #pasta com as variaveis ambientais
caminhosCamadasTemp = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
projectFolder = "/home/anderson/Projetos/Sps artificiais/" #pasta do projeto
AmSulShape = rgdal::readOGR("/home/anderson/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
############################


# ynorm = (dnorm(temp,230,10) - min(dnorm(temp,230,10)))/( max(dnorm(temp,230,10)) - min(dnorm(temp,230,10)) )
# plot(ynorm ~ temp)
# points(betaTemp ~ temp, t='l')
# 
# ynorm = (dnorm(prec,2750,300) - min(dnorm(prec,2750,300)))/( max(dnorm(prec,2750,300)) - min(dnorm(prec,2750,300)) )
# plot(ynorm ~ prec)
# points(betaPrec ~ prec, t='l')
# 
# 
# for (i in 1:length(caminhosCamadasTemp)){
# 
#     ##variaveis preditoras
#     predictors = stack(paste(caminhosCamadasTemp[i],'/bioclim_01.asc',sep=''),paste(caminhosCamadasTemp[i],'/bioclim_12.asc',sep='')) #carregando as variaveis ambientais
#     predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
#     nameScenario = basename(caminhosCamadasTemp[i])
#     
#     ##funcoes especie de clima quente e umido
#     parametersHW <- formatFunctions(bioclim_01=c(fun='betaFun',p1=200,p2=295,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=2000,p2=3500,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais
# 
#     ##funcoes especie de clima quente e seco
#     parametersHD <- formatFunctions(bioclim_01=c(fun='betaFun',p1=200,p2=260,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=50,p2=1800,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais
#  
#     ##funcoes especie de clima frio e seco
#     parametersCD <- formatFunctions(bioclim_01=c(fun='betaFun',p1=50,p2=220,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=50,p2=1800,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais
# 
#     ##criando distribuicao geografica das sps
#     spHW <- generateSpFromFun(predictors, parametersHW) #criando a especie artifical (clima quente e umido)
#     spHD <- generateSpFromFun(predictors, parametersHD) #criando a especie artifical (clima quente e seco)
#     spCD <- generateSpFromFun(predictors, parametersCD) #criando a especie artifical (clima frio e seco)
# 
#     ##empilhando distribuicoes geradas
#     auxVector=stack(c(spHW$suitab.raster,spHD$suitab.raster,spCD$suitab.raster))
#     names(auxVector) = c('spHW', 'spHD', 'spCD')
# 
#     ##ajustando o raster e salvando
#     for(j in 1:dim(auxVector)[3]){
#         projection(auxVector[[j]]) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
#         writeRaster(auxVector[[j]], filename=paste(projectFolder,'NichoReal/',names(auxVector[[j]]),'/',nameScenario,'.asc',sep=""), overwrite=TRUE,prj=TRUE) #salvando o raster do mapa da sp
#         
#     }
# }

# ##definindo prametros e variaveis globais (NOTEBOOK)
envVarFolder = "/home/anderson/gridfiles/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("/home/anderson/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
# maxentFolder = '/home/anderson/R/x86_64-pc-linux-gnu-library/3.4/dismo/java/maxent.jar' #pasta para resultados do maxent
# ## spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
# ## sdmTypes = c('normal','optimized')
# sampleSizes = c(10,20,40,80,160)
# NumRep = 10 #numero de replicas (de cada cenario amostral)
# vies_levels = 5
# ##variaveis preditoras
# ## elevation = raster('/home/anderson/PosDoc/dados_ambientais/DEM/DEM.tif')
predictors = stack(list.files(path=envVarPaths[1],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
predictors = predictors[[c('bioclim_01','bioclim_12')]]
predictors = stack(mask(x=predictors, mask=AmSulShape))
crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
# Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
# statResultsSDMnormal = data.frame() #tabela de estatisticas basicas do modelo
# statResultsSDMimproved = data.frame()

## funcoes autorais :-) #NOTEBOOK
source('/home/anderson/R-Scripts/makeSpeciesSuitability.R')
source('/home/anderson/R-Scripts/makeSpeciesSuitabilityIterator.R')
source('/home/anderson/R-Scripts/rangeByACcontinuous.R')
source('/home/anderson/R-Scripts/makeOutput.R')
source('/home/anderson/R-Scripts/bestModel.R')
source('/home/anderson/R-Scripts/temporalUpdate.R')



#obtaining species distribution using simulation
# spsSuit = makeSpeciesSuitability(predictors)
# SpDistAC = rangeByACcontinuous(spsSuit)
# plot(SpDistAC)
# plot(AmSulShape, add=T)


SAcells = sum(freq(predictors[[2]]>=0, useNA='no')[,2]) #number of cells in South America
prevCorrupted = TRUE
iterLimit = 0
while(any(prevCorrupted)){
  
  #obtendo distribuicao de suitability
  spsSuit = makeSpeciesSuitability(predictors)
  
  #computando prevalencia do suitability (ao longo do tempo, quando houverem conjuntos de dados para cada idade)
  spsPrevalData =  lapply(seq(length(spsSuit)), function(i){
    spsPreval = freq(spsSuit[[i]])[2,2] / SAcells #sps prevalence
    return(spsPreval)
  })
  
  #verificando se a prevalencia esta dentro dos limites estabelecidos (aqui: maior que 5% e menor que 60%)
  prevCorrupted = sapply(seq(length(spsPrevalData)), function(i) {(spsPrevalData[i] < 0.05) | (spsPrevalData[i] > 0.6)})
  iterLimit = iterLimit + 1
  
  #helper para prevenir loop infinito
  if (iterLimit == 11){
    stop("Numero maximo de tentativas realizado sem obter distribuicoes dentro dos parametros. Tente novamente!")
  }
  
}

#updating sps distribution
newSuitabilityMap = rasterSpDistribution
SpDistACup = temporalUpdate(currentSpsRange = SpDistAC, newSuitabilityMap = rasterSpDistribution)
plot(SpDistACup)

## CONTINUAR EM: TEMPORAL UPDATE










#######################################################################

seedCoords = SpDistDF[which(SpDistDF$DistUpdated>0),c('x','y')]

seedCoords = seedCoords[ sample(x=nrow(seedCoords), size=ceiling(0.01*nrow(seedCoords))), ]

spRange_i = SpatialPoints(seedCoords)
crs(spRange_i) = CRS('+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84 ')

###growing up sps' range###
for (i in 1:5){
  
  ##criando Buffer
  rangeBuff = rgeos::gBuffer(spRange_i, width=1) #capStyle="SQUARE"
  
  ##pegando celulas = zero dentro do Buffer e executando ocupacao usando suitability como probabilidade de occ.
  vals = extract(suitT2, rangeBuff, cellnumbers=TRUE)[[1]]
  vals2 = vals[which(vals[,2] > 0),]
  probOcc = runif(n=nrow(vals2))
  vals3 = vals2[vals2[,2] > probOcc, 1]
  spRangeUpdated[vals3] = 1 #; plot(spRangeUpdated)
  
  ##fazendo Buffer sobre Buffer
  #spRange_i = rangeBuff
  
  ##pontos para o novo Buffer
  ptsDistUpadatedRaw = extract(spRangeUpdated, rangeBuff, cellnumbers=TRUE)[[1]]
  ptsDistUpadated = ptsDistUpadatedRaw[which(ptsDistUpadatedRaw[,2]==1),]
  ptsMigrationTotal = xyFromCell(spRangeUpdated,ptsDistUpadated[,1])
  ptsMigrationSelected = ptsMigrationTotal[sample(x=nrow(ptsMigrationTotal), size=ceiling(0.01*nrow(ptsMigrationTotal))),]
  spRange_i = SpatialPoints(ptsMigrationSelected)
  
  #plot(spRangeUpdated)
}
####


r <- raster(nrow=18, ncol=36, xmn=0)
r[251:450] <- 1
plot(r)
plot( boundaries(r, type='inner') )
plot( boundaries(r, type='outer') )
plot( boundaries(r, classes=TRUE) )

bouter = boundaries(r, type='outer')
plot(stack(bouter,r))
smd = sum(stack(bouter,r), na.rm=T)
plot(stack(bouter,r,smd))
plot(smd)
plot(sum(bouter,smd,na.rm=T))




SpDistAC2 = SpDistAC
values(SpDistAC2)[ is.na(values(SpDistAC2)) ] = 0
plot(SpDistAC2)
xx = boundaries(SpDistAC2, asNA=FALSE, classes=TRUE, type='outer')
plot( xx )
xxs = stack(xx,SpDistAC2)
plot( xxs )
plot(sum(xxs, na.rm=T))
#######################################################################




###SEGUNDA PARTE: amostragem de pontos de ocorrencia em diferentes camadas de tempo###

##pacotes
library(raster)

## ##definindo variaveis e parametros (LORIEN)
## envVarFolder = "J:/Pesquisadorxs/Anderson_Eduardo/dados_projeto" #pasta com as variaveis ambientais
## projectFolder = "J:/Pesquisadorxs/Anderson_Eduardo/spsArtificiais" #pasta do projeto
## mainSampleFolder = 'J:/Pesquisadorxs/Anderson_Eduardo/spsArtificiais/Amostras' #caminho para pasta onde a anilha com os pontos amostrados sera salva
## AmSulShape = rgdal::readOGR("J:/Pesquisadorxs/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
## crs(AmSulShape) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
## biasLayer = raster('J:/Pesquisadorxs/Anderson_Eduardo/spsArtificiais/biasLayer.grd')
## #biomodFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/biomod/' #pasta para resultados do maxent
## spsTypes = c('spHW', 'spCD') #c('spHW', 'spHD', 'spCD') #nomes das especies
## sampleSizes = c(10, 50, 100) #c(5,10,20,40,80,160) #tamanhos das amostras
## NumRep = 5 #10 #numero de replicas (de cada cenario amostral)
## Tmax = 22 #idade maxima (no passado)
## bgPoints = 1000 #numero de pontos de background


##definindo variaveis e parametros (NOTEBOOK)
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha com os pontos amostrados sera salva
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
crs(AmSulShape) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
biasLayer = raster('/home/anderson/Documentos/Projetos/Sps artificiais/biasLayer.grd')
spsTypes = c('spHW', 'spCD')  #c('spHW', 'spHD', 'spCD') #nomes das especies
sampleSizes = c(10,50,100) #c(5,10,20,40,80,160) #tamanhos das amostras
NumRep = 5 #numero de replicas (de cada cenario amostral)
Tmax = 22 #idade maxima (no passado)
bgPoints = 5000 #numero de pontos de background


##PARA SDM MULTITEMPORAL E SEM VIES AMOSTRAL


sampleData = data.frame()
sampleDataBG = data.frame()

for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
  
  ##criando uma pasta da especie, se nao exisitir
  if(!file.exists(paste(projectFolder,'/Amostras','/multitemporal/',spsTypes[i],sep=''))){
    dir.create(paste(projectFolder,'/Amostras','/multitemporal/',spsTypes[i],sep=''),recursive=TRUE)}
  
  for (sSize in sampleSizes){ #numero de pontos (registros, dados) na amostra
    
    sampledAges = vector()
    sampledAges = round(runif(sSize,0,Tmax)) #selecionando 'n' camadas de tempo aleatoriamente
    nicheRealFolder = paste(projectFolder,'/NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
    nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da
    
    for (j in 1:NumRep){ #replicas do cenario amostral
      
      for (sAge in unique(sampledAges)){ #amostrando em cada camada de tempo que consta na amostra
        
        ## occ pts
        
        sampleData_i = dismo::randomPoints(mask=raster(nicheRealPath[sAge+1])>0.2,prob=TRUE, n=sum(sAge==sampledAges)) #amostra d ponto
        #sampleData_i = randomPoints(mask=raster(nicheRealPath[sAge+1]), n=1) #amostra d ponto
        scenarioName = basename(nicheRealPath[1:24][sAge+1]) #tempo vinculado ao cenario para variaveis ambientais
        scenarioName = gsub('.asc','',scenarioName) #retirando do nome o '.asc'
        layers_i = extract(
          x=stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='asc', full.names=TRUE)),
          y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
        sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
        
        ## background points
        
        envVarPath = list.files(path=envVarFolder,full.names=TRUE)[sAge+1] #lista com os enderecos das variaveis ambientais no tempo corresposndente a interacao
        envData = list.files(envVarPath,full.names=TRUE)
        sampleDataBG_i = dismo::randomPoints(mask = raster(envData[1], 
                                                           crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')), 
                                             n = sum(sAge==sampledAges)*(round(bgPoints/length(sampledAges)))) #amostra dos pontos
        scenarioName = list.files(path=paste(envVarFolder))[sAge+1] #nome do cenario
        layersBG_i = extract(
          x=stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='asc', full.names=TRUE)),
          y=sampleDataBG_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
        sampleDataBG = rbind(sampleDataBG, data.frame(lon=sampleDataBG_i[,1],lat=sampleDataBG_i[,2],layersBG_i,kyrBP=sAge)) #juntando com os dados das outras camadas de tempo amostradas
        
      }
      
      ## occ pts
      names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
      write.csv(sampleData,paste(projectFolder,'/Amostras/multitemporal/',spsTypes[i],'/occ_',sSize,'pts_multitemporal_', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
      sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
      
      ## background pts
      names(sampleDataBG) = c('lon','lat',names(as.data.frame(layersBG_i)),'kyrBP') #ajustando os nomes
      write.csv(sampleDataBG,paste(projectFolder,'/Amostras/multitemporal/',spsTypes[i],'/bg_',sSize,'pts_multitemporal_', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
      sampleDataBG = data.frame() #devolvendo data.frame vazio para proxima rodada
      
    }
  }
}


##PARA SDM MULTITEMPORAL E COM VIES AMOSTRAL


sampleData = data.frame()
sampleDataBG = data.frame()

for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
  
  ##criando uma pasta da especie, se nao exisitir
  if(!file.exists(paste(projectFolder,'/Amostras','/multitemporal/',spsTypes[i],sep=''))){
    dir.create(paste(projectFolder,'/Amostras','/multitemporal/',spsTypes[i],sep=''), recursive=TRUE)}
  
  for (sSize in sampleSizes){ #numero de pontos (registros, dados) na amostra
    
    sampledAges = vector()
    sampledAges = round(runif(sSize,0,Tmax)) #selecionando 'n' camadas de tempo aleatoriamente
    nicheRealFolder = paste(projectFolder,'/NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
    nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da
    biasLayerAdjusted = projectRaster(biasLayer,raster(nicheRealPath[1])) #alinhando o biasLayer com os rasters do projeto
    
    for (j in 1:NumRep){ #replicas do cenario amostral
      
      for (sAge in unique(sampledAges)){ #amostrando em cada camada de tempo que consta na amostra
        
        ## occ pts
        
        sampleData_i = dismo::randomPoints(mask=raster(nicheRealPath[sAge+1])*biasLayerAdjusted, prob=TRUE, n=sum(sAge==sampledAges)) #amostra d ponto
        #sampleData_i = randomPoints(mask=raster(nicheRealPath[sAge+1])*biasLayerAdjusted, n=1, prob=TRUE) #amostra d ponto
        scenarioName = basename(nicheRealPath[1:24][sAge+1]) #tempo vinculado ao cenario para variaveis ambientais
        scenarioName = gsub('.asc','',scenarioName) #retirando do nome o '.asc'
        layers_i = extract(
          x=stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='asc', full.names=TRUE)),
          y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
        sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
        
        ## background points
        
        envVarPath = list.files(path=envVarFolder,full.names=TRUE)[sAge+1] #lista com os enderecos das variaveis ambientais no tempo corresposndente a interacao
        envData = list.files(envVarPath,full.names=TRUE)
        sampleDataBG_i = dismo::randomPoints(mask = raster(envData[1], 
                                                           crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')), 
                                             n = sum(sAge==sampledAges)*(round(bgPoints/length(sampledAges)))) #amostra dos pontos
        scenarioName = list.files(path=paste(envVarFolder))[sAge+1] #nome do cenario
        layersBG_i = extract(
          x=stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='asc', full.names=TRUE)),
          y=sampleDataBG_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
        sampleDataBG = rbind(sampleDataBG, data.frame(lon=sampleDataBG_i[,1],lat=sampleDataBG_i[,2],layersBG_i,kyrBP=sAge)) #juntando com os dados das outras camadas de tempo amostradas
        
        
      }
      
      ## occ pts
      names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
      write.csv(sampleData,paste(projectFolder,'/Amostras/multitemporal/',spsTypes[i],'/occ_',sSize,'pts_multitemporal_comVIES_', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
      sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
      
      ## background pts
      names(sampleDataBG) = c('lon','lat',names(as.data.frame(layersBG_i)),'kyrBP') #ajustando os nomes
      write.csv(sampleDataBG,paste(projectFolder,'/Amostras/multitemporal/',spsTypes[i],'/bg_',sSize,'pts_multitemporal_', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
      sampleDataBG = data.frame() #devolvendo data.frame vazio para proxima rodada
      
    }
  }
}


##PARA SDM MONOTEMPORAL (pode ser presente ou passado) E SEM VIES AMOSTRAL (obviamente, pq esse vies e para fosseis)


sampleData = data.frame()
sampleDataBg = data.frame()

for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
  
  ##criando uma pasta da especie, se nao exisitir
  if(!file.exists(paste(projectFolder,'/Amostras','/monotemporal/',spsTypes[i],sep=''))){
    dir.create(paste(projectFolder,'Amostras','/monotemporal/',spsTypes[i],sep=''), recursive=TRUE)}
  
  for (sSize in sampleSizes){ #numero de pontos (registros, dados) na amostra
    
    for (j in 1:NumRep){ #replicas do cenario amostral
      
      sampledAge = round(runif(1,0,Tmax)) #selecionando a camada de tempo aleatoriamente
      nicheRealFolder = paste(projectFolder,'/NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
      nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da
      
      sampleData_i = dismo::randomPoints(mask=raster(nicheRealPath[sampledAge])>0.2, prob=TRUE, n=sSize) #amostra do ponto
      scenarioName = basename(nicheRealPath[1:24][sampledAge]) #tempo vinculado ao cenario para variaveis ambientais
      scenarioName = gsub('.asc','',scenarioName) #retirando do nome o '.asc'
      layers_i = extract(
        x=stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='asc', full.names=TRUE)),
        y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
      sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sampledAge)) #juntando com os dados das outras camadas de tempo amostradas
      names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
      write.csv(sampleData,paste(projectFolder,'/Amostras/monotemporal/',spsTypes[i],'/occ_',sSize,'pts_monotemporal_',j,'rep','.csv',sep=''),row.names=FALSE) #salvando
      
      sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
      
      ##background points##
      sampleDataBg_i = dismo::randomPoints(mask=raster(nicheRealPath[sampledAge]),
                                           n=bgPoints) #amostra dos pontos
      #            scenarioName = list.files(path=paste(envVarFolder))[sAge+1] #nome do cenario
      layersBg_i = extract(
        x=stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='asc', full.names=TRUE)),
        y=sampleDataBg_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
      sampleDataBg = rbind(sampleDataBg, data.frame(lon=sampleDataBg_i[,1],lat=sampleDataBg_i[,2],layersBg_i,kyrBP=sampledAge)) #juntando com os dados das outras camadas de tempo amostradas
      names(sampleDataBg) = c('lon','lat',names(as.data.frame(layersBg_i)),'kyrBP') #ajustando os nomes
      write.csv(sampleDataBg,paste(projectFolder,'/Amostras/monotemporal/',spsTypes[i],'/bg_',sSize,'pts_monotemporal_',j,'rep','.csv',sep=''),row.names=FALSE) #salvando
      sampleDataBg = data.frame() #devolvendo data.frame vazio para proxima rodada
      
    }
  }
}


###TERCEIRA PARTE: SDM usando de pontos de ocorrencia em diferentes camadas de tempo (do atual ate 120 kyr BP)###


#######################################################
####################### MAXENT ########################
#######################################################

##pacotes
library(biomod2)

## ##definindo variaveis e parametros (LORIEN)
## options(java.parameters = "-Xmx7g") ###set available memmory to java
## projectFolder =  "J:/Anderson_Eduardo/spsArtificiais" #pasta do projeto
## envVarFolder = "J:/Anderson_Eduardo/dados_projeto" #pasta com as variaveis ambientais
## envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
## AmSulShape = rgdal::readOGR("J:/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
## mainSampleFolder = "J:/Anderson_Eduardo/spsArtificiais/Amostras" #caminho para pasta onde a planilha
## maxentFolder = 'C:/Users/WS/Documents/R/win-library/3.4/dismo/java' #pasta para resultados do maxent
## spsTypes = c('spHW','spCD') #c('spHW', 'spHD', 'spCD') #nomes das especies
## sdmTypes = c('monotemporal') #c('multitemporal','monotemporal')
## #source("/home/anderson/R/R-Scripts/TSSmaxent.R")
## sampleSizes = 50 #c(10,100) #c(5,10,20,40,80,160) #tamanhos das amostras
## NumRep = 3 #10 #numero de replicas (de cada cenario amostral)
## #statResults = data.frame() #tabela de estatisticas basicas do modelo

##definindo variaveis e parametros (NOTEBOOK)
options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais" #pasta do projetoenvVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras' #caminho para pasta onde a planilha com os pontos
maxentFolder = '/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java' #pasta para resultados do maxent
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
spsTypes = c('spHW','spCD') #nomes das especies  #c('spHW', 'spHD', 'spCD') #nomes das especies
sdmTypes = c('multitemporal','monotemporal')
#source("/home/anderson/R/R-Scripts/TSSmaxent.R")
sampleSizes = c(10,50,100) #c(5,10,20,40,80,160) #tamanhos das amostras
NumRep = 5 #numero de replicas (de cada cenario amostral)
timeStart = Sys.time()
#statResults = data.frame() #tabela de estatisticas basicas do modelo


##algoritmo da analise do projeto
for (h in 1:length(sdmTypes)){
  for (i in 1:length(spsTypes)){
    
    statResults = data.frame() #tabela de estatisticas basicas do modelo  
    
    for (j in 1:length(sampleSizes)){
      for (k in 1:NumRep){ #loop sobre o numero de replicas 
        tryCatch({
          
          ##ajustando o diretorio de trabalho (pois o biomod roda e salva tudo simultaneamente)
          if(!file.exists(file.path(projectFolder,'maxent',sdmTypes[h], spsTypes[i],sep=''))){
            dir.create(file.path(projectFolder,'maxent',sdmTypes[h],spsTypes[i],sep=''),recursive=TRUE)
          }
          setwd(file.path(projectFolder,'maxent',sdmTypes[h],spsTypes[i]))
          
          ##definindo variaveis e parametros locais                        
          occPoints = read.csv(paste(mainSampleFolder,sdmTypes[h],'/',spsTypes[i],'/occ_',sampleSizes[j],'pts_',sdmTypes[h],'_',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia    
          backgroundPoints = read.csv(paste(mainSampleFolder,sdmTypes[h],'/',spsTypes[i],'/bg_',sampleSizes[j],'pts_',sdmTypes[h],'_',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de background
          
          
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
              path_to_maxent.jar=maxentFolder,
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
            modelType = sdmTypes[h],
            sp = spsTypes[i],
            sampleSize = sampleSizes[j],
            replicate = k,
            AUC = mean(evaluationScores['ROC','Testing.data',,,]),
            TSS = mean(evaluationScores['TSS','Testing.data',,,]),
            numbOfTimeLayers = length(unique(occPoints$kyrBP)),
            medianKyr = median(occPoints$kyrBP),
            minAge = min(occPoints$kyrBP),
            maxAge = max(occPoints$kyrBP)))
          
          write.csv(statResults,file=paste(projectFolder,'/maxent/',sdmTypes[h],'/',spsTypes[i],'/StatisticalResults-',spsTypes[i],'.csv',sep=''),row.names=FALSE)
          
          ##implementando projecoes do modelo
          for (l in 1:length(envVarPaths[1:24])){
            
            ##definindo variaveis e parametros internos
            predictors = stack(list.files(path=envVarPaths[l],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
            predictors = predictors[[c('bioclim_01','bioclim_12')]]
            ##predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
            crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
            
            ##selecionando o melhor modelo para projecao
            whichModel = names(evaluationScores['TSS','Testing.data',,,][which(evaluationScores['TSS','Testing.data',,,]== max(evaluationScores['TSS','Testing.data',,,]) )])
            modelName = grep(pattern=whichModel, myBiomodModelOut@models.computed, value=TRUE)
            
            ##rodando algortmo de projecao (i.e. rodando a projecao)
            myBiomodProj <- BIOMOD_Projection(
              modeling.output = myBiomodModelOut,
              new.env = predictors,
              proj.name = paste(l-1,'kyr',sep=''),
              selected.models = modelName,
              binary.meth = 'TSS',
              compress = 'TRUE',
              build.clamping.mask = 'TRUE',
              output.format = '.grd')
            
            ##gerando e salvando um mapa binario (threshold 10%)
            ## projStack = get_predictions(myBiomodProj) #extrai as projecoes
            ## projStackBIN = BinaryTransformation(stack(mean(projStack)),'10')
            ## writeRaster(projStackBIN,file=paste(projectFolder,'/maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[j],'.replica',k,'/proj_',l,'kyr/proj_',i,'kyr','.sample',sampleSizes[j],'.replica',k,'_BIN.asc',sep=''),overwrite=TRUE)
            
          }
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
  }
}

##tempo gasto
print(Sys.time() - timeStart)


###QUARTA PARTE: comparando projecao do SDM e a distribuicao espacial real do nicho da sp###


##abrindo pacotes necessarios
library(raster)
library(ecospat)

## ##definindo variaveis e parametros (LORIEN)
## options(java.parameters = "-Xmx7g") ###set available memmory to java
## projectFolder =  "J:/Anderson_Eduardo/spsArtificiais" #pasta do projeto
## envVarFolder = "J:/Anderson_Eduardo/dados_projeto" #pasta com as variaveis ambientais
## envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
## AmSulShape = rgdal::readOGR("J:/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
## mainSampleFolder = "J:/Anderson_Eduardo/spsArtificiais/Amostras" #caminho para pasta onde a planilha
## maxentFolder = 'C:/Users/WS/Documents/R/win-library/3.4/dismo/java' #pasta para resultados do maxent
## spsTypes = c('spHW','spCD') #c('spHW', 'spHD', 'spCD') #nomes das especies
## sdmTypes = c('multitemporal','monotemporal')
## #source("/home/anderson/R/R-Scripts/TSSmaxent.R")
## sampleSizes = c(10,100) #c(5,10,20,40,80,160) #tamanhos das amostras
## NumRep = 3 #10 #numero de replicas (de cada cenario amostral)
## outputData = data.frame()

##definindo variaveis e parametros (NOTEBOOK)
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha com os pontos amostrados sera salva
spsTypes = c('spHW','spCD') #c('spHW', 'spHD', 'spCD') #nomes das especies
sdmTypes = c('multitemporal','monotemporal')
sampleSizes = c(10,50,100) #c(5,10,20,40,80,160) #aqui, deve ser igual ao usasado nas partes anteriores do script
NumRep = 5
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
timeStart = Sys.time()
outputData = data.frame()


##algoritmo da analise do projeto
for (h in 1:length(sdmTypes)){
  for (i in 1:length(spsTypes)){
    
    ##definindo variaveis e parametros locais
    nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
    nicheRealPath = list.files(path=nicheRealFolder,pattern='.asc',full.names=TRUE) #lista com os enderecos dos mapas de distribuicao da sp
    
    ##loop sobre as cadamdas de tempo
    for (l in 1:length(nicheRealPath[1:24])){   #1:length(nicheRealPath[1:24])){ 
      
      ##definindo variaveis e parametros locais
      realNiche = nicheRealPath[l] #nicho real
      
      ##amostrando pontos da distribuicao real para compracao dos SDMs
      binMap = raster(realNiche)>0.2 #mapa binario do real
      realNicheDataOccCoord = dismo::randomPoints(binMap,1000) #amostrando 1000 pontos do binario real
      realNicheDataOccPres = extract(binMap,realNicheDataOccCoord,na.rm=TRUE) #definindo occ e ausencias para os pontos
      realNicheDataOcc = data.frame(longitude=realNicheDataOccCoord[,1], latitude=realNicheDataOccCoord[,2], pres=realNicheDataOccPres) #tabela lon, lat e pres
      predictors = stack(list.files(path=envVarPaths[l],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis
      predictors = predictors[[c('bioclim_01','bioclim_12')]] #selecionando as variaveis usadas
      predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
      projection(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
      realNicheDataPred = extract(x=predictors,y=realNicheDataOcc[,c('longitude','latitude')],na.rm=TRUE) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
      realNicheData = data.frame(realNicheDataOcc, realNicheDataPred) #juntando com os dados das outras camadas de tempo amostradas
      
      ##loop sobre os tamanhos amostrais
      for (m in sampleSizes){ 
        
        ## ##timeSampleData = list.files(path=projectionsFolder, pattern=glob2rx(paste('*Time',l-1,'*Sample',m,'.asc',sep='')),full.names=TRUE)
        ## timeSampleData = list.files(path=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],sep=''), pattern=glob2rx(paste('*',l,'kyr','*sample',m,'*.asc',sep='')), recursive=TRUE, full.names=TRUE) #pasta com as projecoes do cenario
        
        ##loop sobre replicas de cada combinacao de tempo e tamanho amostral
        for(n in 1:NumRep){ 
          tryCatch({
            
            ##definindo variaveis e parametros locais
            ##outputData = data.frame() #tabela de dados de saida
            ## projectionsFolder = paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',m,'.replica',n,sep='') #pasta com as projecoes do cenario
            ## projectionsPath = list.files(path=projectionsFolder, pattern='.asc',recursive=TRUE,full.names=T) #caminhos para os .asc na paste do cenario
            
            sdmNichePath = paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',m,'.replica',n,'/proj_',l-1,'kyr/','proj_',l-1,'kyr_',spsTypes[i],'.sample',m,'.replica',n,'_TSSbin.grd',sep='') #caminho do mapa de suitability gerado por SDM
            sdmNicheStack = stack(sdmNichePath)
            binMapSDM = sdmNicheStack #sum(sdmNicheStack)#>0.5
            #                    binMapSDM = sdmNicheStack[[round(runif(1,0.5,nlayers(sdmNicheStack)))]] #mapa de suitability gerado por SDM
            #                    binMapSDM = biomod2::BinaryTransformation(data=mean(sdmNiche), threshold=10) #fazendo mapa binario, threshold 10%
            
            SDMDataOccCoord = dismo::randomPoints(binMapSDM, 1000)
            SDMDataOccPres = extract(binMapSDM, SDMDataOccCoord, na.rm=TRUE)
            SDMDataOcc = data.frame(longitude=SDMDataOccCoord[,1],latitude=SDMDataOccCoord[,2],pres=as.numeric(SDMDataOccPres))
            SDMDataPred = extract(x=predictors,y=SDMDataOcc[,c('longitude','latitude')],na.rm=TRUE) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
            SDMData = data.frame(SDMDataOcc, SDMDataPred) #juntando com os dados das outras camadas de tempo amostradas
            SDMData = SDMData[complete.cases(SDMData),]
            
            ##The PCA is calibrated on all the sites of the study area
            pca.env <- dudi.pca(rbind(realNicheData,SDMData)[,c('bioclim_01','bioclim_12')],scannf=F,nf=2)
            #ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) #grafico
            
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
            ##OBS: Compares the observed niche overlap between z1 and z2 to overlaps between random niches z1.sim
            ## and z2.sim, which are built from random reallocations of occurences of z1 and z2.
            ##'alternative' argument specifies if you want to test for niche conservatism (alternative = "greater", i.e.  the
            ## niche overlap is more equivalent/similar than random) or for niche divergence (alternative = "lower",
            ## i.e. the niche overlap is less equivalent/similar than random).
            eq.test <- ecospat.niche.equivalency.test(grid.clim.realNiche, grid.clim.SDMniche,rep=100, alternative = "greater")
            
            ##similaridade de nicho
            ##OBS: Compares the observed niche overlap between z1 and z2 to overlaps between z1 and random niches
            ## (z2.sim) as available in the range of z2 (z2$Z). z2.sim has the same pattern as z2 but the center is
            ## randomly translatated in the availabe z2$Z space and weighted by z2$Z densities. If rand.type = 1,
            ## both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly shifted.
            ## 'alternative' specifies if you want to test for niche conservatism (alternative = "greater", i.e.  the
            ## niche overlap is more equivalent/similar than random) or for niche divergence (alternative = "lower",
            ## i.e. the niche overlap is less equivalent/similar than random)
            sim.test <- ecospat.niche.similarity.test(grid.clim.realNiche, grid.clim.SDMniche, rep=100, alternative = "greater")
            
            Dobs_equiv = eq.test$obs$D #indice D observado no teste de equivalencia de nicho
            Iobs_equiv = eq.test$obs$I #indice I observado no teste de equivalencia de nicho
            DpValue_equiv = eq.test$p.D #p-valor indice D no teste de equivalencia de nicho
            IpValue_equiv = eq.test$p.I #p-valor indice I no teste de equivalencia de nicho
            ##
            Dobs_simi = sim.test$obs$D #indice D observado no teste de similaridade de nicho
            Iobs_simi = sim.test$obs$I #indice I observado no teste de similaridade de nicho
            DpValue_simi = sim.test$p.D #p-valor indice D no teste de similaridade de nicho
            IpValue_simi = sim.test$p.I #p-valor indice I no teste de similaridade de nicho
            
            ##abrindo planilha de pontos para extrair dados do cenario
            occPoints = read.csv(paste(mainSampleFolder,'/',sdmTypes[h],'/',spsTypes[i],'/occ_',m,'pts_',sdmTypes[h],'_', n ,'rep.csv',sep=''),header=TRUE) 
            occPoints[occPoints==0] = NA
            occPoints = occPoints[complete.cases(occPoints),]
            occPoints = round(occPoints, digits=2)
            occPoints = occPoints[!duplicated(occPoints),]                 
            
            outputData = rbind(outputData,data.frame(sdmType = sdmTypes[h],
                                                     sp = spsTypes[i],
                                                     kyrBP = l-1,
                                                     sampleSize = m,
                                                     replicate = n,
                                                     numbOfTimeLayers = length(unique(occPoints$kyrBP)),
                                                     medianKyr = median(occPoints$kyrBP),
                                                     minAge = min(occPoints$kyrBP),
                                                     maxAge = max(occPoints$kyrBP),
                                                     Schoeners_D_equiv = Dobs_equiv,
                                                     p_value_equiv = DpValue_equiv,
                                                     Hellinger_I_equiv = Iobs_equiv,
                                                     p_value_equiv = IpValue_equiv,
                                                     Schoeners_D_simi = Dobs_simi,
                                                     p_value_simi = DpValue_simi,
                                                     Hellinger_I_simi = Iobs_simi,
                                                     p_value_simi = IpValue_simi))
            
            write.csv(outputData, file=paste(projectFolder,'/maxent/output.csv',sep=''),row.names=FALSE) #salvando os dados do cenario
            
          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }
      }
    }
  }
}

##registro do tempo
print(Sys.time() - timeStart)



### QUINTA PARTE: construindo graficos dos resultados ###



## ##definindo parametros e variaveis (LORIEN)
## spsTypes = c('spHW', 'spCD') #c('spHW', 'spHD', 'spCD') #nomes das especies
## outputData = list() #tabela de dados de saida
## vetor.nomes = vector()
## projectFolder = "J:/Anderson_Eduardo/spsArtificiais" #pasta do projeto

##definindo parametros e variaveis (NOTEBOOK)
spsTypes = c('spHW', 'spCD') #c('spHW', 'spHD', 'spCD') #nomes das especies
outputData = list() #tabela de dados de saida
vetor.nomes = vector()
projectFolder = "/home/anderson/Projetos/Sps artificiais" #pasta do projeto
#projectFolder = '/media/anderson/PIBi/ANDERSON EDUARDO/Sps artificiais'

### AUC e TSS dos modelos

##multitemporal, spHW

spHWmulti = read.csv(paste(projectFolder,'/maxent/multitemporal/spHW/StatisticalResults-spHW.csv', sep=''), header=TRUE)

##multitemporal, spCD
spCDmulti = read.csv(paste(projectFolder,'/maxent/multitemporal/spCD/StatisticalResults-spCD.csv', sep=''), header=TRUE)

##monotemporal, spHW
spHWmono = read.csv(paste(projectFolder,'/maxent/monotemporal/spHW/StatisticalResults-spHW.csv', sep=''), header=TRUE)

##monotemporal, spCD
spCDmono = read.csv(paste(projectFolder,'/maxent/monotemporal/spCD/StatisticalResults-spCD.csv', sep=''), header=TRUE)


## boxplots modelo X AUC e TSS, dados totais

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotModelos&Acuracia_dadosTotais.jpeg', width=800)
par(mfrow=c(1,2), las=2, mar=c(8,5,1,1), cex=1.3)
boxplot(rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$AUC ~  rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$modelType, ylim=c(0,1), ylab='AUC')
boxplot(rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$TSS ~  rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$modelType, ylim=c(0,1), ylab='TSS')
dev.off()

##testes (ambos deram nao significativos)

kruskal.test( rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$AUC ~  rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$modelType )
kruskal.test( rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$TSS ~  rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$modelType )



## boxplots modelo X AUC e TSS, especies

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotModelos&Acuracia_sps.jpeg', height=600)
par(mfrow=c(2,2), las=2, mar=c(8,5,2,1), cex=1.1)
boxplot(rbind(spHWmulti,spHWmono)$AUC ~  rbind(spHWmulti,spHWmono)$modelType, ylim=c(0,1), ylab=c('AUC'), main='HW species')
boxplot(rbind(spCDmulti,spCDmono)$AUC ~  rbind(spCDmulti,spCDmono)$modelType, ylim=c(0,1), ylab=c('AUC'), main='CD species')
boxplot(rbind(spHWmulti,spHWmono)$TSS ~  rbind(spHWmulti,spHWmono)$modelType, ylim=c(0,1), ylab=c('TSS'), main='HW species')
boxplot(rbind(spCDmulti,spCDmono)$TSS ~  rbind(spCDmulti,spCDmono)$modelType, ylim=c(0,1), ylab=c('TSS'), main='CD species')
dev.off()

##testes (todos deram nao significativos)

kruskal.test(rbind(spHWmulti,spHWmono)$AUC ~  rbind(spHWmulti,spHWmono)$modelType)
kruskal.test(rbind(spCDmulti,spCDmono)$AUC ~  rbind(spCDmulti,spCDmono)$modelType)
kruskal.test(rbind(spHWmulti,spHWmono)$TSS ~  rbind(spHWmulti,spHWmono)$modelType)
kruskal.test(rbind(spCDmulti,spCDmono)$TSS ~  rbind(spCDmulti,spCDmono)$modelType)


## boxplots sp X AUC e TSS, dados totais

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotSps&Acuracia.jpeg')
par(mfrow=c(2,2), mar=c(3,4,5,1),cex=1.1)
boxplot(rbind(spHWmulti,spCDmulti)$AUC ~  rbind(spHWmulti,spCDmulti)$sp, ylim=c(0,1), ylab=c('AUC'), main='Multitemporal')
boxplot(rbind(spHWmono,spCDmono)$AUC ~  rbind(spHWmono,spCDmono)$sp, ylim=c(0,1), ylab=c('AUC'), main='Monotemporal')
boxplot(rbind(spHWmulti,spCDmulti)$TSS ~  rbind(spHWmulti,spCDmulti)$sp, ylim=c(0,1), ylab=c('TSS'), main='Multitemporal')
boxplot(rbind(spHWmono,spCDmono)$TSS ~  rbind(spHWmono,spCDmono)$sp, ylim=c(0,1), ylab=c('TSS'), main='Monotemporal')
dev.off()

##teste (tanto para AUC quanto para TSS, deu diferenca significativa apenas para SDMmulti, e nao para SDMmono)

kruskal.test(rbind(spHWmulti,spCDmulti)$AUC ~  rbind(spHWmulti,spCDmulti)$sp) 
kruskal.test(rbind(spHWmono,spCDmono)$AUC ~  rbind(spHWmono,spCDmono)$sp)
kruskal.test(rbind(spHWmulti,spCDmulti)$TSS ~  rbind(spHWmulti,spCDmulti)$sp)
kruskal.test(rbind(spHWmono,spCDmono)$TSS ~  rbind(spHWmono,spCDmono)$sp)


## boxplots sampleSize X AUC e TSS, dados totais

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotSampleSize&Acuracia_dadosTotais.jpeg')
par(mfrow=c(2,2))
boxplot(rbind(spHWmulti,spCDmulti)$AUC ~  rbind(spHWmulti,spCDmulti)$sampleSize, ylim=c(0,1), ylab='AUC', main='Multitemporal')
boxplot(rbind(spHWmono,spCDmono)$AUC ~  rbind(spHWmono,spCDmono)$sampleSize, ylim=c(0,1), ylab='AUC', main='Monotemporal')
boxplot(rbind(spHWmulti,spCDmulti)$TSS ~  rbind(spHWmulti,spCDmulti)$sampleSize, ylim=c(0,1), ylab='TSS', main='Multitemporal')
boxplot(rbind(spHWmono,spCDmono)$TSS ~  rbind(spHWmono,spCDmono)$sampleSize, ylim=c(0,1), ylab='TSS', main='Monotemporal')
dev.off()

## boxplots sampleSize X AUC, especies

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotSampleSize&Acuracia_spHW.jpeg')
par(mfrow=c(2,2), mar=c(3,4,5,1))
boxplot(rbind(spHWmulti)$AUC ~  rbind(spHWmulti)$sampleSize, ylim=c(0,1), ylab='AUC', main='Multitemporal')
boxplot(rbind(spHWmono)$AUC ~  rbind(spHWmono)$sampleSize, ylim=c(0,1), ylab='AUC', main='Monotemporal')
boxplot(rbind(spHWmulti)$TSS ~  rbind(spHWmulti)$sampleSize, ylim=c(0,1), ylab='TSS', main='Multitemporal')
boxplot(rbind(spHWmono)$TSS ~  rbind(spHWmono)$sampleSize, ylim=c(0,1), ylab='TSS', main='Monotemporal')
title('spHW', outer=TRUE, line=-1)
dev.off()

## boxplots modelo X AUC, especies

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotSampleSize&Acuracia_spCD.jpeg')
par(mfrow=c(2,2), mar=c(3,4,5,1))
boxplot(rbind(spCDmulti)$AUC ~  rbind(spCDmulti)$sampleSize, ylim=c(0,1), ylab='AUC', main='Multitemporal')
boxplot(rbind(spCDmono)$AUC ~  rbind(spCDmono)$sampleSize, ylim=c(0,1), ylab='AUC', main='Monotemporal')
boxplot(rbind(spCDmulti)$TSS ~  rbind(spCDmulti)$sampleSize, ylim=c(0,1), ylab='TSS', main='Multitemporal')
boxplot(rbind(spCDmono)$TSS ~  rbind(spCDmono)$sampleSize, ylim=c(0,1), ylab='TSS', main='Monotemporal')
title('spCD', outer=TRUE, line=-1)
dev.off()


### Sobreposicao de nicho

outputData = read.csv(file=paste(projectFolder,'/maxent/output.csv',sep=''), header=TRUE)
##outputData = read.csv(file=paste(projectFolder,'/Resultados Lorien/output.csv',sep=''),header=TRUE)
##outputData = read.csv(file=paste(projectFolder,'/maxent/output_uni.csv',sep=''), header=TRUE)
#vetor.nomes = append(vetor.nomes,paste(spsTypes[i],sep=''))


## Schoener e Hellinger para resultados totais
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotDadosTotais.jpeg', width=600)
par(mfrow=c(1,2), mar=c(8,3,3,1), cex=1.4, las=2)
boxplot(outputData$Schoeners_D_simi~ outputData$sdmType, ylim=c(0,1), main="Schoeners' D")
boxplot(outputData$Hellinger_I_simi~ outputData$sdmType, ylim=c(0,1), main='Hellinger')
dev.off()

##testes (nao houve diferencas significativas pra nenhum)
kruskal.test(Schoeners_D_simi ~ sdmType, data = outputData)
kruskal.test(Hellinger_I_simi ~ sdmType, data = outputData)


## bosplots das sps
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotSps.jpeg', height=650)
par(mfrow=c(2,2), mar=c(7,4.5,6,1), cex=1.1, las=2)# cex.axis=2.5, cex.lab=3, cex.main=3)
boxplot(outputData[outputData$sp == 'spHW',]$Schoeners_D_simi ~ outputData[outputData$sp == 'spHW',]$sdmType, ylim=c(0,1), ylab="Schoener's D", main='HW species')
boxplot(outputData[outputData$sp == 'spHW',]$Hellinger_I_simi ~ outputData[outputData$sp == 'spHW',]$sdmType, ylim=c(0,1), ylab="Hellinger", main='HW species')
boxplot(outputData[outputData$sp == 'spCD',]$Schoeners_D_simi ~ outputData[outputData$sp == 'spCD',]$sdmType, ylim=c(0,1), ylab="Schoener's D",  main='CD species')
boxplot(outputData[outputData$sp == 'spCD',]$Hellinger_I_simi ~ outputData[outputData$sp == 'spCD',]$sdmType, ylim=c(0,1), ylab="Hellinger", main='CD species')
dev.off()

##testes (nenhuma diferenca significativa)
kruskal.test(outputData[outputData$sp == 'spHW',]$Schoeners_D_simi ~ outputData[outputData$sp == 'spHW',]$sdmType)
kruskal.test(outputData[outputData$sp == 'spHW',]$Hellinger_I_simi ~ outputData[outputData$sp == 'spHW',]$sdmType)
kruskal.test(outputData[outputData$sp == 'spCD',]$Schoeners_D_simi ~ outputData[outputData$sp == 'spCD',]$sdmType)
kruskal.test(outputData[outputData$sp == 'spCD',]$Hellinger_I_simi ~ outputData[outputData$sp == 'spCD',]$sdmType)

## Densidade para dados totais
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/densidadeDadosTotais.jpeg', width=600, height = 400)
par(mfrow=c(1,2), lwd=2, cex=1)
plot(density(outputData[outputData$sdmType == 'multitemporal',]$Schoeners_D_simi),ylim=c(0,5), lwd=2, col='red', main='', xlab="Schoener's D", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal',]$Schoeners_D_simi), lwd=2)
#
plot(density(outputData[outputData$sdmType == 'multitemporal',]$Hellinger_I_simi),ylim=c(0,5), lwd=2, col='red', main='', xlab='Hellinger', ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal',]$Hellinger_I_simi), lwd=2)
##
legend(x='topright', legend=c('Multitemporal calibration', 'Monotemporal calibration'), lty=1, col=c('red','black'), bty='n')
dev.off()

## Densidade para as sps
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/densidade_sps.jpeg')
par(mfrow=c(2,2), mar=c(5,4,3,1), lwd=2, cex=1)
#
plot(density(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi),ylim=c(0,7), lwd=2, col='red', main='HW species', xlab="Schoener's D", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi), lwd=2)
#
plot(density(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi),ylim=c(0,5), lwd=2, col='red', main='HW species', xlab="Hellinger", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi), lwd=2)
legend(x='topright', legend=c('Multitemporal calibration', 'Monotemporal calibration'), lty=1, col=c('red','black'), bty='n', cex=0.8)
#
plot(density(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi),ylim=c(0,5), lwd=2, col='red', main='CD species', xlab="Schoener's D", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi), lwd=2)
#
plot(density(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi),ylim=c(0,5), col='red', lwd=2, main='CD species', xlab="Hellinger", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi), lwd=2)
dev.off()

## Schoener e Hellinger no tempo
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Shoener&HellingerXtempo.jpeg',width=600, height=600)
par(mfrow=c(2,2), mar=c(4,4,4,1), cex=1.2)
plot(outputData[outputData$sdmType == 'multitemporal',]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal',]$kyrBP),type='p',ylab="Schoeners' D", xlab="Time (kyr BP)", ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), main='Multitemporal')
#
plot(outputData[outputData$sdmType == 'multitemporal',]$Hellinger_I_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal',]$kyrBP),type='p',ylab="Hellinger",xlab="Time (kyr BP)",ylim=c(0,1),col=rgb(0,0,0,alpha=0.5), main='Multitemporal')
#
plot(outputData[outputData$sdmType == 'monotemporal',]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal',]$kyrBP),type='p',ylab="Schoeners' D",xlab="Time (kyr BP)",ylim=c(0,1),col=rgb(0,0,0,alpha=0.5), main='Monotemporal')
#
plot(outputData[outputData$sdmType == 'monotemporal',]$Hellinger_I_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal',]$kyrBP),type='p',ylab="Hellinger",xlab="Time (kyr BP)",ylim=c(0,1),col=rgb(0,0,0,alpha=0.5), main='Monotemporal')
dev.off()

## Shoener e Hellinger no tempo - sps
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Shoener&HellingerXtempo_sps.jpeg', width=1200, height=1200)
par(mfrow=c(2,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='HW species', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$kyrBP),type='p',ylab="Hellinger",xlab="Time (kyr BP)", main='HW species',ylim=c(0,1), col=rgb(0,0,0,alpha=0.5))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$kyrBP),type='p',ylab="Schoeners' D",xlab="Time (kyr BP)", main='CD species', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$kyrBP),type='p',ylab="Hellinger",xlab="Time (kyr BP)", main='CD species',ylim=c(0,1), col=rgb(0,0,0,alpha=0.5))
dev.off()

## Schoener X sample size X tempo - spHW
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/SchoenerXtempoXsample_spHW.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
dev.off()

## Schoener X sample size X tempo - spCD
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/SchoenerXtempoXsample_spCD.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
dev.off()


## Hellinger X sample size X tempo - spHW
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/HellingerXtempoXsample_spHW.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
dev.off()

## Hellinger X sample size X tempo - spCD
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/HellingerXtempoXsample_spCD.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
dev.off()



## Tamanho amostral - dados totais
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplot_sampleSize_dadosTotais.jpeg', width=800, height=900)
par(mfrow=c(2,2), cex=1.5)
boxplot(outputData[outputData$sdmType == 'multitemporal',]$Schoeners_D_simi~ outputData[outputData$sdmType == 'multitemporal',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoeners' D", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'monotemporal',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoeners' D", main='Monotemporal')
boxplot(outputData[outputData$sdmType == 'multitemporal',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'multitemporal',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'monotemporal',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Monotemporal')
dev.off()

## Tamanho amostral - spHW
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplot_sampleSize_spHW.jpeg', width=800, height=900)
par(mfrow=c(2,2), cex=1.5)
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoener's D", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoener's D", main='Monotemporal')
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Monotemporal')
dev.off()

## Tamanho amostral - spCD
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplot_sampleSize_spCD.jpeg', width=800, height=900)
par(mfrow=c(2,2), cex=1.5)
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoeners' D", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoener's D", main='Monotemporal')
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Monotemporal')
dev.off()

## Shoener's D e Hellinger X número de camadas no SDMmultitemporal
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplot_NumberOfTimeLayers.jpeg', height=1000, width=600)
par(mfrow=c(3,2), cex=1.3)
boxplot(outputData[outputData$sdmType == 'multitemporal',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'multitemporal',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Schoener's D", main='Full dataset')
##
boxplot(outputData[outputData$sdmType == 'multitemporal',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'multitemporal',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Hellinger", main='Full dataset')
##
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Schoener's D", main='HW species')
##
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Hellinger", main='HW species')
##
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Schoener's D", main='CD species')
##
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Hellinger", main='CD species')
dev.off()

##graficos para clamping

projectFolder = "/home/anderson/Projetos/Sps artificiais/"
sdmTypes = c("multitemporal", "monotemporal")
spsTypes = c("spHW", "spCD")
sampleSizes = c(10, 50, 100)
numRep = 5
clampList = list()
territory = list()


for(h in 1:length(sdmTypes)){
  for(i in 1:length(spsTypes)){
    for(m in 1:length(sampleSizes)){
      for(n in 1:numRep){
        for(l in 1:24){
          ##mapa de clamping
          sdmClampPath = paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[m],'.replica',n,'/proj_',l-1,'kyr/','proj_',l-1,'kyr_ClampingMask.grd',sep='') #caminho do mapa de suitability gerado por SDM
          clampLayer_i = raster(sdmClampPath)
          scenName = paste(sdmTypes[h],'_proj_',l-1,'kyr_',spsTypes[i],'.sample',sampleSizes[m],'.replica',n,sep='')
          clampList[[scenName]] = clampLayer_i
          clamping = (sum(getValues(clampLayer_i)>0, na.rm=TRUE)/ncell(getValues(clampLayer_i))) * 100
          ##mapa de distribuicao da especie
          sdmDistPath = paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[m],'.replica',n,'/proj_',l-1,'kyr/','proj_',l-1,'kyr_',spsTypes[i],'.sample',sampleSizes[m],'.replica',n,'_TSSbin.grd',sep='') #caminho do mapa de suitability gerado por SDM
          distLayer_i = raster(sdmDistPath)
          ##calculo da porporcao de clamping na area de distribuicao modelada
          distUnderClamp = (clampLayer_i + distLayer_i)==2
          distUnderClamp = ( freq(distUnderClamp, value=1)/sum(freq(distUnderClamp)[1:2,2]) ) * 100
          ##tabela de dados final
          outputDF = outputData[ which(outputData$sdmType==sdmTypes[h] & outputData$sp==spsTypes[i] & outputData$sampleSize==sampleSizes[m] & outputData$replicate==n & outputData$kyrBP==l-1 ), ]
          if(nrow(outputDF)>0){
            territory[[scenName]] = data.frame( outputDF,
                                                clamping = clamping,
                                                distUnderClamp = distUnderClamp )
          }
        }
      }
    }
  }
}




##tranformando a lista em stack de gridfiles
clampStack = stack(clampList)


##multitemporal

##todos os mapas de clamping - multitemporal, spHW, sample 10 pts (todos as camadas temporais)
scenNames =  grep(pattern='^multitemporal.*spHW.*sample10.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
scenNames =  grep(pattern='^multitemporal.*spHW.*sample100.*replica1', x=scenNames, value=TRUE, invert=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Multitemporal_SpHW_sample10_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 10',
                     names.attr=c(paste('spHW ',0:22,'kyr BP',sep='')))
dev.off()

##todos os mapas de clamping - multitemporal, spCD, sample 10 pts (todos as camadas temporais)
scenNames =  grep(pattern='^multitemporal.*spCD.*sample10.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
scenNames =  grep(pattern='^multitemporal.*spCD.*sample100.*replica1', x=scenNames, value=TRUE, invert=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Multitemporal_SpCD_sample10_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 10',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()

##todos os mapas de clamping - multitemporal, spHW, sample 50 pts (todos as camadas temporais)
scenNames =  grep(pattern='^multitemporal.*spHW.*sample50.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Multitemporal_SpHW_sample50_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 50',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

##todos os mapas de clamping - multitemporal, spCD, sample 50 pts (todos as camadas temporais)
scenNames =  grep(pattern='^multitemporal.*spCD.*sample50.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Multitemporal_SpCD_sample50_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 50',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()

##todos os mapas de clamping - multitemporal, spHW, sample 100 pts (todos as camadas temporais)
scenNames =  grep(pattern='^multitemporal.*spHW.*sample100.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Multitemporal_SpHW_sample100_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 100',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

##todos os mapas de clamping - multitemporal, spCD, sample 100 pts (todos as camadas temporais)
scenNames =  grep(pattern='^multitemporal.*spCD.*sample100.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Multitemporal_SpCD_sample100_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 100',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()


###monotemporal

##todos os mapas de clamping - monotemporal, spHW, sample 10 pts (todos as camadas temporais)
scenNames =  grep(pattern='^monotemporal.*spHW.*sample10.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
scenNames =  grep(pattern='^monotemporal.*spHW.*sample100.*replica1', x=scenNames, value=TRUE, invert=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Monotemporal_SpHW_sample10_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 10',
                     names.attr=c(paste('HW sp. ',0:23,'kyr BP',sep='')))
dev.off()

##todos os mapas de clamping - monotemporal, spCD, sample 10 pts (todos as camadas temporais)
scenNames =  grep(pattern='^monotemporal.*spCD.*sample10.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Monotemporal_SpCD_sample10_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 10',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()

##todos os mapas de clamping - monotemporal, spHW, sample 50 pts (todos as camadas temporais)
scenNames =  grep(pattern='^monotemporal.*spHW.*sample50.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Monotemporal_SpHW_sample50_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 50',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

##todos os mapas de clamping - monotemporal, spCD, sample 50 pts (todos as camadas temporais)
scenNames =  grep(pattern='^monotemporal.*spCD.*sample50.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Monotemporal_SpCD_sample50_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 50',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()

##todos os mapas de clamping - monotemporal, spHW, sample 100 pts (todos as camadas temporais)
scenNames =  grep(pattern='^monotemporal.*spHW.*sample100.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Monotemporal_SpHW_sample100_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 100',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

##todos os mapas de clamping - monotemporal, spCD, sample 100 pts (todos as camadas temporais)
scenNames =  grep(pattern='^monotemporal.*spCD.*sample100.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Monotemporal_SpCD_sample100_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 100',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()




### graficos de tendencia ###



##dataClamp = data.frame(clamping=as.numeric(territory), scenario=names(territory))
dataClamp = do.call('rbind', territory)
dataClamp$scenario = names(territory)
rownames(dataClamp) = NULL

## ##
## dataClamp$sdm = NA
## dataClamp[ grep('monotemporal', dataClamp$scenario), ]$sdm = 'monotemporal'
## dataClamp[ grep('multitemporal', dataClamp$scenario), ]$sdm = 'multitemporal'
## dataClamp$kyr = c(0:23)
## dataClamp$sp = NA
## dataClamp[grep('spHW', dataClamp$scenario),]$sp = 'spHW'
## dataClamp[grep('spCD', dataClamp$scenario),]$sp = 'spCD'
## dataClamp$sample = NA
## dataClamp[grep('10.replica', dataClamp$scenario),]$sample = 10
## dataClamp[grep('50.replica', dataClamp$scenario),]$sample = 50
## dataClamp[grep('100.replica', dataClamp$scenario),]$sample = 100



## clamping na America do Sul X Schoener's D (D = sobreposicao distribuicao real X distribuicao modelada) ##

jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/clampAmSulXschoenerXsampleXsps.jpeg', width=1200, height=1200)
par(mfrow=c(2,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.5, cex.axis=2, cex.main=2)
##
##spCD monotemporal
##
plot(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("CD species - SDM"["mono"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping), ylim=c(0,1), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping), ylim=c(0,1), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping)), lty=2, col='black')
##
##spCD multitemporal
##
plot(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(South America area (in %))", main=expression("CD species - SDM"["multi"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping), ylim=c(0,1), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping), ylim=c(0,1), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping)), lty=2, col='black')
##
##spHW monotemporal
##
plot(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(South America area (in %))", main=expression("HW species - SDM"["mono"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping)), lty=2, col='black')
##
##spHW multitemporal
##
plot(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(South America area (in %))", main=expression("HW species - SDM"["multi"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping)), lty=2, col='black')
##
dev.off()



## proporcao de clamping na distribuicao modelada X Schoener's D (D = sobreposicao distribuicao real X distribuicao modelada) ##

jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/clampSpsDistXschoenerXsampleXsps.jpeg', width=1200, height=1200)
par(mfrow=c(2,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.5, cex.axis=2, cex.main=2)
##
##spCD monotemporal
##
plot(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-3,1), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("CD species - SDM"["mono"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spCD multitemporal
##
plot(dataClamp[which(dataClamp$distUnderClamp>=0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>=0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-6.5,0), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("CD species - SDM"["multi"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp), pch=20, ylab="Schoener's D", col='black')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spHW monotemporal
##
plot(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-10,5), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("HW species - SDM"["mono"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spHW multitemporal
##
plot(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-10,0), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("HW species - SDM"["multi"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
dev.off()



## Area de clamping total na America do Sul X porporcao de clamping na distribuicao modelada ##
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/clampSpsDistXclampAmSul.jpeg', width=1200, height=1200)
par(mfrow=c(2,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.5, cex.axis=2, cex.main=2)
##
##spCD monotemporal
##
plot(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-3.5,1), ylim=c(-2,3), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="log(South America area (in %))", main=expression("CD species - SDM"["mono"]), col='red')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["mono"]), col='blue')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["mono"]), col='black')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spCD multitemporal
##
plot(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-6.5,-0.5), ylim=c(-3,0.5), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="log(South America area (in %))", main=expression("CD species - SDM"["multi"]), col='red')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["multi"]), col='blue')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["multi"]), col='black')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spHW monotemporal
##
plot(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-3.5,1), ylim=c(-2,3), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="log(South America area (in %))", main=expression("HW species - SDM"["mono"]), col='red')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["mono"]), col='blue')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["mono"]), col='black')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spHW multitemporal
##
plot(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-6.5,-0.5), ylim=c(-3,0.5), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="log(South America area (in %))", main=expression("HW species - SDM"["multi"]), col='red')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["multi"]), col='blue')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["multi"]), col='black')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
dev.off()



## tendencia do clamping no tempo


jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/clampXtempo.jpeg', width=900)
par(mfrow=c(1,2))
plot(dataClamp[which(dataClamp$sdm=='multitemporal'),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal'),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='South América area (in %)', main='Multitemporal')
plot(dataClamp[which(dataClamp$sdm=='monotemporal'),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal'),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='South América area (in %)', main='Monotemporal')
dev.off()




## tendencia do clamping no tempo, considerando sps, sample size e tipo de calibracao

## graficos para spCD
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/clampXtempoXsample_spCD.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.7, cex.axis=2, cex.main=2)
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$kyr), ylim=c(0,15), xlab='',  ylab='', main=expression('SDM'['mono']), col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$kyr), ylim=c(0,15), xlab='',  ylab='', main=expression('SDM'['multi']), col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$kyr), ylim=c(0,15), xlab='',  ylab='South América area (in %)', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$kyr), ylim=c(0,15), xlab='',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
dev.off()


## graficos para spHW
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/clampXtempoXsample_spHW.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.7, cex.axis=2, cex.main=2)
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$kyr), ylim=c(0,15), xlab='',  ylab='', main=expression('SDM'['mono']), col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$kyr), ylim=c(0,15), xlab='',  ylab='', main=expression('SDM'['multi']), col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$kyr), ylim=c(0,15), xlab='',  ylab='South América area (in %)', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$kyr), ylim=c(0,15), xlab='',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
dev.off()




## boxplot clamping (full dataset e sps)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/BoxplotClamp.jpeg', width=1200)
par(mfrow=c(1,3), cex=1.3)
boxplot(log(dataClamp$clamping) ~ dataClamp$sdm, ylab='log(South America area (in %))', main='Full dataset')
##
boxplot(log(dataClamp[which(dataClamp$sp == 'spHW'),]$clamping) ~ dataClamp[which(dataClamp$sp == 'spHW'),]$sdm, ylab='log(South America area (in %))', main='HW species')
##
boxplot(log(dataClamp[which(dataClamp$sp == 'spCD'),]$clamping) ~ dataClamp[which(dataClamp$sp == 'spCD'),]$sdm, ylab='log(South America area (in %))', main='CD species')
dev.off()



##correlacao entre Schoener X % de clamping na distribuicao da especie

##monotemporal - spCD
corSpCDmono10 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$distUnderClamp, method='pearson')

corSpCDmono50 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$distUnderClamp, method='pearson')

corSpCDmono100 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$Schoeners_D_simi,
                          dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$distUnderClamp, method='pearson')

##monotemporal - spHW
corSpHWmono10 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$distUnderClamp, method='pearson')

corSpHWmono50 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$distUnderClamp, method='pearson')

corSpHWmono100 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$Schoeners_D_simi,
                          dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$distUnderClamp, method='pearson')

##multitemporal - spCD
corSpCDmulti10 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$Schoeners_D_simi,
                          dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$distUnderClamp, method='pearson')

corSpCDmulti50 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$Schoeners_D_simi,
                          dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$distUnderClamp, method='pearson')

corSpCDmulti100 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$Schoeners_D_simi,
                           dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$distUnderClamp, method='pearson')

##multitemporal - spHW
corSpHWmulti10 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$Schoeners_D_simi,
                          dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$distUnderClamp, method='pearson')

corSpHWmulti50 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$Schoeners_D_simi,
                          dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$distUnderClamp, method='pearson')

corSpHWmulti100 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$Schoeners_D_simi,
                           dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$distUnderClamp, method='pearson')



##correlacao geral (so por exploracao)

corFullDataset = cor.test(dataClamp$Schoeners_D_simi,
                          dataClamp$distUnderClamp, method='pearson')

corMulti = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal'),]$Schoeners_D_simi,
                    dataClamp[ which(dataClamp$sdm=='multitemporal'),]$distUnderClamp, method='pearson')

corMono = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal'),]$Schoeners_D_simi,
                   dataClamp[ which(dataClamp$sdm=='monotemporal'),]$distUnderClamp, method='pearson')

##tabela
corTable = data.frame( scenario = c('corFullDataset','corMulti','corMono','corSpCDmono10','corSpCDmono50','corSpCDmono100','corSpHWmono10','corSpHWmono50','corSpHWmono100','corSpCDmulti10','corSpCDmulti50','corSpCDmulti100','corSpHWmulti10','corSpHWmulti50','corSpHWmulti100'),
                       correlation = c(as.numeric(corFullDataset$estimate),as.numeric(corMulti$estimate),as.numeric(corMono$estimate),as.numeric(corSpCDmono10$estimate),as.numeric(corSpCDmono50$estimate),as.numeric(corSpCDmono100$estimate),as.numeric(corSpHWmono10$estimate),as.numeric(corSpHWmono50$estimate),as.numeric(corSpHWmono100$estimate),as.numeric(corSpCDmulti10$estimate),as.numeric(corSpCDmulti50$estimate),as.numeric(corSpCDmulti100$estimate),as.numeric(corSpHWmulti10$estimate),as.numeric(corSpHWmulti50$estimate),as.numeric(corSpHWmulti100$estimate)),
                       p.value = c(as.numeric(corFullDataset$p.value),as.numeric(corMulti$p.value),as.numeric(corMono$p.value),as.numeric(corSpCDmono10$p.value),as.numeric(corSpCDmono50$p.value),as.numeric(corSpCDmono100$p.value),as.numeric(corSpHWmono10$p.value),as.numeric(corSpHWmono50$p.value),as.numeric(corSpHWmono100$p.value),as.numeric(corSpCDmulti10$p.value),as.numeric(corSpCDmulti50$p.value),as.numeric(corSpCDmulti100$p.value),as.numeric(corSpHWmulti10$p.value),as.numeric(corSpHWmulti50$p.value),as.numeric(corSpHWmulti100$p.value)) )

corTable[,'p.value'] = round(corTable[,'p.value'], 3)

write.csv(corTable, paste(projectFolder,'correlationTable.csv'), row.names=FALSE)



### mapas para comparar condicoes do presente com 22 kyrBP


temp0kyr = raster('/home/anderson/PosDoc/dados_ambientais/dados_projeto/000/bioclim_01.asc')
temp0kyr = mask(temp0kyr, AmSulShape)
preci0kyr = raster('/home/anderson/PosDoc/dados_ambientais/dados_projeto/000/bioclim_12.asc')
preci0kyr = mask(preci0kyr, AmSulShape)

temp22kyr = raster('/home/anderson/PosDoc/dados_ambientais/dados_projeto/022/bioclim_01.asc')
temp22kyr = mask(temp22kyr, AmSulShape)
preci22kyr = raster('/home/anderson/PosDoc/dados_ambientais/dados_projeto/022/bioclim_12.asc')
preci22kyr = mask(preci22kyr, AmSulShape)

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/temp0&22kyr.jpeg', width=800)
par(mfrow=c(1,2), mar=c(5,5,5,6))
plot(temp0kyr, main='0 kyr BP'); grid()
plot(temp22kyr, main='22 kyr BP'); grid()
dev.off()

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/preci0&22kyr.jpeg', width=800)
par(mfrow=c(1,2), mar=c(5,5,5,6))
plot(preci0kyr, main='0 kyr BP'); grid()
plot(preci22kyr, main='22 kyr BP'); grid()
dev.off()

clamp0 = mask(clampStack$multitemporal_proj_0kyr_spHW.sample50.replica1, AmSulShape)
clamp22 = mask(clampStack$multitemporal_proj_22kyr_spHW.sample50.replica1, AmSulShape)

par(mfrow=c(1,2))
plot(clamp0, main=c(paste('spHW ',0,'kyr BP',sep='')), col=c('lightgrey','red'), legend=FALSE)
plot(clamp22, main=c(paste('spHW ',22,'kyr BP',sep='')), col=c('lightgrey','red'), legend=FALSE)

dev.off()


## ## graficos pres/aus variaveis (sempre com 100 pts)

## ##spHW multitemporal
## arqvs = list.files(path=file.path(projectFolder,'Amostras','multitemporal','spHW'), pattern='100pts', full.names=TRUE)

## spHWmultiOCC = read.csv(arqvs[10], header=TRUE)
## spHWmultiOCC$pres = 1
## spHWmultiBG = read.csv(arqvs[5], header=TRUE)
## spHWmultiBG$pres = 0
## spHWmulti = rbind(spHWmultiOCC, spHWmultiBG)

## mdBio01Mult = glm(pres ~ bioclim_01, family='binomial', data=spHWmulti)
## mdBio12Mult = glm(pres ~ bioclim_12, family='binomial', data=spHWmulti)

## plot(spHWmulti[,'pres'] ~ spHWmulti[,c('bioclim_01')])
## points(predict(mdBio01Mult,newdata=data.frame(bioclim_01=c(-10:310)), type='response'), type='l')

## plot(spHWmulti[,'pres'] ~ spHWmulti[,c('bioclim_12')], xlim=c(-10,20000))
## points(predict(mdBio12Mult,newdata=data.frame(bioclim_12=c(min(spHWmono$bioclim_12):max(spHWmulti$bioclim_12))), type='response'), type='l')


## ##spHW monotemporal
## arqvs = list.files(path=file.path(projectFolder,'Amostras','monotemporal','spHW'), pattern='100pts', full.names=TRUE)

## spHWmonoOCC = read.csv(arqvs[10], header=TRUE)
## spHWmonoOCC$pres = 1
## spHWmonoBG = read.csv(arqvs[5], header=TRUE)
## spHWmonoBG$pres = 0
## spHWmono = rbind(spHWmonoOCC, spHWmonoBG)

## mdBio01Mono = glm(pres ~ bioclim_01, family='binomial', data=spHWmono)
## mdBio12Mono = glm(pres ~ bioclim_12, family='binomial', data=spHWmono)

## plot(spHWmono[,'pres'] ~ spHWmono[,c('bioclim_01')])
## points(predict(mdBio01Mono,newdata=data.frame(bioclim_01=c(-10:310)), type='response'), type='l')

## plot(spHWmono[,'pres'] ~ spHWmono[,c('bioclim_12')])
## points(predict(mdBio12Mono,newdata=data.frame(bioclim_12=c(0:15000)), type='response'), type='l')

## points(spHWmono[,'pres'] ~ spHWmono[,c('bioclim_12')], col='red')
## points(predict(mdBio12Mono,newdata=data.frame(bioclim_12=c(min(spHWmono$bioclim_12):max(spHWmono$bioclim_12))), type='response'), type='l', col='red')


##distribuicao presente, inter e maximo glacial

library(raster)
library(maptools)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp")

### MULTITEMPORAL ###

## spHW ##

##  distribuicao real
HWcurrentReal = raster(paste(projectFolder,'NichoReal/spHW/000.asc',sep='')) > 0.2
HW22Real = raster(paste(projectFolder,'NichoReal/spHW/022.asc',sep='')) > 0.2

## multitemporal
HWModel_0kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/multitemporal/spHW/spHW.sample10.replica1/proj_0kyr/proj_0kyr_spHW.sample10.replica1_TSSbin.grd')
HWModel_0kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/multitemporal/spHW/spHW.sample50.replica1/proj_0kyr/proj_0kyr_spHW.sample50.replica1_TSSbin.grd')
HWModel_0kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/multitemporal/spHW/spHW.sample100.replica1/proj_0kyr/proj_0kyr_spHW.sample100.replica1_TSSbin.grd')
##
HWModel_22kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/multitemporal/spHW/spHW.sample10.replica1/proj_22kyr/proj_22kyr_spHW.sample10.replica1_TSSbin.grd') 
HWModel_22kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/multitemporal/spHW/spHW.sample50.replica1/proj_22kyr/proj_22kyr_spHW.sample50.replica1_TSSbin.grd')
HWModel_22kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/multitemporal/spHW/spHW.sample100.replica1/proj_22kyr/proj_22kyr_spHW.sample100.replica1_TSSbin.grd')

## spCD ##

##  distribuicao real
CDcurrentReal = raster(paste(projectFolder,'NichoReal/spCD/000.asc',sep='')) > 0.2
CD22Real = raster(paste(projectFolder,'NichoReal/spCD/022.asc',sep='')) > 0.2

## SDM multitemporal
CDModel_0kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/multitemporal/spCD/spCD.sample10.replica1/proj_0kyr/proj_0kyr_spCD.sample10.replica1_TSSbin.grd')
CDModel_0kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/multitemporal/spCD/spCD.sample50.replica1/proj_0kyr/proj_0kyr_spCD.sample50.replica1_TSSbin.grd')
CDModel_0kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/multitemporal/spCD/spCD.sample100.replica1/proj_0kyr/proj_0kyr_spCD.sample100.replica1_TSSbin.grd')
##
CDModel_22kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/multitemporal/spCD/spCD.sample10.replica1/proj_22kyr/proj_22kyr_spCD.sample10.replica1_TSSbin.grd') 
CDModel_22kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/multitemporal/spCD/spCD.sample50.replica1/proj_22kyr/proj_22kyr_spCD.sample50.replica1_TSSbin.grd')
CDModel_22kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/multitemporal/spCD/spCD.sample100.replica1/proj_22kyr/proj_22kyr_spCD.sample100.replica1_TSSbin.grd')

##sobreposicoes spHW

jpeg(filename='/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/sobreposicoesHWmulti.jpg', width = 1400 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,20), mar=c(3,3,5,6))
plot(HWcurrentReal*1+HWModel_0kyrSample10*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample50*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample100*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
##legenda
legend("topright",legend=c('Virtual species','Maxent projection','Overlap'),inset=c(-0.7,0),xpd=NA,pch=20,col=c('green','blue','dark green'),cex=2.5)
##
plot(HW22Real*1+HWModel_22kyrSample10*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample50*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample100*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
mtext('HW species',outer=TRUE,cex=4)
dev.off()

##sobreposicoes spCD

jpeg(filename='/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/sobreposicoesCDmulti.jpg', width = 1400 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,20), mar=c(3,3,5,6))
plot(CDcurrentReal*1+CDModel_0kyrSample10*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample50*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample100*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
##legenda
legend("topright",legend=c('Virtual species','Maxent projection','Overlap'),inset=c(-0.7,0),xpd=NA,pch=20,col=c('green','blue','dark green'),cex=2.5)
##
plot(CD22Real*1+CDModel_22kyrSample10*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample50*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample100*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
mtext('CD species',outer=TRUE,cex=4)
dev.off()


### MONOTEMPORAL ###


## spHW ##

##  distribuicao real
HWcurrentReal = raster(paste(projectFolder,'NichoReal/spHW/000.asc',sep='')) > 0.2
HW22Real = raster(paste(projectFolder,'NichoReal/spHW/022.asc',sep='')) > 0.2

## monotemporal
HWModel_0kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample10.replica1/proj_0kyr/proj_0kyr_spHW.sample10.replica1_TSSbin.grd')
HWModel_0kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample50.replica1/proj_0kyr/proj_0kyr_spHW.sample50.replica1_TSSbin.grd')
HWModel_0kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample100.replica1/proj_0kyr/proj_0kyr_spHW.sample100.replica1_TSSbin.grd')
##
HWModel_22kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample10.replica1/proj_22kyr/proj_22kyr_spHW.sample10.replica1_TSSbin.grd') 
HWModel_22kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample50.replica1/proj_22kyr/proj_22kyr_spHW.sample50.replica1_TSSbin.grd')
HWModel_22kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample100.replica1/proj_22kyr/proj_22kyr_spHW.sample100.replica1_TSSbin.grd')

## spCD ##

##  distribuicao real
CDcurrentReal = raster(paste(projectFolder,'NichoReal/spCD/000.asc',sep='')) > 0.2
CD22Real = raster(paste(projectFolder,'NichoReal/spCD/022.asc',sep='')) > 0.2

## SDM monotemporal
CDModel_0kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample10.replica1/proj_0kyr/proj_0kyr_spCD.sample10.replica1_TSSbin.grd')
CDModel_0kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample50.replica1/proj_0kyr/proj_0kyr_spCD.sample50.replica1_TSSbin.grd')
CDModel_0kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample100.replica1/proj_0kyr/proj_0kyr_spCD.sample100.replica1_TSSbin.grd')
##
CDModel_22kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample10.replica1/proj_22kyr/proj_22kyr_spCD.sample10.replica1_TSSbin.grd') 
CDModel_22kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample50.replica1/proj_22kyr/proj_22kyr_spCD.sample50.replica1_TSSbin.grd')
CDModel_22kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample100.replica1/proj_22kyr/proj_22kyr_spCD.sample100.replica1_TSSbin.grd')

##sobreposicoes spHW
jpeg(filename='/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/sobreposicoesHWmono.jpg', width = 1400 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,20), mar=c(3,3,5,6))
plot(HWcurrentReal*1+HWModel_0kyrSample10*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample50*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample100*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
##legenda
legend("topright",legend=c('Virtual species','Maxent projection','Overlap'),inset=c(-0.7,0),xpd=NA,pch=20,col=c('green','blue','dark green'),cex=2.5)
##
plot(HW22Real*1+HWModel_22kyrSample10*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample50*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample100*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
mtext('HW species.',outer=TRUE,cex=4)
dev.off()

##sobreposicoes spCD
jpeg(filename='/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/sobreposicoesCDmono.jpg', width = 1400 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,20), mar=c(3,3,5,6))
plot(CDcurrentReal*1+CDModel_0kyrSample10*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample50*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample100*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
##legenda
legend("topright",legend=c('Virtual species','Maxent projection','Overlap'),inset=c(-0.7,0),xpd=NA,pch=20,col=c('green','blue','dark green'),cex=2.5)
##
plot(CD22Real*1+CDModel_22kyrSample10*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample50*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample100*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
mtext('CD species',outer=TRUE,cex=4)
dev.off()


