##abrindo pacotes e funcoes necessarias##
library(raster)
library(hypervolume)
library(biomod2)
library(sensitivity)
library(ENMeval)
library(dismo)
library(usdm)
##notebook
source('/home/anderson/R-Scripts/paleoextract.R')
source('/home/anderson/R-Scripts/strings2na.R')
source('/home/anderson/R-Scripts/dataInstance.R')
source('/home/anderson/R-Scripts/uplot2.R')
##source('/home/anderson/R-Scripts/uniche.R')
source('/home/anderson/R-Scripts/uniche3.R') #uniche2
source('/home/anderson/R-Scripts/cleanData.R')
source('/home/anderson/R-Scripts/paleobg.R')
source('/home/anderson/R-Scripts/stpv.R')
# ##Yavanna
# source('D:/Anderson_Eduardo/R/paleoextract.R')
# source('D:/Anderson_Eduardo/R/strings2na.R')
# source('D:/Anderson_Eduardo/R/dataInstance.R')
# source('D:/Anderson_Eduardo/R/uplot.R')
# source('D:/Anderson_Eduardo/R/uplot2.R')
# ##source('/home/anderson/R-Scripts/uniche.R')
# source('D:/Anderson_Eduardo/R/uniche2.R') #uniche2
# source('D:/Anderson_Eduardo/R/uniche3.R') #uniche3
# source('D:/Anderson_Eduardo/R/cleanData.R')
# source('D:/Anderson_Eduardo/R/paleobg.R')



##definindo parametros e variaveis globais##
##notebook
projectFolder = "/home/anderson/Projetos/SDM megafauna Sul-Americana" #pasta de trabalho do projeto
envFolder = "/home/anderson/gridfiles/dados_projeto" #pasta das variaveis ambientais
AmSulBorders = rgdal::readOGR('/home/anderson/shapefiles/Am_Sul/borders.shp') #shapefile da Am. do Sul
maxentFolder = "/home/anderson/R/x86_64-pc-linux-gnu-library/3.4/dismo/java" #'/home/anderson/R/x86_64-pc-linux-gnu-library/3.5/dismo/java'
# ##Yavanna
# projectFolder = "D:/Anderson_Eduardo/SDM megafauna Sul-Americana" #pasta de trabalho do projeto
# setwd(projectFolder)
# envFolder = "D:/Anderson_Eduardo/gridfiles/variaveis ambientais AmSul 120kyr" #pasta das variaveis ambientais
# AmSulBorders = rgdal::readOGR('D:/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp') #shapefile da Am. do Sul
# maxentFolder = 'D:/Anderson_Eduardo/maxent'



##analise de incerteza de sensibilidade##



##selecao de variaveis preditoras para analise multitemporal
stpvTable = stpv(path = envFolder, ageMin = 0, ageMax = 22)
preditorsStpv = c('bioclim_04', 'bioclim_10', 'bioclim_15', 'bioclim_16','bioclim_17')

##arquivo do banco de dados
dataSetRaw = read.csv(file='/home/anderson/Projetos/SDM megafauna Sul-Americana/dataset_clean.csv', header=TRUE, dec='.', sep=',')
#dataSetRaw = read.csv(file='D:/Anderson_Eduardo/SDM megafauna Sul-Americana/dataset_clean.csv', header=TRUE, dec='.', sep=',') #abrindo e tratando o banco de dados

##subset do banco de dados
dataSetRaw = dataSetRaw[,c('Species','Longitude','Latitude','Cal..Mean','Min.','Max.')]

##lista dos nomes das especies
sps = as.character( unique(dataSetRaw$Species) )

for (i in seq(length(sps))){
  
  ##dataset da especie atual
  pts = dataSetRaw[which(dataSetRaw$Species == sps[i]),]
  
  if (nrow(pts) <= 5){
    cat("\n ATENÇÃO: A espécie", sps[i], "possui menos que 5 registros, por isso não foi analisada. \n")
    next
  }else{
    cat("\n Rodando para a espécie", sps[i], "...\n")
  }
  
  ##ajustando dados
  pts = strings2na(pts, 'Species') #transformando strings ao longo do dateset em NA
  
  ## ##inspeção visual dos dados
  ## plot(AmSulBorders)
  ## points(pts[[1]][,c('lon','lat')], pch=20, cex=1.5, col='red')
  
  ##analise de incerteza e sensibilidade
  cols = c("Longitude","Latitude", "Cal..Mean","Min.","Max.", preditorsStpv)
  analise.incerteza = uniche3(x=pts, cols=cols, envFolder=envFolder, maxentFolder=maxentFolder) #analise de incerteza e sensibilidade
  
  ##salvando output da analise
  save(analise.incerteza, file=paste(projectFolder,'/Analise de sensibilidade e incerteza/',sps[i],'.R',sep=''))
  
  ##salvando output grafico
  jpeg(paste(projectFolder,'/Analise de sensibilidade e incerteza/',sps[i],'.jpeg',sep=''), 600, 600)
  uplot2(x =  analise.incerteza, shape =  AmSulBorders, legend=FALSE)
  dev.off()
  
  ##na esperanca de dar uma limpada na memoria...
  gc()
  
}



## SDMs - criando as instancias de dados



##working directory
setwd(paste(projectFolder,'/SDMs',sep=''))

##dataset
#dataSetRaw = read.csv(file='/home/anderson/Projetos/SDM megafauna Sul-Americana/dataset_clean.csv', header=TRUE, dec='.', sep=',')
dataSetRaw = read.csv(file= paste(projectFolder,'/dataset_clean.csv',sep=''), header=TRUE, dec='.', sep=',') #abrindo e tratando o banco de dados
maxentFolder = 'D:/Anderson_Eduardo/maxent'

##subset from the dataset
dataSetRaw = dataSetRaw[,c('Species','Longitude','Latitude','Cal..Mean','Min.','Max.')]

##list of species names
sps = as.character( unique(dataSetRaw$Species) )

##pts choosed to be excluded
```
ptsToExclude = 29
``` = list(
  Catonyx_chilensis = NA,
  Catonyx_cuvieri = NA,
  Cuvieronius_hyodon = NA,
  Doedicurus_clavicaudatus = 1,
  Eremotherium_laurillardi = 29,
  Glossotherium_robustum = NA,
  Glyptodon_clavipes = 2,
  Holmesina_occidentalis = 3,
  Holmesina_paulacoutoi = NA,
  Lestodon_armatus = NA,
  Megatherium_americanum = 34,
  Mylodon_darwinii = NA,
  Neosclerocalyptus_paskoensis = NA,
  Notiomastodon_platensis = NA,
  Pampatherium_humboldtii = c(1, 9),
  Pampatherium_typum = NA,
  Panochtus_tubreclatus = NA
)


##loop over species
for (i in seq(length(sps))){
  tryCatch({
    
    ##sps folder
    if (!file.exists(paste(projectFolder,'/SDMs/',sps[i],sep=''))){
      dir.create(paste(projectFolder,'/SDMs/',sps[i],sep=''), recursive=TRUE)
    }
    setwd(paste(projectFolder,'/SDMs/',sps[i],sep=''))
    
    ##abrindo dados
    load(paste(projectFolder,'/Analise de sensibilidade e incerteza/',sps[i],'.R',sep=''))
    
    ##excluding pts
    ptsToExcludeID = ptsToExclude[[agrep(sps[i], names(ptsToExclude))]]
    if ( !is.na(ptsToExcludeID) ){
      occData = analise.incerteza$dataset[ -match(ptsToExcludeID, analise.incerteza$dataset$ID),]
    }else{
      occData = analise.incerteza$dataset
    }
    
    ##cleaning the data
    occData[,c('lon','lat')] = round(occData[,c('lon','lat')], 2)
    occData = occData[!duplicated(occData[,c('lon','lat')]), ]
    
<<<<<<< HEAD
    occDataList = dataInstance(x=occData, col_names=c('ageMean', 'ageMin', 'ageMax'), tempRes=1) #Obs.: put tempres=1 if age already is in the desired temporal resolution.
    
=======
    occDataList = dataInstance(x=occData, col_names=c('ageMean', 'ageMin', 'ageMax'))
  
>>>>>>> adfc8bb7ff6e5a665f292f29750183ce813c84ca
    save(occDataList, file=paste(projectFolder, '/SDMs/', gsub(' ','_',sps[i]), '_dataInstances.R', sep=''))
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



## SDMs - calibracao dos SDMs



##loop over species
for (i in seq(length(sps))){
  tryCatch({
    
    ##diretorio de trabalho
    setwd(paste(projectFolder, '/SDMs/', sps[i], sep=''))
    
    gc() #hoping to clean some computer memory...
    
    rm(occDataList)
    load(occDataList, file=paste(projectFolder, '/SDMs/', gsub(' ','_',sps[i]), '_dataInstances.R', sep=''))
    
    ##cross-temporal background points
    current_bgData = paleobg(x=occDataList[[1]], colNames=c('lon','lat','age'), envFolder=envFolder, n=10000) #sanpling bg points
    current_bgData[,c('lon','lat')] = round(current_bgData[,c('lon','lat')], 2) #round data
    current_bgData = current_bgData[!duplicated(current_bgData[,c('lon','lat','age')]), ] #excluding duplicated points (in a same age)
    
    ##loop over data intances
    for(j in seq(length(occDataList))){ 
      
      names(occDataList[[j]]) = c("lon","lat","originalID","age")
      
      current_occData = paleoextract(x=occDataList[[j]], path=envFolder)
      current_occData = current_occData[complete.cases(current_occData),]
      
      ##cross-temporal background points
      current_bgData$originalID = NA #adjusting to fit occ data.frame
      idxOrdem = match(names(current_occData), names(current_bgData)) #adjusting to fit occ data.frame
      current_bgData = current_bgData[,idxOrdem] #adjusting to fit occ data.frame
      
      ##occ + bg points
      current_occData$type = 'occ' #assigning data type
      current_bgData$type = 'bg' #assigning data type
      current_dataSet = rbind(current_occData, current_bgData) #dataset consolidation
      current_dataSet = current_dataSet[complete.cases(current_dataSet[, grep(pattern='originalID', x=names(current_dataSet), invert=TRUE) ]), ]
      
      ##variaveis e parametros locais especificos para o biomod2
      myRespName <- paste( 'instance_',j, sep='') # nome do cenario atual (para biomod2)
      myResp <- data.frame(type=as.integer(current_dataSet[,c('type')] == 'occ')) # variavel resposta (para biomod2)
      myRespXY <- current_dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
      myExpl = current_dataSet[,grep(pattern="lon|lat|originalID|age|type", x=names(current_dataSet),invert=TRUE)]  #variavel preditora (para biomod2)
      
      ##ajuste de dados de entrada para biomod2
      myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                           expl.var = myExpl,
                                           resp.xy = myRespXY,
                                           resp.name = myRespName)
      
      ##parametrizando os modelos
      myBiomodOption <- BIOMOD_ModelingOptions(
        MAXENT.Phillips=list(
          path_to_maxent.jar=maxentFolder,
          linear=TRUE,
          quadratic=TRUE,
          product=FALSE,
          threshold=FALSE,
          hinge=FALSE,
          doclamp=FALSE,
          maximumiterations=1000,
          convergencethreshold=1.0E-5,
          threads=4))
      
      ##rodando o(s) algoritmo(s) (i.e. SDMs)
      myBiomodModelOut <- BIOMOD_Modeling(
        myBiomodData,
        models = c('MAXENT.Phillips'),
        models.options = myBiomodOption,
        NbRunEval = 50,
        DataSplit = 75,
        VarImport = 10,
        models.eval.meth = c('TSS','ROC'),
        SaveObj = TRUE,
        rescal.all.models = FALSE,
        do.full.models = FALSE,
        modeling.id = myRespName)
      
      ##My output data
      evaluationScores = get_evaluations(myBiomodModelOut)
      
      ##salvando modelo no HD
      save(myBiomodModelOut, file=paste(projectFolder,'/SDMs/',gsub(pattern=' ',x=sps[i], replacement='_'),'_instance_',j,'_SDM','.R',sep=''))
      
      ##variableimportance final model
      varImport = get_variables_importance(myBiomodModelOut)
      varImportMean = rowMeans(varImport[,,1:50,], na.rm=TRUE)
      varImportSD = apply(varImport[,,1:50,],1, sd )
      varImportDF = data.frame(variables=names(varImportMean), importance.mean=varImportMean, impotance.sd=varImportSD)
      write.csv(varImportDF, paste(projectFolder,'/SDMs/',sps[i],'/',gsub(pattern=' ',replacement='_',x=sps[i]),'_VariableImportance.csv',sep=''), row.names=FALSE)
      
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



## SDMs - projecao dos SDMs



##loop over species
for (i in seq(length(sps))){
  tryCatch({
    
    ##diretorio de trabalho
    setwd(paste(projectFolder, '/SDMs/', sps[i], sep=''))
    
    gc() #hoping to clean some computer memory...
    
    ##predicao
    envVarPaths = list.files(envFolder, full.names=TRUE)
    
    for (l in 1:length(envVarPaths[1:24])){
      
      ##definindo variaveis e parametros internos
      predictors = stack(list.files(path=envVarPaths[l], full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
      ##predictors = predictors[[c('bioclim_01','bioclim_12')]]
      ##predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
      crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
      
      myBiomodModelFull <- BIOMOD_Modeling(
        myBiomodData,
        models = c('MAXENT.Phillips'),
        models.options = myBiomodOption,
        NbRunEval = 1,
        DataSplit = 100,
        SaveObj = FALSE,
        rescal.all.models = FALSE,
        do.full.models = TRUE,
        modeling.id = paste(myRespName,'_FULL',sep=''))
      
      ##rodando algortmo de projecao (i.e. rodando a projecao)
      myBiomodProj <- BIOMOD_Projection(
        modeling.output = myBiomodModelFull,
        new.env = predictors,
        proj.name = paste(l-1,'kyr',sep=''),
        selected.models = 'all',
        binary.meth = c('TSS','ROC'),
        compress = 'TRUE',
        build.clamping.mask = 'TRUE',
        output.format = '.grd')
      
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}