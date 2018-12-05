##abrindo pacotes e funcoes necessarias##
library(raster)
library(hypervolume)
library(biomod2)
library(sensitivity)
library(ENMeval)
library(dismo)
# ##notebook
# source('/home/anderson/R-Scripts/paleoextract.R')
# source('/home/anderson/R-Scripts/strings2na.R')
# source('/home/anderson/R-Scripts/dataInstance.R')
# source('/home/anderson/R-Scripts/uplot.R')
# ##source('/home/anderson/R-Scripts/uniche.R')
# source('/home/anderson/R-Scripts/uniche2.R') #uniche2
# source('/home/anderson/R-Scripts/cleanData.R')
##Yavanna
source('D:/Anderson_Eduardo/R/paleoextract.R')
source('D:/Anderson_Eduardo/R/strings2na.R')
source('D:/Anderson_Eduardo/R/dataInstance.R')
source('D:/Anderson_Eduardo/R/uplot.R')
source('D:/Anderson_Eduardo/R/uplot2.R')
##source('/home/anderson/R-Scripts/uniche.R')
source('D:/Anderson_Eduardo/R/uniche2.R') #uniche2
source('D:/Anderson_Eduardo/R/uniche3.R') #uniche3
source('D:/Anderson_Eduardo/R/cleanData.R')


##definindo parametros e variaveis globais##
# ##notebook
# projectFolder = "/home/anderson/Projetos/SDM megafauna Sul-Americana" #pasta de trabalho do projeto
# envFolder = "/home/anderson/gridfiles/dados_projeto" #pasta das variaveis ambientais
# AmSulBorders = rgdal::readOGR('/home/anderson/shapefiles/Am_Sul/borders.shp') #shapefile da Am. do Sul
##Yavanna
projectFolder = "D:/Anderson_Eduardo/SDM megafauna Sul-Americana" #pasta de trabalho do projeto
envFolder = "D:/Anderson_Eduardo/gridfiles/variaveis ambientais AmSul 120kyr" #pasta das variaveis ambientais
AmSulBorders = rgdal::readOGR('D:/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp') #shapefile da Am. do Sul



##analise de incerteza de sensibilidade##



##arquivo do banco de dados
#dataSetRaw = read.csv(file='/home/anderson/Projetos/SDM megafauna Sul-Americana/dataset_clean.csv', header=TRUE, dec='.', sep=',')
dataSetRaw = read.csv(file='D:/Anderson_Eduardo/SDM megafauna Sul-Americana/dataset_clean.csv', header=TRUE, dec='.', sep=',') #abrindo e tratando o banco de dados
maxentFolder = 'D:/Anderson_Eduardo/maxent'

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
  analise.incerteza = uniche3(x=pts, cols=c("Longitude","Latitude", "Cal..Mean","Min.","Max."), envFolder=envFolder, maxentFolder=maxentFolder) #analise de incerteza e sensibilidade
  
  ##salvando output da analise
  save(analise.incerteza, file=paste(projectFolder,'/Analise de sensibilidade e incerteza/',sps[i],'.R',sep=''))
  
  ##salvando output grafico
  jpeg(paste(projectFolder,'/Analise de sensibilidade e incerteza/',sps[i],'.jpeg',sep=''), 600, 600)
  uplot2(analise.incerteza, AmSulBorders, legend=FALSE)
  dev.off()
  
  gc()
  
}



## SDM ##

##pts choosed to be excluded
ptsToExclude = list(
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
  
  occDataList = dataInstance(x=occData, col_names=c('ageMean', 'ageMin', 'ageMax'), tempRes=1) #Obs.: put tempres=1 if age already is in the desired temporal resolution.
  
  iter=1 ##loop over data intances
  for(j in seq(length(occDataList))){
  
    names(occDataList[[iter]]) = c("lon","lat","originalID","age")
    
    current_occData = paleoextract(x=occDataList[[iter]], path=envFolder)
    current_occData = current_occData[complete.cases(current_occData),]
    
    ##cross-temporal background points
    current_bgData = paleobg(x=current_occData, colNames=c('lon','lat','age'), envFolder=envFolder, n=10000) #sanpling bg points
    
    current_bgData[,c('x','y')] = round(current_bgData[,c('x','y')], 2) #round data
    current_bgData = current_bgData[!duplicated(current_bgData[,c('x','y','age')]), ] #excluding duplicated points (in a same age)
    current_bgData$originalID = NA #adjusting to fit occ data.frame
    idxOrdem = match(names(current_occData), names(current_bgData)) #adjusting to fit occ data.frame
    current_bgData = current_bgData[,idxOrdem] #adjusting to fit occ data.frame
    
    ##occ + bg points
    current_occData$type = 'occ' #assigning data type
    current_bgData$type = 'bg' #assigning data type
    current_dataSet = rbind(current_occData, current_bgData) #dataset consolidation
    current_dataSet = current_dataSet[complete.cases(current_dataSet[, grep(pattern='originalID', x=names(current_dataSet), invert=TRUE) ]), ]
    
    gc() #hoping to clean some computer memory...
  
    ##variaveis e parametros locais especificos para o biomod2
    myRespName <- gsub(pattern=' ',x=sps[i], replacement='_') # nome do cenario atual (para biomod2)
    myResp <- data.frame(type=as.integer(current_dataSet[,c('type')] == 'occ')) # variavel resposta (para biomod2)
    myRespXY <- current_dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
    myExpl = current_dataSet[,grep(pattern="originalID|age|type", x=names(current_dataSet),invert=TRUE)]  #variavel preditora (para biomod2)
    
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
      NbRunEval = 100,
      DataSplit = 75,
      VarImport = 10,
      models.eval.meth = c('TSS','ROC'),
      SaveObj = FALSE,
      rescal.all.models = TRUE,
      do.full.models = FALSE,
      modeling.id = paste(myRespName))
    
    ##My output data
    evaluationScores = get_evaluations(myBiomodModelOut)
  
    ##salvando modelo no HD
    save(myBiomodModelOut, file=paste(projectFolder,'/SDMs/',sps[i],'_SDM.R',sep=''))
    
    
  
    ##variableimportance final model
    varImportOutput = rowMeans(SDMmaxent@results[grep("importance", rownames(SDMmaxent@results)), ])
    varImportTable = data.frame(names(varImportOutput), permutation.importance=varImportOutput)
    write.csv(varImportTable, paste(projectFolder,'/maxent/omobranchus_variablesImportance_finalModel.csv',sep=''), row.names=FALSE)
    
    ##predicao AREA NATIVA
    names(predAreaNat) = gsub( pattern='^predAreaNat_', replacement='', x=names(predAreaNat) ) #retirando parte do nome para adequar ao modelo
    predAreaNat = predAreaNat[[ grep(pattern=paste(varNames,collapse='|'), x=names(predAreaNat), value=TRUE) ]] #novo objeto, correto
    SDMpredNativ = dismo::predict(SDMmaxent, predAreaNat, progress='text') #projecao espacial
    crs(SDMpredNativ) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
    
    ##predicao AREA INVADIDA
    names(predAreaInv) = gsub( pattern='^predAreaInv_', replacement='', x=names(predAreaInv) ) #retirando parte do nome para adequar ao modelo
    predAreaInv = predAreaInv[[ grep(pattern=paste(varNames,collapse='|'), x=names(predAreaNat), value=TRUE) ]] #novo objeto, correto
    SDMpredInv = dismo::predict(SDMmaxent, predAreaInv, progress='text') #projecao espacial
    crs(SDMpredInv) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
    
    
    ##salvando gridfiles
    writeRaster(calc(SDMpredNativ,mean), paste(projectFolder,'/maxent/SDMpredNativ_Suitability_Mean.asc',sep=''), overwrite=TRUE) #area nativa
    writeRaster(calc(SDMpredNativ,sd), paste(projectFolder,'/maxent/SDMpredNativ_Suitability_SD.asc',sep=''), overwrite=TRUE) #area nativa
    ##
    writeRaster(calc(SDMpredInv,mean),paste(projectFolder,'/maxent/SDMpredInv_Suitability_Mean.asc',sep=''), overwrite=TRUE) #area nativa
    writeRaster(calc(SDMpredInv,sd),paste(projectFolder,'/maxent/SDMpredInv_Suitability_SD.asc',sep=''), overwrite=TRUE) #area nativa
    
    
    ##calculo do threshold##
    ##valores preditos pelo SDM em cada ponto do dataset na AREA NATIVA
    pred =  extract(x=calc(SDMpredNativ,mean),
                    y=dataSet[dataSet$area!='inv',c('lon','lat')],
                    method='bilinear',
                    na.rm=TRUE)
    
    threshDF = data.frame(occ=dataSet[dataSet$area!='inv',]$occ, pred=pred) #juntando predicoes e observacoes
    threshDF = threshDF[complete.cases(threshDF),] #retirando possiveis NAs
    
    ##OBS.: omissao -> falso negativo -> false negative rate (FNR) = 1 - TPR (True Positive Rate)
    ##TPR = sensitividade
    myRoc = roc(predictor=threshDF$pred, response=threshDF$occ, positive=1) #curva ROC
    rocDF = data.frame( FNR = 1-myRoc$sensitivities, thresholds = myRoc$thresholds ) #data.frame com omissao e thresholds
    thre = max(rocDF[which(rocDF$FNR == max(rocDF[rocDF$FNR <= 0.05 ,]$FNR)),]$thresholds) #thresold de 5% omissao
    
    
    ##calculo do TSS para com o threshold -- OBS: TSS = sensitivity + specificity - 1
    TSSvalue = myRoc$sensitivities[which(myRoc$thresholds == thre)] + myRoc$specificities[which(myRoc$thresholds == thre)] - 1
    
    
    ##salvando output do ENMeval
    SDMevalOutput = as.data.frame(SDMeval@results)
    write.csv(SDMevalOutput, paste(projectFolder,'/maxent/omobranchus_SDMeval.csv',sep=''), row.names=FALSE)
    
    
    ##salvando graficos do ENMeval
    jpeg(paste(projectFolder,'/maxent/omobranchus_SDMeval.jpg',sep=''), width=800)
    par(mfrow=c(1,2), mar=c(5,5,10,5))
    eval.plot(SDMeval@results)
    eval.plot(SDMeval@results, 'Mean.AUC', var='Var.AUC')
    dev.off()
    
    
    ##salvando imagens jpeg no HD
  
    }
    
  }




########################################################




##function for cross-temporal background points

paleobg = function(x, colNames=names(x), envFolder, n=10000){
  
  ### parameters ###
  ##x = data.frame with (at least) occurrences (lon, lat) and age
  ##colNames = names of the columns in x for lon, lat and age
  ##envFolder = path to predictor variables' folder
  ##n = desired number of background points
  ################

  colNames = colNames
  envFolder = envFolder
  occTable = x[,colNames]
  names(occTable) = c('lon','lat','age')
  ages = unique(occTable$age)
  ages = sample(ages, n, replace=TRUE)
  bgData = data.frame()
  
  for (age_i in unique(ages)){
    sampleSize = length(grep(age_i, ages))
    current_occPts = occTable[age_i, c('lon','lat')]
    
    if(age_i %in% as.numeric(list.files(envFolder))){
      envData = list.files(paste(envFolder,'/',age_i,sep=''), full.names=TRUE)
    }
    
    envData = stack(envData)
    crs(envData) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
    
    bgData_i = dismo::randomPoints(mask=envData[[1]], n=sampleSize, p=current_occPts)
    bgData_i = data.frame(bgData_i)
    bgData_i$age = age_i
    names(bgData_i) = c('lon','lat','age')
    bgData_i = data.frame(bgData_i, extract(envData, bgData_i[,c('lon','lat')]))
    bgData = rbind(bgData, bgData_i)
  }
  return(bgData)
}



