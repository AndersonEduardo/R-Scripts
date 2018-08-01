##SCRIPT PARA DISTRIBUICAO DE BARBEIROS (SDMclim) - APENAS VARIAVEIS CLIMATICAS##

library(raster)
library(maptools)
library(adehabitatHR)
library(usdm)
##library(dismo)
library(ENMeval)
library(biomod2)
library(pROC)
#Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8.0_131/bin') # for 64-bit version
#Windows#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_91') # for 64-bit version
library(rJava)
#source("J:/Anderson_Eduardo/TSSmaxent.R") #sempre verificar aqui o caminho para o arquivo TSSmaxent.R
#Sys.setenv(JAVA_HOME = "/usr/lib/jvm/java-7-openjdk-amd64")
options(java.parameters = "Xmx7g")


## Definindo parametros e variaveis globais

# #notebook anderson
# envVarFolder = "/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/Variaveis Climaticas"
# spOccFolder = "/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias"
# projectFolder = "/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos"
# AmSulShape = rgdal::readOGR("/home/anderson/shapefiles/Am_Sul/borders.shp") #abrindo shape da America do Sul
# SAborders = rgdal::readOGR('/home/anderson/shapefiles/ne_50m_land/ne_50m_land.shp') #bordas de continentes
# SOAextent = extent(-81.57551,-34.03384,-57.13385,12.99115)
# SAborders = crop(SAborders,SOAextent)
# biasLayer = raster('/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/reduviidaeBiasLayer.grd')

#yavanna
envVarFolder = "D:/Anderson_Eduardo/Distribuicao de barbeiros com interacao com humanos/Variaveis Climaticas"
spOccFolder = "D:/Anderson_Eduardo/Distribuicao de barbeiros com interacao com humanos/Ocorrencias"
projectFolder = "D:/Anderson_Eduardo/Distribuicao de barbeiros com interacao com humanos"
AmSulShape = rgdal::readOGR("D:/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #abrindo shape da America do Sul
SAborders = rgdal::readOGR('D:/Anderson_Eduardo/shapefiles/ne_50m_land/ne_50m_land.shp') #bordas de continentes
SOAextent = extent(-81.57551,-34.03384,-57.13385,12.99115)
SAborders = crop(SAborders,SOAextent)
biasLayer = raster('D:/Anderson_Eduardo/Distribuicao de barbeiros com interacao com humanos/reduviidaeBiasLayer.grd')




###PRIMEIRA PARTE: limpando dados de occorrencia



##GBIF para Reduviidae##

occDataGBIF = read.csv('/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/GBIF_mar-2018/Reduviidae_CSV/reduviidae_GBIF_occData.csv', header=TRUE, sep=',', dec='.', stringsAsFactors=FALSE, na.strings="") #abrindo dados do computador
occDataGBIFClean = subset(occDataGBIF, coordinateuncertaintyinmeters < 5000 | is.na(coordinateuncertaintyinmeters)) # 377 retirados, com erro > 5 km
occDataGBIFclean[,c('decimallongitude','decimallatitude')] = round(x=occDataGBIFclean[,c('decimallongitude','decimallatitude')], digits=2) #coords com duas casas decimais
occDataGBIFclean = occDataGBIFclean[!duplicated(occDataGBIFclean[,c('decimallongitude','decimallatitude')]),] #retirando pontos duplicados
##definindo dataset limpo para o GBIF
occDataGBIF = occDataGBIFclean
write.csv(occDataGBIF,'/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/GBIF_mar-2018/Reduviidae_CSV/reduviidae_GBIF_occDataClean.csv',row.names=FALSE)


##SpsLink  para Reduviidae##


occDataSpsLink = read.csv('/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/SpsLink_mar-2018/occ_reduviidae-triatominae_speciesLink.csv', header=TRUE, sep=',', dec='.', stringsAsFactors=FALSE, na.strings="", colClasses=c('character',rep('numeric',5))) #abrindo dados
occDataSpsLink[is.na(occDataSpsLink$longitude),'longitude'] = occDataSpsLink[is.na(occDataSpsLink$longitude),'longitude_mun']
occDataSpsLink[is.na(occDataSpsLink$latitude),'latitude'] = occDataSpsLink[is.na(occDataSpsLink$latitude),'latitude_mun']
occDataSpsLinkClean = subset(occDataSpsLink, coordinateprecision < 5000 | is.na(coordinateprecision)) #407 retirados, com erro > 5 km
occDataSpsLinkClean = occDataSpsLinkClean[complete.cases(occDataSpsLinkClean[,c('longitude','latitude')]),] #eliminando dados sem coords
occDataSpsLinkClean[,c('longitude','latitude')] = round(x=occDataSpsLinkClean[,c('longitude','latitude')], digits=2) #coords com duas casas decimais
occDataSpsLinkClean = occDataSpsLinkClean[!duplicated(occDataSpsLinkClean[,c('longitude','latitude')]),]
##definindo dataset limpo para o SpsLink
occDataSpsLink = occDataSpsLinkClean
write.csv(occDataSpsLink, '/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/SpsLink_mar-2018/reduviidae_SpsLink_occDataClean.csv', row.names=FALSE)


##conjunto de dados consolidado para Reduviidae (GBIF + SpsLink)##


occDataGBIF = read.csv('/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/GBIF_mar-2018/Reduviidae_CSV/reduviidae_GBIF_occDataClean.csv', header=TRUE, sep=',', dec='.', stringsAsFactors=FALSE, na.strings="")

occDataSpsLink = read.csv('/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/SpsLink_mar-2018/reduviidae_SpsLink_occDataClean.csv', header=TRUE, sep=',', dec='.', stringsAsFactors=FALSE, na.strings="NA")


reduviidaeDataset = data.frame(sps = c(occDataGBIF$species, occDataSpsLink$scientificname),
                               lon = c(occDataGBIF$decimallongitude, occDataSpsLink$longitude),
                               lat = c(occDataGBIF$decimallatitude, occDataSpsLink$latitude))

write.csv(reduviidaeDataset, '/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/reduviidaeDataset.csv', row.names=FALSE, na="NA")


##KDE para vies amostral a partir do banco de dados de Reduviide##


##abrindo arquivo salvo
reduviidaeDataset = read.csv(file='/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/reduviidaeDataset.csv', header=TRUE, sep=',', dec='.', stringsAsFactors=FALSE, na.strings="")
reduviidaeOcc = reduviidaeDataset[,c('lon','lat')]

##convertendo de data.frame para staialPointsDataFrame
coordinates(reduviidaeOcc) <- ~lon+lat
proj4string(reduviidaeOcc) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')

##abrindo o shapefile que 'cortara' os pontos
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp",p4s=proj4string(reduviidaeOcc))

##ajustando o sistema de coordenadas geograficas entre pontos e shapefile
reduviidaeOcc <- spTransform(reduviidaeOcc, CRS(proj4string(AmSulShape)))
reduviidaeOcc <- reduviidaeOcc[AmSulShape, ]

##sjuatando para toda a area da america do sul e grid file final
#SOAextent = extent(-81.57551,-34.03384,-57.13385,12.99115)

SAbg = predictors[[1]]*0 ##America do Sul como 'pano de fundo'
crs(SAbg) = crs(raster())

##inspecionando pontos
plot(reduviidaeOcc)
plot(AmSulShape, add=TRUE)

##convertendo spatialPoints em data.frame para o KDE2D
reduviidaeOcc = as.data.frame(reduviidaeOcc)

##kernel density estimation com o pacote MASS
dens = MASS::kde2d(x=reduviidaeOcc$lon, y=reduviidaeOcc$lat, n=100)
densRas = raster(dens)
densRas = mask(x=densRas, mask=AmSulShape)

##ajustando projecao e fundindo com pano de fundo
densRas = projectRaster(densRas, crs=proj4string(SAbg), res=res(SAbg), method="bilinear")
densRas = merge(densRas, SAbg, tolerance=0.5)

##salvado raster
writeRaster(x=densRas, file='/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/reduviidaeBiasLayer.grd', overwrite=TRUE)

##inspecionando
plot(AmSulShape)
plot(densRas, add=TRUE)
points(reduviidaeOcc[,c('lon','lat')],cex=0.5)


##Formando o conjunto de dados de cada especie a ser modelada##


##abrindo conjunto de dados
reduviidaeDataset = read.csv(file='/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/reduviidaeDataset.csv', header=TRUE, sep=',', dec='.', stringsAsFactors=FALSE, na.strings="")

##lista com o nome das especies
spsNames = c("Panstrongylus_geniculatus",
             "Panstrongylus_lutzi",
             "Panstrongylus_megistus",
             "Rhodnius_nasutus",
             "Rhodnius_neglectus",
             "Rhodnius_pictipes",
             "Rhodnius_robustus",
             "Triatoma_brasiliensis",
             "Triatoma_infestans",
             "Triatoma_maculata",
             "Triatoma_pseudomaculata",
             "Triatoma_rubrovaria",
             "Triatoma_sordida",
             "Triatoma_vitticeps")

##garantindo ajuste de nomes
reduviidaeDataset$sps =  gsub(pattern=' ', replacement='_', x=reduviidaeDataset$sps)

##salvando dados por especie
for(sp_i in spsNames){
    spLines = agrep(pattern=sp_i, x=reduviidaeDataset$sps, value=FALSE)
    spDataset = reduviidaeDataset[spLines,]
    write.csv(spDataset, paste(spOccFolder,'/sps_occ/',sp_i,'.csv',sep=''), row.names=FALSE)
}



###SEGUNDA PARTE: trabalhando as variaveis climaticas



##passei pra dentro da proxima PARTE - racional para isso esta no paper do Soberon sobre o conjunto M - do BAM###



###TERCEIRA PARTE: rodando SDMclim (sem interacao com humanos) para as especies (e fazendo projecoes)###



##variaveis preditoras

##predictors = stack(list.files(paste(envVarFolder,'/presente/usadas',sep=''), pattern='*bio.*.asc', full.names=TRUE))
##crs(predictors) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
predictorsRaw = stack(list.files(path=paste(envVarFolder,'/presente/bio_2-5m_bil',sep=''), pattern='.bil', full.names=TRUE)) #modificar a extensao .bil de acordo com os arquivos
predictorsRaw = crop(predictorsRaw, SOAextent)

##Criando objeto com a lista dos nomes das especies
occ.sps <- list.files(paste(spOccFolder,'/sps_occ_Lucas',sep=''), pattern="csv")
splist <- unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))

##criando uma tabela vazia para salvador alguns dados
tabRes = data.frame()


##loop para SDM com cada uma das especies##


for (sp_i in splist){
    ##for (sp_i in splist[1:3]){
    tryCatch({

        ##diretorio base de trabalho
        setwd(paste(projectFolder,'/SDM outputs/resultados SDM sem humanos',sep=''))

        ##verifica e cria diretorio para salvar resultados da especie atual
        if (file.exists(sp_i)){
            setwd(sp_i)
        } else {
            dir.create(sp_i)
            setwd(sp_i)
        }
        
        ##dados de ocorrencia
        occPoints = read.csv(paste(spOccFolder,'/sps_occ_Lucas/',sp_i,'.csv',sep=''), header=TRUE, sep=',', dec='.', na.strings='',colClasses=c('character','numeric','numeric')) #abrindo pontos de ocorrencia
        names(occPoints) =  c('sp','lon','lat')
        occPoints = occPoints[,c('lon','lat')]

        ##minimo poligono convexo
        occMat = as.matrix(occPoints[,c('lon','lat')]) #sps occ matrix
                                        #SDMpredBIN = SDMpred > thre
                                        #crs(SDMpredBIN) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
        
        ##espacializacao dos pontos
        coordinates(occPoints) = ~lon+lat
        projection(occPoints) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

        ##centroide da distribuicao espacial
        SpsCentroid = kmeans(occMat, 1)
        SpsCentroid = SpsCentroid$centers

        ##convex hull
        cHull = dismo::convHull(occPoints)

        ##distancia entre pontos e centroide
        euc.dist <- function(x1) sqrt(sum((x1 - SpsCentroid) ^ 2))
        meanDistCentroid = mean( apply(occMat, 1, euc.dist) )

        ##buffer - contendo a hipotese para area de alcance da especie
        SpsBuffer = raster::buffer(polygons(cHull), width=meanDistCentroid)

        ##abrindo e cortando camadas de variaveis ambientais para o presente
        predictors = mask(predictorsRaw, SpsBuffer) #cortando para a hipotese de area de alcance da especie

        ##analisando correlacao das variaveis
        predictorsForVif = predictors

        vif(predictorsForVif)
        predictorsVif1 = vifcor(predictorsForVif, th=0.7)

        ##predictorsVif2 <- vifstep(predictorsForVif, th=10) # identify collinear variables that should be excluded
        ##predictorsVif2

        ##comparando
        ##predictorsVif1@results$Variables
        ##predictorsVif2@results$Variables

        ##definindo variaveis preditoras a serem usadas nos modelos
        predictors = predictorsForVif[[ grep(pattern=paste(predictorsVif1@results$Variables,collapse='|'), x=names(predictors), value=TRUE)  ]]

        ##verifica e cria diretorio para salvar as variaveis preditoras da especie atual
        if (!file.exists(paste(envVarFolder,'/presente/usadas/',sp_i,sep=''))){
            dir.create(paste(envVarFolder,'/presente/usadas/SDMclim/',sp_i,sep=''), recursive=TRUE)
        }

        writeRaster(x=predictors, filename=paste(envVarFolder,'/presente/usadas/SDMclim/',sp_i,'/predictors_',sp_i,'.asc',sep=''), overwrite=TRUE, bylayer=TRUE, suffix=names(predictors))
        
        ##pseudo-ausencia com o mesmo vies dos dados de ocorrencia
        currentBiasLayer = mask(biasLayer, SpsBuffer)
        bgPoints = dismo::randomPoints(mask=currentBiasLayer, n=10000, p=occPoints, prob=TRUE)
        bgPoints = data.frame(lon=bgPoints[,1], lat=bgPoints[,2])
        
        ##consolindando dados de presenca e background
        occPoints = as.data.frame(occPoints)
        dataSet = data.frame(lon=c(occPoints$lon, bgPoints$lon),
                             lat=c(occPoints$lat, bgPoints$lat),
                             occ=c(rep(1,nrow(occPoints)), rep(0,nrow(bgPoints))))

        dataSet[,c('lon','lat')] = round(dataSet[,c('lon','lat')], 2) #arredondando para garantir
        ##variaveis ambientais##
        dataSetVars = extract(x=predictors, y=dataSet[,c('lon','lat')], method='bilinear', na.rm=TRUE) #extraindo variaeis ambientais
        dataSet = data.frame(dataSet, dataSetVars) #juntando ao dataset
        dataSet = dataSet[complete.cases(dataSet),] #retirando dados errados
        dataSet = dataSet[!duplicated(dataSet[,c('lon','lat')]),] #retirando pontos numa mesma celula


        ##ENMeval##

        ENMblock = get.block(occ=dataSet[dataSet$occ==1,c('lon','lat')], bg.coords=dataSet[dataSet$occ==0,c('lon','lat')]) #dividindo em blocos

        ENMblockUser = get.user(occ.grp=ENMblock$occ.grp, bg.grp=ENMblock$bg.grp) #usuario, por causa dos dados proprios de bg-points gerados com vies

        SDMeval <- ENMevaluate(occ = dataSet[dataSet$occ==1,c('lon','lat')],
                               env = predictors,
                               bg.coords = dataSet[dataSet$occ==0,c('lon','lat')],
                               occ.grp = ENMblockUser$occ.grp,
                               bg.grp = ENMblockUser$bg.grp,
                               method = 'user',
                               RMvalues = c(0.5:5.5),
                               fc = c("L", "LQ", "LQH", "LQHP", "LQHPT"), #c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                               clamp = FALSE,
                               parallel = TRUE,
                               numCores = 2)

        
        ##melhor modelo
        bestModel = SDMeval@results[SDMeval@results$AICc==min(SDMeval@results$AICc, na.rm=TRUE),]
        bestModel = bestModel[complete.cases(bestModel),]
        bestModel = bestModel[1,]


        ## RMvalue = as.numeric(unlist(regmatches(bestModel$settings,
        ##                                        gregexpr("[[:digit:]]+\\.*[[:digit:]]*", bestModel$settings))))
        
        ##MAXENT##
        SDMmaxent = maxent(x = dataSet[,grep(pattern='bio',x=names(dataSet),value=TRUE)],
                           p = dataSet[,'occ'],
                           args = c(unlist(make.args(RMvalues=bestModel$rm, fc=bestModel$features, labels=FALSE)),'threads=2','doclamp=FALSE'))

        SDMpred = predict(predictors, SDMmaxent) #projecao espacial
        crs(SDMpred) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

        ## ##salvando um gridfile
        ## writeRaster(SDMpred, paste(sp_i,'Suitability.asc',sep=''))


        ##calculo do threshold##


        ##valores preditos pelo SDM em cada ponto do dataset
        pred =  extract(x=SDMpred,
                        y=dataSet[,c('lon','lat')],
                        method='bilinear',
                        na.rm=TRUE)
        threshDF = data.frame(occ=dataSet$occ, pred=pred) #juntando predicoes e observacoes
        threshDF = threshDF[complete.cases(threshDF),] #retirando possiveis NAs
        
        ##OBS.: omissao -> falso negativo -> false negative rate (FNR) = 1 - TPR (True Positive Rate)
        ##TPR = sensitividade
        
        myRoc = roc(predictor=threshDF$pred, response=threshDF$occ, positive=1) #curva ROC
        rocDF = data.frame( FNR = 1-myRoc$sensitivities, thresholds = myRoc$thresholds ) #data.frame com omissao e thresholds
        thre = max(rocDF[which(rocDF$FNR == max(rocDF[rocDF$FNR <= 0.05 ,]$FNR)),]$thresholds) #thresold de 5% omissao
        
        
        ##calculo do TSS para com o threshold -- OBS: TSS = sensitivity + specificity - 1
        
        
        TSSvalue = myRoc$sensitivities[which(myRoc$thresholds == thre)] + myRoc$specificities[which(myRoc$thresholds == thre)] - 1

        
        ##aplicando bufffer no mapa de suitability
        ##suitability continuo
        SAbg = SDMpred*0
        SDMpredRaw = mask(x=SDMpred, mask=SpsBuffer)
        SDMpred = merge(SDMpred, SAbg)
        writeRaster(SDMpred, paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/',sp_i,'/',sp_i,'Suitability.asc',sep=''), overwrite=TRUE)
        ##binario
        SAbg = SDMpred*0
        SDMpredBIN = SDMpred > thre
        ##SDMpredBINRaw = mask(x=SDMpredBIN, mask=SpsBuffer)
        SDMpredBIN = merge(SDMpredBIN, SAbg)
        writeRaster(SDMpredBIN, paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/',sp_i,'/',sp_i,'SuitabilityBIN.asc',sep=''), overwrite=TRUE)
        
        
        jpeg(paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/',sp_i,'/',sp_i,'_OccData.jpg',sep=''), width=1000, height=600)
        par(mfrow=c(1,2))
        ##occ e background
        plot(SAbg, main=paste(gsub('_',' ',sp_i)), cex=1.3, legend=FALSE)
        plot(AmSulShape, add=TRUE, col='white', lwd=0.3)
        plot(SAborders, add=TRUE)
        points(bgPoints, pch=19, cex=0.4, col=rgb(0.4,0.4,0.4,0.3))
        points(occPoints, pch=20, cex=1.1, col=rgb(1,0,0,0.6))
        plot(cHull, add=TRUE)
        plot(SpsBuffer, ad=TRUE, lty=2)
        box()
        grid()
        legend('bottomright', legend=c('Occurrences', 'Background points', 'Minimum convex polygon', 'Buffer'), pch=c(19,21,NA,NA), col=c('red','grey','black','black'), bg='white', lty=c(NA,NA,1,2))
        ##blocks
        plot(SAbg, main=paste(gsub('_',' ',sp_i)), cex=1.3, legend=FALSE)
        plot(AmSulShape, add=TRUE, col='white', lwd=0.3)
        plot(SAborders, add=TRUE)
        points(occPoints, pch=21, bg=ENMblock$occ.grp)
        box()
        grid()
        dev.off()
        
        
        ##mapa de distribuicao e suitability##
        
        
        jpeg(paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/',sp_i,'/',sp_i,'_Suitability.jpg',sep=''), width=1000, height=500)
        par(mfrow=c(1,2), mar=c(3,3,4,8))
        ##binario
        plot(SDMpred > thre, main=gsub('_', ' ', sp_i), cex=1.3, legend=FALSE)
        plot(AmSulShape, add=TRUE, lwd=0.3)
        plot(SAborders, add=TRUE)
        grid()
        legend('topright',legend=c('non-habitat','habitat'), pch='', fill=c('lightgrey','darkgreen'), bg='white')
        ##continuo
        plot(SDMpred, main=gsub('_', ' ', sp_i), cex=1.3)
        plot(AmSulShape, add=TRUE, lwd=0.3)
        plot(SAborders, add=TRUE)
        grid()
        dev.off()

        
        
        ##salvando output do ENMeval

        
        SDMevalOutput = as.data.frame(SDMeval@results)
        write.csv(SDMevalOutput, paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/',sp_i,'/',sp_i,'_SDMeval.csv',sep=''), row.names=FALSE)

        
        ##salvando graficos do ENMeval

        
        jpeg(paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/',sp_i,'/',sp_i,'_SDMeval.jpg',sep=''), width=800)
        par(mfrow=c(1,2), mar=c(5,5,10,5))
        eval.plot(SDMeval@results); title(paste(sp_i))
        eval.plot(SDMeval@results, 'Mean.AUC', var='Var.AUC'); title(paste(sp_i))
        dev.off()

        
        ##importancia das variaveis

        
        modelPars = SDMeval@models[[bestModel$settings]]
        write.csv(var.importance(modelPars), paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/',sp_i,'/',sp_i,'_variablesImportance.csv',sep=''), row.names=FALSE)

        
        ##salvando tabela de outputs##

        
        tabRes = rbind(tabRes,
                       data.frame(
                           sps = sp_i,
                           ENMsettings = bestModel$settings,
                           threshold = thre,
                           meanAUC = bestModel$Mean.AUC,
                           TSS = TSSvalue,
                           AICc = bestModel$AICc
                       ))
        
        write.csv(tabRes, paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/tabOutputsSDMclim.csv',sep=''), row.names=FALSE)
        gc()
        
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



###calculo do PARTIAL AUC####



#source('/home/anderson/R/R-Scripts/PartialROC.R')
source('D:/Anderson_Eduardo/Distribuicao de barbeiros com interacao com humanos/R scripts/partialROC.R')
library(sqldf)

##Criando objeto com a lista dos nomes das especies
occ.sps <- list.files(paste(spOccFolder,'/sps_occ_Lucas',sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))

##criando uma tabela vazia para salvador alguns dados
pAUCoutput = data.frame()


##loop para SDM com cada uma das especies##


for (sp_i in splist){
    ##for (sp_i in splist[1:3]){

    ##diretorio base de trabalho
    setwd(paste(projectFolder,'/SDM outputs/resultados SDM sem humanos',sep=''))
    
    ##camimhos dos arquivos
    occPointsPath = paste(spOccFolder,'/sps_occ_Lucas/',sp_i,'.csv',sep='') #pontos de ocorrencia
    suitabMapPath = paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/', sp_i, '/', sp_i,'Suitability.asc',sep='')
    
    
    ##PartialROC (PresenceFile, PredictionFile, OmissionVal, RandomPercent, NoOfIteration, OutputFile)
    pAUC = PartialROC(PresenceFile = occPointsPath,
                      PredictionFile = suitabMapPath,
                      OmissionVal = 0.05,
                      RandomPercent = 25,
                      NoOfIteration = 100,
                      OutputFile = paste(sp_i,'pAUC.csv'))

    dev.off()
    
    pAUCoutput = rbind(pAUCoutput,
                       data.frame(sps=sp_i,
                                  AUC_at_Value_0.05 = mean(pAUC$AUC_at_Value_0.05[-1]),
                                  AUC_at_0.5 = mean(pAUC$AUC_at_0.5[-1]),
                                  AUC_ratio = mean(pAUC$AUC_ratio[-1])) )

    write.csv(pAUCoutput, 'pAUCoutputSDMclim.csv', row.names=TRUE)

}



###QUARTA PARTE: gerando mapa de suitability acumulado (i.e., mapa de "riqueza") - SEM impacto humano###




##Criando objeto com a lista dos nomes das especies
occ.sps <- list.files(paste(spOccFolder,'/sps_occ_Lucas',sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))

##criando objeto para armazenas os mapas de suitability de cada especie
spsStackSDMclim = stack()

##loop para pegar os mapas de suitability de cada especie
for(sp_i in splist){
    mapSDMclim_sp_i = raster(paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/',sp_i,'/',sp_i,'Suitability.asc',sep=''))
    spsStackSDMclim = stack(c(spsStackSDMclim,mapSDMclim_sp_i))
}

##calculo do suitability acumylado a partir dos mapas de suitability de cada especie
SDMclimAcumSuit = sum(spsStackSDMclim, na.rm=TRUE)
writeRaster(x=SDMclimAcumSuit, paste(projectFolder,'/SDM outputs/SDMclimAcumSuit.asc', sep=''), overwrite=TRUE)

##graficos - OBS: a figura do mapa de suitabality acumulado para SDMclim sera criado no script para SDMhuman (nesta mesma secao daquele script)
