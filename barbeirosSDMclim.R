##SCRIPT PARA DISTRIBUICAO DE BARBEIROS (PRESENTE E FUTURO) - APENAS VARIAVEIS CLIMATICAS##

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
Sys.setenv(JAVA_HOME = "/usr/lib/jvm/java-7-openjdk-amd64")
options(java.parameters = "Xmx7g")


## Definindo parametros e variaveis globais
## envVarFolder = "J:/Lucas/Modelagem barbeiros/Variaveis Climaticas"
## spOccFolder = "J:/Lucas/Modelagem barbeiros/Ocorrencias"
## projectFolder = "J:/Lucas/Modelagem barbeiros/"
#anderson
envVarFolder = "/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Variaveis Climaticas"
spOccFolder = "/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias"
projectFolder = "/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos"
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #abrindo shape da America do Sul
SAborders = rgdal::readOGR('/home/anderson/PosDoc/shapefiles/continents/continent.shp') #bordas de continentes
SOAextent = extent(-81.57551,-34.03384,-57.13385,12.99115)
SAborders = crop(SAborders,SOAextent)
biasLayer = raster('/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/reduviidaeBiasLayer.grd')




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
reduviidaeDataset = read.csv(file='/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/reduviidaeDataset.csv', header=TRUE, sep=',', dec='.', stringsAsFactors=FALSE, na.strings="")

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



##abrindo e cortando camadas de variaveis ambientais para o presente
predictorsRaw <- stack(list.files(path=paste(envVarFolder,'/presente/bio_2-5m_bil',sep=''), pattern='.bil', full.names=TRUE)) #modificar a extensao .bil de acordo com os arquivos
predictors = crop(x=predictorsRaw, y=SOAextent)
predictors = mask(predictors, AmSulShape) #cortando para Am. do Sul

##abrindo e cortando camadas de variaveis ambientais para projecao
filesProjectionOtimistaRaw <- stack(list.files(path=paste(envVarFolder,"/futuro/cenario_otimista",sep=''), pattern='.asc', full.names=TRUE))
filesProjectionPessimistaRaw <- stack(list.files(path=paste(envVarFolder,"/futuro/cenario_pessimista",sep=''), pattern='.asc', full.names=TRUE)) 
filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul


##analisando correlacao das variaveis##


predictorsForVif = predictors

vif(predictorsForVif)
predictorsVif1 = vifcor(predictorsForVif, th=0.7)
predictorsVif1

predictorsVif2 <- vifstep(predictorsForVif, th=10) # identify collinear variables that should be excluded
predictorsVif2

##comparando
predictorsVif1@results$Variables
predictorsVif2@results$Variables

##definindo variaveis preditoras a serem usadas nos modelos
predictors = predictorsForVif[[ predictorsVif1@results$Variables ]]

writeRaster(x=predictors, filename=paste(envVarFolder,'/presente/usadas/predictors.asc',sep=''), overwrite=TRUE, bylayer=TRUE, suffix=names(predictors))



###TERCEIRA PARTE: rodando SDMs para as especies (e fazendo projecoes)###



##variaveis preditoras
predictors = stack(list.files(paste(envVarFolder,'/presente/usadas',sep=''), pattern='.asc', full.names=TRUE))
crs(predictors) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

##Criando objeto com a lista dos nomes das especies
occ.sps <- list.files(paste(spOccFolder,'/sps_occ_Lucas',sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))

##criando uma tabela vazia para salvador alguns dados
tabRes = data.frame()


##loop para SDM com cada uma das especies##


for (sp_i in splist){
##for (sp_i in splist[1:3]){
    tryCatch({
        

        ##diretorio base de trabalho
        setwd(paste(projectFolder,'/SDM outputs',sep=''))

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

        ##pseudo-ausencia com o mesmo vies dos dados de ocorrencia    
        bgPoints = dismo::randomPoints(mask=biasLayer, n=10000, p=occPoints, prob=TRUE)
        bgPoints = data.frame(lon=bgPoints[,1], lat=bgPoints[,2])

        
        ##consolindando dados de presenca e background
        dataSet = data.frame(lon=c(occPoints$lon, bgPoints$lon),
                             lat=c(occPoints$lat, bgPoints$lat),
                             occ=c(rep(1,nrow(occPoints)),rep(0,nrow(bgPoints))))

        dataSet[,c('lon','lat')] = round(dataSet[,c('lon','lat')], 2) #arredondando para garantir
        ##variaveis ambientais##
        dataSetVars = extract(x=predictors, y=dataSet[,c('lon','lat')], method='bilinear', na.rm=TRUE) #extraindo variaeis ambientais
        dataSet = data.frame(dataSet, dataSetVars) #juntando ao dataset
        dataSet = dataSet[complete.cases(dataSet),] #retirando dados errados
        dataSet = dataSet[!duplicated(dataSet[,c('lon','lat')]),] #retirando pontos numa mesma celula


        ##ENMeval##

        ENMblock = get.block(occ=dataSet[dataSet$occ==1,c('lon','lat')],
                             bg.coords=dataSet[dataSet$occ==0,c('lon','lat')]) #dividindo em blocos

        ENMblockUser = get.user(occ.grp=ENMblock$occ.grp,
                                bg.grp=ENMblock$bg.grp) #usuario, por causa dos dados proprios de bg-points gerados com vies
        

        SDMeval <- ENMevaluate(occ = dataSet[dataSet$occ==1,c('lon','lat')],
                               env = predictors,
                               bg.coords = dataSet[dataSet$occ==0,c('lon','lat')],
                               occ.grp = ENMblockUser$occ.grp,
                               bg.grp = ENMblockUser$bg.grp,
                               method = 'user',
                               RMvalues = c(0.5:5.5),
                               fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
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
                           args = c(unlist(make.args(RMvalues=bestModel$rm, fc=bestModel$features, labels=FALSE)),'threads=2'))

        SDMpred = predict(predictors, SDMmaxent) #projecao espacial

        ##salvando um gridfile
        writeRaster(SDMpred, paste(sp_i,'Suitability.asc',sep=''))


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


        ##minimo poligono convexo##



        for (sp_i in splist){

            ##diretorio base de trabalho
            setwd(paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/',sp_i,sep=''))
            
            ##dados de occ
            occPoints = read.csv(paste(spOccFolder,'/sps_occ_Lucas/',sp_i,'.csv',sep=''), header=TRUE, sep=',', dec='.', na.strings='',colClasses=c('character','numeric','numeric')) #abrindo pontos de ocorrencia
            names(occPoints) =  c('sp','lon','lat')
            occPoints = occPoints[,c('lon','lat')]
            occMat = as.matrix(occPoints[,c('lon','lat')]) #sps occ matrix
            SDMpred = raster(paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/',sp_i,'/',sp_i,'Suitability.asc',sep=''))
            crs(SDMpred) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
            tabOutputsSdm = read.csv('/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/SDM outputs/resultados SDM sem humanos/tabOutputs.csv', header=TRUE)
            thre = tabOutputsSdm[tabOutputsSdm$sps == sp_i, 'threshold']
            SDMpredBIN = SDMpred > thre
            crs(SDMpredBIN) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

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

            ##buffer - hipotese para area de alcance da especie
            SpsBuffer = raster::buffer(polygons(cHull), width=meanDistCentroid)
            
            ##aplicando bufffer no mapa de suitability
            ##suitability continuo
            SAbg = sdmPred*0
            SDMpredRaw = mask(x=SDMpred, mask=SpsBuffer)
            SDMpred = merge(SDMpredRaw, SAbg)
            writeRaster(SDMpred, paste(sp_i,'SuitabilityNOVO.asc',sep=''), overwrite=TRUE)
            ##binario
            SAbg = sdmPredBIN*0
            SDMpredBINRaw = mask(x=SDMpredBIN, mask=SpsBuffer)
            SDMpredBIN = merge(SDMpredBINRaw, SAbg)
            writeRaster(SDMpredBIN, paste(sp_i,'SuitabilityBINNOVO.asc',sep=''), overwrite=TRUE)

               
        
        
        ##mapas com a distribuicao dos pontos##

####helper
        ENMblock = get.block(occ=occMat,bg.coords=bgPoints[,c('lon','lat')]) #dividindo em blocos
###
        
        
        jpeg(paste(sp_i,'_mapsOfPoints.jpg'), width=1000, height=600)
        par(mfrow=c(1,2))
        ##occ e background
        plot(SAbg, main=paste(gsub('_',' ',sp_i)), cex=1.3, legend=FALSE)
        plot(SAborders, add=TRUE, col='white')
        points(bgPoints, pch=19, cex=0.4, col=rgb(0.4,0.4,0.4,0.3))
        points(occPoints, pch=20, cex=1.1, col=rgb(1,0,0,0.6))
        plot(cHull, add=TRUE)
        plot(SpsBuffer, ad=TRUE, lty=2)
        box()
        grid()
        legend('bottomright', legend=c('Occurrences', 'Background points', 'Minimum convex polygon', 'Buffer'), pch=c(19,21,NA,NA), col=c('red','grey','black','black'), bg='white', lty=c(NA,NA,1,2))
        ##blocks
        plot(SAbg, main=paste(gsub('_',' ',sp_i)), cex=1.3, legend=FALSE)
        plot(SAborders, add=TRUE, col='white')
        points(occPoints, pch=21, bg=ENMblock$occ.grp)
        box()
        grid()
        dev.off()
        
        
        ##mapa de distribuicao e suitability##
        
        
        jpeg(paste(sp_i,'_maps.jpg'), width=1000, height=500)
        par(mfrow=c(1,2), mar=c(3,3,4,8),  cex=1.3)
        ##binario
        plot(SDMpred > thre, main=gsub('_', ' ', sp_i), legend=FALSE)
        plot(SAborders, add=T)
        grid()
        legend('topright',legend=c('non-habitat','habitat'), pch='', fill=c('lightgrey','darkgreen'), bg='white')
        ##continuo
        plot(SDMpred,main=gsub('_', ' ', sp_i))
        plot(SAborders, add=TRUE)
        grid()
        dev.off()

            
        }
        



        
        ##salvando output do ENMeval

        
        SDMevalOutput = as.data.frame(SDMeval@results)
        write.csv(SDMevalOutput, paste(sp_i,'_SDMeval.csv',sep=''), row.names=FALSE)

        
        ##salvando graficos do ENMeval

        
        jpeg(paste(sp_i,'_SDMeval.jpg',sep=''), width=800)
        par(mfrow=c(1,2), mar=c(5,5,10,5))
        eval.plot(SDMeval@results); title(paste(sp_i))
        eval.plot(SDMeval@results, 'Mean.AUC', var='Var.AUC'); title(paste(sp_i))
        dev.off()

        
        ##importancia das variaveis

        
        modelPars = SDMeval@models[[bestModel$settings]]
        write.csv(var.importance(modelPars), paste(sp_i, '_variablesImportance.csv',sep=''), row.names=FALSE)

        
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

        write.csv(tabRes, paste(projectFolder,'/SDM outputs/tabOutputs.csv',sep=''), row.names=FALSE)

    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



###calculo do PARTIAL AUC####



source('/home/anderson/R/R-Scripts/PartialROC.R')
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
        setwd(paste(projectFolder,'/SDM outputs',sep=''))
        
        ##camimhos dos arquivos
        occPointsPath = paste(spOccFolder,'/sps_occ_Lucas/',sp_i,'.csv',sep='') #pontos de ocorrencia
        suitabMapPath = paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/', sp_i, '/', sp_i,'SuitabilityNOVO.asc',sep='')
        
        
        ##PartialROC (PresenceFile, PredictionFile, OmissionVal, RandomPercent, NoOfIteration, OutputFile)
        pAUC = PartialROC(PresenceFile = occPointsPath,
                          PredictionFile = suitabMapPath,
                          OmissionVal = 0.05,
                          RandomPercent = 50,
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







###QUARTA PARTE: implementado mapas de riqueza###



library(SSDM)
vignette('SSDM')



##Criando objeto com a lista dos nomes das especies
occ.sps <- list.files(paste(spOccFolder,'/sps_occ_Lucas',sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))

dataSet = data.frame()


for (sp_i in splist){
        
        ##dados de ocorrencia
        occPoints = read.csv(paste(spOccFolder,'/sps_occ_Lucas/',sp_i,'.csv',sep=''), header=FALSE, sep=',', dec='.', na.strings='',colClasses=c('character','numeric','numeric')) #abrindo pontos de ocorrencia
        names(occPoints) =  c('sps','lon','lat')
        occPoints$sps = sp_i
   
        ##pseudo-ausencia com o mesmo vies dos dados de ocorrencia    
        bgPoints = dismo::randomPoints(mask=biasLayer, n=10000, p=occPoints[,c('lon','lat')], prob=TRUE)
        bgPoints = data.frame(sps=sp_i, lon=bgPoints[,1], lat=bgPoints[,2])
        
        ##consolindando dados de presenca e background
        dataSet = rbind(dataSet,
                        data.frame(sps=c(occPoints$sps, bgPoints$sps),
                                   lon=c(occPoints$lon, bgPoints$lon),
                                   lat=c(occPoints$lat, bgPoints$lat),
                                   occ=c(rep(1,nrow(occPoints)),rep(0,nrow(bgPoints)))))

        dataSet[,c('lon','lat')] = round(dataSet[,c('lon','lat')], 2) #arredondando para garantir
        dataSet = dataSet[complete.cases(dataSet),] #retirando dados errados
        dataSet = dataSet[!duplicated(dataSet[,c('lon','lat')]),] #retirando pontos numa mesma celula

}


SSDM <- stack_modelling(c('MAXENT'),
                        dataSet[dataSet$occ==1,c('sps','lon','lat')],
                        predictors,
                        rep = 1,
                        Xcol = 'lon',
                        Ycol = 'lat',
                        Spcol = 'sps',
                        method = "pSSDM",
                        verbose = FALSE)

plot(SSDM@diversity.map, main = 'SSDM\nfor Cryptocarya genus\nwith CTA and SVM algorithms')




###INTERACAO ENTRE ESPECIES COM JSDM



install.packages(c('R2jags', 'MASS', 'MCMCpack', 'abind', 'random', 'mclust'))


coefsDF = data.frame()

for (sp_i in splist){
        
        ##dados de ocorrencia
        occPoints = read.csv(paste(spOccFolder,'/sps_occ_Lucas/',sp_i,'.csv',sep=''), header=FALSE, sep=',', dec='.', na.strings='',colClasses=c('character','numeric','numeric')) #abrindo pontos de ocorrencia
        names(occPoints) =  c('sps','lon','lat')
        occPoints$sps = sp_i
   
        ##pseudo-ausencia com o mesmo vies dos dados de ocorrencia    
        bgPoints = dismo::randomPoints(mask=biasLayer, n=10000, p=occPoints[,c('lon','lat')], prob=TRUE)
        bgPoints = data.frame(sps=sp_i, lon=bgPoints[,1], lat=bgPoints[,2])
        
        ##consolindando dados de presenca e background
        dataSet = data.frame(lon=c(occPoints$lon, bgPoints$lon),
                             lat=c(occPoints$lat, bgPoints$lat),
                             occ=c(rep(1,nrow(occPoints)),rep(0,nrow(bgPoints))))

        dataSet[,c('lon','lat')] = round(dataSet[,c('lon','lat')], 2) #arredondando para garantir
        ##variaveis ambientais##
        dataSetVars = extract(x=predictors, y=dataSet[,c('lon','lat')], method='bilinear', na.rm=TRUE) #extraindo variaeis ambientais
        dataSet = data.frame(dataSet, dataSetVars) #juntando ao dataset
        dataSet = dataSet[complete.cases(dataSet),] #retirando dados errados
        dataSet = dataSet[!duplicated(dataSet[,c('lon','lat')]),] #retirando pontos numa mesma celula


    ##equacao do modelo##
    modelEq = occ ~ predictors_bio10 +
        predictors_bio11 +
        predictors_bio12 +
        predictors_bio13 +
        predictors_bio14 +
        predictors_bio15 +
        predictors_bio16
    

    ##rodando SDM
    GLM <- glm(modelEq, family=binomial(link=logit), data=dataSet)


    ##projecao espacial
    rm(SDMpred)
    SDMpred =  predict(predictors, GLM, type='response')


    ##guardando coeficientes
    modelCoeffs = as.numeric(coef(GLM))
    coefsDF = rbind(coefsDF, modelCoeffs)
    names(coefsDF) = names(coef(GLM))
    

    ##threshold

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


    ##mapa binario
    SDMpredBIN =  SDMpred>thre


    ##empilhando
    SDMpredBINstack = stack(SDMpredBINstack, SDMpredBIN)

}


##arrumando os nomes das especies na pilha de mapas
names(SDMpredBINstack) = splist


##trandormando em dataframe
spsOccMat = as.data.frame(SDMpredBIN, xy=TRUE, na.rm=TRUE) #transformando raster em data.frame
names(spsOccMat) = c('lon', 'lat', splist)
spsOccMat[,splist] = as.integer(spsOccMat[,splist])



##matriz de coeficientes
coefsMat = as.matrix(coefsDF)


##matriz das preditoras
locs = spsOccMat[,c('lon','lat')]
varsMat = extract(predictors, locs, method='bilinear', na.rm=TRUE)
X = as.matrix(varsMat)


##alguns parametros para JSDM
n_env_vars <- nlayers(predictors)
n_sites <- nrow(spsOccMat)
n_species <- length(splist)

Occur <- as.matrix(spsOccMat[,splist])

#Occur <- as.matrix(spsOccMat[,c('AA','BB')])

##parametros mcmc 
n.chains <- 5
n.iter <- 1000 #150000
n.burn <- 0.2*1000  #110000
n.thin <- 40
df <- 1
model_name <- 'teste_teste'

source("/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/scriptSumplementarJSDM.R")

Diagnose(Beta, 'rhat')

Diagnose(Beta, 'effn')

TPLOT(Beta, 2, 2)

TPLOT(Rho, 1, 2)




#################################################################################################################




###TERCEIRA PARTE: gerando mapas de sobreposicao (i.e. mapa de riqueza) - SEM impacto humano###


#sobrepondo distribuicoes para mapa de riqueza
#presente
listaPresente = grep(list.files(paste(projectFolder,"Projecoes",sep=""),full.names=TRUE),pattern='Otimista|Pessimista',inv=T,value=T)
listaPresenteBIN = grep(listaPresente,pattern='BIN.asc',value=TRUE)
camadasPresente = stack(listaPresenteBIN)
mapaRiquezaPresente = sum(camadasPresente)
plot(mapaRiquezaPresente)
writeRaster(x=mapaRiquezaPresente,filename=paste(projectFolder,'Mapas de riqueza/mapaRiquezaPresente.asc',sep=''),overwrite=TRUE)

#futuro otimista
camadasFuturoOtimista = stack(list.files(paste(projectFolder,"Projecoes",sep=""),pattern="OtimistaBIN.asc",full.names=TRUE))
mapaRiquezaFuturoOtimista = sum(camadasFuturoOtimista)
plot(mapaRiquezaFuturoOtimista)
writeRaster(x=mapaRiquezaFuturoOtimista,filename=paste(projectFolder,'Mapas de riqueza/mapaRiquezaFuturoOtimista.asc',sep=''),overwrite=TRUE)

#futuro pessimista
camadasFuturoPessimista = stack(list.files(paste(projectFolder,"Projecoes",sep=""),pattern="PessimistaBIN.asc",full.names=TRUE))
mapaRiquezaFuturoPessimista = sum(camadasFuturoPessimista)
plot(mapaRiquezaFuturoPessimista)
writeRaster(x=mapaRiquezaFuturoPessimista,filename=paste(projectFolder,'Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''),overwrite=TRUE)


###QUARTA PARTE: gerando MAPAS DE RISCO - SEM impacto humano###


##indices para RISCO DE INFECCAO por especie de vetor (dados SUS)
##link: http://portalarquivos.saude.gov.br/images/pdf/2015/agosto/03/2014-020..pdf
tabBarb = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/Taxa de infeccao natural vetores 2007-2011.csv",header=TRUE)
#infecBarb = sort(tabBarb[,3],decreasing=TRUE)
#infecIndOrdered = c(infecBarb[6],infecBarb[3],infecBarb[5],infecBarb[9],infecBarb[10],infecBarb[11],infecBarb[13],infecBarb[1])
#infecInd = sort(1:length(infecBarb),decreasing=TRUE) #idice de infeccao para confeccao dos mapas

#sobrepondo distribuicoes para mapa de risco
#presente
listaPresente = grep(list.files(paste(projectFolder,"Projecoes",sep=""),full.names=TRUE),pattern='Otimista|Pessimista',inv=T,value=T)
listaPresenteBIN = grep(listaPresente,pattern='BIN.asc',value=TRUE)
camadasPresente = stack(listaPresenteBIN)
#
listaNomes = names(camadasPresente)
listaNomes = gsub(pattern='BIN',replacement='',x=listaNomes)
listaNomes = basename(listaNomes)
listaNomes = gsub(pattern='.asc',replacement='',x=listaNomes)
infecIndOrdered = tabBarb$taxaInfeccaonatural[match(listaNomes,tabBarb$sp)]/100 #taxa de infeccao natural na ordem dos rasters (de 0 a 1)
#
mapaRiscoPresente = sum(camadasPresente*infecIndOrdered) 
mapaRiscoPresente = mapaRiscoPresente/length(infecIndOrdered)
plot(mapaRiscoPresente)
writeRaster(x=mapaRiscoPresente,filename=paste(projectFolder,'Mapas de risco/mapaRiscoPresente.asc',sep=''),overwrite=TRUE)

#futuro otimista
camadasFuturoOtimista = stack(list.files(paste(projectFolder,"Projecoes",sep=""),pattern = "OtimistaBIN.asc",full.names=TRUE))
listaNomes = names(camadasFuturoOtimista)
listaNomes = gsub(pattern='OtimistaBIN',replacement='',x=listaNomes)
infecIndOrdered = tabBarb$taxaInfeccaonatural[match(listaNomes,tabBarb$sp)]/100 #taxa de infeccao natural na ordem dos rasters (de 0 a 1)
#
mapaRiscoFuturoOtimista = sum(camadasFuturoOtimista*infecIndOrdered)
mapaRiscoFuturoOtimista = mapaRiscoFuturoOtimista/length(infecIndOrdered)
plot(mapaRiscoFuturoOtimista)
writeRaster(x=mapaRiscoFuturoOtimista,filename=paste(projectFolder,'Mapas de risco/mapaRiscoFuturoOtimista.asc',sep=''),overwrite=TRUE)

#futuro pessimista
camadasFuturoPessimista = stack(list.files(paste(projectFolder,"Projecoes",sep=""),pattern = "PessimistaBIN.asc",full.names=TRUE))
listaNomes = names(camadasFuturoOtimista)
listaNomes = gsub(pattern='OtimistaBIN',replacement='',x=listaNomes)
infecIndOrdered = tabBarb$taxaInfeccaonatural[match(listaNomes,tabBarb$sp)]/100 #taxa de infeccao natural na ordem dos rasters (de 0 a 1)
#
mapaRiscoFuturoPessimista = sum(camadasFuturoPessimista*infecIndOrdered)
mapaRiscoFuturoPessimista = mapaRiscoFuturoPessimista/length(infecIndOrdered)
plot(mapaRiscoFuturoPessimista)
writeRaster(x=mapaRiscoFuturoPessimista,filename=paste(projectFolder,'Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''),overwrite=TRUE)


###QUINTA PARTE: estatisticas sumarias a partir dos mapas - SEM impacto humano###


##definindo a area do Brasil na America do Sul, para os mapas
areaBR = extent(-80.00635,-31.71555,-37.8679,8.156474) #extent do Brasil
AmSulBR = crop(AmSulShape, extent(areaBR))  #america do sul recortada para o BR

#presente 

##abrindo
mapaRiquezaPresente = raster(paste(projectFolder,'Mapas de riqueza/mapaRiquezaPresente.asc',sep=''))
mapaRiscoPresente = raster(paste(projectFolder,'Mapas de risco/mapaRiscoPresente.asc',sep=''))

##cortando para o Brasil
mapaRiquezaPresenteBR = mask(mapaRiquezaPresente, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiscoPresenteBR = mask(mapaRiscoPresente, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))

#tamanho da area do quartil superior (para riqueza e risco), para comparar os cenarios
hightRiqPres= mapaRiquezaPresenteBR > quantile(mapaRiquezaPresenteBR, 0.75,na.rm=TRUE) #raster quartil superior riqueza
hightRiscPres = mapaRiscoPresenteBR > quantile(mapaRiscoPresenteBR, 0.75,na.rm=TRUE) #raster quartil superior risco

percCelRiqPres = freq(hightRiqPres,value=1)/(freq(hightRiqPres,value=0)+freq(hightRiqPres,value=1)) #percentagem celulas no quartil sup.
percCelRiscPres = freq(hightRiscPres,value=1)/(freq(hightRiscPres,value=0)+freq(hightRiscPres,value=1)) #percentagem de celulas no quartil sup.

##correlacao entre riqueza e risco
rm(test)
test <- getValues(stack(hightRiqPres,hightRiscPres))
corPres <- as.data.frame(cor(test, use="complete.obs",method='spearman'))
##write.csv(cor.matrix,'cor_matrix.csv')

#futuro otimista 

##abrindo
mapaRiquezaFuturoOtimista = raster(paste(projectFolder,'Mapas de riqueza/mapaRiquezaFuturoOtimista.asc',sep=''))
mapaRiscoFuturoOtimista = raster(paste(projectFolder,'Mapas de risco/mapaRiscoFuturoOtimista.asc',sep=''))

##cortando para o BR
mapaRiquezaFuturoOtimistaBR = mask(mapaRiquezaFuturoOtimista, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiscoFuturoOtimistaBR = mask(mapaRiscoFuturoOtimista, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))

##tamanho da area do quartil superior (para riqueza e risco), para comparar os cenarios
hightRiqOtim= mapaRiquezaFuturoOtimistaBR > quantile(mapaRiquezaFuturoOtimistaBR, 0.75,na.rm=TRUE) #raster quartil superior riqueza
hightRiscOtim = mapaRiscoFuturoOtimistaBR > quantile(mapaRiscoFuturoOtimistaBR, 0.75,na.rm=TRUE) #raster quartil superior risco

percCelRiqOtim = freq(hightRiqOtim,value=1)/(freq(hightRiqOtim,value=0)+freq(hightRiqOtim,value=1)) #percentagem celulas no quartil sup.
percCelRiscOtim =  freq(hightRiscOtim,value=1)/(freq(hightRiscOtim,value=0)+freq(hightRiscOtim,value=1)) #percentagem de celulas no quartil sup.

##correlacao entre riqueza e risco
rm(test)
test <- getValues(stack(hightRiqOtim,hightRiscOtim))
corOtim <- as.data.frame(cor(test, use="complete.obs",method='spearman'))
##write.csv(cor.matrix,'cor_matrix.csv')

##futuro pessimista 

##abrindo
mapaRiquezaFuturoPessimista = raster(paste(projectFolder,'Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''))
mapaRiscoFuturoPessimista = raster(paste(projectFolder,'Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''))

##cortando para o BR
mapaRiquezaFuturoPessimistaBR = mask(mapaRiquezaFuturoPessimista, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiscoFuturoPessimistaBR = mask(mapaRiscoFuturoPessimista, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))

#tamanho da area do quartil superior (para riqueza e risco), para comparar os cenarios
hightRiqPess = mapaRiquezaFuturoPessimistaBR > quantile(mapaRiquezaFuturoPessimistaBR, 0.75,na.rm=TRUE) #raster quartil superior riqueza
hightRiscPess = mapaRiscoFuturoPessimistaBR > quantile(mapaRiscoFuturoPessimistaBR, 0.75,na.rm=TRUE) #raster quartil superior risco

percCelRiqPess = freq(hightRiqPess,value=1)/(freq(hightRiqPess,value=0)+freq(hightRiqPess,value=1)) #percentagem celulas no quartil sup.
percCelRiscPess = freq(hightRiscPess,value=1)/(freq(hightRiscPess,value=0)+freq(hightRiscPess,value=1))  #percentagem de celulas no quartil sup.

##correlacao entre riqueza e risco
rm(test)
test <- getValues(stack(hightRiqPess,hightRiscPess))
corPess <- as.data.frame(cor(test, use="complete.obs",method='spearman'))
##write.csv(cor.matrix,'cor_matrix.csv')

##salvando resultados em tabelas

rm(tam)

tab = data.frame(scenario=c('pres','fut_otim','fut_pess'),
                 quantile75riq = c(quantile(mapaRiquezaPresenteBR, 0.75,na.rm=TRUE),quantile(mapaRiquezaFuturoOtimistaBR, 0.75,na.rm=TRUE),quantile(mapaRiquezaFuturoPessimistaBR, 0.75,na.rm=TRUE)),
                 quantile75risc = c(quantile(mapaRiscoPresenteBR, 0.75,na.rm=TRUE),quantile(mapaRiscoFuturoOtimistaBR, 0.75,na.rm=TRUE), quantile(mapaRiscoFuturoPessimistaBR, 0.75,na.rm=TRUE)),
                 percCellRiq = c(percCelRiqPres,percCelRiqOtim,percCelRiqPess),
                 percCellRisc = c(percCelRiscPres,percCelRiscOtim,percCelRiscPess),
                 corRiqRisc = c(corPres[1,2],corOtim[1,2],corPess[1,2])
)

write.csv(tab,paste(projectFolder,'statsRes.csv',sep=''),row.names = FALSE)


###SEXTA PARTE: figuras dos mapas - SEM impacto humano###


##abrindo os rasters

mapaRiquezaPresente = raster(paste(projectFolder,'Mapas de riqueza/mapaRiquezaPresente.asc',sep=''))
mapaRiquezaFuturoOtimista = raster(paste(projectFolder,'Mapas de riqueza/mapaRiquezaFuturoOtimista.asc',sep=''))
mapaRiquezaFuturoPessimista = raster(paste(projectFolder,'Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''))
mapaRiscoPresente = raster(paste(projectFolder,'Mapas de risco/mapaRiscoPresente.asc',sep=''))
mapaRiscoFuturoOtimista = raster(paste(projectFolder,'Mapas de risco/mapaRiscoFuturoOtimista.asc',sep=''))
mapaRiscoFuturoPessimista = raster(paste(projectFolder,'Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''))

##definindo a area do Brasil na America do Sul, para os mapas

areaBR = extent(-80.00635,-31.71555,-37.8679,8.156474) #extent do Brasil
AmSulBR = crop(AmSulShape, extent(areaBR))  #america do sul recortada para o BR

##cortando para o Brasil
mapaRiquezaPresenteBR = mask(mapaRiquezaPresente, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiquezaFuturoOtimistaBR = mask(mapaRiquezaFuturoOtimista, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiquezaFuturoPessimistaBR = mask(mapaRiquezaFuturoPessimista, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiscoPresenteBR = mask(mapaRiscoPresente, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiscoFuturoOtimistaBR = mask(mapaRiscoFuturoOtimista, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiscoFuturoPessimistaBR = mask(mapaRiscoFuturoPessimista, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))


##salvado os mapas
jpeg(filename=paste(projectFolder,'mapas.jpeg',sep=''),width=1750,height=850)
par(mfrow=c(2,3), mar=c(5,5,5,20))
##riqueza
plot(crop(mapaRiquezaPresenteBR,areaBR),main='Current climate',legend=FALSE,cex.axis=2,cex.main=4) + plot(AmSulShape[AmSulShape$CNTRY_NAME!='Brazil',],col='lightgray',add=TRUE) + plot(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',],add=TRUE) + box() + grid()
plot(crop(mapaRiquezaFuturoOtimistaBR,areaBR),main='2070 optmistic',legend=FALSE,cex.axis=2,cex.main=4) + plot(AmSulShape[AmSulShape$CNTRY_NAME!='Brazil',],col='lightgray',add=TRUE) + plot(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',],add=TRUE) + box() + grid()
plot(crop(mapaRiquezaFuturoPessimistaBR,areaBR),main='2070 pessimistic',legend=FALSE,cex.axis=2,cex.main=4) + plot(AmSulShape[AmSulShape$CNTRY_NAME!='Brazil',],col='lightgray',add=TRUE) + plot(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',],add=TRUE) + box() + grid()
plot(mapaRiquezaFuturoPessimistaBR,legend.only=TRUE,legend.width=3,axis.args=list(cex.axis=2),legend.args=list(text='Species richness',font=2,side=4,line=4.5,cex=2.2,cex.axis=0.2)) #legenda
##risco
plot(crop(mapaRiscoPresenteBR,areaBR),main='Current climate',legend=FALSE,cex.axis=2,cex.main=4) + plot(AmSulShape[AmSulShape$CNTRY_NAME!='Brazil',],col='lightgray',add=TRUE) + plot(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',],add=TRUE) + box() + grid()
plot(crop(mapaRiscoFuturoOtimistaBR,areaBR),main='2070 optmistic',legend=FALSE,cex.axis=2,cex.main=4) + plot(AmSulShape[AmSulShape$CNTRY_NAME!='Brazil',],col='lightgray',add=TRUE) + plot(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',],add=TRUE) + box() + grid()
plot(crop(mapaRiscoFuturoPessimistaBR,areaBR),main='2070 pessimistic',legend=FALSE,cex.axis=2,cex.main=4) + plot(AmSulShape[AmSulShape$CNTRY_NAME!='Brazil',],col='lightgray',add=TRUE) + plot(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',],add=TRUE) + box() + grid()
plot(mapaRiscoFuturoPessimistaBR,legend.only=TRUE,legend.width=3,axis.args=list(cex.axis=2),legend.args=list(text='Risk of infected vector',font=2,side=4,line=6.5,cex=2.2,cex.axis=0.2)) #legenda
dev.off()
