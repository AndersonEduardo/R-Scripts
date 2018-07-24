##SCRIPT PARA DISTRIBUICAO DE BARBEIROS (SDMhuman) - VARIAVEIS CLIMATICAS + DENSIDADE HUMANA##

library(raster)
library(maptools)
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
envVarFolder = "/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/Variaveis Climaticas"
spOccFolder = "/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias"
projectFolder = "/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos"
##AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #abrindo shape da America do Sul
SAborders = rgdal::readOGR('/home/anderson/shapefiles/ne_50m_land/ne_50m_land.shp') #bordas de continentes
SOAextent = extent(-81.57551,-34.03384,-57.13385,12.99115)
SAborders = crop(SAborders,SOAextent)
biasLayer = raster('/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/reduviidaeBiasLayer.grd')
humDens = raster('/home/anderson/gridfiles/Densidade humana/gpw-v4-population-density-rev10_2015_2pt5_min_asc/gpw_v4_population_density_rev10_2015_2pt5_min.asc')




###SEGUNDA PARTE: trabalhando as variaveis climaticas



##abrindo e cortando camadas de variaveis ambientais para o presente
predictorsRaw <- stack(list.files(path=paste(envVarFolder,'/presente/bio_2-5m_bil',sep=''), pattern='.bil', full.names=TRUE)) #modificar a extensao .bil de acordo com os arquivos
predictors = crop(x=predictorsRaw, y=SOAextent)
predictors = mask(predictors, AmSulShape) #cortando para Am. do Sul

## ##abrindo e cortando camadas de variaveis ambientais para projecao
## filesProjectionOtimistaRaw <- stack(list.files(path=paste(envVarFolder,"/futuro/cenario_otimista",sep=''), pattern='.asc', full.names=TRUE))
## filesProjectionPessimistaRaw <- stack(list.files(path=paste(envVarFolder,"/futuro/cenario_pessimista",sep=''), pattern='.asc', full.names=TRUE)) 
## filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul


##ajustando dados de densidade humana e dados ambietais##


SAbg = predictors[[1]]*0 ##America do Sul como 'pano de fundo'
crs(SAbg) = crs(raster())

humDens = crop(humDens, extent(SAbg))

##ajustando projecao e fundindo com pano de fundo
humDens = projectRaster(humDens, crs=proj4string(SAbg), res=res(SAbg), method="bilinear")
humDens = merge(humDens, SAbg, tolerance=0.3)
humDens = crop(humDens, extent(SAbg))

##salvando raster
writeRaster(x=humDens, file=paste(envVarFolder,'/presente/usadas/humDens.asc',sep=''), overwrite=TRUE)


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
predictors = predictorsForVif[[ as.character(predictorsVif1@results$Variables) ]]

writeRaster(x=predictors, filename=paste(envVarFolder,'/presente/usadas/predictors.asc',sep=''), overwrite=TRUE, bylayer=TRUE, suffix=names(predictors))



###TERCEIRA PARTE: rodando SDMs para as especies (e fazendo projecoes)###



##variaveis preditoras
predictors = stack(list.files(paste(envVarFolder,'/presente/usadas',sep=''), pattern='.asc', full.names=TRUE))
crs(predictors) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

## ##densidade humana
## humDens = raster(paste(envVarFolder,'/presente/usadas/humDens.asc',sep=''))
## crs(humDens) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

## ##adicionando humanos nos preditoras
## varNames = names(predictors)
## predictors = stack(predictors,humDens)
## names(predictors) = c(varNames,'humDens')


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

        ENMblock = get.block(occ=dataSet[dataSet$occ==1,c('lon','lat')], bg.coords=dataSet[dataSet$occ==0,c('lon','lat')]) #dividindo em blocos
        
        ENMblockUser = get.user(occ.grp=ENMblock$occ.grp, bg.grp=ENMblock$bg.grp) #usuario, por causa dos dados proprios de bg-points gerados com vies
                                

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
        SDMmaxent = maxent(x = dataSet[,grep(pattern='bio|humDens',x=names(dataSet),value=TRUE)],
                           p = dataSet[,'occ'],
                           args = c(unlist(make.args(RMvalues=bestModel$rm, fc=bestModel$features, labels=FALSE)),'threads=2'))

        SDMpred = predict(predictors, SDMmaxent) #projecao espacial


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


        occMat = as.matrix(occPoints[,c('lon','lat')]) #sps occ matrix
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
        SAbg = SDMpred*0
        SDMpredRaw = mask(x=SDMpred, mask=SpsBuffer)
        SDMpred = merge(SDMpredRaw, SAbg)
        writeRaster(SDMpred, paste(sp_i,'SuitabilityNOVO.asc',sep=''), overwrite=TRUE)
        ##binario
        SAbg = SDMpredBIN*0
        SDMpredBINRaw = mask(x=SDMpredBIN, mask=SpsBuffer)
        SDMpredBIN = merge(SDMpredBINRaw, SAbg)
        writeRaster(SDMpredBIN, paste(sp_i,'SuitabilityBINNOVO.asc',sep=''), overwrite=TRUE)
        
        
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
    suitabMapPath = paste(projectFolder,'/SDM outputs/resultados SDM com humanos/', sp_i, '/', sp_i,'SuitabilityNOVO.asc',sep='')
    
    
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

    write.csv(pAUCoutput, 'pAUCoutputSDMhuman.csv', row.names=TRUE)

}




##TERCEIRA PARTE: comparando SDMs com e sem interacao com humanos##




##tabela de resultados
tabOutputsSdmH = read.csv('/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/SDM outputs/resultados SDM com humanos/tabOutputsSDMhuman.csv', header=TRUE)
tabOutputsSdm = read.csv('/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/SDM outputs/resultados SDM sem humanos/tabOutputsSDMclim.csv', header=TRUE)

##Criando objeto com a lista dos nomes das especies
occ.sps <- list.files(paste(spOccFolder,'/sps_occ_Lucas',sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))

##criando uma tabela vazia para salvador alguns dados
tabModelComparation = data.frame()


##loop para SDM com cada uma das especies##


for (sp_i in splist){

    ##diretorio base de trabalho
    setwd(paste(projectFolder,'/SDM outputs',sep=''))
    
    ##thresholds
    threshSdm = tabOutputsSdm[tabOutputsSdm$sps == sp_i, 'threshold']
    threshSdmH = tabOutputsSdmH[tabOutputsSdmH$sps == sp_i, 'threshold']

    ##abrindo mapa de projecao
    sdmPred = raster(paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/', sp_i, '/', sp_i,'Suitability.asc',sep=''))
    sdmHPred = raster(paste(projectFolder,'/SDM outputs/resultados SDM com humanos/', sp_i, '/', sp_i,'Suitability.asc',sep=''))

    ##mapas binarios
    sdmPredBIN = sdmPred > threshSdm
    sdmHPredBIN = sdmHPred > threshSdmH

    ##salvando .asc dos mapas de sobreposicao
    ovelapMap = sdmPredBIN + 2*sdmHPredBIN
    writeRaster(x=ovelapMap, filename=paste(projectFolder,'/SDM outputs/',sp_i,'_Overlap.asc',sep=''), overwrite=TRUE)
    
    ##plot da sobreposicao dos mapas binarios
    jpeg(paste(projectFolder,'/SDM outputs/graficos/',sp_i,'mapBinSobrepos.jpeg',sep=''))
    plot(sdmPredBIN + 2*sdmHPredBIN, col=c('white','blue','red','purple'), legend=FALSE, main=gsub('_',' ',sp_i))
    plot(SAborders, add=TRUE)
    grid()
    legend('topright', legend=c(expression(SDM[clim]), expression(SDM[human]), 'Overlapping'), pch=c(21,21,21), pt.bg=c('blue','red','purple'), bty='n', cex=1)
    dev.off()


    ##diferenca pixel-by-pixel entre mapas de suitability continuos
    sdmDiff = sdmHPred - sdmPred
    writeRaster(x=sdmDiff, filename=paste(projectFolder,'/SDM outputs/',sp_i,'_SDMdiff.asc',sep=''), overwrite=TRUE)

    ##plot da diferenca pixel-by-pixel entre mapas de suitability continuos
    minSuit =  sdmDiff@data@min
    maxSuit = sdmDiff@data@max
    maxSuitAbs = 1  #max(abs(c(minSuit,maxSuit)))
    rangeSuit = 2*maxSuitAbs / 21
    
    breaks <- seq(-1*maxSuitAbs,maxSuitAbs, by=rangeSuit)
    cols <- colorRampPalette(c("darkblue", 'blue', "white", 'red', "darkred"))(length(breaks)-1)

    jpeg(paste(projectFolder,'/SDM outputs/graficos/',sp_i,'SDMdiff.jpeg',sep=''))    
    plot(sdmDiff, col=cols, breaks=breaks,
         main=gsub('_',' ',sp_i), legend.mar=5,
         lab.breaks=c(round(min(breaks),3),rep('',length(breaks)-2),round(max(breaks),3)))
    plot(SAborders, add=TRUE)
    grid()
    dev.off()


    ##comparacoes a partir dos mapas binarios
    
    MapBinJoin = sdmPredBIN + 2*sdmHPredBIN
    MapBinJoinTab = as.data.frame(freq(MapBinJoin))

    if(exists('cells1')){rm(cells1)}; if(exists('cells2')){rm(cells2)}; if(exists('cells3')){rm(cells3)}
    cells1 = MapBinJoinTab[MapBinJoinTab$value == 1 & !is.na(MapBinJoinTab$value), 'count'] ##SDMclim
    cells2 = MapBinJoinTab[MapBinJoinTab$value == 2 & !is.na(MapBinJoinTab$value), 'count'] ##SDMhuman
    cells3 = MapBinJoinTab[MapBinJoinTab$value == 3 & !is.na(MapBinJoinTab$value), 'count'] ##sobreposicao

    sdmClim = max(0, cells1/sum(c(cells1,cells2,cells3), na.rm=TRUE), na.rm=TRUE)
    sdmHuman = max(0, cells2/sum(c(cells1,cells2,cells3), na.rm=TRUE), na.rm=TRUE)
    sdmOverlap = max(0, cells3/sum(c(cells1,cells2,cells3), na.rm=TRUE), na.rm=TRUE)

    ##comparacoes a partir da diferencia dos mapas de suitability

    sdmDiffMat = getValues(sdmDiff)
    ##celulas em que SDMhuman > SDMclim
    meanSuitSDMhumanGreater =  mean(sdmDiffMat[sdmDiffMat > 0 & !is.na(sdmDiffMat)], na.rm=TRUE)
    medianSuitSDMhumanGreater = median(sdmDiffMat[sdmDiffMat > 0 & !is.na(sdmDiffMat)], na.rm=TRUE)
    sdSuitSDMhumanGreater = var(sdmDiffMat[sdmDiffMat > 0 & !is.na(sdmDiffMat)], na.rm=TRUE)
    ##celulas em que SDMhuman < SDMclim
    meanSuitSDMclimGreater = mean(sdmDiffMat[sdmDiffMat < 0 & !is.na(sdmDiffMat)], na.rm=TRUE)
    medianSuitSDMclimGreater = median(sdmDiffMat[sdmDiffMat < 0 & !is.na(sdmDiffMat)], na.rm=TRUE)
    sdSuitSDMclimGreater = var(sdmDiffMat[sdmDiffMat < 0 & !is.na(sdmDiffMat)], na.rm=TRUE)

    ##correlacao espacial
    sdmPredStack = stack(sdmHPred, sdmPred)
    sdmPredStackMat = getValues(sdmPredStack)
    dimnames(sdmPredStackMat)[[2]] = c('SDMhuman','SDMclim')
    SDMcor = cor(sdmPredStackMat, use="complete.obs", method='pearson')

    ##proporcao de celulas para diferenca nos mapas de suitabilidade
    if(exists('cells1')){rm(cells1)}; if(exists('cells2')){rm(cells2)}; if(exists('cells3')){rm(cells3)}
    cells1 = as.numeric( freq(sdmDiff > 0)[2,2] )  ##area onde SDMhuman > SDMclim
    cells2 = as.numeric( freq(sdmDiff < 0)[2,2] ) ##area onde SDMhuman < SDMclim

    sdmDiffHuman = max(0, cells1/sum(c(cells1, cells2),na.rm=TRUE),na.rm=TRUE)
    sdmDiffClim = max(0,cells2/sum(c(cells1, cells2),na.rm=TRUE),na.rm=TRUE)


    ##colocando na tabela


    tabModelComparation = rbind(tabModelComparation,
                                data.frame(sps = sp_i,
                                           binaryMaps.proportion.SDMclim = sdmClim,
                                           binaryMaps.proportion.SDMhuman = sdmHuman,
                                           binaryMaps.proportion.Overlap = sdmOverlap,
                                           suitabilityDiff.proportion.SDMhuman.greater = sdmDiffHuman,
                                           suitabilityDiff.proportion.SDMclim.greater = sdmDiffClim,
                                           mean.suitability.SDMhuman.greater = meanSuitSDMhumanGreater,
                                           median.suitability.SDMhuman.greater = medianSuitSDMhumanGreater,
                                           sd.suitability.SDMhuman.greater = sdSuitSDMhumanGreater,
                                           mean.suitability.SDMclim.greater = meanSuitSDMclimGreater,
                                           median.suitability.SDMclim.greater = medianSuitSDMclimGreater,
                                           sd.suitability.SDMclim.greater = sdSuitSDMclimGreater,
                                           SDMcor = SDMcor[1,2]))

    write.csv(tabModelComparation, paste(projectFolder,'/SDM outputs/tabModelComparation.csv',sep=''), row.names=FALSE)
    
}




##QUARTA PARTE: graficos comparando SDMs com e sem interacao com humanos##




##tabela
tabModelComparation = read.csv(paste(projectFolder,'/SDM outputs/tabModelComparation.csv',sep=''), header=TRUE)


##barplot da sobreposicao dos mapas binarios
tabela = tabModelComparation[,c('binaryMaps.proportion.SDMclim', 'binaryMaps.proportion.SDMhuman', 'binaryMaps.proportion.Overlap')]
rownames(tabela) = tabModelComparation$sps
tabela = tabela[order(tabela$binaryMaps.proportion.Overlap, decreasing=TRUE),]
tabela = t(tabela)

jpeg(paste(projectFolder,'/SDM outputs/graficos/BoxplotBinMaps.jpg',sep=''), 600, 600)
par(cex=1.3, mar=c(14,4,4,4))
barplot(as.matrix(tabela), las=2, col=c('blue','red','purple'), names.arg=gsub('_',' ', colnames(tabela)))
abline(h=0)
dev.off()


##barplot das proporcoes (espacialmente) entre da diferenca entre SDmclim e SDMhuman
tabela = tabModelComparation[,c('suitabilityDiff.proportion.SDMhuman.greater', 'suitabilityDiff.proportion.SDMclim.greater')]
rownames(tabela) = tabModelComparation$sps
tabela = tabela[order(tabela$suitabilityDiff.proportion.SDMhuman.greater, decreasing=TRUE),]
tabela = t(tabela)

jpeg(paste(projectFolder,'/SDM outputs/graficos/BoxplotSDMdiffProportions.jpg',sep=''),600,600)
par(cex=1.3, mar=c(14,4,4,4))
barplot(as.matrix(tabela), las=2, col=c('red','blue'), names.arg=gsub('_',' ', colnames(tabela)))
abline(h=0)
dev.off()


##barplot da importancia das vairaveis
for (sp_i in splist){ #loop para coletar os dados salvos
    ##SDMclim
    tabelaSDMclim = read.csv(paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/', sp_i,'/',sp_i,'_variablesImportance.csv',sep=''), header=TRUE)
    ##SDMhuman
    tabelaSDMhuman = read.csv(paste(projectFolder,'/SDM outputs/resultados SDM com humanos/', sp_i,'/',sp_i,'_variablesImportance.csv',sep=''), header=TRUE)
    hvar = colSums(tabelaSDMhuman[c(1:2),-1])
    hvar = data.frame(variable='humDens', percent.contribution = as.numeric(hvar[1]), permutation.importance = as.numeric(hvar[2]))
    tabelaSDMhuman = tabelaSDMhuman[-c(1:2),]
    tabelaSDMhuman = rbind(hvar, tabelaSDMhuman)

    if(match(sp_i, splist)==1){
        varImpSDMclim = data.frame(variable=tabelaSDMclim$variable)
        varImpSDMhuman = data.frame(variable=tabelaSDMhuman$variable)
    }
    
    varImpSDMclim = cbind(varImpSDMclim, sps=tabelaSDMclim[,'permutation.importance'])
    varImpSDMhuman = cbind(varImpSDMhuman, sps=tabelaSDMhuman[,'permutation.importance'])

}

names(varImpSDMclim) = c('variable', splist) #ajustando os nomes
names(varImpSDMhuman) = c('variable', splist) #ajustando os nomes


##produzindo e salvando grafico
jpeg(paste(projectFolder,'/SDM outputs/graficos/varImportanceSDMclim.jpg',sep=''),600,600)
par(mar=c(13,3,1,6))
barplot(as.matrix(varImpSDMclim[,-1]), col=rainbow(nrow(varImpSDMclim)), las=2, names.arg=gsub('_',' ',names(varImpSDMclim)[-1]))
abline(h=0)
par(xpd=TRUE)
legend(18,100, gsub('predictors_bio', 'Bio ',  as.character(varImpSDMclim$variable)), pch=15, col=rainbow(nrow(varImpSDMclim)))
dev.off()

jpeg(paste(projectFolder,'/SDM outputs/graficos/varImportanceSDMhuman.jpg',sep=''),650,600)
par(mar=c(13,3,1,10))
barplot(as.matrix(varImpSDMhuman[,-1]), col=rainbow(nrow(varImpSDMhuman)), las=2, names.arg=gsub('_',' ',names(varImpSDMhuman)[-1]))
abline(h=0)
par(xpd=TRUE)
legend(18,100, c('Human density',gsub('predictors_bio', 'Bio ',  as.character(varImpSDMhuman$variable))[-1]), pch=15, col=rainbow(nrow(varImpSDMhuman)))
dev.off()


##plot para os mapas de suitability acumulado - ambos SDMclim e SDMhuman
spsStackSDMhuman = stack()

for(sp_i in splist){ #loop para coletar dados salvos
    mapSDMhuman_sp_i = raster(paste(projectFolder,'/SDM outputs/resultados SDM com humanos/',sp_i,'/',sp_i,'Suitability.asc',sep=''))
    spsStackSDMhuman = stack(c(spsStackSDMhuman,mapSDMhuman_sp_i))
}

##calculo do suitability acumulado a partir dos mapas de suitability de cada especie
SDMhumanAcumSuit = sum(spsStackSDMhuman, na.rm=TRUE)
writeRaster(x=SDMhumanAcumSuit, paste(projectFolder,'/SDM outputs/SDMhumanAcumSuit.asc', sep=''), overwrite=TRUE)

##abrindo gridfiles de suitability acumulado
SDMclimRich = paste(projectFolder,'/SDM outputs/SDMclimAcumSuit.asc', sep='')
SDMhumanRich = raster(paste(projectFolder,'/SDM outputs/SDMhumanAcumSuit.asc', sep=''))


##graficos
breaks = c(seq(-100,-0.1,length=1), seq(-0.1,0.1,length=1), seq(0.1,100,length=100))
cols <- colorRampPalette(c('white', 'yellow', 'orange', 'red', 'darkred'))(length(breaks)-1)

jpeg(paste(projectFolder,'/SDM outputs/graficos/','accumSuitab.jpeg',sep=''), 700,500)
par(mfrow=c(1,2))
##SDMclim
plot(SDMclimRich, col=cols, main=expression(SDM[clim]))
plot(SAborders, add=TRUE)
grid()
##SDMhuman
plot(SDMhumanRich, col=cols, main=expression(SDM[human]))
plot(SAborders, add=TRUE)
grid()
dev.off()





########################teste regressao espacial###########################


library(nlme)


##Criando objeto com a lista dos nomes das especies
occ.sps <- list.files(paste(spOccFolder,'/sps_occ_Lucas',sep=''),pattern="csv")
splist <- unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))

##human density data (gridfile)
humDens = raster('/home/anderson/gridfiles/Densidade humana/gpw-v4-population-density-rev10_2015_2pt5_min_asc/gpw_v4_population_density_rev10_2015_2pt5_min.asc')

##bias layer
biaslayer = raster('/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/reduviidaeBiasLayer.grd')

##Reduviidae data
reduviidaeDataset = read.csv(file='/home/anderson/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/reduviidaeDataset.csv', header=TRUE, sep=',', dec='.', stringsAsFactors=FALSE, na.strings="")
reduviidaeOcc = reduviidaeDataset[,c('lon','lat')]

##objeto para o output
output = data.frame()

for(sp_i in splist){

    ##dados de ocorrencia
    occPoints = read.csv(paste(spOccFolder,'/sps_occ_Lucas/',sp_i,'.csv',sep=''), header=TRUE, sep=',', dec='.', na.strings='',colClasses=c('character','numeric','numeric')) #abrindo pontos de ocorrencia
    names(occPoints) =  c('sp','lon','lat')
    occPoints = occPoints[,c('lon','lat')]

    ##produzindo background points com e sem correcao de vies amostral
    bgPtsRand = dismo::randomPoints(mask=biaslayer, n=nrow(occPoints), p=occPoints, prob=FALSE)
    bgPtsRand = data.frame(lon=bgPtsRand[,1], lat=bgPtsRand[,2])
    bgPtsCorr = dismo::randomPoints(mask=biaslayer, n=nrow(occPoints), p=occPoints, prob=TRUE)
    bgPtsCorr = data.frame(lon=bgPtsCorr[,1], lat=bgPtsCorr[,2])

    ##consolindando dados de presenca e background (para background SEM CORRECAO de vies amostral)
    dataSetRand = data.frame(lon=c(occPoints$lon, bgPtsRand$lon),
                             lat=c(occPoints$lat, bgPtsRand$lat),
                             occ=c(rep(1,nrow(occPoints)),rep(0,nrow(bgPtsRand))))
    
    dataSetRand[,c('lon','lat')] = round(dataSetRand[,c('lon','lat')], 2) #arredondando para garantir

    ##consolindando dados de presenca e background (para background COM CORRECAO de vies amostral)
    dataSetCorr = data.frame(lon=c(occPoints$lon, bgPtsCorr$lon),
                             lat=c(occPoints$lat, bgPtsCorr$lat),
                             occ=c(rep(1,nrow(occPoints)),rep(0,nrow(bgPtsCorr))))

    dataSetCorr[,c('lon','lat')] = round(dataSetCorr[,c('lon','lat')], 2) #arredondando para garantir

    ##dados de densidade humana nos pontos
    dhRand = extract(x=humDens, y=dataSetRand[,c('lon','lat')], buffer=10000, fun=mean)
    dataSetRand$hd = dhRand
    dhCorr = extract(x=humDens, y=dataSetCorr[,c('lon','lat')], buffer=10000, fun=mean)
    dataSetCorr$hd = dhCorr

    
    ##modelos##

    glmRand = glm(occ ~ log(hd), data=dataSetRand, family= binomial(link = "logit"))
    glmCorr = glm(occ ~ log(hd), data=dataSetCorr, family= binomial(link = "logit"))

    output = rbind(output,
                   data.frame(sps = rep(sp_i,2),
                              Background_points = c('random','biased'),
                              coeficients = c(as.numeric(coef(glmRand)[2]), as.numeric(coef(glmCorr)[2])),
                              Null_deviance = c(glmRand$null.deviance, glmCorr$null.deviance),
                              Res_deviance = c(glmRand$deviance, glmCorr$null.deviance),
                              AIC = c(glmRand$aic, glmCorr$aic)))
    
}

write.csv(output, '/home/anderson/Área de Trabalho/output.csv', row.names=FALSE)
    

################################################################################
##################### FIM DA PARTE NOVA ########################################
################################################################################



    occRaster = biaslayer*0
    tab = table(cellFromXY(occRaster, occPoints))
    occRaster[as.numeric(names(tab))] <- tab

    plot(occRaster)

    nSampled = extract(x=occRaster, y=occPoints)
    hd = extract(x=humDens, y=occPoints)

    plot(nSampled ~ log(hd))

    modelo = glm(nSampled ~ log(hd), family=poisson(link = "log"))
    summary(modelo)

    dad = data.frame(lon=c(1,1,1.7,1.9), lat=c(4.5,4.5,3,2.4), val=1:4)
    xx = aggregate(dad, by=list('lon','lat'), FUN=max)

    library(raster)
    r <- raster(xmn=0, ymn=0, xmx=10, ymx=10, res=1)
    r[] <- 0
    xy <- spsample(as(extent(r), 'SpatialPolygons'), 100, 'random')
    tab <- table(cellFromXY(r,xy))
    r[as.numeric(names(tab))] <- tab


