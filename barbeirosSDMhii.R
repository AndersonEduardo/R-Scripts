##SCRIPT PARA DISTRIBUICAO DE BARBEIROS (PRESENTE E FUTURO) - APENAS VARIAVEIS CLIMATICAS##

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
envVarFolder = "/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Variaveis Climaticas"
spOccFolder = "/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias"
projectFolder = "/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos"
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #abrindo shape da America do Sul
SAborders = rgdal::readOGR('/home/anderson/PosDoc/shapefiles/continents/continent.shp') #bordas de continentes
SOAextent = extent(-81.57551,-34.03384,-57.13385,12.99115)
SAborders = crop(SAborders,SOAextent)
biasLayer = raster('/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/reduviidaeBiasLayer.grd')
humDens = raster('/home/anderson/PosDoc/dados_ambientais/Densidade humana/gpw-v4-population-density-rev10_2015_2pt5_min_asc/gpw_v4_population_density_rev10_2015_2pt5_min.asc')




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

##densidade humana
humDens = raster(paste(envVarFolder,'/presente/usadas/humDens.asc',sep=''))
crs(humDens) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

##adicionando humanos nos preditoras
varNames = names(predictors)
predictors = stack(predictors,humDens)
names(predictors) = c(varNames,'humDens')


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
        occPoints = read.csv(paste(spOccFolder,'/sps_occ_Lucas/',sp_i,'.csv',sep=''), header=FALSE, sep=',', dec='.', na.strings='',colClasses=c('character','numeric','numeric')) #abrindo pontos de ocorrencia
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
        
        
        ##mapas com a distribuicao dos pontos##
        
        
        jpeg(paste(sp_i,'_mapsOfPoints.jpg'), width=1000, height=600)
        par(mfrow=c(1,2), cex=1.3)
        ##occ e background
        plot(predictors[[1]]*0, main=paste(gsub('_',' ',sp_i)), legend=FALSE)
        plot(SAborders, cex=1.3, add=TRUE, col='white')
        points(bgPoints, pch=19, cex=0.4, col=rgb(0.4,0.4,0.4,0.3))
        points(occPoints, pch=20, cex=1.1, col=rgb(1,0,0,0.6))
        box()
        grid()
        legend('topright', legend=c('occurrences', 'background points'), fill=c('red','grey'), bg='white')
        ##blocks
        plot(predictors[[1]]*0, main=paste(gsub('_',' ',sp_i)), legend=FALSE)
        plot(SAborders, cex=1.3, add=TRUE, col='white')
        points(SDMeval@occ.pts, pch=21, bg=SDMeval@occ.grp)
        box()
        grid()
        dev.off()
        
        
        ##mapa de distribuicao e suitability##
        
        
        jpeg(paste(sp_i,'_maps.jpg'), width=1000, height=500)
        par(mfrow=c(1,2), mar=c(3,3,4,8))
        ##binario
        plot(SDMpred > thre, main=gsub('_', ' ', sp_i), cex=1.3, legend=FALSE)
        plot(SAborders, add=T)
        grid()
        legend('topright',legend=c('non-habitat','habitat'), pch='', cex=1.3, fill=c('lightgrey','darkgreen'), bg='white')
        ##continuo
        plot(SDMpred,main=gsub('_', ' ', sp_i), cex=1.3)
        plot(SAborders, add=T)
        grid()
        dev.off()

        ##salvando um gridfile
        writeRaster(SDMpred, paste(sp_i,'Suitability.asc',sep=''))

        
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
        suitabMapPath = paste(projectFolder,'/SDM outputs/resultados SDM com humanos/', sp_i, '/', sp_i,'Suitability.asc',sep='')
        
        
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

}




##TERCEIRA PARTE: comparando SDMs com e sem interacao com humanos##




##tabela de resultados
tabOutputsSdmH = read.csv('/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/SDM outputs/resultados SDM com humanos/tabOutputs.csv', header=TRUE)
tabOutputsSdm = read.csv('/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/SDM outputs/resultados SDM sem humanos/tabOutputs.csv', header=TRUE)

##thresholds
threshSdm = tabOutputsSdm[tabOutputsSdm$sps == sp_i, 'threshold']
threshSdmH = tabOutputsSdmH[tabOutputsSdmH$sps == sp_i, 'threshold']


##Criando objeto com a lista dos nomes das especies
occ.sps <- list.files(paste(spOccFolder,'/sps_occ_Lucas',sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))

##criando uma tabela vazia para salvador alguns dados
tabModelComparation = data.frame()


##loop para SDM com cada uma das especies##


for (sp_i in splist){
    ##for (sp_i in splist[1:3]){

    ##diretorio base de trabalho
    setwd(paste(projectFolder,'/SDM outputs',sep=''))

    ##mapa de projecao
    sdmPred = raster(paste(projectFolder,'/SDM outputs/resultados SDM sem humanos/', sp_i, '/', sp_i,'Suitability.asc',sep=''))
    sdmHPred = raster(paste(projectFolder,'/SDM outputs/resultados SDM com humanos/', sp_i, '/', sp_i,'Suitability.asc',sep=''))

    ##mapas binarios
    sdmPredBIN = sdmPred > threshSdm
    sdmHPredBIN = sdmHPred > threshSdmH

    
    ##plot da sobreposicao dos mapas binarios

    jpeg( ...)
    plot(sdmPredBIN + 2*sdmHPredBIN, col=c('white','blue','red','purple'), legend=FALSE)
    plot(SAborders, add=TRUE)
    grid()
    legend('topright', legend=c(expression(SDM[clim]), expression(SDM[human]), 'Overlapping'), pch=c(21,21,21), pt.bg=c('blue','red','purple'), bty='n', cex=1)


    

    ##plot da diferenca entre mapas de suitability continuos
    
    sdmDiff = sdmHPred-sdmPred
    
    breaks <- seq(-0.1, 0.1, by=0.01)
    cols <- colorRampPalette(c("blue", "white", "red"))(length(breaks)-1)

    plot(sdmDiff, col=cols)
    plot(SAborders, add=TRUE)
    grid()


    ##comparacoes a partir dos mapas binarios
    MapBinJoin = sdmPredBIN + 2*sdmHPredBIN

    MapBinJoinTab = as.data.frame(freq(MapBinJoin))

    cells1 = MapBinJoinTab[MapBinJoinTab$value == 1 & !is.na(MapBinJoinTab$value), 'count'] ##SDMclim
    cells2 = MapBinJoinTab[MapBinJoinTab$value == 2 & !is.na(MapBinJoinTab$value), 'count'] ##SDMhuman
    cells3 = MapBinJoinTab[MapBinJoinTab$value == 3 & !is.na(MapBinJoinTab$value), 'count'] ##sobreposicao

    sdmClim = cells1/(cells1+cells2+cells3)
    sdmHuman = cells2/(cells1+cells2+cells3)
    sdmOverlap = cells3/(cells1+cells2+cells3)

    ##comparacoes a partir da diferencia dos mapas de suitability

    rm(list=c(cells1,cells2,cells3))
    sdmDiff freq(sdmDiff<0)
    
    cells1 = as.numeric( freq(sdmDiff > 0)[2,2] )  ##area onde SDMhuman > SDMclim
    cells2 = as.numeric( freq(sdmDiff < 0)[2,2] ) ##area onde SDMhuman < SDMclim

    sdmDiffHuman = cells1/(cells1 + cells2)
    sdmDiffClim = cells2/(cells1 + cells2)


    ##colocando na tabela


    tabModelComparation = rbind(tabModelComparation,
                                data.frame(sps = sp_i,
                                           binary.proportion.SDMclim = sdmClim,
                                           binary.proportion.SDMhuman = sdmHuman,
                                           binary.proportion.Overlap = sdmOverlap,
                                           proportion.SDMhuman.greater = sdmDiffHuman,
                                           proportion.SDMclim.greater = sdmDiffClim))

    write.csv(paste(projectFolder,'/SDM outputs/tabModelComparation.csv',sep''), row.names=FALSE)
    
}



################################################################################
##################### FIM DA PARTE NOVA ########################################
################################################################################





###PRIMEIRA PARTE: planilha de presencas, backgrownd e variaveis ambientais###

##definindo as pastas de trabalho
## envVarFolder = "J:/Lucas/Modelagem barbeiros/Variaveis Climaticas"
## spOccFolder = "J:/Lucas/Modelagem barbeiros/Ocorrencias"
## projectFolder = "J:/Lucas/Modelagem barbeiros/"
#anderson
envVarFolder = "/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Variaveis Climaticas"
spOccFolder = "/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/"
projectFolder = "/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/resultados nicho climatico + impacto humano/"

##abrindo as variaveis climaticas
##abrindo shape da America do Sul
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp")

##abrindo e cortando camadas de variaveis ambientais para o presente
filesRaw <- stack(list.files(path=envVarFolder, pattern='.asc', full.names=TRUE)) #modificar a extensao .bil de acordo com os arquivos
#files = mask(filesRaw,AmSulShape) #cortando para Am. do Sul

##abrindo e cortando camadas de variaveis ambientais para projecao
filesProjectionOtimistaRaw <- stack(list.files(path=paste(envVarFolder,"/futuro/cenario_otimista",sep=''), pattern='.asc', full.names=TRUE))
filesProjectionPessimistaRaw <- stack(list.files(path=paste(envVarFolder,"/futuro/cenario_pessimista",sep=''), pattern='.asc', full.names=TRUE)) 
##filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul

##remove highly correlated variables Bio1,Bio3,Bio9,Bio13,Bio14
#Lucas
## files.crop.sub = filesRaw[[c('bio7','bio15','bio11','bio16','bio17')]] #choose selected layers
## files.crop.sub.projection.otimista = filesProjectionOtimistaRaw[[c('bio7','bio15','bio11','bio16','bio17')]] #choose selected layers
## files.crop.sub.projection.pessimista = filesProjectionPessimistaRaw[[c('bio7','bio15','bio11','bio16','bio17')]] #choose selected layers
#Anderson
files.crop.sub = filesRaw[[c('bio7','bio15','bio11','bio16','bio17')]] #choose selected layers
files.crop.sub.projection.otimista = filesProjectionOtimistaRaw[[c('bio7','bio15','bio11','bio16','bio17')]] #choose selected layers
files.crop.sub.projection.pessimista = filesProjectionPessimistaRaw[[c('bio7','bio15','bio11','bio16','bio17')]] #choose selected layers

##indice de imapcto humano
hii = raster(x='/home/anderson/PosDoc/dados_ambientais/hii-s-america-geo-grid/res2-5/hii-2-5.asc')
extent(hii) = extent(files.crop.sub)

##definindo os objetos para as variaveis preditoras (COM IMPACTO HUMANO)
predictors <- stack(files.crop.sub,hii)
predictorsProjectionOtimista = stack(files.crop.sub.projection.otimista,hii)
predictorsProjectionPessimista = stack(files.crop.sub.projection.pessimista,hii)

##Criando objeto com a lista de especies
occ.sps <- list.files(paste(spOccFolder,sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))
splist

##criando uma tabela vazia para salvador alguns dados
tabRes = data.frame(c(sp=character(),auc=numeric(),tss=numeric(),threshold=numeric()))


###SEGUNDA PARTE: rodando SDMs para as especies (e fazendo projecoes)###

for (i in 1:length(splist)){

    especie = splist[i] #selecting the species
    sp.file <- read.csv(paste(spOccFolder,"/",especie,".csv",sep=""),header=TRUE) ### read sp occurrence
    sp.occ <- sp.file[,2:3] ## select long-lat
  
    ##CRIANDO E RODANDO O MODELO - Atencao: precisa ter uma pasta para cada sp
    MX <- maxent(x=predictors,p=sp.occ,path=paste(projectFolder,'maxent/',especie,sep=""), 
                 args=c(
                     'responsecurves=true',
                     'jackknife=true',
                     'randomseed=true',
                     'randomtestpoints=25',
                     'betamultiplier=1',
                     'replicates=3',
                     'replicatetype=Subsample',
                     'writebackgroundpredictions=true',
                     'linear=true',
                     'quadratic=true',
                     'product=false',
                     'threshold=false',
                     'hinge=false',
                     'maximumiterations=1000',
                     'convergencethreshold=1.0E-5',
                     'threads=2'))
  
    ##avaliacao do modeo
    mres = read.csv(paste(projectFolder,'maxent/',especie,"/maxentResults.csv",sep=''),header=TRUE)
    threshold = mres$X10.percentile.training.presence.logistic.threshold[nrow(mres)] #threshold
    aucModel = mres$Test.AUC[nrow(mres)]
    tss = TSSmaxent(paste(projectFolder,'maxent/',especie,'/',sep=''))$TSS #Atencao: o script do TSS deve estar na 
  
    ##gravando uma tabela
    tabRes = rbind(tabRes,data.frame(sp=especie,auc=aucModel,tss=tss,threshold=threshold))
    
    ##fazendo as projecoes
    projecaoSuitability = predict(MX,predictors) #projecao para o presente
    projecaoFuturoOtimista = predict(MX,predictorsProjectionOtimista) #projecao para o futuro - cenario otimista
    projecaoFuturoPessimista = predict(MX,predictorsProjectionPessimista) #projecao para o futuro - cenario pessimista
    
    ##gravando arquivo raster dos mapas de projecao gerados pelo modelo
    ##mapas do presente
    writeRaster(mean(projecaoSuitability),filename=paste(projectFolder,'Projecoes/',especie,".asc", sep=""),overwrite=TRUE)
    writeRaster(mean(projecaoSuitability)>threshold,filename=paste(projectFolder,'Projecoes/',especie,"BIN.asc", sep=""),overwrite=TRUE)
    ##mapas do futuro otimista
    writeRaster(mean(projecaoFuturoOtimista),filename=paste(projectFolder,'Projecoes/',especie,"Otimista.asc", sep=""),overwrite=TRUE)
    writeRaster(mean(projecaoFuturoOtimista)>threshold,filename=paste(projectFolder,'Projecoes/',especie,"OtimistaBIN.asc", sep=""),overwrite=TRUE)
    ##mapas do futuro pessimista
    writeRaster(mean(projecaoFuturoPessimista),filename=paste(projectFolder,'Projecoes/',especie,"Pessimista.asc", sep=""),overwrite=TRUE)
    writeRaster(mean(projecaoFuturoPessimista)>threshold,filename=paste(projectFolder,'Projecoes/',especie,"PessimistaBIN.asc", sep=""),overwrite=TRUE)
    
}

##gravando a tabela de resultados completa
write.csv(tabRes,file=paste(projectFolder,'tabela_resultados.csv',sep=''))


###TERCEIRA PARTE: gerando mapas de sobreposicao (i.e. mapas de riqueza de especies) - COM impacto humano###

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


###QUARTA PARTE: gerando MAPAS DE RISCO - COM impacto humano###

##indices para RISCO DE INFECCAO por especie de vetor
##link: http://portalarquivos.saude.gov.br/images/pdf/2015/agosto/03/2014-020..pdf
tabBarb = read.csv("/home/anderson/Documentos/Projetos/macroecologia_de_chagas/Taxa de infeccao natural vetores 2007-2011.csv",header=TRUE)
#infecBarb = sort(tabBarb[,3],decreasing=TRUE)
#infecIndOrdered = c(infecBarb[6],infecBarb[3]) #,infecBarb[5],infecBarb[9],infecBarb[10],infecBarb[11],infecBarb[13],infecBarb[1])
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
listaNomes = names(camadasFuturoPessimista)
listaNomes = gsub(pattern='PessimistaBIN',replacement='',x=listaNomes)
infecIndOrdered = tabBarb$taxaInfeccaonatural[match(listaNomes,tabBarb$sp)]/100 #taxa de infeccao natural na ordem dos rasters (de 0 a 1)
#
mapaRiscoFuturoPessimista = sum(camadasFuturoPessimista*infecIndOrdered)
mapaRiscoFuturoPessimista = mapaRiscoFuturoPessimista/length(infecIndOrdered)
plot(mapaRiscoFuturoPessimista)
writeRaster(x=mapaRiscoFuturoPessimista,filename=paste(projectFolder,'Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''),overwrite=TRUE)


###QUINTA PARTE: estatisticas sumarias a partir dos mapas - COM impacto humano###


##definindo a area do Brasil na America do Sul, para os mapas
areaBR = extent(-80.00635,-31.71555,-37.8679,8.156474) #extent do Brasil
AmSulBR = crop(AmSulShape, extent(areaBR))  #america do sul recortada para o BR

##presente 

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

##futuro otimista 

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

##tamanho da area do quartil superior (para riqueza e risco), para comparar os cenarios
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
rm(tab)
tab = data.frame(scenario=c('pres','fut_otim','fut_pess'),
                 quantile75riq = c(quantile(mapaRiquezaPresenteBR, 0.75,na.rm=TRUE),quantile(mapaRiquezaFuturoOtimistaBR, 0.75,na.rm=TRUE),quantile(mapaRiquezaFuturoPessimistaBR, 0.75,na.rm=TRUE)),
                 quantile75risc = c(quantile(mapaRiscoPresenteBR, 0.75,na.rm=TRUE),quantile(mapaRiscoFuturoOtimistaBR, 0.75,na.rm=TRUE), quantile(mapaRiscoFuturoPessimistaBR, 0.75,na.rm=TRUE)),
                 percCellRiq = c(percCelRiqPres,percCelRiqOtim,percCelRiqPess),
                 percCellRisc = c(percCelRiscPres,percCelRiscOtim,percCelRiscPess),
                 corRiqRisc = c(corPres[1,2],corOtim[1,2],corPess[1,2])
)

write.csv(tab,paste(projectFolder,'statsResHII.csv'),row.names = FALSE)


###SEXTA PARTE: figuras dos mapas - COM impacto humano###

##abrindo
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
jpeg(filename=paste(projectFolder,'mapasHIIbr.jpeg',sep=''),width=1750,height=850)
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
##legenda
plot(mapaRiscoFuturoPessimistaBR,legend.only=TRUE,legend.width=3,axis.args=list(cex.axis=2),legend.args=list(text='Risk of infected vector',font=2,side=4,line=6.5,cex=2.2,cex.axis=0.2))
dev.off()


###SETIMA PARTE: comparando areas de ALTA RIQUEZA entre as projecoes com e sem impacto humano###


##pasta base COMUM ao dois tipos de modelo
projectFolder = '/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/'

##presente

##abrindo
mapaRiquezaPresente = raster(paste(projectFolder,'resultados nicho climatico/Mapas de riqueza/mapaRiquezaPresente.asc',sep=''))
mapaRiquezaPresenteHii = raster(paste(projectFolder,'resultados nicho climatico + impacto humano/Mapas de riqueza/mapaRiquezaPresente.asc',sep=''))

##cortando para BR
mapaRiquezaPresenteBR = mask(mapaRiquezaPresente, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiquezaPresenteHiiBR = mask(mapaRiquezaPresenteHii, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))

##areas de alta riqueza
hightRiqPres = mapaRiquezaPresenteBR > quantile(mapaRiquezaPresenteBR, 0.75, na.rm=TRUE) #raster quartil superior riqueza
hightRiqPresHii = mapaRiquezaPresenteHiiBR > quantile(mapaRiquezaPresenteHiiBR, 0.75, na.rm=TRUE) #raster quartil superior riqueza

##futuro otimista

##abrindo
mapaRiquezaFuturoOtimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de riqueza/mapaRiquezaFuturoOtimista.asc',sep=''))
mapaRiquezaFuturoOtimistaHii = raster(paste(projectFolder,'resultados nicho climatico + impacto humano/Mapas de riqueza/mapaRiquezaFuturoOtimista.asc',sep=''))

##cortanto para o BR
mapaRiquezaFuturoOtimistaBR = mask(mapaRiquezaFuturoOtimista, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiquezaFuturoOtimistaHiiBR = mask(mapaRiquezaFuturoOtimistaHii, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))

##areas de alta riqueza 
hightRiqOtim = mapaRiquezaFuturoOtimistaBR > quantile(mapaRiquezaFuturoOtimistaBR, 0.75,na.rm=TRUE) #raster quartil superior riqueza
hightRiqOtimHii = mapaRiquezaFuturoOtimistaHiiBR > quantile(mapaRiquezaFuturoOtimistaHiiBR, 0.75,na.rm=TRUE) #raster quartil superior riqueza

##futuro pessimista

##abrindo
mapaRiquezaFuturoPessimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''))
mapaRiquezaFuturoPessimistaHii = raster(paste(projectFolder,'resultados nicho climatico + impacto humano/Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''))

##cortando para o Brasil
mapaRiquezaFuturoPessimistaBR = mask(mapaRiquezaFuturoPessimista, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiquezaFuturoPessimistaHiiBR = mask(mapaRiquezaFuturoPessimistaHii, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))

##areas de alta riqueza
hightRiqPess = mapaRiquezaFuturoPessimistaBR > quantile(mapaRiquezaFuturoPessimistaBR, 0.75,na.rm=TRUE) #raster quartil superior riqueza
hightRiqPessHii = mapaRiquezaFuturoPessimistaHiiBR > quantile(mapaRiquezaFuturoPessimistaHiiBR, 0.75,na.rm=TRUE) #raster quartil superior riqueza

##unificando para contabilizar areas
hightRiqPresUni = hightRiqPres*1 + hightRiqPresHii*2
hightRiqOtimUni = hightRiqOtim*1 + hightRiqOtimHii*2
hightRiqPessUni = hightRiqPess*1 + hightRiqPessHii*2

## ##correlacao areas alta riqueza modelos com e sem HII
## ##presente
## testRiqPres <- getValues(stack(hightRiqPres,hightRiqPresHii))
## corRiqPres <- as.data.frame(cor(testRiqPres, use="complete.obs",method='spearman'))
## ##futuro otimista
## testRiqOtim <- getValues(stack(hightRiqOtim,hightRiqOtimHii))
## corRiqOtim <- as.data.frame(cor(testRiqOtim, use="complete.obs",method='spearman'))
## ##futuro pessimista
## testRiqPess <- getValues(stack(hightRiqPess,hightRiqPessHii))
## corRiqPess <- as.data.frame(cor(testRiqPess, use="complete.obs",method='spearman'))

##correlacao mapa riqueza modelos com e sem HII## OBS: ACHEI MELHOR USAR ESSE NO MANUSCRITO
##presente
testRiqPres <- getValues(stack(mapaRiquezaPresenteBR, mapaRiquezaPresenteHiiBR))
corRiqPres <- as.data.frame(cor(testRiqPres, use="complete.obs",method='pearson'))
##futuro otimista
testRiqOtim <- getValues(stack(mapaRiquezaFuturoOtimistaBR, mapaRiquezaFuturoOtimistaHiiBR))
corRiqOtim <- as.data.frame(cor(testRiqOtim, use="complete.obs",method='pearson'))
##futuro pessimista
testRiqPess <- getValues(stack(mapaRiquezaFuturoPessimistaBR, mapaRiquezaFuturoPessimistaHiiBR))
corRiqPess <- as.data.frame(cor(testRiqPess, use="complete.obs",method='pearson'))

rm(tabSobreposicaoRiqueza)

tabSobreposicaoRiqueza = data.frame(
    cenario = c('presente','futuroOtimista','futuroPessimista'),
    total_celulas_nichoClim = c(freq(hightRiqPres,value=1),freq(hightRiqOtim,value=1),freq(hightRiqPess,value=1)),
    percentual_do_Brasil= c(
        (freq(hightRiqPres,value=1)/(ncell(hightRiqPres)-freq(hightRiqPres,value=NA)))*100,
        (freq(hightRiqOtim,value=1)/(ncell(hightRiqOtim)-freq(hightRiqOtim,value=NA)))*100,
        (freq(hightRiqPess,value=1)/(ncell(hightRiqPess)-freq(hightRiqPess,value=NA)))*100 ),
    total_celulas_HII= c(freq(hightRiqPresHii,value=1),freq(hightRiqOtimHii,value=1),freq(hightRiqPessHii,value=1)),
    percentual_do_Brasil= c(
        (freq(hightRiqPresHii,value=1)/(ncell(hightRiqPresHii)-freq(hightRiqPresHii,value=NA)))*100,
        (freq(hightRiqOtimHii,value=1)/(ncell(hightRiqOtimHii)-freq(hightRiqOtimHii,value=NA)))*100,
        (freq(hightRiqPessHii,value=1)/(ncell(hightRiqPessHii)-freq(hightRiqPessHii,value=NA)))*100 ),
    total_cels_sobrepostas = c(freq(hightRiqPresUni,value=3),freq(hightRiqOtimUni,value=3),freq(hightRiqPessUni,value=3)),
    percentual_para_nichoClim=
        c(( freq(hightRiqPresUni,value=3)/freq(hightRiqPres,value=1) ) * 100,
        ( freq(hightRiqOtimUni,value=3)/freq(hightRiqOtim,value=1) ) * 100,
        ( freq(hightRiqPessUni,value=3)/freq(hightRiqPess,value=1) ) * 100),
    percentual_para_HII=
        c(( freq(hightRiqPresUni,value=3)/freq(hightRiqPresHii,value=1) ) * 100,
        ( freq(hightRiqOtimUni,value=3)/freq(hightRiqOtimHii,value=1) ) * 100,
        ( freq(hightRiqPessUni,value=3)/freq(hightRiqPessHii,value=1) ) * 100),
    percentual_geral=
        c(freq(hightRiqPresUni,value=3)/(freq(hightRiqPresUni,value=1)+freq(hightRiqPresUni,value=2)+freq(hightRiqPresUni,value=3)) * 100,
        freq(hightRiqOtimUni,value=3)/(freq(hightRiqOtimUni,value=1)+freq(hightRiqOtimUni,value=2)+freq(hightRiqOtimUni,value=3)) * 100,
        freq(hightRiqPessUni,value=3)/(freq(hightRiqPessUni,value=1)+freq(hightRiqPessUni,value=2)+freq(hightRiqPessUni,value=3)) * 100),
    pearson = c(corRiqPres[2,1],corRiqOtim[2,1],corRiqPess[2,1])
)

write.csv(tabSobreposicaoRiqueza,paste(projectFolder,'tabSobreposicaoRiquezaBR.csv',sep=''),row.names=FALSE)

##mapas de sobreposicao

##definindo a area do Brasil na America do Sul, para os mapas

areaBR = extent(-80.00635,-31.71555,-37.8679,8.156474) #extent do Brasil
AmSulBR = crop(AmSulShape, extent(areaBR))  #america do sul recortada para o BR

##cortando para o BR
hightRiqPresBR = mask(hightRiqPres, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
hightRiqPresHiiBR = mask(hightRiqPresHii, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
hightRiqOtimBR = mask(hightRiqOtim, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
hightRiqOtimHiiBR = mask(hightRiqOtimHii, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
hightRiqPessBR = mask(hightRiqPess, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
hightRiqPessHiiBR = mask(hightRiqPessHii, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
                
##mapas
jpeg(filename=paste(projectFolder,'mapasSobreposicaoRiquezaBR.jpeg'),width=1800,height=600)
par(mfrow=c(1,3),mar=c(5,5,5,5))
plot(crop(hightRiqPres*1 + hightRiqPresHii*2, areaBR),col=c('white','blue','red','purple'),main='Current climate',legend=FALSE,cex.axis=2,cex.main=4) + plot(AmSulShape[AmSulShape$CNTRY_NAME!='Brazil',],col='lightgray',add=TRUE) + plot(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',],add=TRUE) + box() + grid()
plot(crop(hightRiqOtim*1 + hightRiqOtimHii*2, areaBR),col=c('white','blue','red','purple'),main='2070 optimistic',legend=FALSE,cex.axis=2,cex.main=4) + plot(AmSulShape[AmSulShape$CNTRY_NAME!='Brazil',],col='lightgray',add=TRUE) + plot(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',],add=TRUE) + box() + grid()
plot(crop(hightRiqPess*1 + hightRiqPessHii*2, areaBR),col=c('white','blue','red','purple'),main='2070 pessimistic',legend=FALSE,cex.axis=2,cex.main=4)+ plot(AmSulShape[AmSulShape$CNTRY_NAME!='Brazil',],col='lightgray',add=TRUE) + plot(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',],add=TRUE) + box() + grid()
legend('bottomright',legend=c(expression(SDM[clim]), expression(SDM[hii]),'Overlapping'),pch=c(21,21,21),pt.bg=c('blue','red','purple'),bty='n',cex=3)
dev.off()


###OITAVA PARTE: comparando areas de ALTO RISCO entre as projecoes com e sem impacto humano###

##pasta base COMUM ao dois tipos de modelo
projectFolder = '/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/'

##presente

##abrindo
mapaRiscoPresente = raster(paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoPresente.asc',sep=''))
mapaRiscoPresenteHii = raster(paste(projectFolder,'resultados nicho climatico + impacto humano/Mapas de risco/mapaRiscoPresente.asc',sep=''))

##cortando para BR
mapaRiscoPresenteBR = mask(mapaRiscoPresente, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiscoPresenteHiiBR = mask(mapaRiscoPresenteHii, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))

##areas de alto risco
hightRiscPres = mapaRiscoPresenteBR > quantile(mapaRiscoPresenteBR, 0.75,na.rm=TRUE) #raster quartil superior risco
hightRiscPresHii = mapaRiscoPresenteHiiBR > quantile(mapaRiscoPresenteHiiBR, 0.75,na.rm=TRUE) #raster quartil superior risco

##futuro otimista

##abrindo
mapaRiscoFuturoOtimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoFuturoOtimista.asc',sep=''))
mapaRiscoFuturoOtimistaHii = raster(paste(projectFolder,'resultados nicho climatico + impacto humano/Mapas de risco/mapaRiscoFuturoOtimista.asc',sep=''))

##cortando para BR
mapaRiscoFuturoOtimistaBR = mask(mapaRiscoFuturoOtimista, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiscoFuturoOtimistaHiiBR = mask(mapaRiscoFuturoOtimistaHii, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))

##areas de alto risco
hightRiscOtim = mapaRiscoFuturoOtimistaBR > quantile(mapaRiscoFuturoOtimistaBR, 0.75,na.rm=TRUE) #raster quartil superior risco
hightRiscOtimHii = mapaRiscoFuturoOtimistaHiiBR > quantile(mapaRiscoFuturoOtimistaHiiBR, 0.75,na.rm=TRUE) #raster quartil superior risco

##futuro pessimista

##abrindo
mapaRiscoFuturoPessimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''))
mapaRiscoFuturoPessimistaHii = raster(paste(projectFolder,'resultados nicho climatico + impacto humano/Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''))

##cortando para BR
mapaRiscoFuturoPessimistaBR = mask(mapaRiscoFuturoPessimista, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
mapaRiscoFuturoPessimistaHiiBR = mask(mapaRiscoFuturoPessimistaHii, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))

##areas de alto risco
hightRiscPess = mapaRiscoFuturoPessimistaBR > quantile(mapaRiscoFuturoPessimistaBR, 0.75,na.rm=TRUE) #raster quartil superior risco
hightRiscPessHii = mapaRiscoFuturoPessimistaHiiBR > quantile(mapaRiscoFuturoPessimistaHiiBR, 0.75,na.rm=TRUE) #raster quartil superior risco

##unificando para contabilizar areas
hightRiscPresUni = hightRiscPres*1 + hightRiscPresHii*2
hightRiscOtimUni = hightRiscOtim*1 + hightRiscOtimHii*2
hightRiscPessUni = hightRiscPess*1 + hightRiscPessHii*2

## ##correlacao areas alta riqueza modelos com e sem HII
## ##presente
## testRiscPres <- getValues(stack(hightRiscPres,hightRiscPresHii))
## corRiscPres <- as.data.frame(cor(testRiscPres, use="complete.obs",method='spearman'))
## ##futuro otimista
## testRiscOtim <- getValues(stack(hightRiscOtim,hightRiscOtimHii))
## corRiscOtim <- as.data.frame(cor(testRiscOtim, use="complete.obs",method='spearman'))
## ##futuro pessimista
## testRiscPess <- getValues(stack(hightRiscPess,hightRiscPessHii))
## corRiscPess <- as.data.frame(cor(testRiscPess, use="complete.obs",method='spearman'))

##correlacao mapa risco modelos com e sem HII## OBS: ACHEI MELHOR USAR ESSE NO MANUSCRITO
##presente
testRiscPres <- getValues(stack(mapaRiscoPresenteBR,mapaRiscoPresenteHiiBR))
corRiscPres <- as.data.frame(cor(testRiscPres, use="complete.obs",method='pearson'))
##futuro otimista
testRiscOtim <- getValues(stack(mapaRiscoFuturoOtimistaBR,mapaRiscoFuturoOtimistaHiiBR))
corRiscOtim <- as.data.frame(cor(testRiscOtim, use="complete.obs",method='pearson'))
##futuro pessimista
testRiscPess <- getValues(stack(mapaRiscoFuturoPessimistaBR,mapaRiscoFuturoPessimistaHiiBR))
corRiscPess <- as.data.frame(cor(testRiscPess, use="complete.obs",method='pearson'))

rm(tabSobreposicaoRisco)

tabSobreposicaoRisco = data.frame(
    cenario = c('Current','2070 optimistic','2070 pessimistic'),
    total_celulas_nichoClim = c(freq(hightRiscPres,value=1),freq(hightRiscOtim,value=1),freq(hightRiscPess,value=1)),
    percentual_do_Brasil= c(
        (freq(hightRiscPres,value=1)/(ncell(hightRiscPres)-freq(hightRiscPres,value=NA)))*100,
        (freq(hightRiscOtim,value=1)/(ncell(hightRiscOtim)-freq(hightRiscOtim,value=NA)))*100,
        (freq(hightRiscPess,value=1)/(ncell(hightRiscPess)-freq(hightRiscPess,value=NA)))*100 ),
    total_celulas_HII= c(freq(hightRiscPresHii,value=1),freq(hightRiscOtimHii,value=1),freq(hightRiscPessHii,value=1)),
    percentual_do_Brasil= c(
        (freq(hightRiscPresHii,value=1)/(ncell(hightRiscPresHii)-freq(hightRiscPresHii,value=NA)))*100,
        (freq(hightRiscOtimHii,value=1)/(ncell(hightRiscOtimHii)-freq(hightRiscOtimHii,value=NA)))*100,
        (freq(hightRiscPessHii,value=1)/(ncell(hightRiscPessHii)-freq(hightRiscPessHii,value=NA)))*100 ),
    total_cels_sobrepostas = c(freq(hightRiscPresUni,value=3),freq(hightRiscOtimUni,value=3),freq(hightRiscPessUni,value=3)),
    percentual_para_nichoClim=
        c(( freq(hightRiscPresUni,value=3)/freq(hightRiscPres,value=1) )*100,
        ( freq(hightRiscOtimUni,value=3)/freq(hightRiscOtim,value=1) )*100,
        ( freq(hightRiscPessUni,value=3)/freq(hightRiscPess,value=1) )*100),
    percentual_para_HII=
        c(( freq(hightRiscPresUni,value=3)/freq(hightRiscPresHii,value=1) )*100,
        ( freq(hightRiscOtimUni,value=3)/freq(hightRiscOtimHii,value=1) )*100,
        ( freq(hightRiscPessUni,value=3)/freq(hightRiscPessHii,value=1) )*100),
    percentual_geral=
        c((freq(hightRiscPresUni,value=3)/(ncell(hightRiscPresUni)-freq(hightRiscPresUni,value=NA))*100),
        (freq(hightRiscOtimUni,value=3)/(ncell(hightRiscOtimUni)-freq(hightRiscOtimUni,value=NA))*100),
        ( freq(hightRiscPessUni,value=3)/(ncell(hightRiscPessUni)-freq(hightRiscPessUni,value=NA))*100)),
    pearson = c(corRiscPres[2,1],corRiscOtim[2,1],corRiscPess[2,1])
)

write.csv(tabSobreposicaoRisco,paste(projectFolder,'tabSobreposicaoRiscoBR.csv',sep=''),row.names=FALSE)

##mapas de sobreposicao

##definindo a area do Brasil na America do Sul, para os mapas
areaBR = extent(-80.00635,-31.71555,-37.8679,8.156474) #extent do Brasil
AmSulBR = crop(AmSulShape, extent(areaBR))  #america do sul recortada para o BR

##cortando para o BR
hightRiscPresBR = mask(hightRiscPres, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
hightRiscPresHiiBR = mask(hightRiscPresHii, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
hightRiscOtimBR = mask(hightRiscOtim, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
hightRiscOtimHiiBR = mask(hightRiscOtimHii, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
hightRiscPessBR = mask(hightRiscPess, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
hightRiscPessHiiBR = mask(hightRiscPessHii, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))
                
##mapas
jpeg(filename=paste(projectFolder, 'mapasSobreposicaoRiscoBR.jpeg'), width=1800, height=600)
par(mfrow=c(1,3),mar=c(5,5,5,5))
plot(crop(hightRiscPres*1 + hightRiscPresHii*2, areaBR), col=c('white','blue','red','purple'), main='Current climate', legend=FALSE, cex.axis=2, cex.main=4) + plot(AmSulShape[AmSulShape$CNTRY_NAME!='Brazil',], col='lightgray', add=TRUE) + plot(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',], add=TRUE) + box() + grid()
plot(crop(hightRiscOtim*1 + hightRiscOtimHii*2, areaBR), col=c('white','blue','red','purple'), main='2070 optimistic', legend=FALSE, cex.axis=2, cex.main=4) + plot(AmSulShape[AmSulShape$CNTRY_NAME!='Brazil',], col='lightgray',add=TRUE) + plot(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',], add=TRUE) + box() + grid()
plot(crop(hightRiscPess*1 + hightRiscPessHii*2, areaBR), col=c('white','blue','red','purple'), main='2070 pessimistic', legend=FALSE, cex.axis=2, cex.main=4) + plot(AmSulShape[AmSulShape$CNTRY_NAME!='Brazil',], col='lightgray', add=TRUE) + plot(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',], add=TRUE) + box() + grid()
legend('bottomright', legend=c(expression(SDM[clim]), expression(SDM[HII]), 'Overlapping'), pch=c(21,21,21), pt.bg=c('white','blue','red','purple'), bty='n', cex=3)
dev.off()

##graficos das tendencias: percentual de cobertura da america do sul

##abrindo dados
tabSobreposicaoRisco = read.csv(paste(projectFolder,'tabSobreposicaoRiscoBR.csv',sep=''), header=TRUE)
tabSobreposicaoRiqueza = read.csv(paste(projectFolder,'tabSobreposicaoRiquezaBR.csv',sep=''), header=TRUE)

##grafico
jpeg(filename=paste(projectFolder,'tendenciasCoberturaBR.jpeg'),width=500,height=500)
##riqueza SDMclim
plot.default(c(1,2,2),tabSobreposicaoRiqueza$percentual_do_Brasil,axe=FALSE,xlim=c(0.75,2.25), ylim=c(0,40), pch=c(15,16,17), cex=3, lwd=2, col=rgb(0,1,0,0.5), xlab='Scenario', ylab='Extent of Brazil (in %)') + axis(side=1, at=c(1,2,2), labels=c('Current climate','Future climate (2070)','Future climate (2070)')) + axis(side=2, at=c(c(0:5)*10), labels=c(c(0:5)*10)) + box() + grid()
abline(lm(tabSobreposicaoRiqueza$percentual_do_Brasil~c(1,2,2)), col='green',lty=1, lwd=2)
##riqueza SDMhii
points(c(1,2,2),tabSobreposicaoRiqueza$percentual_do_Brasil.1, xlabel='', ylabel='', pch=c(15,16,17), ylim=c(0,50), xlim=c(0.5,2.5), cex=3, col=rgb(0,1,0,0.5)); grid()
abline(lm(tabSobreposicaoRiqueza$percentual_do_Brasil.1~c(1,2,2)), col='green', lty=2, lwd=2)
##risco SDMclim
points(tabSobreposicaoRisco$percentual_do_Brasil~c(1,2,2), xlab='', ylab='', pch=c(15,16,17), ylim=c(0,50), xlim=c(0.5,2.5), cex=3, lwd=2, col=rgb(1,0,0,0.5)); grid()
abline(lm(tabSobreposicaoRisco$percentual_do_Brasil~c(1,2,2)), col='red', lty=1, lwd=2)
##risco SDMhii
points(tabSobreposicaoRisco$percentual_do_Brasil.1~c(1,2,2), xlab='', ylab='', pch=c(15,16,17), ylim=c(0,50), xlim=c(0.5,2.5), cex=3, lwd=2, col=rgb(1,0,0,0.5)); grid()
abline(lm(tabSobreposicaoRisco$percentual_do_Brasil.1~c(1,2,2)), col='red', lty=2, lwd=2)
##legenda
legend('topleft', legend=c('Presente','2070 optimistic','2070 pessimistic',expression(SDM[clim]),expression(SDM[HII]),'Richness','Risk'),pch=c(0,1,2,NA,NA,21,21), lty=c(NA,NA,NA,1,2,NA,NA), col=c(1,1,1,1,1,0,0), pt.bg=c(0,0,0,0,0,'green','red') ,text.col=c('black','black','black','black','black','green','red'))
dev.off()

##graficos das tendencias: correlacao
jpeg(filename=paste(projectFolder,'tendenciasCorrelacaoBR.jpeg'),width=500,height=500)
plot.default(c(1,2,2),tabSobreposicaoRisco$pearson,axe=FALSE,xlim=c(0.75,2.25), ylim=c(0,1), pch=c(15,16,17), cex=3, lwd=2, col=rgb(1,0,0,0.5), xlab='Scenario', ylab='Pearson correlation') + axis(side=1, at=c(1,2,2), labels=c('Current climate','Future climate (2070)','Future climate (2070)')) + axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1.0), labels=c(0,0.2,0.4,0.6,0.8,1.0)) + box() + grid()
abline(lm(tabSobreposicaoRisco$pearson~c(1,2,2)),col='red',lwd=2)
points(tabSobreposicaoRiqueza$pearson~c(1,2,2),pch=c(15,16,17),ylim=c(0,50),xlim=c(0.5,2.5),cex=3,col=rgb(0,1,0,0.5));grid()
abline(lm(tabSobreposicaoRiqueza$pearson~c(1,2,2)),col='green',lwd=2)
legend('topright',legend=c('Current climate','2070 optimistic','2070 pessimistic','Richness','Risk'),pch=c(0,1,2,21,21), pt.bg=c(NA,NA,NA,'green','red'), text.col=c('black','black','black','green','red'))
dev.off()
