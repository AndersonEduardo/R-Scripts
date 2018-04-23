

library("spocc")
library(raster)
library(usdm)
library(ENMeval)
library(pROC)
library(ecospat)



##PARTE 1: variaveis e parametros globais




##caminhos##

projectFolder = '/home/anderson/Projetos/Darlan - macroevolucao caviomorfa'
occDataFolder = '/home/anderson/Projetos/Darlan - macroevolucao caviomorfa/GBIF_Livro'
envDataFolder = '/home/anderson/Projetos/Darlan - macroevolucao caviomorfa/Variaveis selecionadas'


##nomes das especies##

spList = gsub(pattern='.csv', replacement='', x=list.files(path=occDataFolder, pattern='.csv', full.names=FALSE))
spList = gsub(pattern=' ', replacement='_', x=spList)


##array de combinacoes##

namesCombTab = data.frame()

for (sp_i in spList){
    for (sp_j in spList){
        if (sp_i != sp_j){        
            namesCombTab = rbind(namesCombTab,
                                 data.frame(sp1 = sp_i,
                                            sp2 = sp_j)
                                 )
        }
    }
}


##variaveis ambientais para America do Sul##

SAborders = rgdal::readOGR('/home/anderson/shapefiles/ne_50m_land/ne_50m_land.shp') #bordas de continentes
SOAextent = extent(-81.57551,-34.03384,-57.13385,12.99115)
SAborders = crop(SAborders,SOAextent)
predictorsRaw = getData('worldclim', var='bio', res=2.5) #dowload wordclim 
predictorsRaw = crop(predictorsRaw, SAborders)
save(predictorsRaw, file=paste(envDataFolder,'/fullPredictorsSouthAmerica.R',sep=''))


##selecao de variaveis##

load(paste(envDataFolder,'/fullPredictorsSouthAmerica.R',sep=''))
predictorsForVif = predictorsRaw

vif(predictorsForVif)
predictorsVif1 = vifcor(predictorsForVif, th=0.7)
predictorsVif1

predictorsVif2 <- vifstep(predictorsForVif, th=10) # identify collinear variables that should be excluded
predictorsVif2

##comparando
predictorsVif1@results$Variables
predictorsVif2@results$Variables

##definindo variaveis preditoras a serem usadas nos modelos
predictors = predictorsForVif[[ grep(pattern=paste(predictorsVif2@results$Variables,collapse='|'), x=names(predictors), value=TRUE) ]]

save(predictors, file=paste(envDataFolder,'/uncorrelatedPredictorsSouthAmerica.R',sep=''))


##KDE para vies amostral a partir de ocorrencias no GBIF para a familia rodentia##

##abrindo arquivo salvo
rodentiaDataset = read.csv(file='/home/anderson/Projetos/Darlan - macroevolucao caviomorfa/GBIF ordem rodentia Am Sul/0004307-180412121330197.csv', header=TRUE, sep='\t', dec='.', stringsAsFactors=FALSE, na.strings="")
rodentiaDataset = subset(rodentiaDataset, coordinateuncertaintyinmeters < 5000 | is.na(coordinateuncertaintyinmeters))

rodentiaOcc = rodentiaDataset[,c('decimallongitude','decimallatitude')]
names(rodentiaOcc) = c('lon','lat')
rodentiaOcc = rodentiaOcc[which(complete.cases(rodentiaOcc)),]

##convertendo de data.frame para staialPointsDataFrame
coordinates(rodentiaOcc) <- ~lon+lat
proj4string(rodentiaOcc) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')

## ##abrindo o shapefile que 'cortara' os pontos
## wrld = readOGR("/home/anderson/shapefiles/ne_50m_land/ne_50m_land.shp")

##ajustando o sistema de coordenadas geograficas entre pontos e shapefile
rodentiaOcc <- spTransform(rodentiaOcc, crs(SAborders))
rodentiaOcc <- rodentiaOcc[SAborders, ]

## ##sjuatando para toda a area da america do sul e grid file final
## #SOAextent = extent(-81.57551,-34.03384,-57.13385,12.99115)

## SAbg = predictors[[1]]*0 ##America do Sul como 'pano de fundo'
## crs(SAbg) = crs(raster())

## ##inspecionando pontos
## plot(rodentiaOcc)
## plot(wrld, add=TRUE)

## ##convertendo spatialPoints em data.frame para o KDE2D
## reduviidaeOcc = as.data.frame(reduviidaeOcc)


##kernel density estimation com o pacote MASS
dens = MASS::kde2d(x=rodentiaOcc$lon, y=rodentiaOcc$lat, n=100)
densRas = raster(dens)
values(densRas)[values(densRas)<=0] = 0
densRas = mask(x=densRas, mask=SAborders)
densRas = crop(x=densRas, y=SAborders)

##ajustando projecao e fundindo com pano de fundo
predAreaNatBG = predictors[[1]]
                                        #predAreaNatBG = crop(x=predAreaNatBG, y=natAreaExtent)
predAreaNatBG = predAreaNatBG*0
densRas = projectRaster(densRas, crs=proj4string(predAreaNatBG), res=res(predAreaNatBG), method="bilinear")
densRas = merge(densRas, predAreaNatBG, tolerance=0.3)


densRas = crop(x=densRas, y=predAreaNatBG)
extent(densRas) = extent(predAreaNatBG)

##salvado raster
writeRaster(x=densRas, file=paste(projectFolder,'/biasLayerRodentia.asc', sep=''), overwrite=TRUE)




##PARTE 2: modelagem de distribuicao de especies




##preditoras
load(paste(envDataFolder,'/uncorrelatedPredictorsSouthAmerica.R', sep=''))
predictors = crop(x=predictors, y=SAborders)

##bias layer
biasLayer = raster(paste(projectFolder,'/biasLayerRodentia.asc', sep=''))
biasLayer = crop(x=biasLayer, y=SAborders)


##tabela de resultados
tabRes = data.frame()


##loop para SDM com cada uma das especies##

for (sp_i in unique(namesCombTab$sp1)){ 
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
        occPoints = read.csv(paste(occDataFolder,'/',gsub(pattern='_',replacement=' ',x=sp_i),'.csv',sep=''), header=TRUE, sep=',', dec='.', na.strings='',colClasses=c('character','numeric','numeric')) #abrindo pontos de ocorrencia
        names(occPoints) =  c('sp','lon','lat')
        occPoints = occPoints[,c('lon','lat')]

        ##convertendo de data.frame para staialPointsDataFrame
        coordinates(occPoints) <- ~lon+lat
        proj4string(occPoints) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')

        ## ##abrindo o shapefile que 'cortara' os pontos
        ## SAborders = readOGR("/home/anderson/shapefiles/ne_50m_land/ne_50m_land.shp")

        ##ajustando o sistema de coordenadas geograficas entre pontos e shapefile
        occPoints <- spTransform(occPoints, crs(SAborders)) #espacializando
        occPoints <- occPoints[SAborders, ] #lipando pontos fora do mapa
        occPoints = as.data.frame(occPoints) #pontos limpos de volta para data.frame

        ##pseudo-ausencia com o mesmo vies dos dados de ocorrencia    
        bgPoints = randomPoints(mask=biasLayer, n=10000, p=occPoints, prob=TRUE)
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
                               RMvalues = c(1:5),
                               fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                               clamp = FALSE,
                               parallel = TRUE,
                               numCores = 2)

        
        ##melhor modelo
        bestModel = SDMeval@results[SDMeval@results$AICc==min(SDMeval@results$AICc, na.rm=TRUE),]
        bestModel = bestModel[complete.cases(bestModel),]
        bestModel = bestModel[1,]


        ## RMvalue = as.numeric(unlist(regmatches(bestModel$settings,
        ## gregexpr("[[:digit:]]+\\.*[[:digit:]]*", bestModel$settings))))
        
        ##MAXENT##
        SDMmaxent = maxent(x = dataSet[,grep(pattern='bio',x=names(dataSet),value=TRUE)],
                           p = dataSet[,'occ'],
                           path = paste(projectFolder,'/SDM outputs/',sp_i,sep=''),
                           args = c(unlist(make.args(RMvalues=bestModel$rm, fc=bestModel$features, labels=FALSE)),
                                    'randomseed=TRUE',
                                    'replicates=100',
                                    'randomtestpoints=25',
                                    'replicatetype=Subsample',
                                    'jackknife=TRUE',
                                    'responsecurves=TRUE',
                                    'pictures=TRUE',
                                    'extrapolate=FALSE',
                                    'doclamp=FALSE',
                                    'fadebyclamping=FALSE',
                                    'plots=TRUE',
                                    'threads=2'))

        SDMpred = dismo::predict(SDMmaxent, predictors) #projecao espacial
        crs(SDMpred) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

        SDMpredMean = calc(SDMpred,mean) #mapa da média
        SDMpredSD = calc(SDMpred,sd) #mapa do desvio padrao

        ##salvando um gridfile
        writeRaster(SDMpredMean, paste(sp_i,'SuitabilityMean.asc',sep=''))
        writeRaster(SDMpredSD, paste(sp_i,'SuitabilitySD.asc',sep=''))


        ##calculo do threshold##


        ##valores preditos pelo SDM em cada ponto do dataset
        pred =  extract(x=SDMpredMean,
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

        ##matriz de occ
        occMat = as.matrix(dataSet[which(dataSet$occ==1),c('lon','lat')]) #sps occ matrix
        
        ##espacializacao dos pontos
        coordinates(occPoints) = ~lon+lat
        projection(occPoints) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

        ##centroide da distribuicao espacial
        SpsCentroid = kmeans(occMat, 1)
        SpsCentroid = SpsCentroid$centers

        ##convex hull
        cHull = convHull(occPoints)

        ##distancia entre pontos e centroide
        euc.dist <- function(x1) sqrt(sum((x1 - SpsCentroid) ^ 2))
        meanDistCentroid = mean( apply(occMat, 1, euc.dist) )

        ##buffer - hipotese para area de alcance da especie
        SpsBuffer = raster::buffer(polygons(cHull), width=meanDistCentroid)
        
        ##aplicando bufffer no mapa de suitability
        ## ##suitability continuo
        ## SAbg = predictors[[1]]*0
        ## SDMpredRaw = mask(x=SDMpredMean, mask=SpsBuffer)
        ## SDMpredMean = merge(SDMpredRaw, SAbg)
        ## writeRaster(SDMpredMean, paste(sp_i,'Suitability.asc',sep=''), overwrite=TRUE) 
        ##binario
        SDMpredBIN = SDMpredMean > thre
        crs(SDMpredBIN) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
        SAbg = predictors[[1]]*0
        SDMpredBINRaw = mask(x=SDMpredBIN, mask=SpsBuffer)
        SDMpredBIN = merge(SDMpredBINRaw, SAbg)
        writeRaster(SDMpredBIN, paste(sp_i,'SuitabilityBIN.asc',sep=''), overwrite=TRUE)
        
        
        jpeg(paste(projectFolder,'/SDM outputs/',sp_i,'_mapsOfPoints.jpg',sep=''), width=1000, height=600)
        par(mfrow=c(1,2))
        ##occ e background
        plot(SAbg, main=paste(gsub('_',' ',sp_i)), col='white', cex=1.3, legend=FALSE)
        plot(SAborders, add=TRUE, col='white')
        points(bgPoints, pch=19, cex=0.4, col=rgb(0.4,0.4,0.4,0.3))
        points(occPoints, pch=20, cex=1.1, col=rgb(1,0,0,0.6))
        plot(cHull, add=TRUE)
        plot(SpsBuffer, ad=TRUE, lty=2)
        box()
        grid()
        legend('bottomright', legend=c('Occurrences', 'Background points', 'Minimum convex polygon', 'Buffer'), pch=c(19,21,NA,NA), col=c('red','grey','black','black'), bg='white', lty=c(NA,NA,1,2))
        ##blocks
        plot(SAbg, main=paste(gsub('_',' ',sp_i)), col='white', cex=1.3, legend=FALSE)
        plot(SAborders, add=TRUE, col='white')
        points(occPoints, pch=21, bg=ENMblock$occ.grp)
        box()
        grid()
        dev.off()
        
        
        ##mapa de distribuicao e suitability##
        
        
        jpeg(paste(projectFolder,'/SDM outputs/',sp_i,'_maps.jpg',sep=''), width=1000, height=500)
        par(mfrow=c(1,2), mar=c(3,3,4,8),  cex=1.3)
        ##binario
        plot(SDMpredBIN > thre, main=gsub('_', ' ', sp_i), legend=FALSE)
        plot(SAborders, add=T)
        grid()
        legend('topright',legend=c('non-habitat','habitat'), pch='', fill=c('lightgrey','darkgreen'), bg='white')
        ##continuo
        plot(SDMpredMean,main=gsub('_', ' ', sp_i))
        plot(cHull, add=TRUE)
        plot(SpsBuffer, ad=TRUE, lty=2)
        plot(SAborders, add=TRUE)
        grid()
        dev.off()

        
        
        ##salvando output do ENMeval
        
        SDMevalOutput = as.data.frame(SDMeval@results)
        write.csv(SDMevalOutput, paste(projectFolder,'/SDM outputs/',sp_i,'/',sp_i,'_SDMeval.csv',sep=''), row.names=FALSE)

        
        ##salvando graficos do ENMeval
        
        jpeg(paste(projectFolder,'/SDM outputs/',sp_i,'/',sp_i,'_SDMeval.jpg',sep=''), width=800)
        par(mfrow=c(1,2), mar=c(5,5,10,5))
        eval.plot(SDMeval@results); title(paste(sp_i))
        eval.plot(SDMeval@results, 'Mean.AUC', var='Var.AUC'); title(paste(sp_i))
        dev.off()

        
        ##importancia das variaveis

        
        modelPars = SDMeval@models[[bestModel$settings]]
        write.csv(var.importance(modelPars), paste(projectFolder,'/SDM outputs/',sp_i,'/',sp_i, '_variablesImportance.csv',sep=''), row.names=FALSE)

        
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

        write.csv(tabRes, paste(projectFolder,'/SDM outputs/tabOutputsSDMclim.csv',sep=''), row.names=FALSE)

    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}




## PARTE 3: Analise de similaridade de equivalencia de nicho




##abrindo variaveis ambientais
load(paste(envDataFolder,'/uncorrelatedPredictorsSouthAmerica.R', sep=''))
predictors = crop(x=predictors, y=SAborders)

##bias layer
biasLayer = raster(paste(projectFolder,'/biasLayerRodentia.asc', sep=''))
biasLayer = crop(x=biasLayer, y=SAborders)


##loop para SDM com cada uma das especies##

##nomes das especies
spList = gsub(pattern='.csv', replacement='', x=list.files(path=occDataFolder, pattern='.csv', full.names=FALSE))
spList = gsub(pattern=' ', replacement='_', x=spList)

##tabela de output
output = data.frame()


for (sp_i in spList){
    for (sp_j in spList){
        if (sp_i != sp_j){        
            tryCatch({

                ##diretorio base de trabalho
                setwd(paste(projectFolder,'/Niche overlap outputs',sep=''))

                ##verifica e cria diretorio para salvar resultados da especie atual
                if (file.exists(sp_i)){
                    setwd(sp_i)
                } else {
                    dir.create(sp_i)
                    setwd(sp_i)
                }
                
                ##dados de ocorrencia
                sp1Data = read.csv(paste(occDataFolder,'/',gsub(pattern='_',replacement=' ',x=sp_i),'.csv',sep=''), header=TRUE, sep=',', dec='.', na.strings='',colClasses=c('character','numeric','numeric')) #abrindo pontos de ocorrencia
                sp1Data = sp1Data[,2:3];
                names(sp1Data) = c('lon','lat')

                sp2Data = read.csv(paste(occDataFolder,'/',gsub(pattern='_',replacement=' ',x=sp_j),'.csv',sep=''), header=TRUE, sep=',', dec='.', na.strings='',colClasses=c('character','numeric','numeric')) #abrindo pontos de ocorrencia
                sp2Data = sp2Data[,2:3];
                names(sp2Data) = c('lon','lat')

                ##background points
                sp1BG = data.frame( randomPoints(mask=biasLayer, n=1000, p=sp1Data, prob=TRUE) ); names(sp1BG) = c('lon','lat')
                sp2BG = data.frame( randomPoints(mask=biasLayer, n=1000, p=sp2Data, prob=TRUE) ); names(sp2BG) = c('lon','lat')
                
                ##data.frames
                sp1DF = data.frame(rbind(sp1Data,sp1BG), pres=c(rep(1,nrow(sp1Data)), rep(0,nrow(sp1BG)) ))
                sp2DF = data.frame(rbind(sp2Data,sp2BG), pres=c(rep(1,nrow(sp2Data)), rep(0,nrow(sp2BG)) ))

                ##data.frame com variaveis ambientais
                sp1DFenv = extract(predictors, sp1DF[,c('lon','lat')], na.rm=TRUE); sp1DF = data.frame(sp1DF, sp1DFenv)
                sp2DFenv = extract(predictors, sp2DF[,c('lon','lat')], na.rm=TRUE); sp2DF = data.frame(sp2DF, sp2DFenv)

                ##limpando NAs
                sp1DF = sp1DF[complete.cases(sp1DF),]
                sp2DF = sp2DF[complete.cases(sp2DF),]

                ##The PCA is calibrated on all the sites of the study area
                pca.env <- dudi.pca(rbind(sp1DF,sp2DF)[,c(4:ncol(sp1DF))],scannf=F,nf=2)
                ##[,c("Bio13_Precipitation_of_Wettest_Month_","Bio3_Isothermality")]
                ##ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) #grafico

                ##PCA scores for the whole study area
                scores.globclim <- pca.env$li

                ##PCA scores for the species native distribution
                scores.sp.sp1 <- suprow(pca.env,sp1DF[which(sp1DF[,'pres']==1),c(4:ncol(sp1DF))])$li

                ##PCA scores for the species invasive distribution
                scores.sp.sp2 <- suprow(pca.env,sp2DF[which(sp2DF[,'pres']==1),c(4:ncol(sp2DF))])$li

                ##PCA scores for the whole native study area
                scores.clim.sp1 <-suprow(pca.env,sp1DF[,c(4:ncol(sp1DF))])$li

                ##PCA scores for the whole invaded study area
                scores.clim.sp2 <- suprow(pca.env,sp2DF[,c(4:ncol(sp2DF))])$li

                ##gridding the native niche
                grid.clim.sp1 <-ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.sp1, sp=scores.sp.sp1, R=100, th.sp=0)

                ##gridding the invasive niche
                grid.clim.sp2 <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.sp2, sp=scores.sp.sp2, R=100, th.sp=0)

                ##equivalencia de nicho
                ##OBS: Niche equivalency test H1: Is the overlap between the native and invaded niche higher than two random niches
                eq.test.equi <- ecospat.niche.equivalency.test(grid.clim.sp1, grid.clim.sp2, rep=100, alternative="greater")
                eq.test.simi <- ecospat.niche.similarity.test(grid.clim.sp1, grid.clim.sp2, rep=100, alternative="greater") #aleatoriza so uma das sp

                Dobs_equi= eq.test.equi$obs$D #indice D observado
                Iobs_equi = eq.test.equi$obs$I #indice I observado
                DpValue_equi = eq.test.equi$p.D #p-valor indice D
                IpValue_equi = eq.test.equi$p.I #p-valor indice I
                Dobs_simi= eq.test.simi$obs$D #indice D observado
                Iobs_simi = eq.test.simi$obs$I #indice I observado
                DpValue_simi = eq.test.simi$p.D #p-valor indice D
                IpValue_simi = eq.test.simi$p.I #p-valor indice I

                ##montando tabela de resultados e salvando
                output = rbind(output, data.frame(sp1=sp1,
                                                  sp2=sp2,
                                                  D.equi=Dobs_equi,
                                                  p_value.equi=DpValue_equi,
                                                  I.equi=Iobs_equi,
                                                  p_valor.equi=IpValue_equi,
                                                  D.simi=Dobs_simi,
                                                  p_value.simi=DpValue_simi,
                                                  I.simi=Iobs_simi,
                                                  p_valor.simi=IpValue_simi))


                ##salvando resultados no HD
                write.csv(output, paste(projectFolder,'/Niche overlap outputs/NicheOverlapOutputs.csv',sep=''), row.names=FALSE)

            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }
    }
}




## PARTE 4: Sinal Filogenetico




##pacote necessario
library(ape)


##cenarios
Scenario1<-read.tree("Scenario1.txt")
Scenario2<-read.tree("Scenario2.txt")
Filogenia<-read.tree("Arvore_analises.tre")
plot(Filogenia)


##labels
Scenario1$tip.label<-c("Cavia_fulgida","Kerodon_rupestris","Galea_spixii","Hydrochoerus_hydrochaeris","Galea_musteloides","Galea_leucoblephara","Microcavia_australis","Microcavia_niata","Cavia_tschudii","Hydrochoerus_isthmius","Dolichotis_salinicola","Dolichotis_patagonum","Cavia_aperea")
Scenario2$tip.label<-c("Cavia_aperea","Dolichotis_patagonum","Galea_leucoblephara","Microcavia_australis","Galea_musteloides","Cavia_tschudii","Dolichotis_salinicola","Microcavia_niata","Hydrochoerus_isthmius","Kerodon_rupestris","Galea_spixii","Cavia_fulgida","Hydrochoerus_hydrochaeris")


##distancias
dist_scenario1<-dist.topo(Scenario1, Filogenia, method = "PH85")
dist_scenario2<-dist.topo(Scenario2, Filogenia, method = "PH85")

table.nula <- matrix(nrow=1000,ncol=2)
colnames(table.nula) <- c('rodada','distancia')
head(table.nula)


##loop
for(i in 1:1000){ 
    especies <- c("Cavia_aperea","Dolichotis_patagonum","Galea_leucoblephara","Microcavia_australis","Galea_musteloides","Cavia_tschudii","Dolichotis_salinicola","Microcavia_niata","Hydrochoerus_isthmius","Kerodon_rupestris","Galea_spixii","Cavia_fulgida","Hydrochoerus_hydrochaeris")
    x <- sample(especies,size=13, replace=F)
    tree <- pbtree(n=13)
    tree$tip.label <- x
    dist_scenario1 <- dist.topo(tree, Filogenia, method = "PH85")
    table.nula[i,2] <- dist_scenario1
}


##estatisticas
hist(table.nula[,2])
columna_prob<-table.nula[,2]/sum(table.nula[,2])

mean(columna_prob[columna_prob<17])

quantile(table.nula[,2], 0.05)




## PARTE 5: procrustes




##pacotes necessarios
library (geomorph)
load(file="Semilandmark.RData")


##dados
Dentes <- readland.tps("Dentes_inteiros2.TPS", specID = c("ID"))
Peixe_GPA <- gpagen(Dentes, ProcD=F) # ProcD=T->Procrustes distancia criterio; ProcD=F->Binding Energy criterio
Peixe_GPA

slider1 <- define.sliders(c(1,8:12,2),1)
slider1
slider2 <- define.sliders(c(2,13:17,3),1)
slider3 <- define.sliders(c(3,18:27,4),1)
slider4 <- define.sliders(c(4,28:37,5),1)
slider5 <- define.sliders(c(5,38:47,6),1)
slider6 <- define.sliders(c(6,48:57,1),1)
slider7 <- define.sliders(c(1,58:67,7),1)
slider8 <- define.sliders(c(7,68:77,3),1)
curvas <- rbind(slider1,slider2,slider3,slider4,slider5,slider6,slider7,slider8)
curvas

Dentes_GPA2 <- gpagen(Dentes, curves=curvas,ProcD=F) # ProcD=T->Procrustes distancia criterio; ProcD=F->Binding Energy criterio
Dentes_GPA2


##graficos
plot(Dentes_GPA2)




## PARTE 6: PNO e DTT




##pacote necessario
library(geiger)


##caminhos para pastas
setwd("C:/Users/Usuário/Documents/otros proyectos/pablo")
setwd("C:/Users/Pibi Lab/Desktop/Filogenia_Alvarez2017/DTT")


##PNO Analysis
path_bioclim <- list.files('path_bioclim', pattern='.asc',full.names=T)
path_model <- "maxent_asc_raw"
var.names <- sapply(strsplit(basename(path_bioclim),"\\."), function(x) paste(x[1:(length(x)-1)], collapse="."))
sp.names <- sapply(strsplit(basename(list.files("maxent_asc_raw", pattern='.asc',full.names=T)),"\\."), 
                   function(x) paste(x[1:(length(x)-1)], collapse=".")) 


##loop sobre path_bioclim
for (i in 1:length(path_bioclim)){
    test <- pno(path_bioclim[i], path_model, bin_number = 100)
    colnames(test) <- c(var.names[i],sp.names)
    write.csv(test, file=paste(getwd(),"/PNO/",var.names[i],".csv",sep = ""),row.names=F)
    pdf(file=paste(getwd(),"/PNO/PNO Plot/",var.names[i],".pdf",sep = ""), height=5,width=6)
    plotPNO(test,xlab=var.names[i],wm=T)
    dev.off()
}


##plots
Bio2_Mean_Diurnal_Range <- read.csv("Bio2_Mean_Diurnal_Range.csv")
plotPNO(test,xlab="clacita",wm=F)

Bio3_Isothermality <- read.csv("Bio3_Isothermality.csv")
plotPNO(sal, subset = NULL, thinning = NULL, xlab = "salinity", 
        tail_threshold = 0, wm = FALSE, legend.pos = NULL)

tem_min <- read.csv("SSTmin2.csv")
plotPNO(tem_min, subset = NULL, thinning = NULL, xlab = "temp min", 
        tail_threshold = 0, wm = FALSE, legend.pos = NULL)

tem_range <- read.csv("SSTRange2.csv")
plotPNO(tem_range, subset = NULL, thinning = NULL, xlab = "temp range", 
        tail_threshold = 0, wm = FALSE, legend.pos = NULL)


##New tree
tree <- read.tree("Arvore_analises.tre")
plot(tree)

tree$tip.label

bio.list <- list.files('DTT', pattern='.csv', full.names=T)
bio.names <- sapply(strsplit(basename(bio.list),"\\."), function(x) paste(x[1:(length(x)-1)], collapse=".")) 

i=1
for (i in length(bio.list)){
    bio <- read.csv(bio.list[i],h=T)
    bio.anc <- anc.clim(target=tree, pno=bio, n=100, method = 'ML')
    bio.wm <- pno.weighted.mean(bio)
    pdf(paste('DTT/PNO/DTT/', bio.names[i],'.pdf',sep=''),height=10, width=10)
    dtt.bio <- dtt(tree,bio.wm,nsim=1000)
    dev.off()
    pdf(paste('DTT/PNO/PNO_plot/', bio.names[i],'.pdf',sep=''),height=10, width=10)
    plotAncClim(bio.anc, ylab=bio.names[i],cex=1,tipmode = 2)
    dev.off()
}




## PARTE 7: Procrustes - ANCOVA - MANCOVA




##pacotes necessarios
library (geomorph)
load(file="Semilandmark.RData")


##pasta de trabalho
setwd("C:/Users/Darlan/Desktop/Darlan/UFS/PIBI/Morfometria Dentes/Teste3(Denes_inteiros)")


##dados
Dentes<-readland.tps("Dentes.TPS", specID = c("ID"))
head(Dentes)

Dentes_GPA<-gpagen(Dentes, ProcD=F) # ProcD=T->Procrustes distancia criterio; ProcD=F->Binding Energy criterio
head(Dentes_GPA)

slider1 <- define.sliders(c(1,8:17,2),1)
slider2 <- define.sliders(c(2,18:27,3),1)
slider3 <- define.sliders(c(3,28:37,4),1)
slider4 <- define.sliders(c(4,38:47,5),1)
slider5 <- define.sliders(c(5,48:57,6),1)
slider6 <- define.sliders(c(6,58:67,1),1)
slider7 <- define.sliders(c(1,68:77,7),1)
slider8 <- define.sliders(c(7,78:87,3),1)
curvas <- rbind(slider1,slider2,slider3,slider4,slider5,slider6,slider7,slider8)

Dentes_GPA2 <- gpagen(Dentes, curves=curvas,ProcD=F) # ProcD=T->Procrustes distancia criterio; ProcD=F->Binding Energy criterio
plot(Dentes_GPA2)


##PCA

Cs <- Dentes_GPA2$Csize# Cs #centróide
classifier_Dentes <- read.table("Classificatorias.txt", header=T) # para entrar una variable de grupo
##não precisa de tabela pra fazer o PCA, só se quiser diferenciar os pontinhos dos gráficos para cada espécie
names(classifier_Dentes)
Dentes_Sp <- classifier_Dentes[,3] #a coluna 3, onde estão as espécies, são as variáveis 

Dentes2D <- two.d.array(Dentes_GPA2$coords) 

##grafico
plotTangentSpace(Dentes_GPA2$coords, axis1 = 1, label=T, axis2 = 2, mesh=T,warpgrids = TRUE, groups = Dentes_Sp, verbose = FALSE, legend = TRUE)
##faz os componentes principais PCAs
ref <- mshape(Dentes_GPA2$coords) #média da forma, o indivíduo no centro da figura, um modelo que serve de comparação, de referência, indivíduo médio
plotRefToTarget(ref,Dentes_GPA2$coords[,,1],method="TPS") #indivíduo que eu quero comparar com a referência
plotRefToTarget(ref,Dentes_GPA2$coords[,,1],method="vector") #quanto mais comprida a seta, maior a mudança


##indica a direção e a intensidade da deformação: quanto maior o vetor, mais intenso
plotRefToTarget(ref,Dentes_GPA2$coords[,,1],method="point")
##esse pacote não faz variáveis canônicas - as variáveis canônicas maximizam as diferenças

##ANOVA de Procrustes ou MANOVA
Dif_especies<-procD.lm(Dentes2D~Dentes_Sp, iter = 9999,RRPP = F) #Dentes_Sp são os nomes das espécies

##p=1e-04 *** - indica que as espécies são muito diferentes
##Teste de tukey
tukey <- TukeyHSD(Dif_especies)

##MANCOVA
Dif_especiesII<-procD.lm(Dentes2D~Dentes_Sp*Cs, iter = 9999,RRPP = F) #interação entre as espécies
##quer saber se a espécie A e a espécie B variam de forma diferente, se a forma varia em relação ao tamanho
##a covariância da forma nas variadas espécies
summary(Dif_especiesII)




## PARTE 8: Dist. geo. - Dist. Patristica




##pacotes necessarios
library(adephylo)
library(geosphere)
library(ape)


##diretorio de trabalho da analise
setwd("C:/Users/Darlan/Desktop/Darlan/UFS/PIBI/Caviidae/Caviidae_modelo_nicho/Sinal Geografico/Teste 1")


##abrir arvore
arvore <- read.tree("Arvore_Cavi_Alvarez.txt")
plot(arvore)

##analise patristic
patristic <- distTips(arvore, patristic)
patristic


##abrir registros de ocorrencia para cada espécie
Cavia_aperea <- read.csv("Cavia aperea.csv",header=T)
Cavia_fulgida <- read.csv("Cavia fulgida.csv", header=T)
Cavia_tschudii <- read.csv("Cavia tschudii.csv", header=T)
Dolichotis_patagonum <- read.csv("Dolichotis patagonum.csv", header=T)
Dolichotis_salinicola <- read.csv("Dolichotis salinicola.csv", header=T)
Galea_leucoblephara <- read.csv("Galea leucoblephara.csv", header=T)
Galea_musteloides <- read.csv("Galea musteloides.csv", header=T)
Galea_spixii <- read.csv("Galea spixii.csv", header=T)
Hydrochoerus_hydrochaeris <- read.csv("Hydrochoerus hydrochaeris.csv", header=T)
Hydrochoerus_isthmius <- read.csv("Hydrochoerus isthmius.csv", header=T)
Kerodon_rupestris <- read.csv("Kerodon rupestris.csv", header=T)
Microcavia_australis <- read.csv("Microcavia australis.csv", header=T)
Microcavia_niata <- read.csv("Microcavia niata.csv", header=T)
head(Cavia_aperea)


##calcular o centroide geográfico para cada espécie
Cen_C.aperea <- centroid(Cavia_aperea[,2:3])
Cen_C.fulgida <- centroid(Cavia_fulgida[,2:3])
Cen_C.tschudii <- centroid(Cavia_tschudii[,2:3])
Cen_D.patagonum <- centroid(Dolichotis_patagonum[,2:3])
Cen_D.salinicola <- centroid(Dolichotis_salinicola[,2:3])
Cen_G.leucoblephara <- centroid(Galea_leucoblephara[,2:3])
Cen_G.musteloides <- centroid(Galea_musteloides[,2:3])
Cen_G.spixii <- centroid(Galea_spixii[,2:3])
Cen_H.hydrochaeris <- centroid(Hydrochoerus_hydrochaeris[,2:3])
Cen_H.isthmius <- centroid(Hydrochoerus_isthmius[,2:3])
Cen_K.rupestris <- centroid(Kerodon_rupestris[,2:3])
Cen_M.australis <- centroid(Microcavia_australis[,2:3])
Cen_M.niata <- centroid(Microcavia_niata[,2:3])


##salvar tudo em uma tabela
tabela <- rbind(Cen_C.aperea,Cen_C.fulgida,Cen_C.tschudii,Cen_D.patagonum,Cen_D.salinicola,Cen_G.leucoblephara,Cen_G.musteloides,Cen_G.spixii,Cen_H.hydrochaeris,Cen_H.isthmius,Cen_K.rupestris,Cen_M.australis,Cen_M.niata)
tabela


##nomear linhas
rownames(tabela)<-c("Cavia_aperea","Cavia_fulgida","Cavia_tschudii","Dolichotis_patagonum","Dolichotis_salinicola","Galea_leucoblephara","Galea_musteloides","Galea_spixii","Hydrochoerus_hydrochaeris","Hydrochoerus_isthmius","Kerodon_rupestris","Microcavia_australis","Microcavia_niata")
head(tabela)


##analise...
test <- distm(tabela, fun = distHaversine)
head(test)


##nomear linhas
rownames(test) <- c("Cavia_aperea","Cavia_fulgida","Cavia_tschudii","Dolichotis_patagonum","Dolichotis_salinicola","Galea_leucoblephara","Galea_musteloides","Galea_spixii","Hydrochoerus_hydrochaeris","Hydrochoerus_isthmius","Kerodon_rupestris","Microcavia_australis","Microcavia_niata")


##nomear colunas
colnames(test)<-c("Cavia_aperea","Cavia_fulgida","Cavia_tschudii","Dolichotis_patagonum","Dolichotis_salinicola","Galea_leucoblephara","Galea_musteloides","Galea_spixii","Hydrochoerus_hydrochaeris","Hydrochoerus_isthmius","Kerodon_rupestris","Microcavia_australis","Microcavia_niata")
head(test)
plot(test)
write.csv(test, "resultado.csv")
