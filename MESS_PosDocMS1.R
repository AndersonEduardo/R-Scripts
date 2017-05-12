###MESS para os dados do primeiro manuscrito do pos-doc UFS 2016-2019###
library(raster)
library(dismo)
library(maptools)
library(BiodiversityR)

###detectando climas nao analogos com funcao enseble.analogue do pacote BiodiveristyR###
##EXEMPLO DO TUTURIAL##

## get predictor variables
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
predictors <- subset(predictors, subset=c("bio1", "bio5", "bio6", "bio7", "bio8", "bio12", "bio16", "bio17"))
predictors
predictors@title <- "base"

# instead of searching for current analogue of future climate conditions,
# search for analogue in southern hemisphere
future.stack <- stack(crop(predictors, y=extent(-125, -32, 0, 40)))
future.stack@title <- "north"
current.stack <- stack(crop(predictors, y=extent(-125, -32, -56, 0)))
current.stack@title <- "south"

# reference location in Florida
# in this case future.stack and current.stack are both current
ref.loc <- data.frame(t(c(-80.19, 25.76)))
names(ref.loc) <- c("lon", "lat")

# climate analogue analysis based on the Mahalanobis distance
Florida.object.mahal <- ensemble.analogue.object(ref.location=ref.loc, future.stack=future.stack, current.stack=current.stack, name="FloridaMahal", method="mahal", an=10000)
Florida.object.mahal

Florida.analogue.mahal <- ensemble.analogue(x=current.stack, analogue.object=Florida.object.mahal, analogues=50)
Florida.analogue.mahal

# climate analogue analysis based on the Euclidean distance and dividing each variable by the sd
Florida.object.sd <- ensemble.analogue.object(ref.location=ref.loc,future.stack=future.stack, current.stack=current.stack,name="FloridaSD", method="sd", z=2)
Florida.object.sd

Florida.analogue.sd <- ensemble.analogue(x=current.stack,analogue.object=Florida.object.sd, analogues=50)
Florida.analogue.sd

# plot analogues on climatic distance maps
par(mfrow=c(1,2))
analogue.file <- paste(getwd(), "//ensembles//analogue//FloridaMahal_south_analogue.grd", sep="")
plot(raster(analogue.file), main="Mahalanobis climatic distance")
points(Florida.analogue.sd[3:50, "lat"] ~ Florida.analogue.sd[3:50, "lon"], 
    pch=1, col="red", cex=1)
points(Florida.analogue.mahal[3:50, "lat"] ~ Florida.analogue.mahal[3:50, "lon"], 
    pch=3, col="black", cex=1)
points(Florida.analogue.mahal[2, "lat"] ~ Florida.analogue.mahal[2, "lon"], 
    pch=22, col="blue", cex=2)
legend(x="topright", legend=c("closest", "Mahalanobis", "SD"), pch=c(22, 3 , 1), 
    col=c("blue" , "black", "red"))

analogue.file <- paste(getwd(), "//ensembles//analogue//FloridaSD_south_analogue.grd", sep="")
plot(raster(analogue.file), main="Climatic distance normalized by standard deviation")
points(Florida.analogue.mahal[3:50, "lat"] ~ Florida.analogue.mahal[3:50, "lon"], 
    pch=3, col="black", cex=1)
points(Florida.analogue.sd[3:50, "lat"] ~ Florida.analogue.sd[3:50, "lon"], 
    pch=1, col="red", cex=1)
points(Florida.analogue.sd[2, "lat"] ~ Florida.analogue.sd[2, "lon"], 
    pch=22, col="blue", cex=2)
legend(x="topright", legend=c("closest", "Mahalanobis", "SD"), pch=c(22, 3 , 1), 
    col=c("blue" , "black", "red"))
par(mfrow=c(1,1))


##primeiro manuscrito do pos-doc UFS 2016-2019##

## get predictor variables

##definindo as pastas de trabalho
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/"
spOccFolder = "/home/anderson/PosDoc/dados_ocorrencia/PO_unique/"
projectFolder = "/home/anderson/PosDoc/teste/"

##abrindo shape da America do Sul
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp")

##abrindo as variaveis climaticas
##abrindo e cortando camadas de variaveis ambientais para o presente
filesRaw <- stack(list.files(path=paste(envVarFolder,"dados_projeto/000",sep=''), pattern='asc', full.names=TRUE))
files = mask(filesRaw,AmSulShape) #cortando para Am. do Sul
files.crop.sub = files[[c('bioclim_10','bioclim_11','bioclim_16','bioclim_17')]] #choose selected layers

predictors <- stack(files.crop.sub)
crs(predictors) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
predictors@title = 'base'

########## Criando objetos com a lista de especies #############
occ.sps <- list.files(paste(spOccFolder,sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))
##fosseis
occ.sps.fosseis = read.csv(paste(spOccFolder,'/fosseis/',"fosseis.csv",sep=''),header=T)
splist.fosseis = unlist(lapply(occ.sps.fosseis[,1],as.character))

for (i in 1:length(splist)){
    ##abrindo os dados de registros fosseis  para uma especie
    especie = i
    sp.fossil.data = occ.sps.fosseis[occ.sps.fosseis$species==splist[especie],] #ATENCAO: este script nao funciona se houver mais de um registro fossil por camada de tempo usada para projecao

    for(l in 1:nrow(sp.fossil.data)){#loop para cada registro fossil de uma especie

        ##definindoo fossil
        sp.fossil = sp.fossil.data[l,]

        ##abrindo as variaveis ambientais do tempo do fossil
        filesPastRaw <- stack(list.files(path = paste(envVarFolder,"dados_projeto/0",sp.fossil$kyr,sep=""), pattern='asc', full.names=T)) ###abrindo camandas para projecao (passado, futuro, outro local, etc)
        filesPast = mask(filesPastRaw,AmSulShape) #cortando para Am. do Sul
        files.crop.sub.past <- filesPast[[c('bioclim_10','bioclim_11','bioclim_16','bioclim_17')]]
        predictorsPast = files.crop.sub.past #preditoras para o tempo do fossil
        crs(predictorsPast) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

        ## reference location (occ da especie atual, no loop)
        sp.file <- read.csv(paste(spOccFolder, splist[especie],".csv",sep=""),h=T) ### read sp occurrence
        sp.occ <- sp.file[,2:3]
        ref.loc <- sp.occ

        ###CONTINUAR APARTIR DAQUI###
        
        ## climate analogue analysis based on the Mahalanobis distance
        object.mahal <- ensemble.analogue.object(ref.location=ref.loc, future.stack=predictorsPast, current.stack=predictors, name=paste(sp.fossil$species,sp.fossil$kyr,'kyr',sep=''), method="mahal", an=10000)
        xx = object.mahal$target.values[1,]

        analogue.mahal <- ensemble.analogue(x=predictors, analogue.object=object.mahal, analogues=50)

    }
}

        
# instead of searching for current analogue of future climate conditions,
# search for analogue in the past
past.stack.Ccrocodilus11k <- stack(crop(predictors, y=extent(-125, -32, 0, 40)))
future.stack@title <- "north"
current.stack <- stack(crop(predictors, y=extent(-125, -32, -56, 0)))
current.stack@title <- "south"


        



# climate analogue analysis based on the Euclidean distance and dividing each variable by the sd
Florida.object.sd <- ensemble.analogue.object(ref.location=ref.loc,future.stack=future.stack, current.stack=current.stack,name="FloridaSD", method="sd", z=2)
Florida.object.sd

Florida.analogue.sd <- ensemble.analogue(x=current.stack,analogue.object=Florida.object.sd, analogues=50)
Florida.analogue.sd

# plot analogues on climatic distance maps
par(mfrow=c(1,2))
analogue.file <- paste(getwd(), "//ensembles//analogue//FloridaMahal_south_analogue.grd", sep="")
plot(raster(analogue.file), main="Mahalanobis climatic distance")
points(Florida.analogue.sd[3:50, "lat"] ~ Florida.analogue.sd[3:50, "lon"], 
    pch=1, col="red", cex=1)
points(Florida.analogue.mahal[3:50, "lat"] ~ Florida.analogue.mahal[3:50, "lon"], 
    pch=3, col="black", cex=1)
points(Florida.analogue.mahal[2, "lat"] ~ Florida.analogue.mahal[2, "lon"], 
    pch=22, col="blue", cex=2)
legend(x="topright", legend=c("closest", "Mahalanobis", "SD"), pch=c(22, 3 , 1), 
    col=c("blue" , "black", "red"))

analogue.file <- paste(getwd(), "//ensembles//analogue//FloridaSD_south_analogue.grd", sep="")
plot(raster(analogue.file), main="Climatic distance normalized by standard deviation")
points(Florida.analogue.mahal[3:50, "lat"] ~ Florida.analogue.mahal[3:50, "lon"], 
    pch=3, col="black", cex=1)
points(Florida.analogue.sd[3:50, "lat"] ~ Florida.analogue.sd[3:50, "lon"], 
    pch=1, col="red", cex=1)
points(Florida.analogue.sd[2, "lat"] ~ Florida.analogue.sd[2, "lon"], 
    pch=22, col="blue", cex=2)
legend(x="topright", legend=c("closest", "Mahalanobis", "SD"), pch=c(22, 3 , 1), 
    col=c("blue" , "black", "red"))
par(mfrow=c(1,1))
## End(Not run)



#####################################################################
###resultados da analise de climas nao analogos a partir do maxent###
#####################################################################

caminho = '/home/anderson/Área de Trabalho/climas_nao-analogos/ascii_files'
rasters = stack(list.files(path=caminho,pattern='.asc',full.names=TRUE))
rastersCortados = mask(rasters,AmSulShape)
writeRaster(rastersCortados,filename=paste(caminho,'/cortados_AS/',names(rastersCortados),'.asc',sep=''),bylayer=TRUE,overwrite=TRUE)

library(rasterVis)

##11kyr
spaceisLayers11 = stack(list.files(path=paste(caminho,'/cortados_AS',sep=''),pattern='11',full.names=TRUE))
nomesSubgraficos = c('C. crocodilus (clamping)','C. crocodilus (MESS)','C. crocodilus (MoD)',
                     'C. latirostris (clamping)','C. latirostris (MESS)','C. latirostris (MoD)',
                     'C. yacare (clamping)','C. yacare (MESS)','C. yacare (MoD)',
                     'M. niger (clamping)','M. niger (MESS)','M. niger (MoD)')

jpeg(file='/home/anderson/Área de Trabalho/climas_nao-analogos/11kyrClimateAnalysis.jpg', width = 1200, height = 1200)
levelplot(spaceisLayers11,scales=list(x=list(cex=1.5), y=list(cex=1.5)),between=list(x=1.8, y=0.25),par.strip.text=list(cex=1.5),layout=c(3,4), main=list(label='11 kyr BP',cex=3),names.attr=nomesSubgraficos,colorkey=list(space="right",labels=list(cex=2))) + layer(sp.polygons(AmSulShape))
dev.off()

##21kyr
spaceisLayers21 = stack(list.files(path=paste(caminho,'/cortados_AS',sep=''),pattern='21',full.names=TRUE))
nomesSubgraficos = c('C. crocodilus (clamping)','C. crocodilus (MESS)','C. crocodilus (MoD)',
                     'C. latirostris (clamping)','C. latirostris (MESS)','C. latirostris (MoD)',
                     'C. yacare (clamping)','C. yacare (MESS)','C. yacare (MoD)',
                     'M. niger (clamping)','M. niger (MESS)','M. niger (MoD)')

jpeg(file='/home/anderson/Área de Trabalho/climas_nao-analogos/21kyrClimateAnalysis.jpg', width = 1200, height = 1200)
levelplot(spaceisLayers21,scales=list(x=list(cex=1.5), y=list(cex=1.5)),between=list(x=1.8, y=0.25),par.strip.text=list(cex=1.5),layout=c(3,4), main=list(label='21 kyr BP',cex=3),names.attr=nomesSubgraficos,colorkey=list(space="right",labels=list(cex=2))) + layer(sp.polygons(AmSulShape))
dev.off()

##13kyr
spaceisLayers13 = stack(list.files(path=paste(caminho,'/cortados_AS',sep=''),pattern='13',full.names=TRUE))
nomesSubgraficos = c('M. coypus (clamping)','M. coypus (MESS)','M. coypus (MoD)')

jpeg(file='/home/anderson/Área de Trabalho/climas_nao-analogos/13kyrClimateAnalysis.jpg', width = 1200, height = 1200)
levelplot(spaceisLayers13,scales=list(x=list(cex=1.5), y=list(cex=1.5)),between=list(x=1.8, y=0.25),par.strip.text=list(cex=2.5),layout=c(3,1), main=list(label='13 kyr BP',cex=4),names.attr=nomesSubgraficos,colorkey=list(space="right",labels=list(cex=2))) + layer(sp.polygons(AmSulShape))
dev.off()

##14kyr
spaceisLayers14 = stack(list.files(path=paste(caminho,'/cortados_AS',sep=''),pattern='14',full.names=TRUE))
nomesSubgraficos = c('M. coypus (clamping)','M. coypus (MESS)','M. coypus (MoD)')

jpeg(file='/home/anderson/Área de Trabalho/climas_nao-analogos/14kyrClimateAnalysis.jpg', width = 1200, height = 1200)
levelplot(spaceisLayers14,scales=list(x=list(cex=1.5), y=list(cex=1.5)),between=list(x=1.8, y=0.25),par.strip.text=list(cex=2.5),layout=c(3,1), main=list(label='14 kyr BP',cex=4),names.attr=nomesSubgraficos,colorkey=list(space="right",labels=list(cex=2))) + layer(sp.polygons(AmSulShape))
dev.off()

##19kyr
spaceisLayers19 = stack(list.files(path=paste(caminho,'/cortados_AS',sep=''),pattern='19',full.names=TRUE))
nomesSubgraficos = c('L. maximus (clamping)','L. maximus (MESS)','L. maximus (MoD)')

jpeg(file='/home/anderson/Área de Trabalho/climas_nao-analogos/19kyrClimateAnalysis.jpg', width = 1200, height = 1200)
levelplot(spaceisLayers19,scales=list(x=list(cex=1.5), y=list(cex=1.5)),between=list(x=1.8, y=0.25),par.strip.text=list(cex=2.5),layout=c(3,1), main=list(label='19 kyr BP',cex=4),names.attr=nomesSubgraficos,colorkey=list(space="right",labels=list(cex=2))) + layer(sp.polygons(AmSulShape))
dev.off()

##20kyr
spaceisLayers20 = stack(list.files(path=paste(caminho,'/cortados_AS',sep=''),pattern='19',full.names=TRUE))
nomesSubgraficos = c('L. maximus (clamping)','L. maximus (MESS)','L. maximus (MoD)')

jpeg(file='/home/anderson/Área de Trabalho/climas_nao-analogos/20kyrClimateAnalysis.jpg', width = 1200, height = 1200)
levelplot(spaceisLayers20,scales=list(x=list(cex=1.5), y=list(cex=1.5)),between=list(x=1.8, y=0.25),par.strip.text=list(cex=2.5),layout=c(3,1), main=list(label='20 kyr BP',cex=4),names.attr=nomesSubgraficos,colorkey=list(space="right",labels=list(cex=2))) + layer(sp.polygons(AmSulShape))
dev.off()
