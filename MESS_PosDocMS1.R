###MESS para os dados do primeiro manuscrito do pos-doc UFS 2016-2019###
library(raster)
library(dismo)
library(maptools)

###detectando climas nao analogos com funcao enseble.analogue do pacote BiodiveristyR###
##TUTURIAL##
# get predictor variables

library(BiodiversityR)

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
library(dismo)
library(BiodiversityR)

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
predictors@name = 'base'

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
## End(Not run)

