library("spocc")
library(raster)
library("phyloclim")
library(maptools)
library(dismo)
library(ape)
library(geiger)
library(rgdal)
library(rgeos)
library(GISTools)
library(randomForest)
Sys.setenv(JAVA_HOME='/usr/lib/jvm/java-7-openjdk-amd64') # for 64-bit version
#Windows#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_91') # for 64-bit version
library(rJava)

##DEFININDO PASTAS DE TRABALHO##
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/"
spOccFolder = "/home/anderson/PosDoc/dados_ocorrencia/PO_unique/"
projectFolder = "/home/anderson/PosDoc/teste/"

####ABRINDO AS VARIAVEIS CLIMATICAS#####
#abrindo shape da America do Sul
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp")

#abrindo e cortando camads de variaveis ambientais para o presente
filesRaw <- stack(list.files(path=paste(envVarFolder,"dados_projeto/000",sep=''), pattern='asc', full.names=T)) ### stack all rasters in Bioclim folder
#files <- stack(list.files(path = "/home/anderson/R/PosDoc/dados_ambientais/bcmidbi_2-5m _asc/dados_ambientais_para_projeto", pattern='asc', full.names=T))
files = mask(filesRaw,AmSulShape) #cortando para Am. do Sul

##abrindo e cortando camads de variaveis ambientais para o passado
filesProjectionRaw <- stack(list.files(path=paste(envVarFolder,"dados_projeto/021",sep=''), pattern='asc', full.names=T)) ###abrindo camandas para projecao (passado, futuro, outro local, etc)
filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul

##testando correcaloes
##test<-getValues(files)
##cor.matrix <- as.data.frame(cor(test, use="complete.obs"))
##write.csv(cor.matrix,'cor_matrix.csv')

#remove highly correlated variables Bio1,Bio3,Bio9,Bio13,Bio14
files.crop.sub <- dropLayer(files, c(1,2,5,6)) #### remove selected layers
files.crop.sub.projection <- dropLayer(filesProjection, c(1,2,5,6))

#remover as mesmas camadas dos dados para projecao
#test2<-getValues(files.crop.sub)
#cor.matrix2<- cor(test2, use="complete.obs")
#write.csv(cor.matrix2,'cor.matrix2.csv')

#definindo os objetos para as variaveis preditoras
predictors <- files.crop.sub
predictorsProjection = files.crop.sub.projection
rm(files,files.crop.sub)

####DOWLOAD DADOS DE OCORRENCIA DAS ESPECIES####
vetor_de_nomes = c('Lagostomus maximus')

for (i in 1:length(vetor_de_nomes)) {
#for (i in 1:length(unique(Especies$Especie))) {
  #sp<-as.character(Especies[i,1])
  sp <- vetor_de_nomes[i]
  res <- occ(query = sp, from = 'gbif', limit = 10000)
  locs<-occ2df(res)
  #locs = read.csv(file=paste("/home/anderson/PosDoc/dados_ocorrencia/",sp,".csv",sep=""),header=T,stringsAsFactors=FALSE)    
  locs2<-locs[,2:3]
  locs2[locs2 == 0] <- NA
  locs3<-locs2[complete.cases(locs2),]
  locs4<-round(locs3, digits = 4)
  locs5<-locs4[!duplicated(locs4), ]
  locs6<-cbind(sp,locs5) 
  #write.csv(locs6, file=paste("/home/anderson/R/PosDoc/dados_ocorrencia/",unique(Especies$Especie)[i],".csv",sep=""),row.names=F)
  write.csv(locs6, file=paste("/home/anderson/PosDoc/dados_ocorrencia/PO_unique/",sp,".csv",sep=""),row.names=F)
}

########## Criando objetos com a lista de especies #############
occ.sps <- list.files(paste(spOccFolder,sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))
##fosseis
occ.sps.fosseis = read.csv(paste(spOccFolder,'/fosseis/',"fosseis.csv",sep=''),header=T)
splist.fosseis = lapply(occ.sps.fosseis[,1],as.character)

###inspeção visual exploratoria
data(wrld_simpl)

especie = splist[1] #escolher qual especie
sp.file <- read.csv(paste(spOccFolder,especie,".csv",sep=""),h=T) ### read sp occurrence
sp.occ <- sp.file[,2:3] ## select lat long columns
sp.occ = sp.occ[sp.occ$latitude < 0,] #eliminar pontos
sp.occ = sp.occ[sp.occ$longitude > -80,] #eliminar pontos 
plot(wrld_simpl, xlim=c(-80,10), ylim=c(-60,10), axes=TRUE, col= 'light yellow' )
box() # restore the box around the map
# plot points
points(sp.occ, col= ' orange ' , pch=20, cex=0.75)
#plot points again to add a border, for better visibility
points(sp.occ, col= ' red ' , cex=0.75)

sp.occ<-cbind(especie,sp.occ) 
write.csv(sp.occ, file=paste("/home/anderson/R/PosDoc/dados_ocorrencia/PO_unique/",especie,".csv",sep=""),row.names=F)


######################## MAXENT #############################
##################### exploratorio  #########################

options(java.parameters = "-Xmx7g") ###set available memmory to java

for (i in 1:length(splist)){
especie = i #escolher qual especie
#sp.file=read.table("CulpeoBergmann.txt",header=T)
#sp.file <- read.csv(paste(spOccFolder, splist[especie],".csv",sep=""),h=T) ### read sp occurrence
#sp.occ <- sp.file[,2:3] ## select lat long columns
me <- maxent(predictors,sp.occ, args=c("raw","maximumiterations=1000"),path=paste("/home/anderson/PosDoc/pablo/culpeos")) ## run maxent model with raw output
crs <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
r <- predict(me,predictors,crs=crs)
#r <- predict(me,predictorsProjection,crs=crs) #para projetar passado, futuro, etc
names(r) <- splist[especie]
writeRaster(r, filename=paste("/home/anderson/PosDoc/pablo/culpeosPass.asc",sep=""), overwrite=T)
}

######Verificando pontos dos fosseis#########
setwd("/home/anderson/R/PosDoc/teste/Myocastor coypus/raster layers")
mapaAsc = raster("Myocastor coypus - 22mil.asc")
setwd("/home/anderson/R/PosDoc/teste/Myocastor coypus")
jpeg(filename=paste(splist[1],'.jpg',sep=''))
#plot(r)
plot(mapaAsc)
plot(wrld_simpl, add=TRUE, xlim=c(-80,10), ylim=c(-60,10), axes=TRUE)
box() # restore the box around the map
# plot points
points(sp.occ, col= 'blue' , pch=20, cex=0.5)
# plot points again to add a border, for better visibility
points(sp.occ, col= ' red ' , cex=0.5)
#plotando os pontos dos fosseis
points(-41.553056,-12.393417, col='blue', pch=2, cex=1.1)
dev.off()


### GRAFICOS AGRUPADOS ###

library(rasterVis)

#definindo objeto com os nomes
teste = 'Random Forest'

#presente
setwd(paste(projectFolder,teste,'/Raster Layers',sep='')) 
files = list.files(paste(getwd()),full.names=TRUE,pattern='.asc')

##TSS pela minha funcao
source("/home/anderson/R/R-Scripts/TSSfunction.R")
##TSS pela funcao do Pablo Riul
#source("/home/anderson/R/R-Scripts/TSSmaxent.R")
##mapa binario
teste = 'Maxent'
setwd(paste(projectFolder,teste,'/Raster Layers',sep='')) 
files = list.files(paste(getwd()),full.names=TRUE,pattern='.asc')
files = files[c(1,3,5,7,9,11)]
species.layers = stack(files)

rasterNames = gsub("_"," ",names(species.layers))
rasterNames = gsub("BINARIO","",rasterNames)

##planilha de ocorrencia da especie
tssTable=data.frame()
for (i in 1:length(rasterNames)){
    especie = rasterNames[i]
    sp.file <- read.csv(paste(spOccFolder,especie,".csv",sep=""),h=TRUE) ### read sp occurrence
    sp.file = sp.file[,2:3]
    ##TSS
    tssValue = TSSfunction(species.layers[[names(species.layers)[i]]]>0,sp.file)
    tssTable = rbind(tssTable,cbind(teste,especie,tssValue))
}

raster.layer = stack(list.files(paste(projectFolder,teste,'/Passado/Raster Layers',sep=''),full.names=TRUE,pattern='Melanosuchus'))
pontos = data.frame(lon=c(-41.553056,-37.753611), lat=c(-12.393417,-9.926944))

suitability = extract(x=raster.layer,y=pontos[2,])


#caimans + Melnosuchus
files=c(files[1],files[3],files[5],files[9]) 
#lagostomus e myocastor
files=c(files[7],files[11]) 

#passado
#setwd(paste(projectFolder,teste,'/Passado/Raster Layers',sep=''))
filesPass = list.files(paste(projectFolder,teste,'/Passado/Raster Layers',sep=''),full.names=TRUE,pattern='asc')

 #caimans (completo)
filesPass=c(filesPass[1],filesPass[3],filesPass[5],filesPass[7],filesPass[9],filesPass[11],filesPass[13],filesPass[15],filesPass[17],filesPass[19],filesPass[21],filesPass[23],filesPass[29],filesPass[31],filesPass[33],filesPass[35]) 
#caimans (principal)
filesPass=c(filesPass[3],filesPass[5],filesPass[11],filesPass[13],filesPass[19],filesPass[21],filesPass[30],filesPass[31])
#L. maximus e M. coypus
filesPass=c(filesPass[25],filesPass[27],filesPass[37],filesPass[39])

#empilhando os rasters (passado e presente)
species.layers = stack(c(files,filesPass))

#organizando para caimans (completo)
species.layers=species.layers[[ c(1,2,3,4, 5,9,13,17, 6,10,14,18, 7,11,15,19, 8,12,16,20) ]]
#organizando para caimans (principal)
species.layers=species.layers[[ c(1,2,3,4, 5,7,9,11, 6,8,10,12) ]]
#organizando para L. maximus e M. coypus
species.layers=species.layers[[ c(1,2,3,5,4,6) ]]

#cortando e gravando varios rasters com o shapefile da america do sul
for (i in 1:length(names(species.layers))){
    species.layers[[i]] = mask(species.layers[[i]],AmSulShape)
    writeRaster(species.layers[[i]], filename=paste("/home/anderson/R/PosDoc/teste/",teste,'/Raster Layers Cortados/',names(species.layers)[i],".asc",sep=""), overwrite=T)
}

## #ajustando os nomes dos subgraficos
## rasterNames=names(species.layers)
## names(species.layers)
## rasterNames = gsub("Caiman_latirostris_._","",names(species.layers))
## rasterNames = gsub("presente","Presente",rasterNames)

## levelplot ##
setwd(paste("/home/anderson/PosDoc/teste/",teste,sep=''))

## Genero Caiman completo
nomesSubgraficos = c("C. crocodilus","C. latirostris","C. yacare","M. niger","10 kyr BP","10 kyr BP","10 kyr BP","10 kyr BP","11 kyr BP","11 kyr BP","11 kyr BP","11 kyr BP","21 kyr BP","21 kyr BP","21 kyr BP","21 kyr BP","22 kyr BP","22 kyr BP","22 kyr BP","22 kyr BP")
#pdf(file='CaimanCompleto.pdf')
jpeg(file='CaimanCompletoTESTE.jpg', width = 1200, height = 1200)
levelplot(species.layers,scales=list(x=list(cex=1), y=list(cex=1)),between=list(x=1.8, y=0.25),par.strip.text=list(cex=1.8),layout=c(4,5),col.regions=colorRampPalette(c("white","blue","green","yellow","red")), main='',names.attr=nomesSubgraficos,colorkey=list(space="right",labels=list(cex=1.2))) + layer(sp.polygons(AmSulShape)) + layer(panel.xyplot(-41.553056, -12.393417,pch=17,col='red',cex=1.6),rows=c(2,3)) + layer(panel.xyplot(-37.753611,-9.926944,pch=17,col='red',cex=1.6),rows=c(4,5))
dev.off()

## Genero Caiman (figura principal)
nomesSubgraficos = c("C. crocodilus","C. latirostris","C. yacare","M. niger","11 kyr BP","11 kyr BP","11 kyr BP","11 kyr BP","21 kyr BP","21 kyr BP","21 kyr BP","21 kyr BP")
##pdf(file='CaimanPrincipal.pdf')
jpeg(file='CaimanPrincipalTESTE.jpg', width = 1200, height = 1200)
levelplot(species.layers,scales=list(x=list(cex=1.3), y=list(cex=1.3)),between=list(x=1, y=0.25),par.strip.text=list(cex=2.5),layout=c(4,3),col.regions=colorRampPalette(c("white","blue","green","yellow","red")), main='', names.attr=nomesSubgraficos, colorkey=list(space="right",labels=list(cex=1.8))) + layer(sp.polygons(AmSulShape)) + layer(panel.xyplot(-41.553056, -12.393417,pch=17,col="red",cex=2),rows=c(2)) + layer(panel.xyplot(-37.753611,-9.926944,pch=17,col="red",cex=2),rows=c(3)) + layer(panel.xyplot(-53.283333,-33.683333,pch=17,col='red',cex=2),rows=4)
dev.off()

## M. coypus L. maximus
nomesSubgraficos = c("L. maximus","M. coypus","13 kyr BP","19 kyr BP","14 kyr BP","20 kyr BP")
##pdf(file='MyoLago.pdf')
jpeg(file='MyoLago.jpg', width = 1200, height = 1200)
levelplot(species.layers,scales=list(x=list(cex=1.5), y=list(cex=1.5)),between=list(x=1.8, y=0.25),par.strip.text=list(cex=2.1),layout=c(2,3),col.regions=colorRampPalette(c("white","blue","green","yellow","red")), main='', names.attr=nomesSubgraficos, colorkey=list(space="right",labels=list(cex=1.75))) + layer(sp.polygons(AmSulShape)) + layer(panel.xyplot(-55.993283,-34.270064,pch=17,col="red",cex=2),rows=c(2:3),columns=c(1)) + layer(panel.xyplot(-41.553056,-12.393333,pch=17,col="red",cex=2),rows=c(2:3),columns=c(2))
dev.off()


##codigo para graficos a partir dos resultados do Biomod##

spString = 'Melanosuchus.niger'

myOutputFolder = "/home/anderson/PosDoc/teste/biomod/myOutput"
filesPath = list.files(paste(myOutputFolder),full.names=TRUE,pattern=spString)
layers = stack(filesPath)


##Maximum training sensitivity plus specificity logistic threshold
threshold =  0.3824 #Melanosuchus niger 
threshold = 0.02443 # Caiman yacare
threshold = 0.1899 # Caiman latirostris
threshold = 0.3964 #Caiman crocodilus

##media presente (000) para MAXENT, RF e GLM
meanMX000 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_MAXENT.Phillips_000',sep=''))]]) #MX
meanRF000 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_RF_000',sep=''))]]) #RF
meanGLM000 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_GLM_000',sep=''))]]) #GLM

###media para as projecoes para o passado
##10kyrBP
meanMX010 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_MAXENT.Phillips_10kyrBP',sep=''))]]) #MX
meanRF010 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_RF_10kyrBP',sep=''))]]) #RF
meanGLM010 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_GLM_10kyrBP',sep=''))]]) #GLM

###media para as projecoes para o passado
##11kyrBP
meanMX011 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_MAXENT.Phillips_11kyrBP',sep=''))]]) #MX
meanRF011 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_RF_11kyrBP',sep=''))]]) #RF
meanGLM011 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_GLM_11kyrBP',sep=''))]]) #GLM

###media para as projecoes para o passado
##21kyrBP
meanMX021 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_MAXENT.Phillips_21kyrBP',sep=''))]]) #MX
meanRF021 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_RF_21kyrBP',sep=''))]]) #RF
meanGLM021 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_GLM_21kyrBP',sep=''))]]) #GLM

###media para as projecoes para o passado
##22kyrBP
meanMX022 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_MAXENT.Phillips_22kyrBP',sep=''))]]) #MX
meanRF022 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_RF_22kyrBP',sep=''))]]) #RF
meanGLM022 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_GLM_22kyrBP',sep=''))]]) #GLM

###media para as projecoes para o passado
##120kyrBP
meanMX120 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_MAXENT.Phillips_120kyrBP',sep=''))]]) #MX
meanRF120 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_RF_120kyrBP',sep=''))]]) #RF
meanGLM120 = mean(layers[[ c(paste(spString,'_AllData_RUN',1:3,'_GLM_120kyrBP',sep=''))]]) #GLM

##salvando
writeRaster(meanMX000,paste(myOutputFolder,'/',spString,'MX000.asc',sep=''),overwrite=TRUE)
writeRaster(meanRF000,paste(myOutputFolder,'/',spString,'RF000.asc',sep=''),overwrite=TRUE) 
writeRaster(meanGLM000,paste(myOutputFolder,'/',spString,'GLM000.asc',sep=''),overwrite=TRUE) 

writeRaster(meanMX010>threshold,paste(myOutputFolder,'/',spString,'MX010-BINARIO.asc',sep=''),overwrite=TRUE)
writeRaster(meanRF010>threshold,paste(myOutputFolder,'/',spString,'RF010-BINARIO.asc',sep=''),overwrite=TRUE) 
writeRaster(meanGLM010>threshold,paste(myOutputFolder,'/',spString,'GLM010-BINARIO.asc',sep=''),overwrite=TRUE) 

writeRaster(meanMX011>threshold,paste(myOutputFolder,'/',spString,'MX011-BINARIO.asc',sep=''),overwrite=TRUE)
writeRaster(meanRF011>threshold,paste(myOutputFolder,'/',spString,'RF011-BINARIO.asc',sep=''),overwrite=TRUE) 
writeRaster(meanGLM011>threshold,paste(myOutputFolder,'/',spString,'GLM011-BINARIO.asc',sep=''),overwrite=TRUE) 

writeRaster(meanMX021>threshold,paste(myOutputFolder,'/',spString,'MX021-BINARIO.asc',sep=''),overwrite=TRUE)
writeRaster(meanRF021>threshold,paste(myOutputFolder,'/',spString,'RF021-BINARIO.asc',sep=''),overwrite=TRUE) 
writeRaster(meanGLM021>threshold,paste(myOutputFolder,'/',spString,'GLM021-BINARIO.asc',sep=''),overwrite=TRUE) 

writeRaster(meanMX022>threshold,paste(myOutputFolder,'/',spString,'MX022-BINARIO.asc',sep=''),overwrite=TRUE)
writeRaster(meanRF022>threshold,paste(myOutputFolder,'/',spString,'RF022-BINARIO.asc',sep=''),overwrite=TRUE) 
writeRaster(meanGLM022>threshold,paste(myOutputFolder,'/',spString,'GLM022-BINARIO.asc',sep=''),overwrite=TRUE) 

writeRaster(meanMX120>threshold,paste(myOutputFolder,'/',spString,'MX120-BINARIO.asc',sep=''),overwrite=TRUE)
writeRaster(meanRF120>threshold,paste(myOutputFolder,'/',spString,'RF120-BINARIO.asc',sep=''),overwrite=TRUE) 
writeRaster(meanGLM120>threshold,paste(myOutputFolder,'/',spString,'GLM120-BINARIO.asc',sep=''),overwrite=TRUE) 



### CORRELACAO ENTRE OS DIRERENTES MODELOS ###

algoritmos = c('Maxent','GLM','Random Forest')
cor.list = list()

for (i in 1:16){ #OBS.: 16 e o numero de projecoes a serem comparadas (ver na pasta 'raster layers')
    layersToCompare=list()
    for (j in 1:length(algoritmos)){
        setwd(paste(projectFolder,algoritmos[j],'/Passado/Raster Layers',sep='')) #apenas passado
        filesPass = list.files(paste(getwd()),full.names=TRUE,pattern='.asc')
        filesPass=c(filesPass[1],filesPass[3],filesPass[5],filesPass[7],filesPass[9],filesPass[11],filesPass[13],filesPass[15],filesPass[17],filesPass[19],filesPass[21],filesPass[23],filesPass[25],filesPass[27],filesPass[29],filesPass[31]) #selecionando camadas e consertando a ordem das camadas
        layersToCompare = append(layersToCompare, filesPass[i])
    }
    stackToTest = stack(layersToCompare)
    valuesToTest = getValues(stackToTest)
    cor.matrix = as.data.frame(cor(valuesToTest, use="complete.obs"))
    cor.list[[i]] = data.frame(cor.matrix)
    write.csv(cor.list[[i]], file=paste("/home/anderson/PosDoc/teste/Correlacoes/","correlacao_",names(cor.matrix)[1],".csv",sep=""),row.names=T)
}


################################################################
#################### RANDOM FOREST #############################
################################################################

#registrando o tempo de processamento
ptm <- proc.time()

#abrindo um data.frame para armazenar os resultados de AUC
resultados.evaluacion.RF<-data.frame(Species=character(), auc=numeric(), stringsAsFactors=FALSE)
fila=0 #auxiliara na criacao do data.frame durante o loop
fossilPointsSuitability = NULL #objetvo em que serao gravadas as projcoes de suitability especificamente para cada ponto de registro fossil

for (i in 1:length(splist)){
    #especie = 1 #para fazer na mao bruta
    especie = i
     print(paste('Rodando Random Forest para a especie',splist[especie]))

    #presencas
    sp.file <- read.csv(paste(spOccFolder, splist[especie],".csv",sep=""),h=T) ### read sp occurrence
    sp.occ <- sp.file[,2:3]
    #coordinates(sp.occ)<- c("longitude","latitude")

    #extraindo dados da variavel climatica nos pontos de ocorrencia
    presclim <- extract(predictors, sp.occ, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)

    #criando um vetor de presenca para usar em uma coluna de presenca/ausencia na tabela final
    pres = rep(1, nrow(presclim))

    #juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia
    presclim = cbind(presclim,pres,sp.occ)
    presencias<-data.frame(presclim)
    presencias = presencias[complete.cases(presencias),]

    #criando ausencias para o background
    #rand = round(0.75*runif(nrow(presencias)))
    #presenciasTrain = presencias[rand==0,]
    #ausencias
    pseudoausencia1 <- randomPoints(mask=predictors[[1]], n=nrow(presencias), p=presencias[ , c("latitude","longitude")], excludep=TRUE) #este sera usado no loop para gerar ausencias de teste, la embaixo
    pseudoausencia2 <- round(pseudoausencia1[,1:2], digits=4)
    pseudoausencia3<-pseudoausencia2[!duplicated(pseudoausencia2),]
    pseudoausencia4<-pseudoausencia3[complete.cases(pseudoausencia3),]
    pseudoausencia<-data.frame(pseudoausencia4)
    colnames(pseudoausencia) <- c("longitude", "latitude")

    #extraindo dados da variavel climatica nos pontos de background
    ausclim <- extract(predictors, pseudoausencia, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)

    #criando um vetor de ausencias para usar em uma coluna de presenca/ausencia na tabela final
    pres = rep(0, nrow(ausclim))

    #juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia    
    ausclim = data.frame(ausclim, pres)
    ausencias <- cbind(ausclim,pseudoausencia)

    #ausenciasTrain = ausencias[rand==0,]
    
    #juntando presenca e ausencia
    #presausTrainRaw<-rbind(presenciasTrain, ausenciasTrain)
    #presausTrainRaw = data.frame(presausTrainRaw)
    #presausTrainRaw = presausTrainRaw[!duplicated(presausTrainRaw[,7:8]),] #selecionar colunas de longitude e latitude
    #presausTrainRaw<-presausTrainRaw[complete.cases(presausTrainRaw),]
    #presausTrain = presausTrainRaw
    
    #write.table(presausTrain, file=paste("/home/anderson/R/PosDoc/teste/Random Forest/",splist[i],'/',splist[i],"_presausTrain.csv", sep=""))

    ##CRIANDO E RODANDO O MODELO##    
    #model <- pres ~ bioclim_01+bioclim_04+bioclim_10+bioclim_11+bioclim_12+bioclim_15+bioclim_16+bioclim_17

    #RF <- randomForest(model, data=presausTrain, ntree=500)

    #avaliando o modelo
    V <- numeric()#abrir un vector vacío 
    
    for (j in 1:10){
        tryCatch({# bootstrapping with 10 replications

            #reparando uma porcao dos dados de presenca e ausencia (background) para calibrar (treinar) o modelo
            rand = round(0.75*runif(nrow(presencias)))
            presenciasTrain = presencias[rand==0,]
            ausenciasTrain = ausencias[rand==0,]

            #juntando presencas e ausencias da calibracao
            presausTrainRaw <- rbind(presenciasTrain, ausenciasTrain)
            presausTrainRaw = data.frame(presausTrainRaw)
            presausTrainRaw = presausTrainRaw[!duplicated(presausTrainRaw[,7:8]),] #selecionar colunas de longitude e latitude
            presausTrainRaw<-presausTrainRaw[complete.cases(presausTrainRaw),]
            presausTrain = presausTrainRaw

            ##CRIANDO E RODANDO O MODELO##    
            ##model <- pres ~ bioclim_10+bioclim_11+bioclim_16+bioclim_17
            model <- pres ~ bioclim_10 + I(bioclim_10^2) + bioclim_11 + I(bioclim_11^2) + bioclim_16 + I(bioclim_16^2) + bioclim_17 + I(bioclim_17^2)
            RF <- randomForest(model, data=presausTrain, ntree=500)

            #porcentajepres = round(0.25*nrow(presencias)) #seleccionar un porcentajes de filas de un data.frame
            #presencias.evaluacion<-presencias[sample(nrow(presencias), porcentajepres), ] #seleccionar ese porcentaje de filas aleatorias.

            #pegando a porcao dos dados separados para a avaliacao (validacao) do modelo
            presencias.evaluacion = presencias[rand==1,]
            presencias.evaluacion <- cbind(presencias.evaluacion$longitude,presencias.evaluacion$latitude)
            #porcentajeabs = round(0.25*nrow(pseudoausencia)) # choose pseudoabsences points evaluation
            #pseudoausencias.evaluacion <- pseudoausencia[sample(nrow(pseudoausencia), porcentajeabs), ]
            pseudoausencias.evaluacion = ausencias[rand==1,]
            pseudoausencias.evaluacion = cbind(pseudoausencias.evaluacion$longitude,pseudoausencias.evaluacion$latitude)

            ##RODANDO A AVALIACAO DO MODELO##
            evaluacion=evaluate(presencias.evaluacion, pseudoausencias.evaluacion, RF, predictors)

            #registrando o valor de AUC em um objeto
            V[j]<-evaluacion@"auc" #sacamos el valor de auc (fíjate que es una @ en lugar de $ para mirar dentro de los slots)y guardamos en vector
        }, error=function(e){cat("ERROR :",conditionMessage(e),"\n")}
   )}  

    #registrando a media de AUC
    auc<-mean(V, na.rm=T)#media de los vectores de las iternaciones de j

    #registrando o valor de AUC medio em uma tabela
    fila=fila+1
    resultados.evaluacion.RF[fila, "Species"]<-splist[especie]
    resultados.evaluacion.RF[fila, "auc"]<-evaluacion@auc
    
    #gravando um PDF com a AUC do modelo
    pdf(file=paste(projectFolder,"Random Forest/",splist[especie],"/",splist[especie],'_ROC',".pdf",sep=""))
    plot(evaluacion, "ROC", cex=0.3)
    dev.off()

    ##PROJETANDO o nicho no espaco atraves do modelo ajustado##
    projecaoSuitability <- predict(predictors, RF)

    #gravando um raster com o mapa de projecao gerado pelo modelo
    writeRaster(projecaoSuitability,filename=paste(projectFolder,"Random Forest/",splist[especie],"/",splist[especie],".asc", sep=""),overwrite=T)

    #gravando um PDF com o mapa gerado pelo modelo
    pdf(file=paste(projectFolder,"Random Forest/",splist[i],"/",splist[i],".pdf",sep=""))
    plot(projecaoSuitability, main=paste(splist[i]))
    plot(AmSulShape,add=T)
    points(presencias$longitude,presencias$latitude,col='orange',pch=20,cex=0.75)
    points(presencias$longitude,presencias$latitude,col='red',cex=0.7)
    dev.off()

    #criando um mapa binario
    threshold <- threshold(evaluacion,'spec_sens')
    bin <- projecaoSuitability > threshold #apply threshold to transform logistic output into binary maps

    #salvando um raster com o mapa binario
    writeRaster(bin,filename=paste(projectFolder,"Random Forest/",splist[i],"/",splist[especie],"BINARIO.asc",sep=""),overwrite=T)

    #gravando um PDF com o mapa gerado pelo modelo
    pdf(file=paste(projectFolder,"Random Forest/",splist[especie],"/",splist[i],"-BINARIO.pdf",sep=""))
    plot(bin, main= paste(splist[especie]))
    plot(AmSulShape,add=T)
    points(presencias$longitude,presencias$latitude,col='orange',pch=20,cex=0.75)
    points(presencias$longitude,presencias$latitude,col='red',cex=0.7)
    dev.off()

    ###PROJECAO PARA O PASSADO###

    #abrindo os dados de registros fosseis  para uma especie
    sp.fossil.data = occ.sps.fosseis[occ.sps.fosseis$species==splist[especie],] #ATENCAO: este script nao funciona se houver mais de um registro fossil por camada de tempo usada para projecao

    for(l in 1:nrow(sp.fossil.data)){
        #loop para cada registro fossil de uma especie

        #definindoo fossil
        sp.fossil = sp.fossil.data[l,]

        #abrindo as variaveis ambientais do tempo do fossil
        filesProjectionRaw <- stack(list.files(path = paste(envVarFolder,"dados_projeto/0",sp.fossil$kyr,sep=""), pattern='asc', full.names=T)) ###abrindo camandas para projecao (passado, futuro, outro local, etc)
        filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul
        files.crop.sub.projection <- dropLayer(filesProjection, c(1,2,5,6)) #removendo as camadas que mostraram correlacao
        predictorsProjection = files.crop.sub.projection #preditoras para o tempo do fossil

        ##PROJETANDO o nicho no espaco atraves do modelo ajustado##
        projecaoSuitabilityPassado <- predict(predictorsProjection, RF) #PASSADO

        #criando um objeto com as coordenadas do registro fossil
        fossilPoints = sp.fossil
        fossilPoints = cbind(fossilPoints$longitude, fossilPoints$latitude)

        #obtendo a projecao de qualidade de habitat especificamente para o ponto do fossil
        fossilPointsVars = extract(predictorsProjection,fossilPoints)
        fossilPoints.RF = predict(RF, fossilPointsVars)
        fossilPointsSuitability = rbind(fossilPointsSuitability,data.frame(splist[especie],sp.fossil$kyr,fossilPoints.RF))

        #salvando um raster com a projecao do modelo para o tempo do fossil
        writeRaster(projecaoSuitabilityPassado,filename=paste(projectFolder,"Random Forest/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr," K years BP.asc", sep=""),overwrite=T)

        #salvando um PDF com a projecao do modelo para o tempo do fossil
        pdf(file=paste(projectFolder,"Random Forest/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr," K years BP.pdf", sep=""))
        plot(projecaoSuitabilityPassado, main=paste(splist[especie],'-',sp.fossil$kyr,'K years BP',sep=''))
        plot(AmSulShape,add=T)
        #plotando as coordenadas do refistro fossil
        points(fossilPoints,col='orange',pch=20,cex=0.75)
        points(fossilPoints,col='red',cex=0.8)
        dev.off()

        #criando um mapa binario para a projecao do modelo (empregando o threshold que ja foi criado apos a avaliacao do modelo)
        bin <- projecaoSuitabilityPassado > threshold#apply threshold to transform logistic output into binary maps

        #salvando um raster com a projecao do modelo para o tempo do fossil
        writeRaster(bin,filename=paste(projectFolder,"Random Forest/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr,"K years BP - BINARIO.asc",sep=""),overwrite=T)

        #salvando um PDF com a projecao do modelo para o tempo do fossil
        pdf(file=paste(projectFolder,"Random Forest/Passado/",splist[especie],"/",splist[especie],"-",sp.fossil$kyr," K years BP - BINARIO.pdf",sep="")) 
        plot(bin, main= paste(splist[especie],"-",sp.fossil$kyr," K years BP.pdf",sep=""))
        plot(AmSulShape,add=T)
        #plotando as coordenadas do refistro fossil
        points(fossilPoints,col='orange',pch=20,cex=0.75)
        points(fossilPoints,col='red',cex=0.8)
        dev.off()

    }
}

#salvando a tabela de dados da avaliacao dos modelos
write.table(resultados.evaluacion.RF,file=paste(projectFolder,"Random Forest/","AUCmodelos.csv",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")

write.table(fossilPointsSuitability,file=paste(projectFolder,"Random Forest/","suitabilityNoPontoFossil.csv",sep=""))

#fechando e informando o tempo de processamento
msgm= proc.time() - ptm
print(paste('Tempo gasto para rodar RANDOM FOREST: ', msgm[3]/60,' minutos',sep=''))

################################################################
########################## GLM #################################
################################################################

#registrando o tempo de processamento
ptm = proc.time()

#abrindo um data.frame para armazenar os resultados de AUC
resultados.evaluacion.GLM<-data.frame(Species=character(), auc=numeric(), stringsAsFactors=FALSE)
fila=0 #auxiliara na criacao do data.frame durante o loop
fossilPointsSuitability = NULL #objetvo em que serao gravadas as projcoes de suitability especificamente para cada ponto de registro fossil

for (i in 1:length(splist)){
    #especie = 1 #para fazer na mao bruta
    especie = i
    print(paste('Rodando GLM para a especie',splist[especie]))
    
    #presencas
    sp.file <- read.csv(paste(spOccFolder, splist[especie],".csv",sep=""),h=T) ### read sp occurrence
    sp.occ <- sp.file[,2:3]
    #coordinates(sp.occ)<- c("longitude","latitude")

    #extraindo dados da variavel climatica nos pontos de ocorrencia
    presclim <- extract(predictors, sp.occ, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)

    #criando um vetor de presenca para usar em uma coluna de presenca/ausencia na tabela final
    pres = rep(1, nrow(presclim))

    #juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia
    presclim = cbind(presclim,pres,sp.occ)
    presencias<-data.frame(presclim)
    presencias = presencias[complete.cases(presencias),]

    #criando ausencias para o background
    #rand = round(0.75*runif(nrow(presencias)))
    #presenciasTrain = presencias[rand==0,]
    #ausencias
    pseudoausencia1 <- randomPoints(mask=predictors[[1]], n=nrow(presencias), p=presencias[ , c("latitude","longitude")], excludep=TRUE) #este sera usado no loop para gerar ausencias de teste, la embaixo
    pseudoausencia2 <- round(pseudoausencia1[,1:2], digits=4)
    pseudoausencia3<-pseudoausencia2[!duplicated(pseudoausencia2),]
    pseudoausencia4<-pseudoausencia3[complete.cases(pseudoausencia3),]
    pseudoausencia<-data.frame(pseudoausencia4)
    colnames(pseudoausencia) <- c("longitude", "latitude")

    #extraindo dados da variavel climatica nos pontos de background
    ausclim <- extract(predictors, pseudoausencia, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)

    #criando um vetor de ausencias para usar em uma coluna de presenca/ausencia na tabela final
    pres = rep(0, nrow(ausclim))

    #juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia    
    ausclim = data.frame(ausclim, pres)
    ausencias <- cbind(ausclim,pseudoausencia)

    #ausenciasTrain = ausencias[rand==0,]
    
    #juntando presenca e ausencia
    #presausTrainRaw<-rbind(presenciasTrain, ausenciasTrain)
    #presausTrainRaw = data.frame(presausTrainRaw)
    #presausTrainRaw = presausTrainRaw[!duplicated(presausTrainRaw[,7:8]),] #selecionar colunas de longitude e latitude
    #presausTrainRaw<-presausTrainRaw[complete.cases(presausTrainRaw),]
    #presausTrain = presausTrainRaw
    
    #write.table(presausTrain, file=paste("/home/anderson/R/PosDoc/teste/GLM/",splist[i],'/',splist[i],"_presausTrain.csv", sep=""))

    ##CRIANDO E RODANDO O MODELO##    
    #model <- pres ~ bioclim_01+bioclim_04+bioclim_10+bioclim_11+bioclim_12+bioclim_15+bioclim_16+bioclim_17

    #GLM <- randomForest(model, data=presausTrain, ntree=500)

    #avaliando o modelo
    V <- numeric()#abrir un vector vacío 
    
    for (j in 1:10){
        tryCatch({# bootstrapping with 10 replications

            #reparando uma porcao dos dados de presenca e ausencia (background) para calibrar (treinar) o modelo
            rand = round(0.75*runif(nrow(presencias)))
            presenciasTrain = presencias[rand==0,]
            ausenciasTrain = ausencias[rand==0,]

            #juntando presencas e ausencias da calibracao
            presausTrainRaw <- rbind(presenciasTrain, ausenciasTrain)
            presausTrainRaw = data.frame(presausTrainRaw)
            presausTrainRaw = presausTrainRaw[!duplicated(presausTrainRaw[,7:8]),] #selecionar colunas de longitude e latitude
            presausTrainRaw<-presausTrainRaw[complete.cases(presausTrainRaw),]
            presausTrain = presausTrainRaw

            ##CRIANDO E RODANDO O MODELO##                
            #model <- pres ~ bioclim_10 + I(bioclim_10^2) + bioclim_11 + I(bioclim_11^2) + bioclim_16 + I(bioclim_16^2) + bioclim_17 + I(bioclim_17^2)
            model <- pres ~ bioclim_11 + bioclim_17#melhor modelo
            GLM <- glm(model, family=binomial(link=logit), data=presausTrainRaw)
            
            #porcentajepres = round(0.25*nrow(presencias)) #seleccionar un porcentajes de filas de un data.frame
            #presencias.evaluacion<-presencias[sample(nrow(presencias), porcentajepres), ] #seleccionar ese porcentaje de filas aleatorias.

            #pegando a porcao dos dados separados para a avaliacao (validacao) do modelo
            presencias.evaluacion = presencias[rand==1,]
            presencias.evaluacion <- cbind(presencias.evaluacion$longitude,presencias.evaluacion$latitude)
            #porcentajeabs = round(0.25*nrow(pseudoausencia)) # choose pseudoabsences points evaluation
            #pseudoausencias.evaluacion <- pseudoausencia[sample(nrow(pseudoausencia), porcentajeabs), ]
            pseudoausencias.evaluacion = ausencias[rand==1,]
            pseudoausencias.evaluacion = cbind(pseudoausencias.evaluacion$longitude,pseudoausencias.evaluacion$latitude)

            ##RODANDO A AVALIACAO DO MODELO##
            evaluacion=evaluate(presencias.evaluacion, pseudoausencias.evaluacion, GLM, predictors,type='response')

            #registrando o valor de AUC em um objeto
            V[j]<-evaluacion@"auc" #sacamos el valor de auc (fíjate que es una @ en lugar de $ para mirar dentro de los slots)y guardamos en vector
        }, error=function(e){cat("ERROR :",conditionMessage(e),"\n")}
   )}  

    #registrando a media de AUC
    auc<-mean(V, na.rm=T)#media de los vectores de las iternaciones de j

    #registrando o valor de AUC medio em uma tabela
    fila=fila+1
    resultados.evaluacion.GLM[fila, "Species"]<-splist[i]
    resultados.evaluacion.GLM[fila, "auc"]<-evaluacion@auc
    
    #gravando um PDF com a AUC do modelo
    pdf(file=paste(projectFolder,"GLM/",splist[i],"/",splist[i],'_ROC',".pdf",sep=""))
    plot(evaluacion, "ROC", cex=0.3)
    dev.off()

    ##PROJETANDO o nicho no espaco atraves do modelo ajustado##
    projecaoSuitability <- predict(predictors, GLM, type='response')

    #gravando um raster com o mapa de projecao gerado pelo modelo
    writeRaster(projecaoSuitability,filename=paste(projectFolder,"GLM/",splist[i],"/",splist[i],".asc", sep=""),overwrite=TRUE)

    #gravando um PDF com o mapa gerado pelo modelo
    pdf(file=paste(projectFolder,"GLM/",splist[i],"/",splist[i],".pdf",sep=""))
    plot(projecaoSuitability, main=paste(splist[i]))
    plot(AmSulShape,add=TRUE)
    points(presencias$longitude,presencias$latitude,col='orange',pch=20,cex=0.75)
    points(presencias$longitude,presencias$latitude,col='red',cex=0.7)
    dev.off()

    #criando um mapa binario
    threshold = threshold(evaluacion,'spec_sens')
    bin <- projecaoSuitability > threshold #apply threshold to transform logistic output into binary maps

    #salvando um raster com o mapa binario
    writeRaster(bin,filename=paste(projectFolder,"GLM/",splist[i],"/",splist[i],"BINARIO.asc",sep=""),overwrite=T)

    #gravando um PDF com o mapa gerado pelo modelo
    pdf(file=paste(projectFolder,"GLM/",splist[i],"/",splist[i],"-BINARIO.pdf",sep=""))
    plot(bin, main= paste(splist[i]))
    plot(AmSulShape,add=T)
    points(presencias$longitude,presencias$latitude,col='orange',pch=20,cex=0.75)
    points(presencias$longitude,presencias$latitude,col='red',cex=0.7)
    dev.off()

    ###PROJECAO PARA O PASSADO###
        
    #abrindo os dados de registros fosseis  para uma especie
    sp.fossil.data = occ.sps.fosseis[occ.sps.fosseis$species==splist[especie],] #ATENCAO: este script nao funciona se houver mais de um registro fossil por camada de tempo usada para projecao
        
    for(l in 1:nrow(sp.fossil.data)){
        #loop para cada registro fossil de uma especie

        #definindoo fossil
        sp.fossil = sp.fossil.data[l,]
    
        #abrindo as variaveis ambientais do tempo do fossil
        filesProjectionRaw <- stack(list.files(path = paste(envVarFolder,"dados_projeto/0",sp.fossil$kyr,sep=""), pattern='asc', full.names=T)) ###abrindo camandas para projecao (passado, futuro, outro local, etc)
        filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul
        files.crop.sub.projection <- dropLayer(filesProjection, c(1,2,5,6)) #removendo as camadas que mostraram correlacao
        predictorsProjection = files.crop.sub.projection #preditoras para o tempo do fossil

        ##PROJETANDO o nicho no espaco atraves do modelo ajustado##
        projecaoSuitabilityPassado <- predict(predictorsProjection,GLM,type='response') #PASSADO

        #criando um objeto com as coordenadas do registro fossil
        fossilPoints = sp.fossil
        fossilPoints = cbind(fossilPoints$longitude,fossilPoints$latitude)

        #obtendo a projecao de qualidade de habitat especificamente para o ponto do fossil
        fossilPointsVars = extract(predictorsProjection,fossilPoints)
        fossilPoints.GLM = predict(GLM, as.data.frame(fossilPointsVars),type='response')
        fossilPointsSuitability = rbind(fossilPointsSuitability,data.frame(splist[especie],sp.fossil$kyr,fossilPoints.GLM))

        
        #salvando um raster com a projecao do modelo para o tempo do fossil
        writeRaster(projecaoSuitabilityPassado,filename=paste(projectFolder,"GLM/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr," K years BP.asc", sep=""),overwrite=T)
            
        #salvando um PDF com a projecao do modelo para o tempo do fossil
        pdf(file=paste(projectFolder,"GLM/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr," K years BP.pdf", sep=""))
        plot(projecaoSuitabilityPassado, main=paste(splist[especie],'-',sp.fossil$kyr,'K years BP',sep=''))
        plot(AmSulShape,add=T)
        #plotando as coordenadas do refistro fossil
        points(fossilPoints,col='orange',pch=20,cex=0.75)
        points(fossilPoints,col='red',cex=0.8)
        dev.off()

        #criando um mapa binario para a projecao do modelo (empregando o threshold que ja foi criado apos a avaliacao do modelo)
        bin <- projecaoSuitabilityPassado > threshold#apply threshold to transform logistic output into binary maps

        #salvando um raster com a projecao do modelo para o tempo do fossil
        writeRaster(bin,filename=paste(projectFolder,"GLM/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr,"K years BP - BINARIO.asc",sep=""),overwrite=T)

        #salvando um PDF com a projecao do modelo para o tempo do fossil
        pdf(file=paste(projectFolder,"GLM/Passado/",splist[especie],"/",splist[especie],"-",sp.fossil$kyr," K years BP - BINARIO.pdf",sep="")) 
        plot(bin, main= paste(splist[especie],"-",sp.fossil$kyr," K years BP.pdf",sep=""))
        plot(AmSulShape,add=T)
        #plotando as coordenadas do refistro fossil
        points(fossilPoints,col='orange',pch=20,cex=0.75)
        points(fossilPoints,col='red',cex=0.8)
        dev.off()
    }
}

#salvando a tabela de dados da avaliacao dos modelos 
write.table(resultados.evaluacion.GLM,file=paste(projectFolder,"GLM/","AUCmodelos.csv",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")

write.csv(fossilPointsSuitability,file=paste(projectFolder,"GLM/","suitabilityNoPontoFossil.csv",sep=""),row.names=F)

#registrando e informando o tempo de processamento
msgm= proc.time() - ptm
print(paste('Tempo gasto para rodar GLM: ', msgm[3]/60,' minutos',sep=''))


################################################################
########################## BIOCLIM #############################
################################################################

#registrando o tempo de processamento
ptm = proc.time()

#abrindo um data.frame para armazenar os resultados de AUC
resultados.evaluacion.BIOC<-data.frame(Species=character(), auc=numeric(), stringsAsFactors=FALSE)
fila=0 #auxiliara na criacao do data.frame durante o loop
fossilPointsSuitability = NULL #objetvo em que serao gravadas as projcoes de suitability especificamente para cada ponto de registro fossil


for (i in 1:length(splist)){
    #especie = 1 #para fazer na mao bruta
    especie = i
    print(paste('Rodando BIOCLIM para a especie',splist[especie]))


    #presencas
    sp.file <- read.csv(paste(spOccFolder, splist[especie],".csv",sep=""),h=T) ### read sp occurrence
    sp.occ <- sp.file[,2:3]
    #coordinates(sp.occ)<- c("longitude","latitude")

    #extraindo dados da variavel climatica nos pontos de ocorrencia
    presclim <- extract(predictors, sp.occ, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)

    #criando um vetor de presenca para usar em uma coluna de presenca/ausencia na tabela final
    pres = rep(1, nrow(presclim))

    #juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia
    presclim = cbind(presclim,pres,sp.occ)
    presencias<-data.frame(presclim)
    presencias = presencias[complete.cases(presencias),]

    #criando ausencias para o background
    #rand = round(0.75*runif(nrow(presencias)))
    #presenciasTrain = presencias[rand==0,]
    #ausencias
    pseudoausencia1 <- randomPoints(mask=predictors[[1]], n=nrow(presencias), p=presencias[ , c("latitude","longitude")], excludep=TRUE) #este sera usado no loop para gerar ausencias de teste, la embaixo
    pseudoausencia2 <- round(pseudoausencia1[,1:2], digits=4)
    pseudoausencia3<-pseudoausencia2[!duplicated(pseudoausencia2),]
    pseudoausencia4<-pseudoausencia3[complete.cases(pseudoausencia3),]
    pseudoausencia<-data.frame(pseudoausencia4)
    colnames(pseudoausencia) <- c("longitude", "latitude")

    #extraindo dados da variavel climatica nos pontos de background
    ausclim <- extract(predictors, pseudoausencia, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)

    #criando um vetor de ausencias para usar em uma coluna de presenca/ausencia na tabela final
    pres = rep(0, nrow(ausclim))

    #juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia    
    ausclim = data.frame(ausclim, pres)
    ausencias <- cbind(ausclim,pseudoausencia)

    #ausenciasTrain = ausencias[rand==0,]
    
    #juntando presenca e ausencia
    #presausTrainRaw<-rbind(presenciasTrain, ausenciasTrain)
    #presausTrainRaw = data.frame(presausTrainRaw)
    #presausTrainRaw = presausTrainRaw[!duplicated(presausTrainRaw[,7:8]),] #selecionar colunas de longitude e latitude
    #presausTrainRaw<-presausTrainRaw[complete.cases(presausTrainRaw),]
    #presausTrain = presausTrainRaw
    
    #write.table(presausTrain, file=paste("/home/anderson/R/PosDoc/teste/BIOC/",splist[i],'/',splist[i],"_presausTrain.csv", sep=""))

    ##CRIANDO E RODANDO O MODELO##    
    #model <- pres ~ bioclim_01+bioclim_04+bioclim_10+bioclim_11+bioclim_12+bioclim_15+bioclim_16+bioclim_17

    #BIOC <- randomForest(model, data=presausTrain, ntree=500)

    #avaliando o modelo
    V <- numeric()#abrir un vector vacío 
    
    for (j in 1:100){
        tryCatch({# bootstrapping with 10 replications

            #reparando uma porcao dos dados de presenca e ausencia (background) para calibrar (treinar) o modelo
            rand = round(0.75*runif(nrow(presencias)))
            presenciasTrain = presencias[rand==0,]
            ausenciasTrain = ausencias[rand==0,]

            #juntando presencas e ausencias da calibracao
            presausTrainRaw <- rbind(presenciasTrain, ausenciasTrain)
            presausTrainRaw = data.frame(presausTrainRaw)
            presausTrainRaw = presausTrainRaw[!duplicated(presausTrainRaw[,7:8]),] #selecionar colunas de longitude e latitude
            presausTrainRaw<-presausTrainRaw[complete.cases(presausTrainRaw),]
            presausTrain = presausTrainRaw

            ##CRIANDO E RODANDO O MODELO##    
            BIOC <- bioclim(predictors, cbind(presenciasTrain$longitude, presenciasTrain$latitude))

            #porcentajepres = round(0.25*nrow(presencias)) #seleccionar un porcentajes de filas de un data.frame
            #presencias.evaluacion<-presencias[sample(nrow(presencias), porcentajepres), ] #seleccionar ese porcentaje de filas aleatorias.

            #pegando a porcao dos dados separados para a avaliacao (validacao) do modelo
            presencias.evaluacion = presencias[rand==1,]
            presencias.evaluacion <- cbind(presencias.evaluacion$longitude,presencias.evaluacion$latitude)
            #porcentajeabs = round(0.25*nrow(pseudoausencia)) # choose pseudoabsences points evaluation
            #pseudoausencias.evaluacion <- pseudoausencia[sample(nrow(pseudoausencia), porcentajeabs), ]
            pseudoausencias.evaluacion = ausencias[rand==1,]
            pseudoausencias.evaluacion = cbind(pseudoausencias.evaluacion$longitude,pseudoausencias.evaluacion$latitude)

            ##RODANDO A AVALIACAO DO MODELO##
            evaluacion=evaluate(presencias.evaluacion, pseudoausencias.evaluacion, BIOC, predictors)
            
            #registrando o valor de AUC em um objeto
            V[j]<-evaluacion@"auc" #sacamos el valor de auc (fíjate que es una @ en lugar de $ para mirar dentro de los slots)y guardamos en vector
        }, error=function(e){cat("ERROR :",conditionMessage(e),"\n")}
   )}  

    #registrando a media de AUC
    auc<-mean(V, na.rm=T)#media de los vectores de las iternaciones de j

    #registrando o valor de AUC medio em uma tabela
    fila=fila+1
    resultados.evaluacion.BIOC[fila, "Species"]<-splist[i]
    resultados.evaluacion.BIOC[fila, "auc"]<-evaluacion@auc
    
    #gravando um PDF com a AUC do modelo
    pdf(file=paste(projectFolder,"BIOC/",splist[i],"/",splist[i],'_ROC',".pdf",sep=""))
    plot(evaluacion, "ROC", cex=0.3)
    dev.off()

    ##PROJETANDO o nicho no espaco atraves do modelo ajustado##
    projecaoSuitability <- predict(predictors, BIOC)

    #gravando um raster com o mapa de projecao gerado pelo modelo
    writeRaster(projecaoSuitability,filename=paste(projectFolder,"BIOC/",splist[i],"/",splist[i],".grd", sep=""),overwrite=T)

    #gravando um PDF com o mapa gerado pelo modelo
    pdf(file=paste(projectFolder,"BIOC/",splist[i],"/",splist[i],".pdf",sep=""))
    plot(projecaoSuitability, main=paste(splist[i]))
    plot(AmSulShape,add=T)
    points(presencias$longitude,presencias$latitude,col='orange',pch=20,cex=0.75)
    points(presencias$longitude,presencias$latitude,col='red',cex=0.7)
    dev.off()

    #criando um mapa binario
    threshold <- evaluacion@t[which.max(evaluacion@TPR + evaluacion@TNR)]
    bin <- projecaoSuitability > threshold #apply threshold to transform logistic output into binary maps

    #salvando um raster com o mapa binario
    writeRaster(bin,filename=paste(projectFolder,"BIOC/",splist[i],"/",splist[i],"BINARIO.asc",sep=""),overwrite=T)

    #gravando um PDF com o mapa gerado pelo modelo
    pdf(file=paste(projectFolder,"BIOC/",splist[i],"/",splist[i],"-BINARIO.pdf",sep=""))
    plot(bin, main= paste(splist[i]))
    plot(AmSulShape,add=T)
    points(presencias$longitude,presencias$latitude,col='orange',pch=20,cex=0.75)
    points(presencias$longitude,presencias$latitude,col='red',cex=0.7)
    dev.off()

    ###PROJECAO PARA O PASSADO###
        
    #abrindo os dados de registros fosseis  para uma especie
    sp.fossil.data = occ.sps.fosseis[occ.sps.fosseis$species==splist[especie],] #ATENCAO: este script nao funciona se houver mais de um registro fossil por camada de tempo usada para projecao
        
    for(l in 1:nrow(sp.fossil.data)){
        #loop para cada registro fossil de uma especie

        #definindoo fossil
        sp.fossil = sp.fossil.data[l,]
    
        #abrindo as variaveis ambientais do tempo do fossil
        filesProjectionRaw <- stack(list.files(path = paste(envVarFolder,"dados_projeto/0",sp.fossil$kyr,sep=""), pattern='asc', full.names=T)) ###abrindo camandas para projecao (passado, futuro, outro local, etc)
        filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul
        files.crop.sub.projection <- dropLayer(filesProjection, c(1,2,5,6)) #removendo as camadas que mostraram correlacao
        predictorsProjection = files.crop.sub.projection #preditoras para o tempo do fossil

        ##PROJETANDO o nicho no espaco atraves do modelo ajustado##
        projecaoSuitabilityPassado <- predict(predictorsProjection, BIOC) #PASSADO

        #criando um objeto com as coordenadas do registro fossil
        fossilPoints = sp.fossil
        fossilPoints = cbind(fossilPoints$longitude, fossilPoints$latitude)

        #obtendo a projecao de qualidade de habitat especificamente para o ponto do fossil
        fossilPointsVars = extract(predictorsProjection,fossilPoints)
        fossilPoints.BIOC = predict(BIOC, fossilPointsVars)
        fossilPointsSuitability = rbind(fossilPointsSuitability,data.frame(splist[especie],sp.fossil$kyr,fossilPoints.BIOC))

        
        #salvando um raster com a projecao do modelo para o tempo do fossil
        writeRaster(projecaoSuitabilityPassado,filename=paste(projectFolder,"BIOC/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr," K years BP.asc", sep=""),overwrite=T)
            
        #salvando um PDF com a projecao do modelo para o tempo do fossil
        pdf(file=paste(projectFolder,"BIOC/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr," K years BP.pdf", sep=""))
        plot(projecaoSuitabilityPassado, main=paste(splist[especie],'-',sp.fossil$kyr,'K years BP',sep=''))
        plot(AmSulShape,add=T)
        #plotando as coordenadas do refistro fossil
        points(fossilPoints,col='orange',pch=20,cex=0.75)
        points(fossilPoints,col='red',cex=0.8)
        dev.off()

        #criando um mapa binario para a projecao do modelo (empregando o threshold que ja foi criado apos a avaliacao do modelo)
        bin <- projecaoSuitabilityPassado > threshold#apply threshold to transform logistic output into binary maps

        #salvando um raster com a projecao do modelo para o tempo do fossil
        writeRaster(bin,filename=paste(projectFolder,"BIOC/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr,"K years BP - BINARIO.asc",sep=""),overwrite=T)

        #salvando um PDF com a projecao do modelo para o tempo do fossil
        pdf(file=paste(projectFolder,"BIOC/Passado/",splist[especie],"/",splist[especie],"-",sp.fossil$kyr," K years BP - BINARIO.pdf",sep="")) 
        plot(bin, main= paste(splist[especie],"-",sp.fossil$kyr," K years BP.pdf",sep=""))
        plot(AmSulShape,add=T)
        #plotando as coordenadas do refistro fossil
        points(fossilPoints,col='orange',pch=20,cex=0.75)
        points(fossilPoints,col='red',cex=0.8)
        dev.off()
            
    }
}

#salvando a tabela de dados da avaliacao dos modelos 
write.table(resultados.evaluacion.BIOC,file=paste(projectFolder,"BIOC/","AUCmodelos.csv",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")

write.csv(fossilPointsSuitability,file=paste(projectFolder,"BIOC/","suitabilityNoPontoFossil.csv",sep=""),row.names=F)

#registrando e informando o tempo de processamento
msgm= proc.time() - ptm
print(paste('Tempo gasto para rodar BIOCLIM: ', msgm[3]/60,' minutos',sep=''))


################################################################
########################### MAXENT #############################
################################################################


#registrando o tempo de processamento
ptm = proc.time()

#abrindo um data.frame para armazenar os resultados de AUC
resultados.evaluacion.MX<-data.frame(Species=character(), auc=numeric(), stringsAsFactors=FALSE)
fila=0 #auxiliara na criacao do data.frame durante o loop
fossilPointsSuitability = NULL #objeto em que serao gravadas as projcoes de suitability especificamente para cada ponto de registro fossil

for (i in 1:length(splist)){
    #especie = 1 #para fazer na mao bruta
    especie = i
    print(paste('Rodando Maxent para a especie',splist[especie]))

    #presencas
    sp.file <- read.csv(paste(spOccFolder, splist[especie],".csv",sep=""),h=T) ### read sp occurrence
    sp.occ <- sp.file[,2:3]
    #coordinates(sp.occ)<- c("longitude","latitude")

    #extraindo dados da variavel climatica nos pontos de ocorrencia
    presclim <- extract(predictors, sp.occ, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)

    #criando um vetor de presenca para usar em uma coluna de presenca/ausencia na tabela final
    #pres = rep(1, nrow(presclim))
    ## pres = rep(1, nrow(sp.occ))

    #juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia
    presclim = cbind(presclim,sp.occ)
    presencias<-data.frame(presclim)
    presencias = presencias[complete.cases(presencias),]

    #criando ausencias para o background
    #rand = round(0.75*runif(nrow(presencias)))
    #presenciasTrain = presencias[rand==0,]

    #ausencias
    pseudoausencia1 <- randomPoints(mask=predictors[[1]], n=nrow(presencias), p=presencias, excludep=TRUE) #este sera usado no loop para gerar ausencias de teste, la embaixo
    pseudoausencia2 <- round(pseudoausencia1[,1:2], digits=4)
    pseudoausencia3 <- pseudoausencia2[!duplicated(pseudoausencia2),]
    pseudoausencia4 <- pseudoausencia3[complete.cases(pseudoausencia3),]
    pseudoausencia<-data.frame(pseudoausencia4)
    colnames(pseudoausencia) <- c("longitude", "latitude")

    #extraindo dados da variavel climatica nos pontos de background
    ausclim <- extract(predictors, pseudoausencia, method='bilinear', buffer=NULL, fun=NULL, df=TRUE)

    #criando um vetor de ausencias para usar em uma coluna de presenca/ausencia na tabela final
    ## pres = rep(0, nrow(pseudoausencia)) 

    #juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia    
    #ausclim = data.frame(ausclim, pres)
    ausencias <- cbind(ausclim,pseudoausencia)
        
    #juntando presenca e ausencia
    #presausTrainRaw<-rbind(presenciasTrain, ausenciasTrain)
    #presausTrainRaw = data.frame(presausTrainRaw)
    #presausTrainRaw = presausTrainRaw[!duplicated(presausTrainRaw[,7:8]),] #selecionar colunas de longitude e latitude
    #presausTrainRaw<-presausTrainRaw[complete.cases(presausTrainRaw),]
    #presausTrain = presausTrainRaw
    
    #write.table(presausTrain, file=paste("/home/anderson/R/PosDoc/teste/Maxent/",splist[i],'/',splist[i],"_presausTrain.csv", sep=""))

    ##CRIANDO E RODANDO O MODELO##    
    #model <- pres ~ bioclim_01+bioclim_04+bioclim_10+bioclim_11+bioclim_12+bioclim_15+bioclim_16+bioclim_17

    #MX <- randomForest(model, data=presausTrain, ntree=500)

    #avaliando o modelo
    V <- numeric()#abrir un vector vacío 
    
    for (j in 1:3){
        tryCatch({# bootstrapping with 10 replications

            #reparando uma porcao dos dados de presenca e ausencia (background) para calibrar (treinar) o modelo
            rand = round(0.75*runif(nrow(presencias)))
            presenciasTrain = presencias[rand==0,]
            presenciasTrain = cbind(presenciasTrain$longitude, presenciasTrain$latitude)
            ## ausenciasTrain = ausencias[rand==0,]

            #juntando presencas e ausencias da calibracao
            ## presausTrainRaw <- rbind(presenciasTrain, ausenciasTrain)
            ## presausTrainRaw = data.frame(presausTrainRaw)
            ## presausTrainRaw = presausTrainRaw[!duplicated(presausTrainRaw),] #selecionar colunas de longitude e latitude
            ## presausTrainRaw<-presausTrainRaw[complete.cases(presausTrainRaw),]
            ## presausTrain = presausTrainRaw

            ##CRIANDO E RODANDO O MODELO##    
            MX <- maxent(predictors, presenciasTrain)

            #porcentajepres = round(0.25*nrow(presencias)) #seleccionar un porcentajes de filas de un data.frame
            #presencias.evaluacion<-presencias[sample(nrow(presencias), porcentajepres), ] #seleccionar ese porcentaje de filas aleatorias.

            #pegando a porcao dos dados separados para a avaliacao (validacao) do modelo
            presencias.evaluacion = presencias[rand==1,]
            presencias.evaluacion <- cbind(presencias.evaluacion$longitude,presencias.evaluacion$latitude)
            #porcentajeabs = round(0.25*nrow(pseudoausencia)) # choose pseudoabsences points evaluation
            #pseudoausencias.evaluacion <- pseudoausencia[sample(nrow(pseudoausencia), porcentajeabs), ]
            pseudoausencias.evaluacion = ausencias[rand==1,]
            pseudoausencias.evaluacion = cbind(pseudoausencias.evaluacion$longitude,pseudoausencias.evaluacion$latitude)

            ##RODANDO A AVALIACAO DO MODELO##
            evaluacion=evaluate(presencias.evaluacion, pseudoausencias.evaluacion, MX, predictors)

            #registrando o valor de AUC em um objeto
            V[j]<-evaluacion@"auc" #sacamos el valor de auc (fíjate que es una @ en lugar de $ para mirar dentro de los slots)y guardamos en vector
        }, error=function(e){cat("ERROR :",conditionMessage(e),"\n")}
   )}  

    #registrando a media de AUC
    auc<-mean(V, na.rm=T)#media de los vectores de las iternaciones de j

    #registrando o valor de AUC medio em uma tabela
    fila=fila+1
    resultados.evaluacion.MX[fila, "Species"]<-splist[especie]
    resultados.evaluacion.MX[fila, "auc"]<-evaluacion@auc
    
    #gravando um PDF com a AUC do modelo
    pdf(file=paste(projectFolder,"Maxent/",splist[especie],"/",splist[especie],'_ROC',".pdf",sep=""))
    plot(evaluacion, "ROC", cex=0.3)
    dev.off()

    ##PROJETANDO o nicho no espaco atraves do modelo ajustado##
    projecaoSuitability <- predict(predictors, MX)

    #gravando um raster com o mapa de projecao gerado pelo modelo
    writeRaster(projecaoSuitability,filename=paste(projectFolder,"Maxent/",splist[especie],"/",splist[especie],".grd", sep=""),overwrite=T)

    #gravando um PDF com o mapa gerado pelo modelo
    pdf(file=paste(projectFolder,"Maxent/",splist[especie],"/",splist[especie],".pdf",sep=""))
    plot(projecaoSuitability, main=paste(splist[especie]))
    plot(AmSulShape,add=T)
    points(presencias$longitude,presencias$latitude,col='orange',pch=20,cex=0.75)
    points(presencias$longitude,presencias$latitude,col='red',cex=0.7)
    dev.off()

    #criando um mapa binario
    threshold = threshold(evaluacion,'spec_sens')
    bin <- projecaoSuitability > threshold #apply threshold to transform logistic output into binary maps

    #salvando um raster com o mapa binario
    writeRaster(bin,filename=paste(projectFolder,"Maxent/",splist[especie],"/",splist[especie],"BINARIO.asc",sep=""),overwrite=T)

    #gravando um PDF com o mapa gerado pelo modelo
    pdf(file=paste(projectFolder,"Maxent/",splist[especie],"/",splist[especie],"-BINARIO.pdf",sep=""))
    plot(bin, main= paste(splist[especie]))
    plot(AmSulShape,add=T)
    points(presencias$longitude,presencias$latitude,col='orange',pch=20,cex=0.75)
    points(presencias$longitude,presencias$latitude,col='red',cex=0.7)
    dev.off()

    ###PROJECAO PARA O PASSADO###
        
    #abrindo os dados de registros fosseis  para uma especie
    sp.fossil.data = occ.sps.fosseis[occ.sps.fosseis$species==splist[especie],] #ATENCAO: este script nao funciona se houver mais de um registro fossil por camada de tempo usada para projecao
        
    for(l in 1:nrow(sp.fossil.data)){
        #loop para cada registro fossil de uma especie

        #definindoo fossil
        sp.fossil = sp.fossil.data[l,]
    
        #abrindo as variaveis ambientais do tempo do fossil
        filesProjectionRaw <- stack(list.files(path = paste(envVarFolder,"dados_projeto/0",sp.fossil$kyr,sep=""), pattern='asc', full.names=T)) ###abrindo camandas para projecao (passado, futuro, outro local, etc)
        filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul
        files.crop.sub.projection <- dropLayer(filesProjection, c(1,2,5,6)) #removendo as camadas que mostraram correlacao
        predictorsProjection = files.crop.sub.projection #preditoras para o tempo do fossil

        ##PROJETANDO o nicho no espaco atraves do modelo ajustado##
        projecaoSuitabilityPassado <- predict(predictorsProjection, MX) #PASSADO

        #criando um objeto com as coordenadas do registro fossil
        fossilPoints = sp.fossil
        fossilPoints = cbind(fossilPoints$longitude, fossilPoints$latitude)

        #obtendo a projecao de qualidade de habitat especificamente para o ponto do fossil
        fossilPointsVars = extract(predictorsProjection,fossilPoints)
        fossilPoints.MX = predict(MX, fossilPointsVars)
        fossilPointsSuitability = rbind(fossilPointsSuitability,data.frame(splist[especie],sp.fossil$kyr,fossilPoints.MX))

        #salvando um raster com a projecao do modelo para o tempo do fossil
        writeRaster(projecaoSuitabilityPassado,filename=paste(projectFolder,"Maxent/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr," K years BP.asc", sep=""),overwrite=T)
            
        #salvando um PDF com a projecao do modelo para o tempo do fossil
        pdf(file=paste(projectFolder,"Maxent/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr," K years BP.pdf", sep=""))
        plot(projecaoSuitabilityPassado, main=paste(splist[especie],'-',sp.fossil$kyr,'K years BP',sep=''))
        plot(AmSulShape,add=T)
        #plotando as coordenadas do refistro fossil
        points(fossilPoints,col='orange',pch=20,cex=0.75)
        points(fossilPoints,col='red',cex=0.8)
        dev.off()

        #criando um mapa binario para a projecao do modelo (empregando o threshold que ja foi criado apos a avaliacao do modelo)
        bin <- projecaoSuitabilityPassado > threshold#apply threshold to transform logistic output into binary maps

        #salvando um raster com a projecao do modelo para o tempo do fossil
        writeRaster(bin,filename=paste(projectFolder,"Maxent/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr,"K years BP - BINARIO.asc",sep=""),overwrite=T)

        #salvando um PDF com a projecao do modelo para o tempo do fossil
        pdf(file=paste(projectFolder,"Maxent/Passado/",splist[especie],"/",splist[especie],"-",sp.fossil$kyr," K years BP - BINARIO.pdf",sep="")) 
        plot(bin, main= paste(splist[especie],"-",sp.fossil$kyr," K years BP.pdf",sep=""))
        plot(AmSulShape,add=T)
        #plotando as coordenadas do refistro fossil
        points(fossilPoints,col='orange',pch=20,cex=0.75)
        points(fossilPoints,col='red',cex=0.8)
        dev.off()

    }
}

#salvando a tabela de dados da avaliacao dos modelos
write.table(resultados.evaluacion.MX,file=paste(projectFolder,"Maxent/","AUCmodelos.csv",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")

write.csv(fossilPointsSuitability,file=paste(projectFolder,"Maxent/","suitabilityNoPontoFossil.csv",sep=""),row.names=F)
#registrando e informando o tempo de processamento
msgm = proc.time() - ptm
print(paste('Tempo gasto para rodar MAXENT: ', msgm[3]/60,' minutos'sep=''))
