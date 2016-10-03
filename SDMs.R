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

#abrindo e cortando camads de variaveis ambientais para o passado
filesProjectionRaw <- stack(list.files(path=paste(envVarFolder,"dados_projeto/011",sep=''), pattern='asc', full.names=T)) ###abrindo camandas para projecao (passado, futuro, outro local, etc)
filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul

#testando correcaloes
## test<-getValues(files)
## cor.matrix <- as.data.frame(cor(test, use="complete.obs"))
#write.csv(cor.matrix,'cor_matrix.csv')

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
  #res <- occ(query = sp, from = 'gbif', limit = 10000)
  #locs<-occ2df(res)
  locs = read.csv(file=paste("/home/anderson/R/PosDoc/dados_ocorrencia/",sp,".csv",sep=""),header=T,stringsAsFactors=FALSE)    
  locs2<-locs[,2:3]
  locs2[locs2 == 0] <- NA
  locs3<-locs2[complete.cases(locs2),]
  locs4<-round(locs3, digits = 4)
  locs5<-locs4[!duplicated(locs4), ]
  locs6<-cbind(sp,locs5) 
  #write.csv(locs6, file=paste("/home/anderson/R/PosDoc/dados_ocorrencia/",unique(Especies$Especie)[i],".csv",sep=""),row.names=F)
  write.csv(locs6, file=paste("/home/anderson/R/PosDoc/dados_ocorrencia/PO_unique/",sp,".csv",sep=""),row.names=F)
}

########## Criando objetos com a lista de especies #############
occ.sps <- list.files(paste(spOccFolder,sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))
##fosseis
occ.sps.fosseis = read.csv(paste(spOccFolder,"fosseis.csv",sep=''),header=T)
splist.fosseis = lapply(occ.sps.fosseis[,1],as.character)

###inspeção visual exploratoria
data(wrld_simpl)

especie = splist[1] #escolher qual especie
sp.file <- read.csv(paste(spOccFolder,especie,".csv",sep=""),h=T) ### read sp occurrence
sp.occ <- sp.file[,2:3] ## select lat long columns
sp.occ = sp.occ[sp.occ$lat > -50,] #eliminar pontos 
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
sp.file <- read.csv(paste(spOccFolder, splist[especie],".csv",sep=""),h=T) ### read sp occurrence
sp.occ <- sp.file[,2:3] ## select lat long columns
me <- maxent(predictors,sp.occ, args=c("raw","maximumiterations=1000"),path=paste("/home/anderson/R/PosDoc/teste/",splist[especie],sep="")) ## run maxent model with raw output
crs <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
r <- predict(me,predictors,crs=crs)
#r <- predict(me,predictorsProjection,crs=crs) #para projetar passado, futuro, etc
names(r) <- splist[especie]
writeRaster(r, filename=paste("/home/anderson/R/PosDoc/teste/",splist[especie],"/",splist[especie],".asc",sep=""), overwrite=T)
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


##GRAFICOS AGRUPADOS##

library(rasterVis)

#definindo objeto com os nomes
teste = 'Maxent'
especie = 'Caiman crocodilus'

setwd(paste(projectFolder,teste,'/Raster Layers',sep='')) #presente
files = list.files(paste(getwd()),full.names=TRUE,pattern='.asc')
files=c(files[2],files[4],files[6],files[8],files[10]) #selecionando camadas e consertando a ordem das camadas

setwd(paste(projectFolder,teste,'/Passado/Raster Layers',sep='')) #passado
filesPass = list.files(paste(getwd()),full.names=TRUE,pattern='.asc')
filesPass=c(filesPass[2],filesPass[4],filesPass[6],filesPass[8],filesPass[10],filesPass[12],filesPass[14],filesPass[16],filesPass[18],filesPass[20],filesPass[22],filesPass[24],filesPass[26],filesPass[28],filesPass[30],filesPass[32],filesPass[34],filesPass[36],filesPass[38],filesPass[40]) #selecionando camadas e consertando a ordem das camadas

species.layers = stack(c(files,filesPass))

#cortando e gravando varios rasters com o shapefile da america do sul
for (i in 1:length(names(species.layers))){
    species.layers[[i]] = mask(species.layers[[i]],AmSulShape)
    writeRaster(species.layers[[i]], filename=paste("/home/anderson/R/PosDoc/teste/",teste,'/Raster Layers Cortados/',names(species.layers)[i],".asc",sep=""), overwrite=T)
}

#ajustando os nomes dos subgraficos
names(species.layers)
rasterNames = gsub("Caiman_latirostris_._","",names(species.layers))
rasterNames = gsub("presente","Presente",rasterNames)

cols <- colorRampPalette(brewer.pal(9,"YlGn")) #escala de cores amarelo-verde

setwd(paste("/home/anderson/R/PosDoc/teste/",teste,sep=''))
jpeg(filename='Caiman.jpg')
levelplot(species.layers,main='',col.regions=cols,names.attr=rasterNames) + layer(sp.polygons(AmSul)) + layer(panel.xyplot(-41.553056, -12.393417,pch=17,col='red',cex=1),columns=1) + layer(panel.xyplot(-37.753611,-9.926944,pch=17,col='red',cex=1),columns=2)
dev.off()


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
            model <- pres ~ bioclim_10+bioclim_11+bioclim_16+bioclim_17
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
    writeRaster(projecaoSuitability,filename=paste(projectFolder,"Random Forest/",splist[especie],"/",splist[especie],".grd", sep=""),overwrite=T)

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
            model <- pres ~ bioclim_10 + I(bioclim_10^2) + bioclim_11 + I(bioclim_11^2) + bioclim_16 + I(bioclim_16^2) + bioclim_17 + I(bioclim_17^2)
            GLM <- glm(model, family=binomial(link=logit), data=presausTrain)

            
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
            evaluacion=evaluate(presencias.evaluacion, pseudoausencias.evaluacion, GLM, predictors)

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
    writeRaster(projecaoSuitability,filename=paste(projectFolder,"GLM/",splist[i],"/",splist[i],".grd", sep=""),overwrite=T)

    #gravando um PDF com o mapa gerado pelo modelo
    pdf(file=paste(projectFolder,"GLM/",splist[i],"/",splist[i],".pdf",sep=""))
    plot(projecaoSuitability, main=paste(splist[i]))
    plot(AmSulShape,add=T)
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
        projecaoSuitabilityPassado <- predict(predictorsProjection, GLM,type='response') #PASSADO

        #criando um objeto com as coordenadas do registro fossil
        fossilPoints = sp.fossil
        fossilPoints = cbind(fossilPoints$longitude, fossilPoints$latitude)

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
fossilPointsSuitability = NULL #objetvo em que serao gravadas as projcoes de suitability especificamente para cada ponto de registro fossil

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
    
    for (j in 1:10){
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
print(paste('Tempo gasto para rodar MAXENT: ', msgm[3]/60,' minutos',sep=''))


