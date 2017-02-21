##########################################################################
############################## SDM WITH MAXENT  ##########################
###########################Anderson A. Eduardo############################
##########################################################################
###names()
##help()
library(raster)
library(maptools)
library(dismo)
Sys.setenv(JAVA_HOME='/usr/lib/jvm/java-7-openjdk-amd64') # for 64-bit version
#Windows#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_91') # for 64-bit version
library(rJava)

###PRIMEIRA PARTE: planilha de presencas, backgrownd e variaveis ambientais

##definindo as pastas de trabalho
envVarFolder = "C:/Users/Lucas/Documents/Projeto Chagas/Variáveis Climáticas"
spOccFolder = "C:/Users/Lucas/Documents/Projeto Chagas/Tabelas Maxent"
projectFolder = "C:/Users/Lucas/Documents/Projeto Chagas/Barbeiros"
#anderosn
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/bcmidbi_2-5m _asc_America_Sul"
spOccFolder = "/home/anderson/PosDoc/dados_ocorrencia/PO_unique/"
projectFolder = "/home/anderson/PosDoc/teste/"

##abrindo as variaveis climaticas
##abrindo shape da America do Sul
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp")

##abrindo e cortando camadas de variaveis ambientais para o presente
filesRaw <- stack(list.files(path=envVarFolder, pattern='grd', full.names=TRUE)) 
#files = mask(filesRaw,AmSulShape) #cortando para Am. do Sul

##abrindo e cortando camads de variaveis ambientais para projecao
##filesProjectionRaw <- stack(list.files(path=paste(envVarFolder,"dados_projeto/021",sep=''), pattern='asc', full.names=T)) ###abrindo camandas para projecao (passado, futuro, outro local, etc)
##filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul


##remove highly correlated variables Bio1,Bio3,Bio9,Bio13,Bio14
files.crop.sub = filesRaw[[c('bcmidbi15','bcmidbi11','bcmidbi16','bcmidbi17','bcmidbi7')]] #choose selected layers
##files.crop.sub.projection <- dropLayer(filesProjection, c(1,2,5,6))

##definindo os objetos para as variaveis preditoras
predictors <- stack(files.crop.sub)
##predictorsProjection = files.crop.sub.projection

##Criando objeto com a lista de especies
occ.sps <- list.files(paste(spOccFolder,sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))


###SEGUNDA PARTE: rodando SDMs para as especies (e fazendo projecoes)

for (i in 1:length(splist)){
    
    especie = splist[i] #selecting the species
    sp.file <- read.csv(paste(spOccFolder,especie,".csv",sep=""),header=TRUE) ### read sp occurrence
    sp.occ <- sp.file[,2:3] ## select lat-long 
    
            ##CRIANDO E RODANDO O MODELO (esquema do SWD - sample with data)##    
            MX <- maxent(x=predictors,p=sp.occ,args=c('responsecurves=true','jackknife=true','randomseed=true','randomtestpoints=0','betamultiplier=1','replicates=1','replicatetype=Subsample','writebackgroundpredictions=true','linear=true','quadratic=true','product=false','threshold=false','hinge=false','maximumiterations=1000','convergencethreshold=1.0E-5','threads=2'))

            ##fazendo a projecao
            projecaoSuitability = predict(predictors,MX)
                        
    ##gravando um raster com o mapa de projecao gerado pelo modelo
    writeRaster(projecaoSuitability,filename=paste(projectFolder,especie,especie,".asc", sep=""),overwrite=TRUE)
}
    
