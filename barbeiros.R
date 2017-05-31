##SCRIPT PARA MODELAGEM DE DISTRIBUICAO DE ESPECIES (PRESENTE E FUTURO) USANDO MAXENT##

library(raster)
library(maptools)
library(dismo)
Sys.setenv(JAVA_HOME='C:/Program Files/Java/jdk1.8.0_73/jre') # for 64-bit version
#Windows#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_91') # for 64-bit version
library(rJava)
#source("C:/Users/Lucas/Documents/Projeto Chagas/Barbeiros/TSSmaxent.R") #sempre verificar aqui o caminho para o arquivo TSSmaxent.R
source("/home/anderson/R/R-Scripts/TSSmaxent.R")

###PRIMEIRA PARTE: planilha de presencas, backgrownd e variaveis ambientais###

##definindo as pastas de trabalho
#envVarFolder = "C:/Users/Lucas/Documents/Projeto Chagas/Variaveis Climaticas"
#spOccFolder = "C:/Users/Lucas/Documents/Projeto Chagas/Ocorrencias"
#projectFolder = "C:/Users/Lucas/Documents/Projeto Chagas/Barbeiros/"
#anderson
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/bcmidbi_2-5m _asc_America_Sul"
spOccFolder = "/home/anderson/PosDoc/dados_ocorrencia/PO_unique/"
projectFolder = "/home/anderson/Documentos/Projetos/Barbeiros_Lucas/"

##abrindo as variaveis climaticas
##abrindo shape da America do Sul
##AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp")

##abrindo e cortando camadas de variaveis ambientais para o presente
filesRaw <- stack(list.files(path=envVarFolder, pattern='.grd', full.names=TRUE)) #modificar a extensao .bil de acordo com os arquivos
#files = mask(filesRaw,AmSulShape) #cortando para Am. do Sul

##abrindo e cortando camadas de variaveis ambientais para projecao
filesProjectionOtimistaRaw <- stack(list.files(path=paste(envVarFolder,"/futuro/cenario_otimista",sep=''), pattern='.asc', full.names=TRUE))
filesProjectionPessimistaRaw <- stack(list.files(path=paste(envVarFolder,"/futuro/cenario_pessimista",sep=''), pattern='.asc', full.names=TRUE)) 
##filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul

##remove highly correlated variables Bio1,Bio3,Bio9,Bio13,Bio14
files.crop.sub = filesRaw[[c('bcmidbi7','bcmidbi15','bcmidbi11','bcmidbi16','bcmidbi17')]] #choose selected layers
files.crop.sub.projection.otimista = filesProjectionOtimistaRaw[[c('bio7','bio15','bio11','bio16','bio17')]] #choose selected layers
files.crop.sub.projection.pessimista = filesProjectionPessimistaRaw[[c('bio7','bio15','bio11','bio16','bio17')]] #choose selected layers

##definindo os objetos para as variaveis preditoras
predictors <- stack(files.crop.sub)
predictorsProjectionOtimista = files.crop.sub.projection.otimista
predictorsProjectionPessimista = files.crop.sub.projection.pessimista

##adicionando impacto humano as variaveis ambientais
##abrindo
hii = raster(x='/home/anderson/Downloads/hii-s-america-geo-grid/hii_s_america_grid/hii_s_amer')

##crop hii
areaExt = extent(predictors[[1]])
hiiSAraw = crop(hii,areaExt)

##hii ajustado
raster1 = res(predictors[[1]])
raster2 = res(hiiSAraw)
fator = round(raster1[1]/raster2[1])  ##raster maior / raster menor
hiiSA = aggregate(hiiSAraw,fator)
##verificando
res(hiiSA)
res(dado)

predictorsHII <- stack(c(predictors,hiiSA))
predictorsProjectionOtimista = files.crop.sub.projection.otimista 
predictorsProjectionPessimista = files.crop.sub.projection.pessimista

##Criando objeto com a lista de especies
occ.sps <- list.files(paste(spOccFolder,sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))

##criando uma tabela vazia para salvador alguns dados
tabRes = data.frame(c(sp=character(),auc=numeric(),tss=numeric(),threshold=numeric()))


###SEGUNDA PARTE: rodando SDMs para as especies (e fazendo projecoes) - SEM impacto humano###

for (i in 1:length(splist)){
    especie = splist[i] #selecting the species
    sp.file <- read.csv(paste(spOccFolder,"/",especie,".csv",sep=""),header=TRUE) ### read sp occurrence
    sp.occ <- sp.file[,2:3] ## select long-lat
    
	##CRIANDO E RODANDO O MODELO - Atencao: precisa ter uma pasta para cada sp
	MX <- maxent(x=predictors,p=sp.occ,path=paste(projectFolder,especie,sep=""), 
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
	mres = read.csv(paste(projectFolder,especie,"/maxentResults.csv",sep=''),header=TRUE)
	threshold = mres$X10.percentile.training.presence.logistic.threshold[nrow(mres)] #threshold
	aucModel = mres$Test.AUC[nrow(mres)]
	tss = TSSmaxent(paste(projectFolder,especie,'/',sep=''))$TSS #Atencao: o script do TSS deve estar na 

	##gravando uma tabela
	tabRes = rbind(tabRes,data.frame(sp=especie,auc=aucModel,tss=tss,threshold=threshold))

	##fazendo as projecoes
    projecaoSuitability = predict(MX,predictors) #projecao para o presente
    projecaoFuturoOtimista = predict(MX,predictorsProjectionOtimista) #projecao para o futuro - cenario otimista
    projecaoFuturoPessimista = predict(MX,predictorsProjectionPessimista) #projecao para o futuro - cenario pessimista
                        
    ##gravando arquivo raster dos mapas de projecao gerados pelo modelo
	##mapas do presente
    writeRaster(projecaoSuitability,filename=paste(projectFolder,especie,"/",especie,".asc", sep=""),overwrite=TRUE)
    writeRaster(projecaoSuitability>threshold,filename=paste(projectFolder,especie,"/",especie,"BIN.asc", sep=""),overwrite=TRUE)
	##mapas do futuro otimista
    writeRaster(projecaoFuturoOtimista,filename=paste(projectFolder,especie,"/",especie,"Otimista.asc", sep=""),overwrite=TRUE)
    writeRaster(projecaoFuturoOtimista>threshold,filename=paste(projectFolder,especie,"/",especie,"OtimistaBIN.asc", sep=""),overwrite=TRUE)
	##mapas do futuro pessimista
    writeRaster(projecaoFuturoPessimista,filename=paste(projectFolder,especie,"/",especie,"Pessimista.asc", sep=""),overwrite=TRUE)
    writeRaster(projecaoFuturoPessimista>threshold,filename=paste(projectFolder,especie,"/",especie,"PessimistaBIN.asc", sep=""),overwrite=TRUE)

}

##gravando a tabela de resultados completa
write.csv(tabRes,path=paste(projectFolder,tabela_resultados.csv,sep=''))


###TERCEIRA PARTE: rodando SDMs para as especies (e fazendo projecoes) - COM impacto humano###

##criando uma tabela vazia para salvador alguns dados
tabResHII = data.frame(c(sp=character(),auc=numeric(),tss=numeric(),threshold=numeric()))


for (i in 1:length(splist)){
    especie = splist[i] #selecting the species
    sp.file <- read.csv(paste(spOccFolder,"/",especie,".csv",sep=""),header=TRUE) ### read sp occurrence
    sp.occ <- sp.file[,2:3] ## select long-lat
    
    ##CRIANDO E RODANDO O MODELO - Atencao: precisa ter uma pasta para cada sp
    MX <- maxent(x=predictorsHII,p=sp.occ,path=paste(projectFolder,"impacto_humano",especie,sep=""), 
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
    mres = read.csv(paste(projectFolder,"impacto_humano",especie,"/maxentResults.csv",sep=''),header=TRUE)
    threshold = mres$X10.percentile.training.presence.logistic.threshold[nrow(mres)] #threshold
    aucModel = mres$Test.AUC[nrow(mres)]
    tss = TSSmaxent(paste(projectFolder,especie,'/',sep=''))$TSS #Atencao: o script do TSS deve estar na 
    
    ##gravando uma tabela
    tabResHII = rbind(tabRes,data.frame(sp=especie,auc=aucModel,tss=tss,threshold=threshold))
    
    ##fazendo as projecoes
    projecaoSuitability = predict(MX,predictors) #projecao para o presente
    projecaoFuturoOtimista = predict(MX,predictorsProjectionOtimista) #projecao para o futuro - cenario otimista
    projecaoFuturoPessimista = predict(MX,predictorsProjectionPessimista) #projecao para o futuro - cenario pessimista
    
    ##gravando arquivo raster dos mapas de projecao gerados pelo modelo
    ##mapas do presente
    writeRaster(mean(projecaoSuitability),filename=paste(projectFolder,"impacto_humano",especie,"/",especie,".asc", sep=""),overwrite=TRUE)
    writeRaster(projecaoSuitability>threshold,filename=paste(projectFolder,"impacto_humano",especie,"/",especie,"BIN.asc", sep=""),overwrite=TRUE)
    ##mapas do futuro otimista
    writeRaster(projecaoFuturoOtimista,filename=paste(projectFolder,"impacto_humano",especie,"/",especie,"Otimista.asc", sep=""),overwrite=TRUE)
    writeRaster(projecaoFuturoOtimista>threshold,filename=paste(projectFolder,"impacto_humano",especie,"/",especie,"OtimistaBIN.asc", sep=""),overwrite=TRUE)
    ##mapas do futuro pessimista
    writeRaster(projecaoFuturoPessimista,filename=paste(projectFolder,"impacto_humano",especie,"/",especie,"Pessimista.asc", sep=""),overwrite=TRUE)
    writeRaster(projecaoFuturoPessimista>threshold,filename=paste(projectFolder,"impacto_humano",especie,"/",especie,"PessimistaBIN.asc", sep=""),overwrite=TRUE)
    
}

##gravando a tabela de resultados completa
write.csv(tabResHII,path=paste(projectFolder,"impacto_humano",tabela_resultados.csv,sep=''))


##indices para RISCO DE INFECCAO por especie de vetor
tabBarb = read.csv("/home/anderson/Documentos/Projetos/Vetor x IDH para Doenca de Chagas/Taxa de infeccao natural vetores 2007-2011.csv",header=TRUE)
infecBarb = sort(tabBarb[,3],decreasing=TRUE)
infecInd = sort(1:length(infecBarb),decreasing=TRUE) #idice de infeccao para confeccao dos mapas


