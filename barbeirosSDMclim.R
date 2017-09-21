##SCRIPT PARA DISTRIBUICAO DE BARBEIROS (PRESENTE E FUTURO) - APENAS VARIAVEIS CLIMATICAS##

library(raster)
library(maptools)
library(dismo)
#Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8.0_131/bin') # for 64-bit version
#Windows#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_91') # for 64-bit version
library(rJava)
#source("J:/Anderson_Eduardo/TSSmaxent.R") #sempre verificar aqui o caminho para o arquivo TSSmaxent.R
source("/home/anderson/R/R-Scripts/TSSmaxent.R")


###PRIMEIRA PARTE: planilha de presencas, backgrownd e variaveis ambientais###


##definindo as pastas de trabalho
## envVarFolder = "J:/Lucas/Modelagem barbeiros/Variaveis Climaticas"
## spOccFolder = "J:/Lucas/Modelagem barbeiros/Ocorrencias"
## projectFolder = "J:/Lucas/Modelagem barbeiros/"
#anderson
envVarFolder = "/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Variaveis Climaticas"
spOccFolder = "/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/Ocorrencias/"
projectFolder = "/home/anderson/Documentos/Projetos/Distribuicao de barbeiros com interacao com humanos/"

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


##definindo os objetos para as variaveis preditoras (SEM IMPACTO HUMANO)
predictors <- stack(files.crop.sub)
predictorsProjectionOtimista = files.crop.sub.projection.otimista 
predictorsProjectionPessimista = files.crop.sub.projection.pessimista

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


###TERCEIRA PARTE: gerando mapas de sobreposicao (riqueza) - SEM impacto humano###


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
listaPresente = grep(list.files(paste(projectFolder,"resultados nicho climatico/Projecoes",sep=""),full.names=TRUE),pattern='Otimista|Pessimista',inv=T,value=T)
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
writeRaster(x=mapaRiscoPresente,filename=paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoPresente.asc',sep=''),overwrite=TRUE)

#futuro otimista
camadasFuturoOtimista = stack(list.files(paste(projectFolder,"resultados nicho climatico/Projecoes",sep=""),pattern = "OtimistaBIN.asc",full.names=TRUE))
listaNomes = names(camadasFuturoOtimista)
listaNomes = gsub(pattern='OtimistaBIN',replacement='',x=listaNomes)
infecIndOrdered = tabBarb$taxaInfeccaonatural[match(listaNomes,tabBarb$sp)]/100 #taxa de infeccao natural na ordem dos rasters (de 0 a 1)
#
mapaRiscoFuturoOtimista = sum(camadasFuturoOtimista*infecIndOrdered)
mapaRiscoFuturoOtimista = mapaRiscoFuturoOtimista/length(infecIndOrdered)
plot(mapaRiscoFuturoOtimista)
writeRaster(x=mapaRiscoFuturoOtimista,filename=paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoFuturoOtimista.asc',sep=''),overwrite=TRUE)

#futuro pessimista
camadasFuturoPessimista = stack(list.files(paste(projectFolder,"resultados nicho climatico/Projecoes",sep=""),pattern = "PessimistaBIN.asc",full.names=TRUE))
listaNomes = names(camadasFuturoOtimista)
listaNomes = gsub(pattern='OtimistaBIN',replacement='',x=listaNomes)
infecIndOrdered = tabBarb$taxaInfeccaonatural[match(listaNomes,tabBarb$sp)]/100 #taxa de infeccao natural na ordem dos rasters (de 0 a 1)
#
mapaRiscoFuturoPessimista = sum(camadasFuturoPessimista*infecIndOrdered)
mapaRiscoFuturoPessimista = mapaRiscoFuturoPessimista/length(infecIndOrdered)
plot(mapaRiscoFuturoPessimista)
writeRaster(x=mapaRiscoFuturoPessimista,filename=paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''),overwrite=TRUE)


###QUINTA PARTE: estatisticas sumarias a partir dos mapas - SEM impacto humano###


#presente 

#tamanho da area do quartil superior (para riqueza e risco), para comparar os cenarios
mapaRiquezaPresente = raster(paste(projectFolder,'resultados nicho climatico/Mapas de riqueza/mapaRiquezaPresente.asc',sep=''))
mapaRiscoPresente = raster(paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoPresente.asc',sep=''))
hightRiqPres= mapaRiquezaPresente > quantile(mapaRiquezaPresente, 0.75,na.rm=TRUE) #raster quartil superior riqueza
hightRiscPres = mapaRiscoPresente > quantile(mapaRiscoPresente, 0.75,na.rm=TRUE) #raster quartil superior risco

percCelRiqPres = freq(hightRiqPres,value=1)/(freq(hightRiqPres,value=0)+freq(hightRiqPres,value=1)) #percentagem celulas no quartil sup.
percCelRiscPres = freq(hightRiscPres,value=1)/(freq(hightRiscPres,value=0)+freq(hightRiscPres,value=1)) #percentagem de celulas no quartil sup.

#correlacao entre riqueza e risco
test <- getValues(stack(hightRiqPres,hightRiscPres))
corPres <- as.data.frame(cor(test, use="complete.obs",method='spearman'))
##write.csv(cor.matrix,'cor_matrix.csv')

#futuro otimista 

#tamanho da area do quartil superior (para riqueza e risco), para comparar os cenarios
mapaRiquezaFuturoOtimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de riqueza/mapaRiquezaFuturoOtimista.asc',sep=''))
mapaRiscoFuturoOtimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoFuturoOtimista.asc',sep=''))
hightRiqOtim= mapaRiquezaFuturoOtimista > quantile(mapaRiquezaFuturoOtimista, 0.75,na.rm=TRUE) #raster quartil superior riqueza
hightRiscOtim = mapaRiscoFuturoOtimista > quantile(mapaRiscoFuturoOtimista, 0.75,na.rm=TRUE) #raster quartil superior risco

percCelRiqOtim = freq(hightRiqOtim,value=1)/(freq(hightRiqOtim,value=0)+freq(hightRiqOtim,value=1)) #percentagem celulas no quartil sup.
percCelRiscOtim =  freq(hightRiscOtim,value=1)/(freq(hightRiscOtim,value=0)+freq(hightRiscOtim,value=1)) #percentagem de celulas no quartil sup.

#correlacao entre riqueza e risco
test <- getValues(stack(hightRiqOtim,hightRiscOtim))
corOtim <- as.data.frame(cor(test, use="complete.obs",method='spearman'))
##write.csv(cor.matrix,'cor_matrix.csv')

#futuro pessimista 

#tamanho da area do quartil superior (para riqueza e risco), para comparar os cenarios
mapaRiquezaFuturoPessimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''))
mapaRiscoFuturoPessimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''))
hightRiqPess= mapaRiquezaFuturoPessimista > quantile(mapaRiquezaFuturoPessimista, 0.75,na.rm=TRUE) #raster quartil superior riqueza
hightRiscPess = mapaRiscoFuturoPessimista > quantile(mapaRiscoFuturoPessimista, 0.75,na.rm=TRUE) #raster quartil superior risco

percCelRiqPess = freq(hightRiqPess,value=1)/(freq(hightRiqPess,value=0)+freq(hightRiqPess,value=1)) #percentagem celulas no quartil sup.
percCelRiscPess = freq(hightRiscPess,value=1)/(freq(hightRiscPess,value=0)+freq(hightRiscPess,value=1))  #percentagem de celulas no quartil sup.

#correlacao entre riqueza e risco
test <- getValues(stack(hightRiqPess,hightRiscPess))
corPess <- as.data.frame(cor(test, use="complete.obs",method='spearman'))
##write.csv(cor.matrix,'cor_matrix.csv')

#salvando resultados em tabelas
tab = data.frame(scenario=c('pres','fut_otim','fut_pess'),
                 quantile75riq = c(quantile(mapaRiquezaPresente, 0.75,na.rm=TRUE),quantile(mapaRiquezaFuturoOtimista, 0.75,na.rm=TRUE),quantile(mapaRiquezaFuturoPessimista, 0.75,na.rm=TRUE)),
                 quantile75risc = c(quantile(mapaRiscoPresente, 0.75,na.rm=TRUE),quantile(mapaRiscoFuturoOtimista, 0.75,na.rm=TRUE), quantile(mapaRiscoFuturoPessimista, 0.75,na.rm=TRUE)),
                 percCellRiq = c(percCelRiqPres,percCelRiqOtim,percCelRiqPess),
                 percCellRisc = c(percCelRiscPres,percCelRiscOtim,percCelRiscPess),
                 corRiqRisc = c(corPres[1,2],corOtim[1,2],corPess[1,2])
)

write.csv(tab,paste(projectFolder,'resultados nicho climatico/statsRes.csv',sep=''),row.names = FALSE)


###SEXTA PARTE: figuras dos mapas - SEM impacto humano###


mapaRiquezaPresente = raster(paste(projectFolder,'resultados nicho climatico/Mapas de riqueza/mapaRiquezaPresente.asc',sep=''))
mapaRiquezaFuturoOtimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''))
mapaRiquezaFuturoPessimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''))
mapaRiscoPresente = raster(paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoPresente.asc',sep=''))
mapaRiscoFuturoOtimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoFuturoOtimista.asc',sep=''))
mapaRiscoFuturoPessimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''))






############################
###RECORTAR PARA O BRASIL####
############################






mapaRiquezaPresenteBR = mask(mapaRiquezaPresente, mask=(AmSulShape[AmSulShape$CNTRY_NAME=='Brazil',]))

areaBR = extent(-78.00635,-31.71555,-37.8679,8.156474)

extent(mapaRiquezaPresenteBR) = areaBR

plot(mapaRiquezaPresenteBR)
plot(AmSulShape,add=T)


jpeg(filename=paste(projectFolder,'resultados nicho climatico/mapas.jpeg',sep=''),width=1700,height=1200)
par(mfrow=c(2,3),mar=c(5,5,5,14))
##riqueza
plot(mapaRiquezaPresente,main='Presente',legend=FALSE,cex.axis=2,cex.main=4);plot(AmSulShape,add=TRUE); grid()
plot(mapaRiquezaFuturoOtimista,main='2070 otimista',legend=FALSE,cex.axis=2,cex.main=4); plot(AmSulShape,add=TRUE); grid()
plot(mapaRiquezaFuturoPessimista,main='2070 pessimista',legend=FALSE,cex.axis=2,cex.main=4); plot(AmSulShape,add=TRUE); grid()
plot(mapaRiquezaFuturoPessimista,legend.only=TRUE,legend.width=3,axis.args=list(cex.axis=2),legend.args=list(text='Riqueza',font=2,side=4,line=4.5,cex=2.2,cex.axis=0.2)) #legenda
##risco
plot(mapaRiscoPresente,legend=FALSE,cex.axis=2); plot(AmSulShape,add=TRUE); grid()
plot(mapaRiscoFuturoOtimista,legend=FALSE,cex.axis=2,cex.lab=2); plot(AmSulShape,add=TRUE); grid()
plot(mapaRiscoFuturoPessimista,legend=FALSE,cex.axis=2); plot(AmSulShape,add=TRUE); grid()
plot(mapaRiscoFuturoPessimista,legend.only=TRUE,legend.width=3,axis.args=list(cex.axis=2),legend.args=list(text='Risco de vetor infectado',font=2,side=4,line=5.7,cex=2.2,cex.axis=0.2)) #legenda
dev.off()
