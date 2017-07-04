##SCRIPT PARA DISTRIBUICAO DE BARBEIROS (PRESENTE E FUTURO) - VARIAVEIS CLIMATICAS + INDICE DE IMPACTO HUMANO##

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
envVarFolder = "/home/anderson/Documentos/Projetos/Barbeiros_Lucas/Variaveis Climaticas"
spOccFolder = "/home/anderson/Documentos/Projetos/Barbeiros_Lucas/Ocorrencias/"
projectFolder = "/home/anderson/Documentos/Projetos/Barbeiros_Lucas/resultados nicho climatico + impacto humano/"

##abrindo as variaveis climaticas
##abrindo shape da America do Sul
AmSulShape = readShapePoly("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp")

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

##definindo os objetos para as variaveis preditoras (SEM IMPACTO HUMANO)
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


###QUINTA PARTE: estatisticas sumarias a partir dos mapas - SEM impacto humano###

#presente 

#tamanho da area do quartil superior (para riqueza e risco), para comparar os cenarios
mapaRiquezaPresente = raster(paste(projectFolder,'Mapas de riqueza/mapaRiquezaPresente.asc',sep=''))
mapaRiscoPresente = raster(paste(projectFolder,'Mapas de risco/mapaRiscoPresente.asc',sep=''))
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
mapaRiquezaFuturoOtimista = raster(paste(projectFolder,'Mapas de riqueza/mapaRiquezaFuturoOtimista.asc',sep=''))
mapaRiscoFuturoOtimista = raster(paste(projectFolder,'Mapas de risco/mapaRiscoFuturoOtimista.asc',sep=''))
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
mapaRiquezaFuturoPessimista = raster(paste(projectFolder,'Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''))
mapaRiscoFuturoPessimista = raster(paste(projectFolder,'Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''))
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

write.csv(tab,paste(projectFolder,'statsRes.csv'),row.names = FALSE)


###SEXTA PARTE: figuras dos mapas - SEM impacto humano###
mapaRiquezaPresente = raster(paste(projectFolder,'Mapas de riqueza/mapaRiquezaPresente.asc',sep=''))
mapaRiquezaFuturoOtimista = raster(paste(projectFolder,'Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''))
mapaRiquezaFuturoPessimista = raster(paste(projectFolder,'Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''))
mapaRiscoPresente = raster(paste(projectFolder,'Mapas de risco/mapaRiscoPresente.asc',sep=''))
mapaRiscoFuturoOtimista = raster(paste(projectFolder,'Mapas de risco/mapaRiscoFuturoOtimista.asc',sep=''))
mapaRiscoFuturoPessimista = raster(paste(projectFolder,'Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''))

jpeg(filename=paste(projectFolder,'mapas.jpeg'),,width=1700,height=1200)
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


###SETIMA PARTE: comparando areas de ALTA RIQUEZA entre as projecoes com e sem impacto humano###

projectFolder = '/home/anderson/Documentos/Projetos/Barbeiros_Lucas/'

##presente
mapaRiquezaPresente = raster(paste(projectFolder,'resultados nicho climatico/Mapas de riqueza/mapaRiquezaPresente.asc',sep=''))
mapaRiquezaPresenteHii = raster(paste(projectFolder,'resultados nicho climatico + impacto humano/Mapas de riqueza/mapaRiquezaPresente.asc',sep=''))
#
hightRiqPres= mapaRiquezaPresente > quantile(mapaRiquezaPresente, 0.75,na.rm=TRUE) #raster quartil superior riqueza
hightRiqPresHii= mapaRiquezaPresenteHii > quantile(mapaRiquezaPresenteHii, 0.75,na.rm=TRUE) #raster quartil superior riqueza

##futuro otimista
mapaRiquezaFuturoOtimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de riqueza/mapaRiquezaFuturoOtimista.asc',sep=''))
mapaRiquezaFuturoOtimistaHii = raster(paste(projectFolder,'resultados nicho climatico + impacto humano/Mapas de riqueza/mapaRiquezaFuturoOtimista.asc',sep=''))
#
hightRiqOtim= mapaRiquezaFuturoOtimista > quantile(mapaRiquezaFuturoOtimista, 0.75,na.rm=TRUE) #raster quartil superior riqueza
hightRiqOtimHii= mapaRiquezaFuturoOtimistaHii > quantile(mapaRiquezaFuturoOtimistaHii, 0.75,na.rm=TRUE) #raster quartil superior riqueza

##futuro pessimista
mapaRiquezaFuturoPessimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''))
mapaRiquezaFuturoPessimistaHii = raster(paste(projectFolder,'resultados nicho climatico + impacto humano/Mapas de riqueza/mapaRiquezaFuturoPessimista.asc',sep=''))
#
hightRiqPess= mapaRiquezaFuturoPessimista > quantile(mapaRiquezaFuturoPessimista, 0.75,na.rm=TRUE) #raster quartil superior riqueza
hightRiqPessHii= mapaRiquezaFuturoPessimistaHii > quantile(mapaRiquezaFuturoPessimistaHii, 0.75,na.rm=TRUE) #raster quartil superior riqueza

##unificando para contabilizar areas
hightRiqPresUni = hightRiqPres*1 + hightRiqPresHii*2
hightRiqOtimUni = hightRiqOtim*1 + hightRiqOtimHii*2
hightRiqPessUni = hightRiqPess*1 + hightRiqPessHii*2

##correlacao areas alta riqueza modelos com e sem HII
##presente
testRiqPres <- getValues(stack(hightRiqPres,hightRiqPresHii))
corRiqPres <- as.data.frame(cor(testRiqPres, use="complete.obs",method='spearman'))
##futuro otimista
testRiqOtim <- getValues(stack(hightRiqOtim,hightRiqOtimHii))
corRiqOtim <- as.data.frame(cor(testRiqOtim, use="complete.obs",method='spearman'))
##futuro pessimista
testRiqPess <- getValues(stack(hightRiqPess,hightRiqPessHii))
corRiqPess <- as.data.frame(cor(testRiqPess, use="complete.obs",method='spearman'))

##correlacao mapa riqueza modelos com e sem HII## OBS: ACHEI MELHOR USAR ESSE NO MANUSCRITO
##presente
testRiqPres <- getValues(stack(mapaRiquezaPresente,mapaRiquezaPresenteHii))
corRiqPres <- as.data.frame(cor(testRiqPres, use="complete.obs",method='pearson'))
##futuro otimista
testRiqOtim <- getValues(stack(mapaRiquezaFuturoOtimista,mapaRiquezaFuturoOtimistaHii))
corRiqOtim <- as.data.frame(cor(testRiqOtim, use="complete.obs",method='pearson'))
##futuro pessimista
testRiqPess <- getValues(stack(mapaRiquezaFuturoPessimista,mapaRiquezaFuturoPessimistaHii))
corRiqPess <- as.data.frame(cor(testRiqPess, use="complete.obs",method='pearson'))

tabSobreposicaoRiqueza = data.frame(
    cenario = c('presente','futuroOtimista','futuroPessimista'),
    total_celulas_nichoClim = c(freq(hightRiqPres,value=1),freq(hightRiqOtim,value=1),freq(hightRiqPess,value=1)),
    percentual_da_AmSul = c(
        (freq(hightRiqPres,value=1)/(ncell(hightRiqPres)-freq(hightRiqPres,value=NA)))*100,
        (freq(hightRiqOtim,value=1)/(ncell(hightRiqPres)-freq(hightRiqPres,value=NA)))*100,
        (freq(hightRiqPess,value=1)/(ncell(hightRiqPres)-freq(hightRiqPres,value=NA)))*100 ),
    total_celulas_HII= c(freq(hightRiqPresHii,value=1),freq(hightRiqOtimHii,value=1),freq(hightRiqPessHii,value=1)),
    percentual_da_AmSul = c(
        (freq(hightRiqPresHii,value=1)/(ncell(hightRiqPresHii)-freq(hightRiqPresHii,value=NA)))*100,
        (freq(hightRiqOtimHii,value=1)/(ncell(hightRiqPresHii)-freq(hightRiqPresHii,value=NA)))*100,
        (freq(hightRiqPessHii,value=1)/(ncell(hightRiqPresHii)-freq(hightRiqPresHii,value=NA)))*100 ),
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
    spearman = c(corRiqPres[2,1],corRiqOtim[2,1],corRiqPess[2,1])
)

write.csv(tabSobreposicaoRiqueza,paste(projectFolder,'tabSobreposicaoRiqueza.csv',sep=''),row.names=FALSE)

##mapas de sobreposicao

jpeg(filename=paste(projectFolder,'mapasSobreposicaoRiqueza.jpeg'),width=1800,height=900)
par(mfrow=c(1,3),mar=c(5,5,5,14))
plot(hightRiqPres*1 + hightRiqPresHii*2,col=c('white','green','red','gray'),main='Presente',legend=FALSE,cex.axis=2,cex.main=4);plot(AmSulShape,add=TRUE); grid()
plot(hightRiqOtim*1 + hightRiqOtimHii*2,col=c('white','green','red','gray'),main='2070 otimista',legend=FALSE,cex.axis=2,cex.main=4);plot(AmSulShape,add=TRUE); grid()
plot(hightRiqPess*1 + hightRiqPessHii*2,col=c('white','green','red','gray'),main='2017 pessimista',legend=FALSE,cex.axis=2,cex.main=4);plot(AmSulShape,add=TRUE); grid()
legend('bottomright',legend=c('Variáveis climáticas','Var. clim. & impacto humano','Sobreposição'),pch=c(21,21,21),pt.bg=c('green','red','gray'),bty='n',cex=2.5)
dev.off()


###OITAVA PARTE: comparando areas de ALTO RISCO entre as projecoes com e sem impacto humano###

projectFolder = '/home/anderson/Documentos/Projetos/Barbeiros_Lucas/'

##presente
mapaRiscoPresente = raster(paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoPresente.asc',sep=''))
mapaRiscoPresenteHii = raster(paste(projectFolder,'resultados nicho climatico + impacto humano/Mapas de risco/mapaRiscoPresente.asc',sep=''))
#
hightRiscPres = mapaRiscoPresente > quantile(mapaRiscoPresente, 0.75,na.rm=TRUE) #raster quartil superior risco
hightRiscPresHii = mapaRiscoPresenteHii > quantile(mapaRiscoPresenteHii, 0.75,na.rm=TRUE) #raster quartil superior risco

##futuro otimista
mapaRiscoFuturoOtimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoFuturoOtimista.asc',sep=''))
mapaRiscoFuturoOtimistaHii = raster(paste(projectFolder,'resultados nicho climatico + impacto humano/Mapas de risco/mapaRiscoFuturoOtimista.asc',sep=''))
#
hightRiscOtim = mapaRiscoFuturoOtimista > quantile(mapaRiscoFuturoOtimista, 0.75,na.rm=TRUE) #raster quartil superior risco
hightRiscOtimHii = mapaRiscoFuturoOtimistaHii > quantile(mapaRiscoFuturoOtimistaHii, 0.75,na.rm=TRUE) #raster quartil superior risco

##futuro pessimista
mapaRiscoFuturoPessimista = raster(paste(projectFolder,'resultados nicho climatico/Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''))
mapaRiscoFuturoPessimistaHii = raster(paste(projectFolder,'resultados nicho climatico + impacto humano/Mapas de risco/mapaRiscoFuturoPessimista.asc',sep=''))
#
hightRiscPess = mapaRiscoFuturoPessimista > quantile(mapaRiscoFuturoPessimista, 0.75,na.rm=TRUE) #raster quartil superior risco
hightRiscPessHii = mapaRiscoFuturoPessimistaHii > quantile(mapaRiscoFuturoPessimistaHii, 0.75,na.rm=TRUE) #raster quartil superior risco

##unificando para contabilizar areas
hightRiscPresUni = hightRiscPres*1 + hightRiscPresHii*2
hightRiscOtimUni = hightRiscOtim*1 + hightRiscOtimHii*2
hightRiscPessUni = hightRiscPess*1 + hightRiscPessHii*2

##correlacao areas alta riqueza modelos com e sem HII
##presente
testRiscPres <- getValues(stack(hightRiscPres,hightRiscPresHii))
corRiscPres <- as.data.frame(cor(testRiscPres, use="complete.obs",method='spearman'))
##futuro otimista
testRiscOtim <- getValues(stack(hightRiscOtim,hightRiscOtimHii))
corRiscOtim <- as.data.frame(cor(testRiscOtim, use="complete.obs",method='spearman'))
##futuro pessimista
testRiscPess <- getValues(stack(hightRiscPess,hightRiscPessHii))
corRiscPess <- as.data.frame(cor(testRiscPess, use="complete.obs",method='spearman'))

##correlacao mapa risco modelos com e sem HII## OBS: ACHEI MELHOR USAR ESSE NO MANUSCRITO
##presente
testRiscPres <- getValues(stack(mapaRiscoPresente,mapaRiscoPresenteHii))
corRiscPres <- as.data.frame(cor(testRiqPres, use="complete.obs",method='pearson'))
##futuro otimista
testRiscOtim <- getValues(stack(mapaRiscoFuturoOtimista,mapaRiscoFuturoOtimistaHii))
corRiscOtim <- as.data.frame(cor(testRiqOtim, use="complete.obs",method='pearson'))
##futuro pessimista
testRiscPess <- getValues(stack(mapaRiscoFuturoPessimista,mapaRiscoFuturoPessimistaHii))
corRiscPess <- as.data.frame(cor(testRiqPess, use="complete.obs",method='pearson'))

tabSobreposicaoRisco = data.frame(
    cenario = c('presente','futuroOtimista','futuroPessimista'),
    total_celulas_nichoClim = c(freq(hightRiscPres,value=1),freq(hightRiscOtim,value=1),freq(hightRiscPess,value=1)),
    percentual_da_AmSul = c(
        (freq(hightRiscPres,value=1)/(ncell(hightRiscPres)-freq(hightRiscPres,value=NA)))*100,
        (freq(hightRiscOtim,value=1)/(ncell(hightRiscPres)-freq(hightRiscPres,value=NA)))*100,
        (freq(hightRiscPess,value=1)/(ncell(hightRiscPres)-freq(hightRiscPres,value=NA)))*100 ),
    total_celulas_HII= c(freq(hightRiscPresHii,value=1),freq(hightRiscOtimHii,value=1),freq(hightRiscPessHii,value=1)),
    percentual_da_AmSul = c(
        (freq(hightRiscPresHii,value=1)/(ncell(hightRiscPres)-freq(hightRiscPresHii,value=NA)))*100,
        (freq(hightRiscOtimHii,value=1)/(ncell(hightRiscPres)-freq(hightRiscPresHii,value=NA)))*100,
        (freq(hightRiscPessHii,value=1)/(ncell(hightRiscPres)-freq(hightRiscPresHii,value=NA)))*100 ),
    total_cels_sobrepostas = c(freq(hightRiscPresUni,value=3),freq(hightRiscOtimUni,value=3),freq(hightRiscPessUni,value=3)),
    percentual_para_nichoClim=
        c(( freq(hightRiscPresUni,value=3)/freq(hightRiscPres,value=1) ) * 100,
        ( freq(hightRiscOtimUni,value=3)/freq(hightRiscOtim,value=1) ) * 100,
        ( freq(hightRiscPessUni,value=3)/freq(hightRiscPess,value=1) ) * 100),
    percentual_para_HII=
        c(( freq(hightRiscPresUni,value=3)/freq(hightRiscPresHii,value=1) ) * 100,
        ( freq(hightRiscOtimUni,value=3)/freq(hightRiscOtimHii,value=1) ) * 100,
        ( freq(hightRiscPessUni,value=3)/freq(hightRiscPessHii,value=1) ) * 100),
    percentual_geral=
        c((freq(hightRiscPresUni,value=3)/(ncell(hightRiscPresUni)-freq(hightRiscPresUni,value=NA)) * 100),
        (freq(hightRiscOtimUni,value=3)/(ncell(hightRiscOtimUni)-freq(hightRiscOtimUni,value=NA)) * 100),
        ( freq(hightRiscPessUni,value=3)/(ncell(hightRiscPessUni)-freq(hightRiscPessUni,value=NA)) * 100)),
    spearman = c(corRiscPres[2,1],corRiscOtim[2,1],corRiscPess[2,1])
)

write.csv(tabSobreposicaoRisco,paste(projectFolder,'tabSobreposicaoRico.csv',sep=''),row.names=FALSE)

##mapas de sobreposicao

jpeg(filename=paste(projectFolder,'mapasSobreposicaoRisco.jpeg'),width=1800,height=900)
par(mfrow=c(1,3),mar=c(5,5,5,14))
plot(hightRiscPres*1 + hightRiscPresHii*2,col=c('white','green','red','gray'),main='Presente',legend=FALSE,cex.axis=2,cex.main=4);plot(AmSulShape,add=TRUE); grid()
plot(hightRiscOtim*1 + hightRiscOtimHii*2,col=c('white','green','red','gray'),main='2070 otimista',legend=FALSE,cex.axis=2,cex.main=4);plot(AmSulShape,add=TRUE); grid()
plot(hightRiscPess*1 + hightRiscPessHii*2,col=c('white','green','red','gray'),main='2017 pessimista',legend=FALSE,cex.axis=2,cex.main=4);plot(AmSulShape,add=TRUE); grid()
legend('bottomright',legend=c('Variáveis climáticas','Var. clim. & impacto humano','Sobreposição'),pch=c(21,21,21),pt.bg=c('green','red','gray'),bty='n',cex=2.5)
dev.off()

##graficos das tendencias: percentual de cobertura da america do sul
jpeg(filename=paste(projectFolder,'tendenciasCobertura.jpeg'))
plot(tabSobreposicaoRisco$percentual_da_AmSul~c(1,2,2),main='Percentual de cobertura da Am. do Sul',xlab='Cenário',ylab='Percentual',pch=c(15,16,17),ylim=c(0,50),xlim=c(0.5,2.5),cex=3,col=rgb(1,0,0,0.5));grid()
abline(lm(tabSobreposicaoRisco$percentual_da_AmSul~c(1,2,2)),col='red')
points(tabSobreposicaoRiqueza$percentual_da_AmSul~c(1,2,2),xlabel='',ylabel='',pch=c(15,16,17),ylim=c(0,50),xlim=c(0.5,2.5),cex=3,col=rgb(0,1,0,0.5));grid()
abline(lm(tabSobreposicaoRiqueza$percentual_da_AmSul~c(1,2,2)),col='green')
#
points(tabSobreposicaoRisco$percentual_da_AmSul.1~c(1,2,2),xlab='',ylab='',pch=c(15,16,17),ylim=c(0,50),xlim=c(0.5,2.5),cex=3,lwd=2,col=rgb(1,0,0,0.5));grid()
abline(lm(tabSobreposicaoRisco$percentual_da_AmSul.1~c(1,2,2)),col='red',lty=2,lwd=2)
points(tabSobreposicaoRiqueza$percentual_da_AmSul.1~c(1,2,2),xlab='',ylab='',pch=c(15,16,17),ylim=c(0,50),xlim=c(0.5,2.5),cex=3,lwd=2,col=rgb(0,1,0,0.5));grid()
abline(lm(tabSobreposicaoRiqueza$percentual_da_AmSul.1~c(1,2,2)),col='green',lty=2,lwd=2)
legend('topright',legend=c('Presente','2070 otimista','2070 pessimista','Apenas variáveis climáticas','Var. clim. & impacto humano','Riqueza','Risco'),pch=c(0,1,2,NA,NA,NA,NA),lty=c(NA,NA,NA,1,2,NA,NA),text.col=c('black','black','black','black','black','green','red'))
dev.off()

##graficos das tendencias: correlacao
jpeg(filename=paste(projectFolder,'tendenciasCorrelacao.jpeg'))
plot(tabSobreposicaoRisco$spearman~c(1,2,2),main='Percentual de cobertura da Am. do Sul',xlab='Cenário',ylab='Percentual',pch=c(15,16,17),ylim=c(0,1),xlim=c(0.5,2.5),cex=3,col=rgb(1,0,0,0.5));grid()
abline(lm(tabSobreposicaoRisco$spearman~c(1,2,2)),col='red')
points(tabSobreposicaoRiqueza$spearman~c(1,2,2),xlabel='',ylabel='',pch=c(15,16,17),ylim=c(0,50),xlim=c(0.5,2.5),cex=3,col=rgb(0,1,0,0.5));grid()
abline(lm(tabSobreposicaoRiqueza$spearman~c(1,2,2)),col='green')
legend('topright',legend=c('Presente','2070 otimista','2070 pessimista','Riqueza','Risco'),pch=c(0,1,2,NA,NA),text.col=c('black','black','black','green','red'))
dev.off()
