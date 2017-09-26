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
