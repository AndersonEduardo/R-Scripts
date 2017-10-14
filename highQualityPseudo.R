## Este script contem algoritmos para gerar 'n' sps artificiais aleatoriamente (Parte 1), modelar ausencias (primeiro SDM) e modelar distribuicao com pseudoausencias melhoradas
## Anderson A. Eduardo, outubro/2014


##PARTE 1: criando as especies artificiais

##abrindo pacotes necessarios
library(raster)

##definindo prametros e variaveis globais
projectFolder = "/home/anderson/Documentos/Projetos/Improved pseudo-absences" #pasta do projeto

predictors = stack(paste(caminhosCamadasTemp[1], '/bioclim_01.asc', sep=''), paste(caminhosCamadasTemp[1], '/bioclim_12.asc', sep='')) #carregando as variaveis ambientais
predictors = mask(predictors, AmSulShape) #recortando as variaveis ambientais

datMat = as.data.frame(predictors, xy=TRUE, na.rm=TRUE) #transformando raster em data.frame
names(datMat) = c('lon', 'lat', 'bio1', 'bio12') #ajustando os nomes das colunas do data.frame

Nsp = 100 #numero de especies a serem criadas

##equacoes para as dimensoes do nicho das especies
betaBio1 = runif(n=Nsp, min=0.01, max=2) #parametro para cada equacao de cada especie
betaBio12 = runif(n=Nsp, min=0.01, max=2) #parametro para cada equacao de cada especie
alphaBio1 = runif(n=Nsp, min=quantile(x=varBio1,probs=0.1,na.rm=TRUE),max=quantile(x=varBio1,probs=0.9,na.rm=TRUE)) #parametro para cada equacao de cada especie
alphaBio12 = runif(n=Nsp, min=quantile(x=varBio12,probs=0.1,na.rm=TRUE), max=quantile(x=varBio12,probs=0.9,na.rm=TRUE)) #parametro para cada equacao de cada especie

## betaBio1 = abs(rnorm(n=Nsp,mean=0.1,sd=0.1)) #parametro para cada equacao de cada especie
## betaBio12 = abs(rnorm(n=Nsp,mean=0.001,sd=0.1)) #parametro para cada equacao de cada especie
## alphaBio1 = abs(rnorm(n=Nsp,mean=quantile(x=varBio1,probs=0.5,na.rm=TRUE))) #parametro para cada equacao de cada especie
## alphaBio12 = abs(rnorm(n=Nsp,mean=quantile(x=varBio12,probs=0.5,na.rm=TRUE))) #parametro para cada equacao de cada especie
varBio1 = datMat$bio1 #variavel ambiental bioclim01
varBio12 = datMat$bio12 #variavel ambiental bioclim12

##solucao numerica para a equacao de cada especie
for (i in 1:Nsp){
    ##equacoes do nicho da especie
    fBio1Sp_i = as.integer( 1/(1+exp(-betaBio1[i]*(varBio1-alphaBio1[i]))) > 0.1 ) #solucao da equacao com output binario ("suitability")
    fBio12Sp_i = as.integer( 1/(1+exp(-betaBio12[i]*(varBio12-alphaBio12[i]))) > 0.1 ) #solucao da equacao com output binario ("suitability")
    #fBio1Sp_i = 1/(1+exp(-betaBio1[i]*(varBio1-alphaBio1[i]))) #solucao da equacao com output continuo ("suitability")
    #fBio12Sp_i = 1/(1+exp(-betaBio12[i]*(varBio12-alphaBio12[i]))) #solucao da equacao com output continuo ("suitability")
    datMat = data.frame(cbind(datMat,fSp=fBio1Sp_i*fBio12Sp_i)) #adicionando ao data.frame
    names(datMat)[ncol(datMat)] = paste('sp',i,sep='') #ajustando os nomes das especies no data.farme
    ##salvando graficos das equacoes de cada especie
#    jpeg(filename=paste('/home/anderson/Documentos/Projetos/divSpsSid/','functions_sp',i,'.jpeg',sep=''))
#    par(mfrow=c(1,2))
#    plot(fBio1Sp_i~varBio1,xlab='Bioclim 01',ylab='Suitability',ylim=c(0,1))
#    plot(fBio12Sp_i~varBio12,xlab='Bioclim 12',ylab='Suitability',ylim=c(0,1))
#    dev.off()
}

##criando a coluna de riqueza de especies##
datMat = cbind(datMat, Richness=rowSums(datMat[,grep('sp',names(datMat))]))

##raster da distribuicao modelada
for(i in 1:Nsp){
    SpDist = datMat[,c('lon','lat',paste('sp',i,sep=''))] #extraindo lon/lat e suitability (ou pres-aus) de cada especie
    coordinates(SpDist) = ~lon+lat #definindo colunas das coordenadas
    gridded(SpDist) = TRUE #definindo gridded
    proj4string(SpDist) = '+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84' #definindo proj
    rasterSpDist = raster(SpDist) #criando objeto raster
    ##criando imagem da distribuicao de cada especie
    jpeg(filename=paste(projectFolder,'/sp',i,'.jpeg',sep=''))
    plot(rasterSpDist)
    dev.off()
    writeRaster(x=rasterSpDist,filename=paste(projectFolder,'/sp',i,'.jpeg',sep=''))
}





##########################################
#CONTINUAR DAQUI: ajustar o codigo para o arranjo das virtual species implementadas na parte 1, acima
##########################################




##PARTE 2: modelando ausencias e SDM com pseudoausencias melhoradas


##pacotes necessarios
library(biomod2)

##definindo variaveis e parametros
options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder = "/home/anderson/Documentos/Projetos/Improved pseudo-absences" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = "/home/anderson/Documentos/Projetos/Improved pseudo-absences/Amostras" #caminho para pasta onde a planilha 
maxentFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/maxent' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
sdmTypes = c('normal','optimized')
source("/home/anderson/R/R-Scripts/TSSmaxent.R")
sampleSizes = c(10,20,40,80)
NumRep = 10 #numero de replicas (de cada cenario amostral)
##variaveis preditoras
predictors = stack(list.files(path=envVarPaths[1],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
predictors = predictors[[c('bioclim_01','bioclim_12')]]
crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS

##algoritmo da analise do projeto
for(h in 1:length(sdmTypes)){
    for(i in 1:length(spsTypes)){
        for(j in 1:length(sampleSizes)){

            ##definindo variaveis e parametros locais
            statResultsSDMnormal = data.frame() #tabela de estatisticas basicas do modelo
            occPoints = read.csv(paste(mainSampleFolder,sdmTypes[h],'/',spsTypes[i],'/occ',sampleSizes[j],'.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia

            myResp <- occPoints[,c("lon","lat")]
            coordinates(myResp) <- ~ lon + lat #transformando em spatialPoints
            crs(myResp) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #transformando em spatialPoints

            ##variaveis e parametros locais especificos para o biomod2
            myRespName <- paste(spsTypes[i],'_sample',sampleSizes[j],sep='') # nome do cenario atual (para biomod2)
#            myResp <- dataSet[,c('pres')] # variavel resposta (para biomod2)
#            myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
            myExpl = predictors  #variavel preditora (para biomod2)

            ##ajuste de dados de entrada para biomod2
            myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                                 expl.var = myExpl,
                                                 resp.name = paste(myRespName,'_absencesModel',sep=''),
                                                 PA.nb.rep = 1
                                                 )

            ## ##inspecionando o objeto gerado pela funcao do biomod2
            ## myBiomodData
            ## plot(myBiomodData)
            
            ##parametrizando os modelos
            myBiomodOption <- BIOMOD_ModelingOptions(
                MAXENT.Phillips=list(
                    path_to_maxent.jar="/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java",
                    maximumiterations=1000,
                    linear=TRUE,
                    quadratic=TRUE,
                    product=FALSE,
                    threshold=FALSE,
                    hinge=FALSE,
                    maximumiterations=1000,
                    convergencethreshold=1.0E-5,
                    threads=2))

            ##rodando o(s) algoritmo(s) (i.e. SDMs)
            myBiomodModelOut <- BIOMOD_Modeling(
                myBiomodData,
                models = c('MAXENT.Phillips'),
                models.options = myBiomodOption,
                NbRunEval = 10,
                DataSplit = 75,
                VarImport = 5,
                models.eval.meth = c('TSS','ROC'),
                SaveObj = TRUE,
                rescal.all.models = TRUE,
                do.full.models = FALSE,
                modeling.id = paste(sdmTypes[h], sep=''))
            
            ##My output data
            evaluationScores = get_evaluations(myBiomodModelOut)

            ##gravando estatistcas basicas do modelo
            statResultsSDMnormal = rbind(statResultsSDMnormal,cbind(
                                                         sdm = sdmTypes[h],
                                                         sp = spsTypes[i],
                                                         sampleSize = sampleSizes[j],
                                                         AUC = mean(evaluationScores['ROC','Testing.data',,,]),
                                                         TSS = mean(evaluationScores['TSS','Testing.data',,,])
                                                     )
                                         )
            
            write.csv(statResultsSDMnormal,file=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/StatisticalResults-',spsTypes[i],'.csv',sep=''),row.names=FALSE)

            ##rodando algortmo de projecao (i.e. rodando a projecao)
            myBiomodProj <- BIOMOD_Projection(
                modeling.output = myBiomodModelOut,
                new.env = predictors,
                proj.name = paste(sdmTypes[h],'_',spsTypes[i],'_',sampleSizes[j],'pts',sep=''),
                selected.models = 'all',
                compress = 'FALSE',
                build.clamping.mask = 'FALSE',
                output.format = '.grd')
            
            ##gerando e salvando um mapa binario (threshold 10%)
            projStack = get_predictions(myBiomodProj) #extrai as projecoes
            projStackBIN = BinaryTransformation(projStack,10)
            
            ##writeRaster(projStackBIN,file=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[j],'.replica',k,'/proj_',l,'kyr/proj_',i,'kyr','.sample',sampleSizes[j],'.replica',k,'_BIN.asc',sep=''),row.names=FALSE)

            ##projStackBIN = projStack>0.5  #BinaryTransformation(projStack,"10")
            
            for (k in 1:nlayers(projStackBIN)){

                ##definindo variaveis e parametros locais
                betterPseudo = list()
                betterPseudoVar = list()
                statResults = data.frame() 

                ## >> AMOSTRANDO PSEUDOAUSENCIAS MELHORADAS <<
                betterPseudoPoints = dismo::randomPoints(mask=projStackBIN[[k]], n=1000) #sorteando pontos da distribuicao modelada
                betterPseudoDF = extract(projStackBIN[[k]], betterPseudoPoints) #distinguindo entre occ e ausencia
                betterPseudo[[k]] =  data.frame(lon=betterPseudoPoints[,1], lat=betterPseudoPoints[,2], occ=betterPseudoDF) #data.frame
                betterPseudo[[k]] = betterPseudo[[k]][which(betterPseudo[[k]]$occ==0),] #excluindo as presencas
                betterPseudoVar[[k]] = extract(predictors,betterPseudo[[k]][,c('lon','lat')]) #obtendo as variaveis preditoras nos pontos
                betterPseudo[[k]] = data.frame(betterPseudo[[k]],betterPseudoVar[[k]]) #motando dataset
                
                ## plot(projStackBIN[[k]])
                ## points(betterPseudo[[k]][,c('lon','lat')])

                ##definindo variaveis e parametros locais para o biomod2 (que rodara a seguir)
                backgroundPoints = betterPseudo[[k]]
                occPoints = data.frame(occPoints[,c('lon','lat')], occ=1, occPoints[,names(predictors)]) #assumindo que o stack predictors contem apenas as variaveis empregadas no projeto atual
                
                ##agrupando ocorrencias e pseudo-ausencias melhoradas
                dataSet = data.frame(rbind(occPoints,backgroundPoints)) #planilha de dados no formato SWD
                ##variaveis e parametros locais especificos para o biomod2
                myRespName <- paste(spsTypes[i],'_sample',sampleSizes[j],'_replica',k,sep='') # nome do cenario atual (para biomod2)
                myResp <- dataSet[,c('occ')] # variavel resposta (para biomod2)
                myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
                myExpl = dataSet[,c('bioclim_01','bioclim_12')]  #variavel preditora (para biomod2)
                
                ##ajuste de dados de entrada para biomod2
                myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                                     expl.var = myExpl,
                                                     resp.xy = myRespXY,
                                                     resp.name = myRespName)
                
                ## ##inspecionando o objeto gerado pela funcao do biomod2
                ## myBiomodData
                ## plot(myBiomodData)
                
                ##parametrizando os modelos
                myBiomodOption <- BIOMOD_ModelingOptions(
                    MAXENT.Phillips=list(
                        path_to_maxent.jar="/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java",
                        maximumiterations=1000,
                        linear=TRUE,
                        quadratic=TRUE,
                        product=FALSE,
                        threshold=FALSE,
                        hinge=FALSE,
                        maximumiterations=1000,
                        convergencethreshold=1.0E-5,
                        threads=2))

                ##rodando o(s) algoritmo(s) (i.e. SDMs)
                myBiomodModelOut <- BIOMOD_Modeling(
                    myBiomodData,
                    models = c('MAXENT.Phillips'),
                    models.options = myBiomodOption,
                    NbRunEval = 3,
                    DataSplit = 75,
                    VarImport = 5,
                    models.eval.meth = c('TSS','ROC'),
                    SaveObj = TRUE,
                    rescal.all.models = TRUE,
                    do.full.models = FALSE,
                    modeling.id = paste(myRespName))
                
                ##My output data
                evaluationScores = get_evaluations(myBiomodModelOut)
                
                ##gravando estatistcas basicas do modelo
                statResults = rbind(statResults,cbind(
                                                    sp = spsTypes[i],
                                                    sampleSize = sampleSizes[j],
                                                    AUC = mean(evaluationScores['ROC','Testing.data',,,]),
                                                    TSS = mean(evaluationScores['TSS','Testing.data',,,])
                                                )
                                    )

#write.csv(statResults,file=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/StatisticalResults-',spsTypes[i],'.csv',sep=''),row.names=FALSE)

                ##implementando projecoes do modelo
                
                ##rodando algortmo de projecao (i.e. rodando a projecao)
                myBiomodProj <- BIOMOD_Projection(
                    modeling.output = myBiomodModelOut,
                    new.env = predictors,
                    proj.name = paste(myRespName),
                    selected.models = 'all',
                    compress = 'FALSE',
                    build.clamping.mask = 'FALSE',
                    output.format = '.grd')
                
                ##gerando e salvando um mapa binario (threshold 10%)
                projStack = get_predictions(myBiomodProj) #extrai as projecoes
                projStackBIN = BinaryTransformation(stack(mean(projStack)),'10')
                
#writeRaster(projStackBIN,file=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[j],'.replica',k,'/proj_',l,'kyr/proj_',i,'kyr','.sample',sampleSizes[j],'.replica',k,'_BIN.asc',sep=''),row.names=FALSE)

            }
        }
    }
}
