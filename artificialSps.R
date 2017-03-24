#########################################################################################
####SCRIPT PARA FAZER A DISTRIBUICAO REAL DAS ESPECIES EM VARIOS MOMENTOS DO TEMPO####
#########################################################################################
library(virtualspecies)
library(maptools)
library(dismo)
library(raster)
library(phyloclim) #para funcao niche.overlap()
source("/home/anderson/R/R-Scripts/TSSmaxent.R")

###PRIMEIRA PARTE: criando sps virtuais###

###Parametros necessarios###
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
caminhosCamadasTemp = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
############################

for (i in 1:length(caminhosCamadasTemp)){
    predictors = stack(paste(caminhosCamadasTemp[i],'/bioclim_01.asc',sep=''),paste(caminhosCamadasTemp[i],'/bioclim_12.asc',sep='')) #carregando as variaveis ambientais
    predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
    nameScenario = basename(caminhosCamadasTemp[i])
    
    #Especie de clima quente e umido
    parametersHW <- formatFunctions(bioclim_01=c(fun='betaFun',p1=200,p2=295,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=2000,p2=3500,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais

    #Especie de clima quente e seco
    parametersHD <- formatFunctions(bioclim_01=c(fun='betaFun',p1=200,p2=260,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=50,p2=1800,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais
 
    #Especie de clima frio e seco
    parametersCD <- formatFunctions(bioclim_01=c(fun='betaFun',p1=50,p2=220,alpha=1,gamma=1),bioclim_12=c(fun='betaFun',p1=50,p2=1800,alpha=1,gamma=1)) #criando as respostas da especie às variaveis ambientais

    spHW <- generateSpFromFun(predictors, parametersHW) #criando a especie artifical (clima quente e umido)
    spHD <- generateSpFromFun(predictors, parametersHD) #criando a especie artifical (clima quente e seco)
    spCD <- generateSpFromFun(predictors, parametersCD) #criando a especie artifical (clima frio e seco)

    auxVector=stack(c(spHW$suitab.raster,spHD$suitab.raster,spCD$suitab.raster))
    names(auxVector) = c('spHW', 'spHD', 'spCD')
    
    for(j in 1:dim(auxVector)[3]){
        projection(auxVector[[j]]) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
        writeRaster(auxVector[[j]], filename=paste(projectFolder,'NichoReal/',names(auxVector[[j]]),'/',nameScenario,'.asc',sep=""), overwrite=TRUE,prj=TRUE) #salvando o raster do mapa da sp
    }
}

###SEGUNDA PARTE: amostragem de pontos de ocorrencia em diferentes camadas de tempo###

envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto/" #pasta com as variaveis ambientais
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha com os pontos amostrados sera salva
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
biomodFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/biomod/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
sampleSizes = c(5,15,25,35,45,55,65,75,85,95) #tamanhos das amostras
NumRep = 10 #numero de replicas (de cada cenario amostral)
Tmax = 22 #idade maxima (no passado)
bgPoints = 1000 #numero de pontos de background
sampleData = data.frame()


#ocorrencias
for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
    for (sSize in sampleSizes){ #numero de pontos (registros, dados) na amostra
        sampledAges = vector()
        sampledAges = round(runif(sSize,0,Tmax)) #selecionando 'n' camadas de tempo aleatoriamente
        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #lista com os enderecos dos mapas de distribuicao da
        for (j in 1:NumRep){ #replicas do cenario amostral
            for (sAge in sampledAges){ #amostrando em cada camada de tempo que consta na amostra
                sampleData_i = sampleOccurrences(x=raster(nicheRealPath[sAge+1]),n=1,plot=FALSE)$sample.points[,1:2] #amostra d ponto
                scenarioName = basename(nicheRealPath[1:24][sAge+1]) #tempo vinculado ao cenario para variaveis ambientais
                scenarioName = gsub('.asc','',scenarioName) #retirando do nome o '.asc'
                layers_i = extract(
                    x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                    y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
                sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
            }
            names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
            write.csv(sampleData,paste(projectFolder,'Amostras/',spsTypes[i],'/occ',sSize,'pts', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
            sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
        }
    }
}

#background points
for (i in 1:length(spsTypes)){ #loop sobre os 'tipos de especies'
    for (sSize in sampleSizes){ #numero de pontos de ocoorrencia (que parearah com a amostra de background dessa  iteracao)
        sampledAges = vector()
        sampledAges = round(runif(bgPoints,0,Tmax)) #selecionando 'n' camadas de tempo aleatoriamente ('n' = bgPoints)
#        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        for (j in 1:NumRep){ #replicas do cenario amostral
            for (sAge in sampledAges){ #amostrando em cada camada de tempo que consta na amostra

                envVarPath = list.files(path=paste(envVarFolder,list.files(path=paste(envVarFolder))[sAge+1],sep=''), full.names=TRUE, pattern='.asc') #lista com os enderecos das variaveis ambientais no tempo corresposndente a interacao 
                
                sampleData_i = randomPoints(mask=raster(envVarPath[1],crs=CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),n=1) #amostra do ponto
                scenarioName = list.files(path=paste(envVarFolder))[sAge+1] #nome do cenario
                layers_i = extract(
                    x=stack(list.files(path=paste(envVarFolder,scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                    y=sampleData_i) #extraindo variaveis ambientais do ponto, em sua respectiva camada de tempo
                sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #juntando com os dados das outras camadas de tempo amostradas
            }
            names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #ajustando os nomes
            write.csv(sampleData,paste(projectFolder,'Amostras/',spsTypes[i],'/bg',sSize,'pts', j ,'rep.csv',sep=''),row.names=FALSE) #salvando
            sampleData = data.frame() #devolvendo data.frame vazio para proxima rodada
        }
    }
}




###TERCEIRA PARTE: SDM usando de pontos de ocorrencia em diferentes camadas de tempo (do atual ate 120 kyr BP)###

#######################################################
####################### MAXENT ########################
#######################################################

options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=T) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Amostras/' #caminho para pasta onde a planilha 
maxentFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/' #pasta para resultados do maxent
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
source("/home/anderson/R/R-Scripts/TSSmaxent.R")
evaluation = list()
TSSvector = data.frame()
sampleSizes = c(5,15,25,35,45,55,65,75,85,95)

for (i in 1:length(spsTypes)){
    for (j in sampleSizes){
        for (k in 1:NumRep){ #loop sobre o numero de replicas 
        
            occPoints = read.csv(paste(mainSampleFolder,spsTypes[i],'/occ',j,'pts',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia
            backgroundPoints = read.csv(paste(mainSampleFolder,spsTypes[i],'/bg',j,'pts',k,'rep.csv',sep=''),header=TRUE) #abrindo pontos de background
            names(backgroundPoints) = names(occPoints) #certificando que os nomes das colunas estão iguais (cuidado aqui...)
            dataSet = data.frame(cbind(rbind(occPoints,backgroundPoints),pres=c(rep(1,nrow(occPoints)),rep(0,nrow(backgroundPoints))))) #planilha de dados no formato SWD
            
            me = maxent(
                ##x=dataSet[,c("bioclim_01","bioclim_04","bioclim_10","bioclim_11","bioclim_12","bioclim_15","bioclim_16","bioclim_17")],
                x=dataSet[,c("bioclim_01","bioclim_12")],
                p=dataSet$pres,
                path=paste(maxentFolder,spsTypes[i],sep=''),
                args=c('responsecurves=TRUE',
                       'jackknife=TRUE',
                       'randomseed=TRUE',
                       'randomtestpoints=25',
                       'maximumbackground=5000',
                       'replicates=10',
                       'replicatetype=subsample',
                       'writebackgroundpredictions=TRUE',
                       'linear=TRUE',
                       'quadratic=TRUE',
                       'product=FALSE',
                       'threshold=FALSE',
                       'hinge=FALSE',
                       'maximumiterations=1000',
                       'convergencethreshold=1.0E-5',
                       'threads=2'
                       ))
            
            ##rodando a avaliacao do modelo
            TSSvector = rbind(TSSvector, TSSmaxent(paste(maxentFolder,spsTypes[i],'/',sep='')))
            threshold = as.data.frame(read.csv(paste(maxentFolder,spsTypes[i],'/maxentResults.csv',sep=''),header=TRUE))$X10.percentile.training.presence.logistic.threshold[11] #threshold 10 percentile training occ
            output = rbind(sp=spsTypes[i],sampleSize=j,replicate=k,AUC=TSSvector$testAUC,TSS=TSSvector$TSS,threshold=threshold,numbrTimeSlices=length(unique(occPoints$kyrBP)),medianSampledAges=median(unique(occPoints$kyrBP)),smallerAgeSampled=min(unique(occPoints$kyrBP)),largerAgeSampled=max(unique(occPoints$kyrBP)))
            write.csv(output,file=paste(projectFolder,'Maxent/',spsTypes[i],'/StatisticsResults-',spsTypes[i],'.csv',sep=''),row.names=FALSE)


            ##CONTINAR VERIFICANDO E ARRUMANDO DAQUI##
            ##acrescentar maxent com ausencias reais mesmo??
            

            
            for (l in 1:length(envVarPaths[1:24])){
                predictors = stack(list.files(path=envVarPaths[l],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
                predictors = predictors[[c('bioclim_01','bioclim_12')]]
                predictors = mask(predictors,AmSulShape) #recortando as variaveis ambientais
                crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
                proj = predict(me,predictors,crs=crs) #realizando projetacoes (para cada replica)
                writeRaster(mean(proj),paste(maxentFolder,spsTypes[i],'/projections/projection-Time',l-1,'kyrBP','-Replica',k,'-Sample',j,'.asc',sep=''),overwrite=TRUE) #salvando a projecao media
            }
        }

        ##criando um mapa binario
        ## thresholdValues = NULL
        ## aucValues = NULL
        ## for (m in 1:length(evaluation)){
        ##     thresholdValues <- append(thresholdValues, threshold(evaluation[[m]],'spec_sens'))
        ##     aucValues = append(aucValues, evaluation[[m]]@auc)
        ## }
        ## aucMean = mean(aucValues)
        ## thresholdMean = mean(thresholdValues)
        ## TSSmean = mean(TSSvector$TSS)
        ## statResults = data.frame(sp=spsTypes[i],AUCmean=aucMean,TSSmean=TSSmean,ThresholdMean=thresholdMean)
        ## write.csv(statResults,file=paste(projectFolder,'Maxent/',spsTypes[i],'/StatisticsResults-',spsTypes[i],'.csv',sep=''),row.names=FALSE)
        
        ##esvaziando o vetor para a proxima especie
        evaluation = list()
        TSSvector = data.frame()
    }
}
#######################################################
#######################################################
#######################################################

###QUARTA PARTE: comparando projecao do SDM e a distribuicao espacial real do nicho da sp###

projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto
spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
scenarioModel = c('8varLinearModel','8varQuadModel','2varLinearModel','2varQuadModel')
sampleSizes = c(5,15,25,35,45,55,65,75,85,95) #aqui, deve ser igual ao usasado nas partes anteriores do script

model1 = pres ~ bioclim_01 + bioclim_04 + bioclim_10 + bioclim_11 + bioclim_12 + bioclim_15 + bioclim_16 + bioclim_17 
model2 = pres ~ bioclim_01 + I(bioclim_01^2) + bioclim_04 + I(bioclim_04^2) + bioclim_10 + I(bioclim_10^2) + bioclim_11 + I(bioclim_11^2) + bioclim_12 + I(bioclim_12^2) + bioclim_15 + I(bioclim_15^2) + bioclim_16 + I(bioclim_16^2) + bioclim_17 + I(bioclim_17^2)
model3 = pres ~ bioclim_01 + bioclim_12
model4 = pres ~ bioclim_01 + I(bioclim_01^2) + bioclim_12 + I(bioclim_12^2)
model = c(model1,model2,model3,model4)

for (i in 1:length(spsTypes)){
    
    nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
    nicheRealPath = list.files(path=nicheRealFolder,pattern='asc',full.names=TRUE) #lista com os enderecos dos mapas de distribuicao da sp
    projectionsFolder = paste(projectFolder,'Maxent/',spsTypes[i],'/projections',sep='') #pasta com as projecoes do cenario
    projectionsPath = list.files(path=projectionsFolder, pattern='asc',full.names=T) #caminhos para os .asc na paste do cenario
    outputData = data.frame(kyrBP=numeric(),sampleSize=numeric(),Schoeners_D=numeric(),Warren_I=numeric())
    
    for (l in 1:length(nicheRealPath[1:24])){ #loop sobre as cadamdas de tempo

        realNiche = raster(nicheRealPath[l]) #nicho real
        realNiche.spgrid = as(realNiche,'SpatialGridDataFrame')
        
        for (m in sampleSizes){ #loop sobre os tamanhos amostrais

            timeSampleData =  stack(list.files(path=projectionsFolder, pattern=glob2rx(paste('*Time',l-1,'*Sample',m,'.asc',sep='')),full.names=TRUE))
            output_i = data.frame(kyrBP=numeric(),sapleSize=numeric(),Schoeners_D=numeric(),Warren_I=numeric()) #dataframe vazio para o loop abaixo
            
            for(n in 1:5){ #loop sobre replicas de cada combinacao de tempo e tamanho amostral
              
                sdmNiche = timeSampleData[[n]] #mapa de suitability gerado por SDM
                sdmNiche.spgrid = as(sdmNiche,'SpatialGridDataFrame')
                nicheOverlap = niche.overlap(c(realNiche.spgrid,sdmNiche.spgrid))
                output_i= rbind(output_i,cbind(kyrBP=l-1,sampleSize=m,Schoeners_D=nicheOverlap[1,2],Warren_I=nicheOverlap[2,1]))
                
            }
        
            outputMeans = colMeans(output_i) #media das iteracoes (para a camada de tempo atual)
            outputData = data.frame(rbind(outputData,outputMeans))
            
        }
    }
    
    names(outputData) = c('kyrBP','sampleSize','Schoeners_D','Warren_I')  
    write.csv(outputData, file=paste(projectionsFolder,"/NO.csv",sep=""),row.names=FALSE) #salvando os dados do cenario
    
}


### QUINTA PARTE: construindo graficos dos resultados ###

spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
###abrindo as planilhas de dados
outputData = list()
vetor.nomes = vector()
projectFolder = "/home/anderson/Documentos/Projetos/Sps artificiais/" #pasta do projeto

for (i in 1:length(spsTypes)){
    outputData[[i]] = read.csv(file=paste(projectFolder,'Maxent/',spsTypes[i],'/projections/NO.csv',sep=''),header=TRUE)
    vetor.nomes = append(vetor.nomes,paste(spsTypes[i],sep=''))
}
names(outputData) = vetor.nomes

###graficos

##Boxplots 

##HW
HWdataD = data.frame(kyrBP=outputData$spHW$kyrBP,indexD=outputData$spHW$Schoeners_D)
HWdataH = data.frame(kyrBP=outputData$spHW$kyrBP,indexH=outputData$spHW$Hellinger_distances)

##HD
HDdataD = data.frame(kyrBP=outputData$spHW$kyrBP,indexD=outputData$spHD$Schoeners_D)
HDdataH = data.frame(kyrBP=outputData$spHW$kyrBP,indexH=outputData$spHD$Hellinger_distances)

##CD
CDdataD = data.frame(kyrBP=outputData$spHW$kyrBP,indexD=outputData$spCD$Schoeners_D)
CDdataH = data.frame(kyrBP=outputData$spHW$kyrBP,indexH=outputData$spCD$Hellinger_distances)

jpeg(file='/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/graficos/boxplots.jpg',width=900,height=600)
par(mfrow=c(1,2))
#D
Ddata = rbind(data.frame(kyrBP=HWdataD$kyrBP,indexD=HWdataD$indexD,spsType='Hot and Wet'),data.frame(kyrBP=HDdataD$kyrBP,indexD=HDdataD$indexD,spsType='Hot and Dry'),data.frame(kyrBP=CDdataD$kyrBP,indexD=CDdataD$indexD,spsType='Cold and Dry'))
boxplot(indexD~spsType,data=Ddata,ylim=c(0.2,1),names=c('H&W','H&D','C&D'),main="Schoeners' D",cex.main=2,cex.lab=1.5,cex.axis=1.5,lwd=3)
stripchart(indexD~spsType,data=Ddata,vertical=TRUE,method="jitter",pch=20,cex=1.5,col=rgb(0.5,0.5,0.5,0.2),add=TRUE) 
#H
Hdata = rbind(data.frame(kyrBP=HWdataH$kyrBP,indexH=HWdataH$indexH,spsType='Hot and Wet'),data.frame(kyrBP=HDdataH$kyrBP,indexH=HDdataH$indexH,spsType='Hot and Dry'),data.frame(kyrBP=CDdataH$kyrBP,indexH=CDdataH$indexH,spsType='Cold and Dry'))
boxplot(indexH~spsType,data=Hdata,ylim=c(0.2,1),names=c('H&W','H&D','C&D'),main="Warren's I",cex.main=2,cex.lab=1.5,cex.axis=1.5,lwd=3)
stripchart(indexH~spsType,data=Hdata,vertical=TRUE,method="jitter",pch=20,cex=1.5,col=rgb(0.5,0.5,0.5,0.2),add=TRUE) 
dev.off()

##Schoeners' D ao longo do tempo##

jpeg(file='/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/graficos/metricsXtime.jpg',width=1100,height=750)
par(mfrow=c(1,2))
#D
plot(indexD~kyrBP,data=HWdataD,type='p',ylim=c(0,1),pch=1,col='black',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(indexD~kyrBP,data=HDdataD,type='p',ylim=c(0,1),pch=2,col='blue',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(indexD~kyrBP,data=CDdataD,type='p',ylim=c(0,1),pch=3,col='red',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
legend(x=19,y=1,legend=c('H&W','H&D','C&D'),pch=20,col=c('black','blue','red'),cex=1.5)
#H
plot(indexH~kyrBP,data=HWdataH,type='p',ylim=c(0,1),pch=1,col='black',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(indexH~kyrBP,data=HDdataH,type='p',ylim=c(0,1),pch=2,col='blue',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(indexH~kyrBP,data=CDdataH,type='p',ylim=c(0,1),pch=3,col='red',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
legend(x=19,y=1,legend=c('H&W','H&D','C&D'),pch=20,col=c('black','blue','red'),cex=1.5)
dev.off()

##metricas X tamanho amostral

#maxent
jpeg(file='/home/anderson/Documentos/Projetos/Sps artificiais/Maxent/graficos/metricsXsampleSize.jpg',width=1100,height=750)
par(mfrow=c(1,2),oma=c(1,2,1,1))
#D
plot(outputData$spHW$Schoeners_D~outputData$spHW$sampleSize,type='p',ylim=c(0,1),xlab='Sample size',ylab="Shoeners' D",pch=1,col='black',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(outputData$spHD$Schoeners_D~outputData$spHD$sampleSize,data=HDdataD,type='p',ylim=c(0,1),pch=2,col='blue',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(outputData$spCD$Schoeners_D~outputData$spCD$sampleSize,data=CDdataD,type='p',ylim=c(0,1),pch=3,col='red',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
legend(x=6,y=1,legend=c('H&W','H&D','C&D'),pch=20,col=c('black','blue','red'),cex=1.5)
#H
plot(outputData$spHW$Hellinger_distances~outputData$spHW$sampleSize,type='p',ylim=c(0,1),xlab='Sample size',ylab="Warren's I",pch=1,col='black',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(outputData$spHD$Hellinger_distances~outputData$spHD$sampleSize,data=HDdataD,type='p',ylim=c(0,1),pch=2,col='blue',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
points(outputData$spCD$Hellinger_distances~outputData$spCD$sampleSize,data=CDdataD,type='p',ylim=c(0,1),pch=3,col='red',cex.lab=1.9,cex.axis=1.5,lwd=2,cex=2)
legend(x=6,y=1,legend=c('H&W','H&D','C&D'),pch=20,col=c('black','blue','red'),cex=1.5)
dev.off()


##distribuicao presente, inter e maximo glacial

threHW = read.csv(paste(projectFolder,'Maxent/spHW/StatisticsResults-spHW.csv',sep=''),header=TRUE)$ThresholdMean
HWcurrentReal = raster(paste(projectFolder,'NichoReal/spHW/000.asc',sep='')) > 0.1
HW22Real = raster(paste(projectFolder,'NichoReal/spHW/022.asc',sep='')) > 0.1
HWModel_0kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threHW
HWModel_0kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threHW
HWModel_0kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threHW
HWModel_22kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threHW
HWModel_22kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threHW
HWModel_22kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHW/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threHW

threHD = read.csv(paste(projectFolder,'Maxent/spHD/StatisticsResults-spHD.csv',sep=''),header=TRUE)$ThresholdMean
HDcurrentReal = raster(paste(projectFolder,'NichoReal/spHD/000.asc',sep='')) > 0.1
HD22Real = raster(paste(projectFolder,'NichoReal/spHD/022.asc',sep='')) > 0.1
HDModel_0kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threHD
HDModel_0kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threHD
HDModel_0kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threHD
HDModel_22kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threHD
HDModel_22kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threHD
HDModel_22kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spHD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threHD

threCD = read.csv(paste(projectFolder,'Maxent/spCD/StatisticsResults-spCD.csv',sep=''),header=TRUE)$ThresholdMean
CDcurrentReal = raster(paste(projectFolder,'NichoReal/spCD/000.asc',sep='')) > 0.1
CD22Real = raster(paste(projectFolder,'NichoReal/spCD/022.asc',sep='')) > 0.1
CDModel_0kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threCD
CDModel_0kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threCD
CDModel_0kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time0*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threCD
CDModel_22kyrSample5 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',5,'.asc',sep='')),full.names=TRUE)) ) > threCD
CDModel_22kyrSample45 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',45,'.asc',sep='')),full.names=TRUE)) ) > threCD
CDModel_22kyrSample95 =  mean( stack(list.files(path=paste(projectFolder,'Maxent/spCD/projections/',sep=''), pattern=glob2rx(paste('*Time22*Sample',95,'.asc',sep='')),full.names=TRUE)) ) > threCD

library(rasterVis)
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp")

##sobreposicoes spHW
jpeg(filename=paste(projectFolder,'Maxent/graficos/sobreposicoesHW.jpg',sep=''), width = 1100 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,0))
plot(HWcurrentReal*1+HWModel_0kyrSample5*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample45*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample95*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample5*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample45*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample95*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
grid()
mtext('Hot & Wet sp.',outer=TRUE,cex=4)
dev.off()

##sobreposicoes spHD
jpeg(filename=paste(projectFolder,'Maxent/graficos/sobreposicoesHD.jpg',sep=''), width = 1100 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,0))
plot(HDcurrentReal*1+HDModel_0kyrSample5*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(HDcurrentReal*1+HDModel_0kyrSample45*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(HDcurrentReal*1+HDModel_0kyrSample95*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
grid()
plot(HD22Real*1+HDModel_22kyrSample5*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(HD22Real*1+HDModel_22kyrSample45*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(HD22Real*1+HDModel_22kyrSample95*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
grid()
mtext('Hot & Dry sp.',outer=TRUE,cex=4)
dev.off()

##sobreposicoes spCD
jpeg(filename=paste(projectFolder,'Maxent/graficos/sobreposicoesCD.jpg',sep=''), width = 1100 , height = 1100) 
par(mfrow=c(2,3)) 
plot(CDcurrentReal*1+CDModel_0kyrSample5*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample45*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample95*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample5*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 5 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample45*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 45 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample95*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 95 records',cex=2)
legend("bottomright", inset=c(-1,0),legend=c('Virtual species','Maxent projection','Overlap'),pch=20,col=c('blue','green','dark green'),cex=2)
grid()
mtext('Cold & Dry sp.',outer=TRUE,cex=4)
dev.off()
