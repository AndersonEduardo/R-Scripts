## Este script contem algoritmos para gerar 'n' sps artificiais aleatoriamente (Parte 1), modelar ausencias (Parte 2), modelar distribuicao com pseudoausencias melhoradas (Parte 3) e analise dos resultados atraves de graficos simples (Parte 4)
## Anderson A. Eduardo, outubro/2014

##registrando hora do incio
timeOne = Sys.time()

##abrindo pacotes necessarios
library(raster)
library(biomod2)

##definindo prametros e variaveis globais
projectFolder = "/home/anderson/Documentos/Projetos/Improved pseudo-absences" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
maxentFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/maxent' #pasta para resultados do maxent
## spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
## sdmTypes = c('normal','optimized')
sampleSizes = c(10,20,40,80)
NumRep = 10 #numero de replicas (de cada cenario amostral)
##variaveis preditoras
predictors = stack(list.files(path=envVarPaths[1],full.names=TRUE, pattern='.asc')) #predictors com todas as variaveis (presente)
predictors = predictors[[c('bioclim_01','bioclim_12')]]
predictors = stack(mask(x=predictors, mask=AmSulShape))
crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #ajustando CRS
Nsp = NumRep #numero de especies a serem criadas e trabalhadas igual ao numero de replicas
statResultsSDMnormal = data.frame() #tabela de estatisticas basicas do modelo
statResultsSDMimproved = data.frame() 

##definindo diretorio de trabalho (importante porque o biomod2 salva tudo automaticamente)
setwd(projectFolder)

for(i in 1:Nsp){
    for(j in 1:length(sampleSizes)){
        tryCatch({


            ##PARTE 1: criando as especies artificiais

            
            ##definindo parametros e variaveis locais
            ##matriz##
            datMat = as.data.frame(predictors, xy=TRUE, na.rm=TRUE) #transformando raster em data.frame
            datMat = data.frame(datMat, fSp=0)
            names(datMat) = c('lon', 'lat', 'bio1', 'bio12', paste('sp',i,sep='')) #ajustando os nomes das colunas do data.frame
            
            ## x = 1:100
            ## a=1 #altura do pico
            ## b=10 #posicao do centro
            ## c=1 #largura
            ## #
            ## fx = exp(-((x-b)^2/(2*c^2)))
            ## plot(fx~x,ylim=c(0,1.2))

            ##condicao para nao permitir distribuicoes vazias (i.e. inexistente) ou tbm sobre a Am. Sul toda. Condicao: distribuicao > 1% ou <95% da america do sul
            while( (sum(datMat[,5]) < 0.01*(nrow(datMat))) | (sum(datMat[,5]) > 0.95*(nrow(datMat))) ){
                
                ##equacoes para as dimensoes do nicho das especies
                betaBio1 = runif(n=1, min=0.001, max=1)*sample(x=c(-1,1), size=1) #parametro para cada equacao de cada especie
                betaBio12 = runif(n=1, min=0.001, max=1)*sample(x=c(-1,1), size=1) #parametro para cada equacao de cada especie
                alphaBio1 = runif(n=1, min=0.1*max(datMat$bio1), max=0.9*max(datMat$bio1)) #parametro para cada equacao de cada especie
                alphaBio12 = runif(n=1, min=0.1*max(datMat$bio12), max=0.9*max(datMat$bio12)) #parametro para cada equacao de cada especie

                ## betaBio1 = abs(rnorm(n=Nsp,mean=0.1,sd=0.1)) #parametro para cada equacao de cada especie
                ## betaBio12 = abs(rnorm(n=Nsp,mean=0.001,sd=0.1)) #parametro para cada equacao de cada especie
                ## alphaBio1 = abs(rnorm(n=Nsp,mean=quantile(x=varBio1,probs=0.5,na.rm=TRUE))) #parametro para cada equacao de cada especie
                ## alphaBio12 = abs(rnorm(n=Nsp,mean=quantile(x=varBio12,probs=0.5,na.rm=TRUE))) #parametro para cada equacao de cada especie
                varBio1 = datMat$bio1 #variavel ambiental bioclim01
                varBio12 = datMat$bio12 #variavel ambiental bioclim12

                ##solucao numerica para a equacoes do nicho de cada especie
                fBio1Sp_i = as.integer( 1/(1+exp(betaBio1*(varBio1-alphaBio1))) > 0.1 ) #solucao da equacao com output binario ("suitability")
                fBio12Sp_i = as.integer( 1/(1+exp(-betaBio12*(varBio12-alphaBio12))) > 0.1 ) #solucao da equacao com output binario ("suitability")
                ##fBio1Sp_i = 1/(1+exp(-betaBio1[i]*(varBio1-alphaBio1[i]))) #solucao da equacao com output continuo ("suitability")
                ##fBio12Sp_i = 1/(1+exp(-betaBio12[i]*(varBio12-alphaBio12[i]))) #solucao da equacao com output continuo ("suitability")
                ##fBio1Sp_i = as.integer( exp(-((varBio1-alphaBio1)^2/(2*betaBio1^2))) > 0.1 ) #solucao da equacao com output binario ("suitability")
                ##fBio12Sp_i = as.integer( exp(-((varBio12-alphaBio12)^2/(2*betaBio12^2))) > 0.1 ) #solucao da equacao com output binario ("suitability")
                ##datMat = data.frame(cbind(datMat,fSp=fBio1Sp_i*fBio12Sp_i)) #adicionando ao data.frame
                ##names(datMat)[ncol(datMat)] = paste('sp',i,sep='') #ajustando os nomes das especies no data.farme
                datMat[,5] = fBio1Sp_i*fBio12Sp_i
                ##salvando graficos das equacoes de cada especie
                ##jpeg(filename=paste('/home/anderson/Documentos/Projetos/divSpsSid/','functions_sp',i,'.jpeg',sep=''))
                ##par(mfrow=c(1,2))
                ## plot(fBio1Sp_i~varBio1,xlab='Bioclim 01',ylab='Suitability',ylim=c(0,1))
                ## plot(fBio12Sp_i~varBio12,xlab='Bioclim 12',ylab='Suitability',ylim=c(0,1))
                ##dev.off()
            }            
            
            ##raster da distribuicao modelada
            SpDist = datMat[,c('lon','lat',paste('sp',i,sep=''))] #extraindo lon/lat e suitability (ou pres-aus) de cada especie
            coordinates(SpDist) = ~lon+lat #definindo colunas das coordenadas
            gridded(SpDist) = TRUE #definindo gridded
            proj4string(SpDist) = '+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84' #definindo proj
            rasterSpDist = raster(SpDist) #criando objeto raster
            ##criando imagem da distribuicao de cada especie
            jpeg(filename=paste(projectFolder,'/virtual species/sp',i,'.jpeg',sep=''))
            plot(rasterSpDist)
            dev.off()
            writeRaster(x=rasterSpDist,filename=paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''),overwrite=TRUE)

            
            ##PARTE 2: modelando ausencias 


            ##diretorio para o biomod2 salvar resultados para SDMnormal
            setwd(file.path(projectFolder,'SDMnormal'))
            
            ##definindo variaveis e parametros locais
            ##occPoints = read.csv(paste(mainSampleFolder,sdmTypes[h],'/',spsTypes[i],'/occ',sampleSizes[j],'.csv',sep=''),header=TRUE) #abrindo pontos de ocorrencia
            values(rasterSpDist)[values(rasterSpDist)==0] = NA
            
            ##amostra de pontos
            occPoints = dismo::randomPoints(mask=rasterSpDist, n=sampleSizes[j]) #sorteando pontos da distribuicao modelada
            occPoints = data.frame(lon=occPoints[,1],lat=occPoints[,2])

            myResp <- data.frame(lon=occPoints[,1], lat=occPoints[,2])
            coordinates(myResp) <- ~ lon + lat #transformando em spatialPoints
            crs(myResp) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #transformando em spatialPoints

            ##variaveis e parametros locais especificos para o biomod2
            myRespName <- paste('sp',i,sep='') # nome do cenario atual (para biomod2)
            ##myResp <- dataSet[,c('pres')] # variavel resposta (para biomod2)
            ##myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
            myExpl = predictors  #variavel preditora (para biomod2)

            ##ajuste de dados de entrada para biomod2
            myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                                 expl.var = myExpl,
                                                 resp.name = paste(myRespName,'_sample',sampleSizes[j],'_SDMnormal',sep=''),
                                                 PA.nb.rep = 1)

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
                models.eval.meth = c('TSS','ROC'),
                do.full.models = FALSE,
                modeling.id = paste(myRespName,'_sample',sampleSizes[j],'_SDMnormal',sep=''))
            
            ##My output data
            evaluationScores = get_evaluations(myBiomodModelOut)

            ##gravando estatistcas basicas do modelo
            statResultsSDMnormal = rbind(statResultsSDMnormal,
                                         cbind(sp = i,
                                               sampleSize = sampleSizes[j],
                                               AUC = mean(evaluationScores['ROC','Testing.data',,,]),
                                               TSS = mean(evaluationScores['TSS','Testing.data',,,])))
            
            write.csv(statResultsSDMnormal, file=paste(projectFolder,'/StatisticalResults_SDMnormal','.csv',sep=''), row.names=FALSE)

            ##rodando algortmo de projecao (i.e. rodando a projecao)
            myBiomodProj <- BIOMOD_Projection(
                modeling.output = myBiomodModelOut,
                new.env = predictors,
                proj.name = paste('sp',i,'_sample',sampleSizes[j],'_SDMnormal',sep=''),
                selected.models = 'all',
                compress = 'FALSE',
                build.clamping.mask = 'FALSE',
                output.format = '.grd')
            
            ##gerando e salvando um mapa binario (threshold 10%)
            projStack = get_predictions(myBiomodProj) #extrai as projecoes
            projStackBIN = BinaryTransformation(projStack, 10)
            
            ##writeRaster(projStackBIN,file=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[j],'.replica',k,'/proj_',l,'kyr/proj_',i,'kyr','.sample',sampleSizes[j],'.replica',k,'_BIN.asc',sep=''),row.names=FALSE)

            ##projStackBIN = projStack>0.5  #BinaryTransformation(projStack,"10")


            ##PARTE 3: SDM com pseudoausencias melhoradas

            
            ##diretorio para biomod salvar os resultados
            setwd(file.path(projectFolder,'SDMimproved'))

            ##projecoes de ausencias do SDM (rodado na etapa 2, acima)
            projAbs = sum(projStackBIN)
            
            ##for (k in 1:nlayers(projStackBIN)){

            ##definindo variaveis e parametros locais
            betterPseudo = list()
            betterPseudoVar = list()

            ## >> AMOSTRANDO PSEUDOAUSENCIAS MELHORADAS <<
            values(projAbs)[values(projAbs) != 0] = NA            
            betterPseudoPoints = dismo::randomPoints(mask=projAbs, n=1000) #sorteando pontos da distribuicao modelada
            ##betterPseudoDF = extract(projStackBIN[[k]], be tterPseudoPoints) #distinguindo entre occ e ausencia
            betterPseudoDF = data.frame(lon=betterPseudoPoints[,1], lat=betterPseudoPoints[,2], occ=0) #distinguindo entre occ e ausencia
            ##betterPseudo[[k]] =  data.frame(lon=betterPseudoPoints[,1], lat=betterPseudoPoints[,2], occ=betterPseudoDF) #data.frame
            ##betterPseudo[[k]] = betterPseudo[[k]][which(betterPseudo[[k]]$occ==0),] #excluindo as presencas
            betterPseudoVar = extract(predictors,betterPseudoDF[,c('lon','lat')]) #obtendo as variaveis preditoras nos pontos
            betterPseudoDF = data.frame(betterPseudoDF,betterPseudoVar) #motando dataset
            
            ## plot(projStackBIN[[k]])
            ## points(betterPseudo[[k]][,c('lon','lat')])

            ##definindo variaveis e parametros locais para o biomod2 (que rodara a seguir)
            occPointsDF = data.frame(occPoints, occ=1, extract(predictors,occPoints)) #assumindo que o stack predictors contem apenas as variaveis empregadas no projeto atual
            
            ##agrupando ocorrencias e pseudo-ausencias melhoradas
            dataSet = data.frame(rbind(occPointsDF,betterPseudoDF)) #planilha de dados no formato SWD
            ##variaveis e parametros locais especificos para o biomod2
            myRespName <- paste('sp',i,sep='')  # nome do cenario atual (para biomod2)
            myResp <- dataSet[,c('occ')] # variavel resposta (para biomod2)
            myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
            myExpl = dataSet[,c('bioclim_01','bioclim_12')]  #variavel preditora (para biomod2)
            
            ##ajuste de dados de entrada para biomod2
            myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                                 expl.var = myExpl,
                                                 resp.xy = myRespXY,
                                                 resp.name = paste(myRespName,'_sample',sampleSizes[j],'_SDMimproved',sep=''))
            
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
                models.eval.meth = c('TSS','ROC'),
                do.full.models = FALSE,
                modeling.id = paste(myRespName,'_sample',sampleSizes[j],'_SDMimproved',sep=''))
            
            ##My output data
            evaluationScores = get_evaluations(myBiomodModelOut)
            
            ##gravando estatistcas basicas do modelo
            statResultsSDMimproved = rbind(statResultsSDMimproved,
                                           cbind(sp = i,
                                                 sampleSize = sampleSizes[j],
                                                 AUC = mean(evaluationScores['ROC','Testing.data',,,]),
                                                 TSS = mean(evaluationScores['TSS','Testing.data',,,])))
            
            write.csv(statResultsSDMimproved,file=paste(projectFolder,'/StatisticalResults_SDMimproved','.csv',sep=''),row.names=FALSE)

            ##implementando projecoes do modelo
            
            ##rodando algortmo de projecao (i.e. rodando a projecao)
            myBiomodProj <- BIOMOD_Projection(
                modeling.output = myBiomodModelOut,
                new.env = predictors,
                proj.name = paste(myRespName,'_sample',sampleSizes[j],'_SDMimproved',sep=''),
                selected.models = 'all',
                compress = 'FALSE',
                build.clamping.mask = 'FALSE',
                output.format = '.grd')
            
            ##gerando e salvando um mapa binario (threshold 10%)
            projStack = get_predictions(myBiomodProj) #extrai as projecoes
            projStackBIN = BinaryTransformation(stack(mean(projStack)),'10')
            
            ##writeRaster(projStackBIN,file=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[j],'.replica',k,'/proj_',l,'kyr/proj_',i,'kyr','.sample',sampleSizes[j],'.replica',k,'_BIN.asc',sep=''),row.names=FALSE)

        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
}

##calculando tempo gasto para processamento
timeOne - Sys.time()


##PARTE 4: graficos e analise dos resultados


##ambiente para salvar os graficos
setwd(projectFolder)

##planilhas de dados
#statResultsSDMimproved = read.csv(paste(projectFolder,'/StatisticalResults_SDMimproved','.csv',sep=''), header=TRUE)
#statResultsSDMnormal = read.csv(paste(projectFolder,'/StatisticalResults_SDMnormal','.csv',sep=''), header=TRUE)

statResultsSDMimproved = read.csv(paste(projectFolder,'/improved','.csv',sep=''), header=TRUE)
statResultsSDMnormal = read.csv(paste(projectFolder,'/normal','.csv',sep=''), header=TRUE)

##boxplot AUC
jpeg(filename='boxplotAUC.jpeg')
boxplot(statResultsSDMnormal$AUC, statResultsSDMimproved$AUC, ylim=c(0,1), names=c('SDM normal','SDM improved'), ylab='AUC')
dev.off()

##boxplot TSS
jpeg(filename='boxplotTSS.jpeg')
boxplot(statResultsSDMnormal$TSS, statResultsSDMimproved$TSS, ylim=c(0,1), names=c('SDM normal','SDM improved'), ylab='TSS')
dev.off()

##AUC x sample size
jpeg(filename='AUC_&_sampleSize.jpeg')
plot(statResultsSDMnormal$AUC~statResultsSDMnormal$sampleSize, ylim=c(0,1), xlim=c(0,90), cex=2, pch=19, col=rgb(0,0,0,0.5), xlab='Sample size', ylab='AUC')
tendenciaSDMnormalAUC = lm(statResultsSDMnormal$AUC~statResultsSDMnormal$sampleSize)
abline(tendenciaSDMnormalAUC)
points(statResultsSDMimproved$AUC~statResultsSDMimproved$sampleSize, ylim=c(0,1), xlim=c(0,90), cex=1.5, pch=20, col=rgb(0,0,1,0.5))
tendenciaSDMimprovedAUC = lm(statResultsSDMimproved$AUC~statResultsSDMimproved$sampleSize)
abline(tendenciaSDMimprovedAUC, col='blue')
legend('bottomleft',legend=c('SDM normal','SDM improved'), pch=c(19,20), col=c(rgb(0,0,0,0.5),rgb(0,0,1,0.5)))
dev.off()

##TSS x sample size
jpeg(filename='TSS_&_sampleSize.jpeg')
plot(statResultsSDMnormal$TSS~statResultsSDMnormal$sampleSize, ylim=c(0,1), xlim=c(0,90), cex=2, pch=19, col=rgb(0,0,0,0.5), xlab='Sample size', ylab='TSS')
tendenciaSDMnormalTSS = lm(statResultsSDMnormal$TSS~statResultsSDMnormal$sampleSize)
abline(tendenciaSDMnormalTSS)
points(statResultsSDMimproved$TSS~statResultsSDMimproved$sampleSize, ylim=c(0,1), xlim=c(0,90), cex=1.5, pch=20, col=rgb(0,0,1,0.5))
tendenciaSDMimprovedTSS = lm(statResultsSDMimproved$TSS~statResultsSDMimproved$sampleSize)
abline(tendenciaSDMimprovedTSS, col='blue')
legend('bottomleft',legend=c('SDM normal','SDM improved'), pch=c(19,20), col=c(rgb(0,0,0,0.5),rgb(0,0,1,0.5)))
dev.off()


##teste de significancia

aucTest =  wilcox.test(statResultsSDMimproved$AUC,statResultsSDMnormal$AUC)
aucTest

tssTest =  wilcox.test(statResultsSDMimproved$TSS,statResultsSDMnormal$TSS)
tssTest

##comparacao em termos de porcentagem

median(statResultsSDMimproved$AUC)/median(statResultsSDMnormal$AUC) 
