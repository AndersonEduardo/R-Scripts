## Este script contem algoritmos para gerar 'n' sps artificiais aleatoriamente (Parte 1), modelar ausencias (Parte 2), modelar distribuicao com pseudoausencias melhoradas (Parte 3) e analise dos resultados atraves de graficos simples (Parte 4)
## Anderson A. Eduardo, outubro/2014

##registrando hora do incio
timeOne = Sys.time()

##abrindo pacotes necessarios
library(raster)
library(biomod2)

##definindo prametros e variaveis globais
projectFolder = "/home/anderson/Documentos/Projetos/Improved pseudo-absences_TESTE" #pasta do projeto
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
maxentFolder = '/home/anderson/Documentos/Projetos/Sps artificiais/maxent' #pasta para resultados do maxent
## spsTypes = c('spHW', 'spHD', 'spCD') #nomes das especies
## sdmTypes = c('normal','optimized')
sampleSizes = 100  #c(10,20,40,80,160)
NumRep = 2 #10 #numero de replicas (de cada cenario amostral)
##variaveis preditoras
elevation = raster('/home/anderson/PosDoc/dados_ambientais/DEM/DEM.tif')
predictors = stack(list.files(path=envVarPaths[1],full.names=TRUE, pattern='.asc'),elevation) #predictors com todas as variaveis (presente)
predictors = predictors[[c('bioclim_01','bioclim_12','DEM')]]
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
            names(datMat) = c('lon', 'lat', 'bio1', 'bio12', 'elevation', paste('sp',i,sep='')) #ajustando os nomes das colunas do data.frame
            
            ## x = 1:100
            ## a=1 #altura do pico
            ## b=10 #posicao do centro
            ## c=1 #largura
            ## #
            ## fx = exp(-((x-b)^2/(2*c^2)))
            ## plot(fx~x,ylim=c(0,1.2))

            ##condicao para nao permitir distribuicoes vazias (i.e. inexistente) ou tbm sobre a Am. Sul toda. Condicao: distribuicao > 1% ou <95% da america do sul
            while( (sum(datMat[,paste('sp',i,sep='')]) < 0.05*(nrow(datMat))) | (sum(datMat[,paste('sp',i,sep='')]) > 0.5*(nrow(datMat))) ){

                ##equacoes para as dimensoes do nicho das especies
                betaBio1 = runif(n=1, min=0.001, max=1)*sample(x=c(-1,1), size=1) #parametro para cada equacao de cada especie
                betaBio12 = runif(n=1, min=0.001, max=1)*sample(x=c(-1,1), size=1) #parametro para cada equacao de cada especie
                betaElev = runif(n=1, min=0.001, max=1)*sample(x=c(-1,1), size=1) #parametro para cada equacao de cada especie
                ##
                alphaBio1 = runif(n=1, min=quantile(datMat$bio1, probs=0.25, na.rm=TRUE), max=quantile(datMat$bio1, probs=0.75, na.rm=TRUE)) #parametro para cada equacao de cada especie
                alphaBio12 = runif(n=1, min=quantile(datMat$bio12, probs=0.25, na.rm=TRUE), max=quantile(datMat$bio12, probs=0.75, na.rm=TRUE)) #parametro para cada equacao de cada especie
                alphaElev = runif(n=1, min=quantile(datMat$elevation, probs=0.25, na.rm=TRUE), max=quantile(datMat$elevation, probs=0.75, na.rm=TRUE)) #parametro para cada equacao de cada especie

                
                ## betaBio1 = abs(rnorm(n=Nsp,mean=0.1,sd=0.1)) #parametro para cada equacao de cada especie
                ## betaBio12 = abs(rnorm(n=Nsp,mean=0.001,sd=0.1)) #parametro para cada equacao de cada especie
                ## alphaBio1 = abs(rnorm(n=Nsp,mean=quantile(x=varBio1,probs=0.5,na.rm=TRUE))) #parametro para cada equacao de cada especie
                ## alphaBio12 = abs(rnorm(n=Nsp,mean=quantile(x=varBio12,probs=0.5,na.rm=TRUE))) #parametro para cada equacao de cada especie
                varBio1 = datMat$bio1 #variavel ambiental bioclim01
                varBio12 = datMat$bio12 #variavel ambiental bioclim12
                varElev = datMat$elevation
                
                ##solucao numerica para a equacoes do nicho de cada especie
                fBio1Sp_i = as.integer( 1/(1+exp(betaBio1*(varBio1-alphaBio1))) > 0.1 ) #solucao da equacao com output binario ("suitability")
                fBio12Sp_i = as.integer( 1/(1+exp(-betaBio12*(varBio12-alphaBio12))) > 0.1 ) #solucao da equacao com output binario ("suitability")
                fElevSp_i = as.integer( 1/(1+exp(-betaElev*(varElev-alphaElev))) > 0.1 ) #solucao da equacao com output binario ("suitability")
                
                ## fBio1Sp_i = 1/(1+exp(-betaBio1[i]*(varBio1-alphaBio1[i]))) #solucao da equacao com output continuo ("suitability")
                ## fBio12Sp_i = 1/(1+exp(-betaBio12[i]*(varBio12-alphaBio12[i]))) #solucao da equacao com output continuo ("suitability")
                ## fBio1Sp_i = as.integer( exp(-((varBio1-alphaBio1)^2/(2*betaBio1^2))) > 0.1 ) #solucao da equacao com output binario ("suitability")
                ## fBio12Sp_i = as.integer( exp(-((varBio12-alphaBio12)^2/(2*betaBio12^2))) > 0.1 ) #solucao da equacao com output binario ("suitability")
                
                ##datMat = data.frame(cbind(datMat,fSp=fBio1Sp_i*fBio12Sp_i)) #adicionando ao data.frame
                ##names(datMat)[ncol(datMat)] = paste('sp',i,sep='') #ajustando os nomes das especies no data.farme
                datMat[,paste('sp',i,sep='')] = fBio1Sp_i*fBio12Sp_i *fElevSp_i
                
                ##salvando graficos das equacoes de cada especie
                ##jpeg(filename=paste('/home/anderson/Documentos/Projetos/divSpsSid/','functions_sp',i,'.jpeg',sep=''))
                ##par(mfrow=c(1,2))
                ## plot(fBio1Sp_i~varBio1,xlab='Bioclim 01',ylab='Suitability',ylim=c(0,1))
                ## plot(fBio12Sp_i~varBio12,xlab='Bioclim 12',ylab='Suitability',ylim=c(0,1))
                ## plot(fElevSp_i~varElev,xlab='Elevation',ylab='Suitability',ylim=c(0,1)) 
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
            writeRaster(x=rasterSpDist, filename=paste(projectFolder,'/virtual species/sp',i,'.asc',sep=''), overwrite=TRUE)

            
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
                MAXENT.Phillips = list(path_to_maxent.jar="/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java",
                                       linear=TRUE,
                                       quadratic=TRUE,
                                       product=FALSE,
                                       threshold=FALSE,
                                       hinge=FALSE,
                                       maximumiterations=1000,
                                       convergencethreshold=1.0E-5,
                                       threads=2),
                GLM = list( type = 'quadratic',
                           interaction.level = 0,
                           myFormula = NULL,
                           test = 'AIC',
                           family = binomial(link = 'logit'),
                           mustart = 0.5,
                           control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE)),
                GAM = list( algo = 'GAM_mgcv',
                           type = 's_smoother',
                           k = -1,
                           interaction.level = 0,
                           myFormula = NULL,
                           family = binomial(link = 'logit'),
                           method = 'GCV.Cp',
                           optimizer = c('outer','newton'),
                           select = FALSE,
                           knots = NULL,
                           paraPen = NULL,
                           control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07,
                                          maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15,
                                          rank.tol = 1.49011611938477e-08,
                                          nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0),
                                          optim = list(factr=1e+07),
                                          newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0),
                                          outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE,
                                          scale.est = 'fletcher', edge.correct = FALSE)),
                MARS = list( type = 'simple',
                            interaction.level = 0,
                            myFormula = NULL,
                            nk = NULL,
                            penalty = 2,
                            thresh = 0.001,
                            nprune = NULL,
                            pmethod = 'backward'),
                CTA = list( method = 'class',
                           parms = 'default',
                           cost = NULL,
                           control = list(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 25)),
                GBM = list( distribution = 'bernoulli',
                           n.trees = 2500,
                           interaction.depth = 7,
                           n.minobsinnode = 5,
                           shrinkage = 0.001,
                           bag.fraction = 0.5,
                           train.fraction = 1,
                           cv.folds = 3,
                           keep.data = FALSE,
                           verbose = FALSE,
                           perf.method = 'cv'),
                RF = list( do.classif = TRUE,
                          ntree = 500,
                          mtry = 'default',
                          nodesize = 5,
                          maxnodes = NULL)
            )            
            
            ##rodando o(s) algoritmo(s) (i.e. SDMs)
            myBiomodModelOut <- BIOMOD_Modeling(
                myBiomodData,
                models = c('MAXENT.Phillips','GLM', 'GAM', 'MARS', 'CTA', 'GBM', 'RF'),
                models.options = myBiomodOption,
                NbRunEval = 10,
                DataSplit = 75,
                models.eval.meth = c('TSS','ROC'),
                do.full.models = FALSE,
                modeling.id = paste(myRespName,'_sample',sampleSizes[j],'_SDMnormal',sep=''))
            
            ##My output data
            evaluationScores = get_evaluations(myBiomodModelOut)

            ##maior valore de especificidade de cada algoritmo implementado (tanto para TSS qto para AUC)
            TSSspec = as.numeric(apply(evaluationScores['TSS','Specificity',,,], 1, max))
            AUCspec = as.numeric(apply(evaluationScores['ROC','Specificity',,,], 1, max))

            ##tabela auxiliar para obtencao das informacoes do melhor modelo
            tabBestScoresTSS = data.frame(evaluationScores['TSS','Specificity',,,], bestvalue=TSSspec)
            tabBestScoresAUC = data.frame(evaluationScores['ROC','Specificity',,,], bestvalue=AUCspec)

            ##vetores vazios (para os nomes dos melhores modelos)
            bestModelRunTSS = vector()
            bestModelRunAUC = vector()
            TSSvalues = vector()
            AUCvalues = vector()
            
            for (l in 1:nrow(tabBestScoresTSS)){
                bestRunNameTSS = names(tabBestScoresTSS)[match(tabBestScoresTSS$bestvalue[l], tabBestScoresTSS[l, 1:ncol(tabBestScoresTSS)-1])] #pegando o nome do melhor modelo (e.g. RUN2)
                TSSvalueBestModel = evaluationScores['TSS','Testing.data',l,bestRunNameTSS,] #pegando o TSS do melhor modelo
                bestModelRunTSS = append(bestModelRunTSS, bestRunNameTSS) #empilhando os nomes dos melhores modelos
                TSSvalues = append(TSSvalues, TSSvalueBestModel)
                ##
                bestRunNameAUC = names(tabBestScoresAUC)[match(tabBestScoresAUC$bestvalue[l], tabBestScoresAUC[l, 1:ncol(tabBestScoresAUC)-1])] #pegando o nome do melhor modelo (e.g. RUN2)
                AUCvalueBestModel = evaluationScores['ROC','Testing.data',l,bestRunNameTSS,] #pegando o AUC do melhor modelo
                bestModelRunAUC = append(bestModelRunAUC, bestRunNameAUC)
                AUCvalues = append(AUCvalues, AUCvalueBestModel)
            }
            
            
            ##gravando estatisticas basicas do melhor modelo de cada algoritmo
            statResultsSDMnormal = rbind(statResultsSDMnormal,
                                         cbind(sp = i,
                                               ##tss
                                               model = colnames(evaluationScores[,'Testing.data',,,]),
                                               meanTSSValue = as.numeric(apply(evaluationScores['TSS','Testing.data',,,], 1, mean)),
                                               meanTSSspecificity = as.numeric(apply(evaluationScores['TSS','Specificity',,,], 1, mean)),
                                               maxTSSspecificity = TSSspec,
                                               bestModelTSS = bestModelRunTSS,
                                               TSSvalue_bestModel= TSSvalues,
                                               ##auc
                                               meanAUCValue = as.numeric(apply(evaluationScores['ROC','Testing.data',,,], 1, mean)),
                                               meanAUCspecificity = as.numeric(apply(evaluationScores['ROC','Specificity',,,], 1, mean)),
                                               maxAUCspecificity = AUCspec,
                                               bestModelAUC = bestModelRunAUC,
                                               AUCvalue_bestModel= AUCvalues),
                                         stringsAsFactors = FALSE)
            
            write.csv(statResultsSDMnormal, file=paste(projectFolder,'/StatisticalResults_SDMnormal','.csv',sep=''), row.names=FALSE)
            
            ##selecao do modelo de maior sensibilidade

            ## tssMax = max(as.numeric(statResultsSDMnormal[which(statResultsSDMnormal$sp == i),]$TSSspec))
            ## aucMax = max(as.numeric(statResultsSDMnormal[which(statResultsSDMnormal$sp == i),]$AUCspec))
            ## bestAlgorithmTSS = statResultsSDMnormal[which(statResultsSDMnormal$TSSspec==tssMax),]$model
            ## bestAlgorithmAUC = statResultsSDMnormal[which(statResultsSDMnormal$TSSspec==tssMax),]$model
            ## ind = match(bestAlgorithmTSS, bestAlgorithmAUC)
            ## modelToProj = bestAlgorithmAUC[ind]
            ## modelNames = grep(pattern=paste(modelToProj,collapse='|'), x=myBiomodModelOut@models.computed, value=TRUE)

            ##maior TSS e AUC para especie da iteracao atual
            tssMax = max(TSSspec)
            aucMax = max(AUCspec)

            ##formacao de vetor com nome e RUN do melhor modelo (TSS)
            bestAlgorithmTSS = statResultsSDMnormal[which(statResultsSDMnormal$maxTSSspecificity==tssMax),]$model
            bestRunTSS = statResultsSDMnormal[which(statResultsSDMnormal$maxTSSspecificity==tssMax),]$bestModelTSS
            patternsTSS = paste(bestRunTSS,bestAlgorithmTSS,sep='_')

            ##formacao de vetor com nome e RUN do melhor modelo (AUC)
            bestAlgorithmAUC = statResultsSDMnormal[which(statResultsSDMnormal$maxAUCspecificity==aucMax),]$model
            bestRunAUC = statResultsSDMnormal[which(statResultsSDMnormal$maxAUCspecificity==aucMax),]$bestModelAUC
            patternsAUC = paste(bestRunAUC,bestAlgorithmAUC,sep='_')

            ##nomes dos melhores modelos
            modelNames = grep(pattern=paste(c(patternsTSS,patternsAUC),collapse='|'), x=myBiomodModelOut@models.computed, value=TRUE) #
                       
            ## if (bestAlgorithmTSS == bestAlgorithmAUC){
            ##     modelToProj = bestAlgorithmTSS
            ## }else{
            ##     modelToProj = c(bestAlgorithmAUC, bestAlgorithmTSS)
            ## }
            ## ind = grep(pattern=modelToProj, x=myBiomodModelOut@models.computed) ##pegando os indices
            ## modelNames = myBiomodModelOut@models.computed[ind] ##pegando os nomes dos modelos para projecao (aqueles com maior especificidade)
            
            ##rodando algortmo de projecao (i.e. rodando a projecao)
            myBiomodProj <- BIOMOD_Projection(
                modeling.output = myBiomodModelOut,
                new.env = predictors,
                proj.name = paste('sp',i,'_sample',sampleSizes[j],'_SDMnormal',sep=''),
                selected.models = modelNames,
                compress = 'TRUE',
                build.clamping.mask = 'TRUE')

            ##rodando o algoritmo de consenso dos modelos (i.e. ensemble model)
            myBiomodEM = BIOMOD_EnsembleModeling(
                modeling.output = myBiomodModelOut,
                chosen.models = modelNames)
            
            ##forecasting com o consenso dos algoritmos (i.e. ensemble projection)
            myBiomodEF = BIOMOD_EnsembleForecasting(
                EM.output = myBiomodEM,
                binary.meth = c('TSS','ROC'),
                projection.output = myBiomodProj)
            
            ##writeRaster(projStackBIN,file=paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[j],'.replica',k,'/proj_',l,'kyr/proj_',i,'kyr','.sample',sampleSizes[j],'.replica',k,'_BIN.asc',sep=''),row.names=FALSE)

            ##projStackBIN = projStack>0.5  #BinaryTransformation(projStack,"10")


            ##PARTE 3: SDM com pseudoausencias melhoradas

            
            ##diretorio para biomod salvar os resultados
            setwd(file.path(projectFolder,'SDMimproved'))

            ##definindo variaveis e parametros locais
            betterPseudo = list()
            betterPseudoVar = list()
            
            ##projecoes de ausencias do SDM (rodado na etapa 2, acima)
            binTSS = raster(paste(projectFolder,'/SDMnormal/','sp',i,'.sample',sampleSizes[j],'.SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal_sp',i,'.sample',sampleSizes[j],'.SDMnormal_ensemble_TSSbin.grd' ,sep=''))
            binAUC = raster(paste(projectFolder,'/SDMnormal/','sp',i,'.sample',sampleSizes[j],'.SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal','/proj_sp',i,'_sample',sampleSizes[j],'_SDMnormal_sp',i,'.sample',sampleSizes[j],'.SDMnormal_ensemble_ROCbin.grd' ,sep=''))
            projStackBIN = stack(binTSS,binAUC) #empilhando mapas binarios (feitos com threshold a partir do AUC e TSS)
            projAbs = sum(projStackBIN) #somando (para depois pegar areas ausencia que ambos os thresolds concordam)
                        
            ## amostrando pontos diretamente das areas de ausencia (abaixo do threshold) obtidas na etapa 1 ##
            values(projAbs)[values(projAbs) != 0] = NA  #tranformando areas diferentes de zero em NA (retando somente os dados de ausencia)          
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
                MAXENT.Phillips = list(path_to_maxent.jar="/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java",
                                       maximumiterations=1000,
                                       linear=TRUE,
                                       quadratic=TRUE,
                                       product=FALSE,
                                       threshold=FALSE,
                                       hinge=FALSE,
                                       maximumiterations=1000,
                                       convergencethreshold=1.0E-5,
                                       threads=2),
                GLM = list( type = 'quadratic',
                           interaction.level = 0,
                           myFormula = NULL,
                           test = 'AIC',
                           family = binomial(link = 'logit'),
                           mustart = 0.5,
                           control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE)),
                GAM = list( algo = 'GAM_mgcv',
                           type = 's_smoother',
                           k = -1,
                           interaction.level = 0,
                           myFormula = NULL,
                           family = binomial(link = 'logit'),
                           method = 'GCV.Cp',
                           optimizer = c('outer','newton'),
                           select = FALSE,
                           knots = NULL,
                           paraPen = NULL,
                           control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07,
                                          maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15,
                                          rank.tol = 1.49011611938477e-08,
                                          nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0),
                                          optim = list(factr=1e+07),
                                          newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0),
                                          outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE,
                                          scale.est = 'fletcher', edge.correct = FALSE)),
                MARS = list( type = 'simple',
                            interaction.level = 0,
                            myFormula = NULL,
                            nk = NULL,
                            penalty = 2,
                            thresh = 0.001,
                            nprune = NULL,
                            pmethod = 'backward'),
                CTA = list( method = 'class',
                           parms = 'default',
                           cost = NULL,
                           control = list(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 25)),
                GBM = list( distribution = 'bernoulli',
                           n.trees = 2500,
                           interaction.depth = 7,
                           n.minobsinnode = 5,
                           shrinkage = 0.001,
                           bag.fraction = 0.5,
                           train.fraction = 1,
                           cv.folds = 3,
                           keep.data = FALSE,
                           verbose = FALSE,
                           perf.method = 'cv'),
                RF = list( do.classif = TRUE,
                          ntree = 500,
                          mtry = 'default',
                          nodesize = 5,
                          maxnodes = NULL)
            )            

            ##rodando o(s) algoritmo(s) (i.e. SDMs)
            myBiomodModelOut <- BIOMOD_Modeling(
                myBiomodData,
                models = c('MAXENT.Phillips','GLM', 'GAM', 'MARS', 'CTA', 'GBM', 'RF'),
                models.options = myBiomodOption,
                NbRunEval = 10,
                DataSplit = 75,
                models.eval.meth = c('TSS','ROC'),
                do.full.models = FALSE,
                modeling.id = paste(myRespName,'_sample',sampleSizes[j],'_SDMimproved',sep=''))

            
            ##My output data
            rm(evaluationScores) ##apagando a inforacao da etapa 1
            evaluationScores = get_evaluations(myBiomodModelOut)

            ##maior valore de especificidade de cada algoritmo implementado (tanto para TSS qto para AUC)
            TSSspec = as.numeric(apply(evaluationScores['TSS','Specificity',,,], 1, max))
            AUCspec = as.numeric(apply(evaluationScores['ROC','Specificity',,,], 1, max))

            ##tabela auxiliar para obtencao das informacoes do melhor modelo
            tabBestScoresTSS = data.frame(evaluationScores['TSS','Specificity',,,], bestvalue=TSSspec)
            tabBestScoresAUC = data.frame(evaluationScores['ROC','Specificity',,,], bestvalue=AUCspec)

            ##vetores vazios (para os nomes dos melhores modelos)
            bestModelRunTSS = vector()
            bestModelRunAUC = vector()
            TSSvalues = vector()
            AUCvalues = vector()
            
            for (l in 1:nrow(tabBestScoresTSS)){
                bestRunNameTSS = names(tabBestScoresTSS)[match(tabBestScoresTSS$bestvalue[l], tabBestScoresTSS[l, 1:ncol(tabBestScoresTSS)-1])] #pegando o nome do melhor modelo (e.g. RUN2)
                TSSvalueBestModel = evaluationScores['TSS','Testing.data',l,bestRunNameTSS,] #pegando o TSS do melhor modelo
                bestModelRunTSS = append(bestModelRunTSS, bestRunNameTSS) #empilhando os nomes dos melhores modelos
                TSSvalues = append(TSSvalues, TSSvalueBestModel)
                ##
                bestRunNameAUC = names(tabBestScoresAUC)[match(tabBestScoresAUC$bestvalue[l], tabBestScoresAUC[l, 1:ncol(tabBestScoresAUC)-1])] #pegando o nome do melhor modelo (e.g. RUN2)
                AUCvalueBestModel = evaluationScores['ROC','Testing.data',l,bestRunNameTSS,] #pegando o AUC do melhor modelo
                bestModelRunAUC = append(bestModelRunAUC, bestRunNameAUC)
                AUCvalues = append(AUCvalues, AUCvalueBestModel)
            }
            
            
            ##gravando estatisticas basicas do melhor modelo de cada algoritmo
            statResultsSDMimproved = rbind(statResultsSDMimproved,
                                           cbind(sp = i,
                                                 ##tss
                                                 model = colnames(evaluationScores[,'Testing.data',,,]),
                                                 meanTSSValue = as.numeric(apply(evaluationScores['TSS','Testing.data',,,], 1, mean)),
                                                 meanTSSspecificity = as.numeric(apply(evaluationScores['TSS','Specificity',,,], 1, mean)),
                                                 maxTSSspecificity = TSSspec,
                                                 bestModelTSS = bestModelRunTSS,
                                                 TSSvalue_bestModel= TSSvalues,
                                                 ##auc
                                                 meanAUCValue = as.numeric(apply(evaluationScores['ROC','Testing.data',,,], 1, mean)),
                                                 meanAUCspecificity = as.numeric(apply(evaluationScores['ROC','Specificity',,,], 1, mean)),
                                                 maxAUCspecificity = AUCspec,
                                                 bestModelAUC = bestModelRunAUC,
                                                 AUCvalue_bestModel= AUCvalues),
                                           stringsAsFactors = FALSE)
            
            write.csv(statResultsSDMimproved, file=paste(projectFolder,'/StatisticalResults_SDMimproved','.csv',sep=''), row.names=FALSE)

            ##implementando projecoes do modelo
            
            ##rodando algortmo de projecao (i.e. rodando a projecao)
            myBiomodProj <- BIOMOD_Projection(
                modeling.output = myBiomodModelOut,
                new.env = predictors,
                proj.name = paste(myRespName,'_sample',sampleSizes[j],'_SDMimproved',sep=''),
                selected.models = 'all',
                binary.meth = 'TSS',
                compress = 'TRUE',
                build.clamping.mask = 'TRUE',
                output.format = '.grd')

        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
}

##calculando tempo gasto para processamento
timeOne - Sys.time()


##PARTE 4: graficos e analise dos resultados


##ambiente para salvar os graficos
setwd(projectFolder)

##planilhas de dados
statResultsSDMimproved = read.csv(paste(projectFolder,'/StatisticalResults_SDMimproved','.csv',sep=''), header=TRUE)
statResultsSDMnormal = read.csv(paste(projectFolder,'/StatisticalResults_SDMnormal','.csv',sep=''), header=TRUE)

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
