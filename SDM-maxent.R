##########################################################################
############################## SDM WITH MAXENT  ##########################
###########################Anderson A. Eduardo############################
##########################################################################
library(raster)
library(maptools)
library(dismo)
Sys.setenv(JAVA_HOME='/usr/lib/jvm/java-7-openjdk-amd64') # for 64-bit version
#Windows#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_91') # for 64-bit version
library(rJava)
source("/home/anderson/R/R-Scripts/TSSfunction.R")

###PRIMEIRA PARTE: planilha de presencas, backgrownd e variaveis ambientais

##definindo as pastas de trabalho
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/"
spOccFolder = "/home/anderson/PosDoc/dados_ocorrencia/PO_unique/"
projectFolder = "/home/anderson/PosDoc/teste/"

##abrindo as variaveis climaticas
##abrindo shape da America do Sul
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp")

##abrindo e cortando camadas de variaveis ambientais para o presente
filesRaw <- stack(list.files(path=paste(envVarFolder,"dados_projeto/000",sep=''), pattern='asc', full.names=TRUE)) 
files = mask(filesRaw,AmSulShape) #cortando para Am. do Sul

##abrindo e cortando camads de variaveis ambientais para projecao
##filesProjectionRaw <- stack(list.files(path=paste(envVarFolder,"dados_projeto/021",sep=''), pattern='asc', full.names=T)) ###abrindo camandas para projecao (passado, futuro, outro local, etc)
##filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul

##testando correcaloes
## test<-getValues(files)
## cor.matrix <- as.data.frame(cor(test, use="complete.obs"))
#write.csv(cor.matrix,'cor_matrix.csv')

##remove highly correlated variables Bio1,Bio3,Bio9,Bio13,Bio14
#files.crop.sub <- dropLayer(files, c(1,2,5,6)) #remove selected layers
files.crop.sub = files[[c('bioclim_10','bioclim_11','bioclim_16','bioclim_17')]] #choose selected layers
##files.crop.sub.projection <- dropLayer(filesProjection, c(1,2,5,6))

##definindo os objetos para as variaveis preditoras
predictors <- stack(files.crop.sub)
##predictorsProjection = files.crop.sub.projection

##Criando objeto com a lista de especies
occ.sps <- list.files(paste(spOccFolder,sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))
##fosseis
occ.sps.fosseis = read.csv(paste(spOccFolder,"fosseis/fosseis.csv",sep=''),header=TRUE)
splist.fosseis = lapply(occ.sps.fosseis[,1],as.character)

###SEGUNDA PARTE: rodando SDMs para as especies (e fazendo projecoes)

##registrando o tempo de processamento
ptm <- proc.time()

##abrindo um data.frame para armazenar os resultados de AUC
resultsEvaluationMX<-data.frame(algorithm=character(), species=character(), auc=numeric(),tss=numeric(),stringsAsFactors=FALSE)
index=0 #auxiliara na criacao do data.frame durante o loop
fossilPointsSuitability = data.frame() #objetvo em que serao gravadas as projcoes de suitability especificamente para cada ponto de registro fossil

for (i in 1:length(splist)){
    
    especie = splist[i] #selecting the species
    sp.file <- read.csv(paste(spOccFolder,especie,".csv",sep=""),header=TRUE) ### read sp occurrence
    sp.occ <- sp.file[,2:3] ## select lat-long columns
    
    ##extraindo dados da variavel climatica nos pontos de ocorrencia
    presencesVars <- extract(predictors, sp.occ, method='bilinear', buffer=NULL, fun=NULL)
    
    ##criando um vetor de presenca para usar em uma coluna de presenca/ausencia na tabela final
    pres = rep(1, nrow(presencesVars))
    
    ##juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia
    presencesData = data.frame(cbind(sp.occ,pres,presencesVars))
    presencesData = presencesData[complete.cases(presencesData),]
    presencesData = presencesData[complete.cases(presencesData),]
    
    ##criando pontos de background
    background1 <- randomPoints(mask=predictors[[1]], n=1000, p=presencesData[,c("latitude","longitude")], excludep=TRUE) #10*nrow(presencesData)
    background2 <- round(background1, digits=4)
    background3 <- background2[!duplicated(background2),]
    background4 <- background3[complete.cases(background3),]
    background <- data.frame(background4)
    colnames(background) <- c("longitude", "latitude")

    ##extraindo dados da variavel climatica nos pontos de background
    backgroundVars <- extract(predictors, background, method='bilinear', buffer=NULL, fun=NULL)

    ##criando um vetor de ausencias para usar em uma coluna de presenca/ausencia na tabela final
    pres = rep(0, nrow(backgroundVars))
    
    ##juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia    
    backgroundData = data.frame(cbind(background,pres,backgroundVars))
    backgroundData = backgroundData[complete.cases(backgroundData),]    
    
    ##planilha de dados final
    dataSet = data.frame(rbind(presencesData,backgroundData))

    ##avaliando o modelo
    projecaoSuitability = list()
    evaluation = list()
    
    for (j in 1:50){
        tryCatch({# bootstrapping with 10 replications
            
            ##reparando uma porcao dos dados de presenca e ausencia (background) para calibrar (treinar) o modelo
            rand = round(0.75*runif(nrow(presencesData)))
            presencesTrain = presencesData[rand==0,]
            backgroundTrain = backgroundData[rand==0,]
            
            ##juntando presencas e ausencias da calibracao
            PresBackTrainRaw <- data.frame(rbind(presencesTrain, backgroundTrain))
            PresBackTrainRaw = PresBackTrainRaw[!duplicated(PresBackTrainRaw[,c('longitude','latitude')]),] #selecionar colunas de longitude e latitude
            PresBackTrainRaw <- PresBackTrainRaw[complete.cases(PresBackTrainRaw),]
            PresBackTrain = PresBackTrainRaw

            ##CRIANDO E RODANDO O MODELO (esquema do SWD - sample with data)##    
            MX <- maxent(x=PresBackTrain[c(-1,-2,-3)],p=PresBackTrain$pres,args=c('responsecurves=true','jackknife=true','randomseed=true','randomtestpoints=0','betamultiplier=1','replicates=1','replicatetype=Subsample','writebackgroundpredictions=true','linear=true','quadratic=true','product=false','threshold=false','hinge=false','maximumiterations=1000','convergencethreshold=1.0E-5','threads=2'))

            ##fazendo a projecao
            MX.projection = maxent(x=dataSet[c(-1,-2,-3)],p=dataSet$pres,args=c('responsecurves=true','jackknife=true','randomseed=true','randomtestpoints=0','betamultiplier=1','replicates=1','replicatetype=Subsample','writebackgroundpredictions=true','linear=true','quadratic=true','product=false','threshold=false','hinge=false','maximumiterations=1000','convergencethreshold=1.0E-5','threads=2')) #calibrando o modelo com todos os dados de ocorrência
            projecaoSuitability = append(projecaoSuitability, predict(predictors,MX.projection))
            
            ##pegando a porcao dos dados separados para a avaliacao (validacao) do modelo
            presencesTest = presencesData[rand==1,]
            presencesTest <- presencesTest[,c('longitude','latitude')]
            backgroundTest = backgroundData[rand==1,]
            backgroundTest = backgroundTest[,c('longitude','latitude')]
            
            ##rodando a avaliacao do modelo
            evaluation = append(evaluation, evaluate(p=presencesTest,a=backgroundTest,m=MX,x=predictors,type='response'))
            
        }, error=function(e){cat("ERROR :",conditionMessage(e),"\n")}
        )}  

    ##PROJETANDO o nicho no espaco atraves do modelo ajustado##
    projecaoSuitabilityStack = stack(projecaoSuitability)
    projecaoSuitabilityMean = mean(projecaoSuitabilityStack)
    ##projecaoSuitability <- predict(predictors, MX, type='response')
    
    ##gravando um raster com o mapa de projecao gerado pelo modelo
    writeRaster(projecaoSuitabilityMean,filename=paste(projectFolder,"Maxent/",splist[i],"/",splist[i],".asc", sep=""),overwrite=TRUE)

    ##criando um mapa binario
    thresholdValues = NULL
    aucValues = NULL
    for (i in 1:length(evaluation)){
        thresholdValues <- append(thresholdValues, threshold(evaluation[[i]],'spec_sens'))
        aucValues = append(aucValues, evaluation[[i]]@auc)
    }
    aucMean = mean(aucValues)
    thresholdMean = mean(thresholdValues)
    bin <- projecaoSuitabilityMean > thresholdMean #apply threshold to transform logistic output into binary maps

    #salvando um raster com o mapa binario
    writeRaster(bin,filename=paste(projectFolder,"Maxent/",splist[i],"/",splist[i],".asc",sep=""),overwrite=TRUE)

    ##calculando o TSS a partir do mapa binario (usando minha propria funcao TSS)
    tss = TSSfunction(binaryMap=bin>0, spOccPoints=presencesData[,c('longitude','latitude')])
    
    #registrando o valor de AUC medio em uma tabela
    index=index+1
    resultsEvaluationMX[index, "algorithm"] = 'MX'
    resultsEvaluationMX[index, "species"] = especie
    resultsEvaluationMX[index, "auc"] = aucMean
    resultsEvaluationMX[index, "tss"] = tss

    ##PROJECAO PARA O PASSADO##

    ##abrindo os dados de registros fosseis  para uma especie
    sp.fossil.data = occ.sps.fosseis[occ.sps.fosseis$species==especie,] #ATENCAO: este script nao funciona se houver mais de um registro fossil por camada de tempo usada para projecao

    for(l in 1:nrow(sp.fossil.data)){#loop para cada registro fossil de uma especie

        ##definindoo fossil
        sp.fossil = sp.fossil.data[l,]

        ##abrindo as variaveis ambientais do tempo do fossil
        filesProjectionRaw <- stack(list.files(path = paste(envVarFolder,"dados_projeto/0",sp.fossil$kyr,sep=""), pattern='asc', full.names=T)) ###abrindo camandas para projecao (passado, futuro, outro local, etc)
        filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul
        files.crop.sub.projection <- dropLayer(filesProjection, c(1,2,5,6)) #removendo as camadas que mostraram correlacao
        predictorsProjection = files.crop.sub.projection #preditoras para o tempo do fossil

        ##PROJETANDO o nicho no espaco atraves do modelo ajustado##
        projecaoSuitabilityPassado <- predict(predictorsProjection, MX.projection)

        ##salvando um raster com a projecao do modelo para o tempo do fossil
        writeRaster(projecaoSuitabilityPassado,filename=paste(projectFolder,"Maxent/Passado/",splist[i],"/",splist[i],'-',sp.fossil$kyr," K years BP.asc", sep=""),overwrite=TRUE)

        ##criando um mapa binario para a projecao do modelo (empregando o threshold que ja foi criado apos a avaliacao do modelo)
        bin <- projecaoSuitabilityPassado > thresholdMean#apply threshold to transform logistic output into binary maps
        
        ##salvando um raster com a projecao do modelo para o tempo do fossil
        writeRaster(bin,filename=paste(projectFolder,"Maxent/Passado/",splist[i],"/",splist[i],'-',sp.fossil$kyr,"K years BP - BINARIO.asc",sep=""),overwrite=TRUE)

        ##criando um objeto com as coordenadas do registro fossil
        fossilPoints = sp.fossil
        fossilPoints = fossilPoints[,c('longitude','latitude')]

        ##obtendo a projecao de qualidade de habitat especificamente para o ponto do fossil
        ##fossilPointsVars = extract(predictorsProjection,fossilPoints)
        ##predict(MX, fossilPointsVars)
        fossilPoints.MX = extract(projecaoSuitabilityPassado,fossilPoints,method='bilinear') 
        fossilPointsSuitability = rbind(fossilPointsSuitability,data.frame(algorithm='Maxent',species=especie,kyr=sp.fossil$kyr,suitability=fossilPoints.MX))
        
    }
}

#salvando a tabela de dados da avaliacao dos modelos
write.table(resultsEvaluationMX,file=paste(projectFolder,"Maxent/","AUC&TSS.csv",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")

write.table(fossilPointsSuitability,file=paste(projectFolder,"Maxent/","suitabilityNoPontoFossil.csv",sep=""))

#fechando e informando o tempo de processamento
msgm= proc.time() - ptm
print(paste('Tempo gasto para rodar MX: ', msgm[3]/60,' minutos',sep=''))
