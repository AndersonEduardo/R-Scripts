##########################################################################
######################SDM WITH RANDOM FOREST##############################
###########################Anderson A. Eduardo############################
##########################################################################
library(randomForest)
library(raster)
library(maptools)
library(dismo)

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
files.crop.sub <- dropLayer(files, c(1,2,5,6)) #### remove selected layers
##files.crop.sub.projection <- dropLayer(filesProjection, c(1,2,5,6))

##definindo os objetos para as variaveis preditoras
predictors <- files.crop.sub
##predictorsProjection = files.crop.sub.projection

##Criando objeto com a lista de especies
occ.sps <- list.files(paste(spOccFolder,sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))
##fosseis
occ.sps.fosseis = read.csv(paste(spOccFolder,"fosseis/fosseis.csv",sep=''),header=TRUE)
splist.fosseis = lapply(occ.sps.fosseis[,1],as.character)

###SEGUNDA PARTE: rodando SDMs para as especies (e fazendo projecoes)

#registrando o tempo de processamento
ptm <- proc.time()

#abrindo um data.frame para armazenar os resultados de AUC
resultsEvaluationRF<-data.frame(algorithm=character(), species=character(), auc=numeric(),tss=numeric(),stringsAsFactors=FALSE)
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

    ##criando pontos de background
    background1 <- randomPoints(mask=predictors[[1]], n=10*nrow(presencesData), p=presencesData[,c("latitude","longitude")], excludep=TRUE)
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
    
    ##planilha de dados final
    dataSet = data.frame(rbind(presencesData,backgroundData))

    ##avaliando o modelo
    V <- numeric()#abrir un vector vacío 
    
    for (j in 1:10){
        tryCatch({# bootstrapping with 10 replications

            ##reparando uma porcao dos dados de presenca e ausencia (background) para calibrar (treinar) o modelo
            rand = round(0.75*runif(nrow(presencesData)))
            presencesTrain = presencesData[rand==0,]
            backgroundTrain = backgroundData[rand==0,]

            ##juntando presencas e ausencias da calibracao
            PresBackTrainRaw <- data.frame(rbind(presencesTrain, backgroundTrain))
            PresBackTrainRaw = PresBackTrainRaw[!duplicated(PresBackTrainRaw[,c('longituude','latitude')]),] #selecionar colunas de longitude e latitude
            PresBackTrainRaw <- PresBackTrainRaw[complete.cases(PresBackTrainRaw),]
            PresBackTrain = PresBackTrainRaw


            ##criando e rodando o modelo
            RF <- randomForest(x=PresBackTrain[,c('bioclim_10','bioclim_11','bioclim_16','bioclim_17')],y=PresBackTrain[,c('pres')], data=PresBackTrain, ntree=500)            
            
            ##pegando a porcao dos dados separados para a avaliacao (validacao) do modelo
            
            presencesTest = presencesData[rand==1,]
            presencesTest <- presencesTest[,c('longitude','latitude')]
            backgroundTest = backgroundData[rand==1,]
            backgroundTest = backgroundTest[,c('longitude','latitude')]
            
            ##rodando a avaliacao do modelo
            evaluation = evaluate(p=presencesTest,a=backgroundTest,m=RF,x=predictors)

            #registrando o valor de AUC em um objeto
            V[j]<-evaluation@"auc" #sacamos el valor de auc (fíjate que es una @ en lugar de $ para mirar dentro de los slots)y guardamos en vector
        }, error=function(e){cat("ERROR :",conditionMessage(e),"\n")}
   )}  

    #registrando a media de AUC
    auc = mean(V, na.rm=TRUE)#media de los vectores de las iternaciones de j

    ##PROJETANDO o nicho no espaco atraves do modelo ajustado##
    projecaoSuitability <- predict(predictors, RF)

    ##gravando um raster com o mapa de projecao gerado pelo modelo
    writeRaster(projecaoSuitability,filename=paste(projectFolder,"Random Forest/",splist[especie],"/",splist[especie],".asc", sep=""),overwrite=TRUE)

    ##criando um mapa binario
    threshold <- threshold(evaluation,'spec_sens')
    bin <- projecaoSuitability > threshold #apply threshold to transform logistic output into binary maps

    #salvando um raster com o mapa binario
    writeRaster(bin,filename=paste(projectFolder,"Random Forest/",splist[i],"/",splist[especie],"BINARIO.asc",sep=""),overwrite=T)

    ##calculando o TSS a partir do mapa binario (usando minha propria funcao TSS)
    tss = TSSfunction(binaryMap=bin, spOccPoints=presencesData[,c('longitude','latitude')])
    
    #registrando o valor de AUC medio em uma tabela
    index=index+1
    resultsEvaluationRF[index, "algorithm"] = 'Random Forest'
    resultsEvaluationRF[index, "species"] = especie
    resultsEvaluationRF[index, "auc"] = evaluation@auc
    resultsEvaluationRF[index, "tss"] = tss

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
        projecaoSuitabilityPassado <- predict(predictorsProjection, RF)

        ##salvando um raster com a projecao do modelo para o tempo do fossil
        writeRaster(projecaoSuitabilityPassado,filename=paste(projectFolder,"Random Forest/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr," K years BP.asc", sep=""),overwrite=TRUE)

        ##criando um mapa binario para a projecao do modelo (empregando o threshold que ja foi criado apos a avaliacao do modelo)
        bin <- projecaoSuitabilityPassado > threshold#apply threshold to transform logistic output into binary maps
        
        ##salvando um raster com a projecao do modelo para o tempo do fossil
        writeRaster(bin,filename=paste(projectFolder,"Random Forest/Passado/",splist[especie],"/",splist[especie],'-',sp.fossil$kyr,"K years BP - BINARIO.asc",sep=""),overwrite=TRUE)

        ##criando um objeto com as coordenadas do registro fossil
        fossilPoints = sp.fossil
        fossilPoints = fossilPoints[,c('longitude','latitude')]

        ##obtendo a projecao de qualidade de habitat especificamente para o ponto do fossil
        ##fossilPointsVars = extract(predictorsProjection,fossilPoints)
        ##predict(RF, fossilPointsVars)
        fossilPoints.RF = extract(projecaoSuitabilityPassado,fossilPoints,method='bilinear') 
        fossilPointsSuitability = rbind(fossilPointsSuitability,data.frame(algorithm='Random Forest',species=especie,kyr=sp.fossil$kyr,suitability=fossilPoints.RF))
    }
}

#salvando a tabela de dados da avaliacao dos modelos
write.table(resultsEvaluationRF,file=paste(projectFolder,"Random Forest/","AUC&TSS.csv",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")

write.table(fossilPointsSuitability,file=paste(projectFolder,"Random Forest/","suitabilityNoPontoFossil.csv",sep=","))

#fechando e informando o tempo de processamento
msgm= proc.time() - ptm
print(paste('Tempo gasto para rodar RANDOM FOREST: ', msgm[3]/60,' minutos',sep=''))
