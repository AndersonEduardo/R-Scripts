library(biomod2)
library(maptools)
library(dismo)

##########################################################################
#########################TESTANDO BIOMOD##################################

##planilha de presencas, backgrownd e variaveis ambientais

##DEFININDO PASTAS DE TRABALHO##
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/"
spOccFolder = "/home/anderson/PosDoc/dados_ocorrencia/PO_unique/"
projectFolder = "/home/anderson/PosDoc/teste/"

####ABRINDO AS VARIAVEIS CLIMATICAS#####
#abrindo shape da America do Sul
AmSulShape = readShapePoly("/home/anderson/PosDoc/Am_Sul/borders.shp")

#abrindo e cortando camads de variaveis ambientais para o presente
filesRaw <- stack(list.files(path=paste(envVarFolder,"dados_projeto/000",sep=''), pattern='asc', full.names=T)) ### stack all rasters in Bioclim folder
#files <- stack(list.files(path = "/home/anderson/R/PosDoc/dados_ambientais/bcmidbi_2-5m _asc/dados_ambientais_para_projeto", pattern='asc', full.names=T))
files = mask(filesRaw,AmSulShape) #cortando para Am. do Sul

#abrindo e cortando camads de variaveis ambientais para o passado
filesProjectionRaw <- stack(list.files(path=paste(envVarFolder,"dados_projeto/021",sep=''), pattern='asc', full.names=T)) ###abrindo camandas para projecao (passado, futuro, outro local, etc)
filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul

#testando correcaloes
## test<-getValues(files)
## cor.matrix <- as.data.frame(cor(test, use="complete.obs"))
#write.csv(cor.matrix,'cor_matrix.csv')

#remove highly correlated variables Bio1,Bio3,Bio9,Bio13,Bio14
files.crop.sub <- dropLayer(files, c(1,2,5,6)) #### remove selected layers
files.crop.sub.projection <- dropLayer(filesProjection, c(1,2,5,6))

#remover as mesmas camadas dos dados para projecao
#test2<-getValues(files.crop.sub)
#cor.matrix2<- cor(test2, use="complete.obs")
#write.csv(cor.matrix2,'cor.matrix2.csv')

#definindo os objetos para as variaveis preditoras
predictors <- files.crop.sub
predictorsProjection = files.crop.sub.projection

########## Criando objetos com a lista de especies #############
occ.sps <- list.files(paste(spOccFolder,sep=''),pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\.csv")))
##fosseis
occ.sps.fosseis = read.csv(paste(spOccFolder,"fosseis/fosseis.csv",sep=''),header=T)
splist.fosseis = lapply(occ.sps.fosseis[,1],as.character)


for (i in 1:length(splist)){
    especie = splist[i] #escolher qual especie
    sp.file <- read.csv(paste(spOccFolder,especie,".csv",sep=""),h=T) ### read sp occurrence
    sp.occ <- sp.file[,2:3] ## select lat long columns
    
    ##extraindo dados da variavel climatica nos pontos de ocorrencia
    presencesVars <- extract(predictors, sp.occ, method='bilinear', buffer=NULL, fun=NULL)

    ##criando um vetor de presenca para usar em uma coluna de presenca/ausencia na tabela final
    pres = rep(1, nrow(presencesVars))

    ##juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia
    presencesData = data.frame(cbind(presencesVars,pres,sp.occ))
    presencesData = presencesData[complete.cases(presencesData),]

    ##criando ausencias para o background
    background1 <- randomPoints(mask=predictors[[1]], n=5000, p=presencesData[,c("latitude","longitude")], excludep=TRUE)
    background2 <- round(background1, digits=4)
    background3 <- background2[!duplicated(background2),]
    background4 <- background3[complete.cases(background3),]
    background <- data.frame(background4)
    colnames(background) <- c("longitude", "latitude")

    ##extraindo dados da variavel climatica nos pontos de background
    ausencesVars <- extract(predictors, background, method='bilinear', buffer=NULL, fun=NULL)

    ##criando um vetor de ausencias para usar em uma coluna de presenca/ausencia na tabela final
    pres = rep(0, nrow(ausencesVars))

    ##juntando dados das variaveis climaticas nos pontos de ocorrencia, coordenadas de ocorrencia e o vetor (coluna na tabela) para presenca/ausencia    
    ausencesData = data.frame(cbind(ausencesVars,pres,background))
    
    ##planilha de dados final
    dataSet = data.frame(rbind(presencesData,ausencesData))

    ###DADOS DE ENTRADA PARA O BIOMOD2###

    setwd(paste(projectFolder,'biomod',sep=''))
    ##myResp = rep(1,nrow(myRespXY))
    myResp = dataSet[,'pres']
    predictors = stack(predictors)
    ##myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
    myRespXY = dataSet[,c('longitude','latitude')]
    myRespName = splist[i]
    

    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = predictors,
                                         resp.xy = myRespXY,
                                         resp.name = myRespName)
    
    myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips=list(path_to_maxent.jar="/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java",maximumiterations=2000,memory_allocated=NULL))
    
    myBiomodModelOut <- BIOMOD_Modeling(
        myBiomodData,
        models = c('GLM','RF','MAXENT.Phillips'),
        models.options = myBiomodOption,
        NbRunEval = 3,
        DataSplit = 75,
        VarImport = 3,
        models.eval.meth = c('TSS','ROC'),
        SaveObj = TRUE,
        rescal.all.models = TRUE,
        do.full.models = FALSE,
        modeling.id = paste(myRespName,"FirstModeling",sep=""))
    
###PROJECAO PARA O PRESENTE###
    myBiomodProj <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = predictors,
        proj.name = '000',
        selected.models = paste(myBiomodModelOut@models.computed,sep=''),
        binary.meth = 'TSS',
        compress = FALSE,
        clamping.mask = TRUE,
        output.format = '.grd',
        on_0_1000 = FALSE)

    ##My output data
    projStack = get_predictions(myBiomodProj)
    varImportance = get_variables_importance(myBiomodModelOut)
    evaluationScores = get_evaluations(myBiomodModelOut)
    ##
    writeRaster(projStack,filename=paste(projectFolder,'biomod/myOutput/',names(projStack),'_000',sep=''),bylayer=TRUE,format='ascii',overwrite=TRUE)
    write.csv(data.frame(varImportance),paste(projectFolder,'biomod/myOutput/varImportance/varImportance_',myRespName,'_000.csv',sep=''),row.names=TRUE)
    write.csv(data.frame(evaluationScores),paste(projectFolder,'biomod/myOutput/evaluationScores/evaluationScores_',myRespName,'_000.csv',sep=''),row.names=TRUE)
        
###PROJECAO PARA O PASSADO###
    
    ##abrindo os dados de registros fosseis  para uma especie
    sp.fossil.data = occ.sps.fosseis[occ.sps.fosseis$species==especie,] #ATENCAO: este script nao funciona se houver mais de um registro fossil por camada de tempo usada para projecao
    
    for(l in 1:nrow(sp.fossil.data)){#loop para cada registro fossil de uma especie
        
        ##definindoo fossil
        sp.fossil = sp.fossil.data[l,]
        
        ##abrindo as variaveis ambientais do tempo do fossil
        filesProjectionRaw <- stack(list.files(path = paste(envVarFolder,"dados_projeto/0",sp.fossil$kyr,sep=""), pattern='asc', full.names=T)) ###abrindo camandas para projecao (passado, futuro, outro local, etc)
        filesProjection = mask(filesProjectionRaw,AmSulShape) #cortando para Am. do Sul
        files.crop.sub.projection <- dropLayer(filesProjection, c(1,2,5,6)) #removendo as camadas que mostraram correlacao
        predictorsProjection = stack(files.crop.sub.projection) #preditoras para o tempo do fossil
        
        ##PROJETANDO o nicho no espaco atraves do modelo ajustado##
        
        myBiomodProj <- BIOMOD_Projection(
            modeling.output = myBiomodModelOut,
            new.env = predictorsProjection,
            proj.name = paste(sp.fossil$kyr,'kyrBP',sep=''),
            selected.models = paste(myBiomodModelOut@models.computed,sep=''),
            binary.meth = 'TSS',
            compress = TRUE,
            clamping.mask = TRUE,
            output.format = '.grd',
            on_0_1000 = FALSE)        
        
        ##My outputs
        projStackPass = get_predictions(myBiomodProj)
        varImportancePass = get_variables_importance(myBiomodModelOut)
        evaluationScoresPass = get_evaluations(myBiomodModelOut)
        ##
        writeRaster(projStackPass,filename=paste(projectFolder,'biomod/myOutput/',names(projStack),'_',sp.fossil$kyr,'kyrBP',sep=''),bylayer=TRUE,format='ascii',overwrite=TRUE)
        ## write.csv(data.frame(varImportancePass),paste(projectFolder,'biomod/myOutput/varImportance/varImportance_',sp.fossil$kyr,'kyrBP.csv',sep=''),row.names=TRUE)
        ## write.csv(data.frame(evaluationScoresPass),paste(projectFolder,'biomod/myOutput/evaluationScores/evaluationScores_',sp.fossil$kyr,'kyrBP.csv',sep=''),row.names=TRUE)

        ##suitability no ponto fossil:
        ##criando um objeto com as coordenadas do registro fossil
        fossilPoints = sp.fossil
        fossilPoints = cbind(fossilPoints$longitude, fossilPoints$latitude)
        ##extratindo valor do suitability nas coordenadas do registro fossil
        suitabNoPontoFossil = extract(projStackPass,fossilPoints)
        write.csv(suitabNoPontoFossil,paste(projectFolder,'biomod/myOutput/suitabilityNoPontoFossil/',sp.fossil$species,sp.fossil$kyr,'kyrBP',sep=''))
    }
}   

    

