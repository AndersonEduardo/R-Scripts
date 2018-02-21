## script para implementacao de SDM e analise de reducao de incerteza taxonomica de registros fosseis

##pacotes
library(biomod2)
library(raster)

##parametros e variaveis globais
rm(list=ls()) #limpando area de trabalho
maxentFolder = '/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java' #pasta para resultados do maxent
projectFolder = "/home/anderson/PosDoc/teste/biomod2"
spOccFolder = "/home/anderson/PosDoc/dados_ocorrencia/PO_unique"
envVarFolder = "/home/anderson/PosDoc/dados_ambientais/dados_projeto"
AmSulShape = rgdal::readOGR("/home/anderson/PosDoc/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
##variaveis ambientais para o presente
predictors = stack(list.files(paste(envVarFolder,'/000',sep=''), full.names=TRUE))[[c('bioclim_10','bioclim_11','bioclim_16','bioclim_17')]] #carregando as variaveis ambientais
predictors = mask(x=predictors, mask=AmSulShape) #recortando as variaveis ambientais
##lista das especies
splist = list.files(path=spOccFolder, pattern='.csv')
splist = gsub(pattern='.csv', replacement='', x=splist)
##lista dos dados fosseis
occ.sps.fosseis = read.csv(paste(spOccFolder,'/fosseis/',"fosseis.csv",sep=''),header=T)
##data.frame para registrar e salvar output teste de significancia para todas as especies
outputFossilPoints = data.frame()
evaluationScoresOutput = data.frame()

### Funcao para teste de significancia ###

fossilSuitabTestFunc = function(dataSet, varNames, spId, biomodOptions, SDMmodel, predictorsProjection, fossilCoords, numbIter){
    ##(dataSet, sps, varNames, SDMmodel, biomodOptions, fossilCoords, numbIter){
    
    ##dados de entrada
    dataSetForTest = as.data.frame(dataSet)
    spId = as.character(spId)
    coords = as.character(varNames[1:2])
    occ = as.character(varNames[3])
    envVars = varNames[4:length(varNames)]
    myBiomodOptionTest = biomodOptions
    SDMmodel = as.character(SDMmodel)
    predictorsProjectionTest = stack(predictorsProjection)
    fossilPoints = fossilCoords[,c('longitude','latitude')]
    numbIter = as.numeric(numbIter)
    suitNullDist = vector()

    ##iteracoes
    for (iter in 1:numbIter){
        
        ##aleatorizando presencas e pseudo-ausencias
        occRandom = sample(dataSetForTest[,occ])
        dataSetForTest[,occ] = occRandom
        
        ##variaveis e parametros locais especificos para o biomod2
        myRespNameTest <- spId # nome do cenario atual (para biomod2)
        myRespTest <- dataSetForTest[,occ] # variavel resposta (para biomod2)
        myRespXYTest <- dataSetForTest[,coords] # coordenadas associadas a variavel resposta (para biomod2)
        myExplTest = dataSetForTest[,envVars]  #variavel preditora (para biomod2)
        
        ##ajuste de dados de entrada para biomod2
        myBiomodDataTest <- BIOMOD_FormatingData(resp.var = myRespTest,
                                                 expl.var = myExplTest,
                                                 resp.xy = myRespXYTest,
                                                 resp.name = myRespNameTest)
        
        ##rodando o(s) algoritmo(s) (i.e. SDMs)
        myBiomodModelOutTest <- BIOMOD_Modeling(
            data = myBiomodDataTest,
            models = SDMmodel,
            models.options = myBiomodOptionTest,
            NbRunEval = 1,
            DataSplit = 100,
            VarImport = 0,
            models.eval.meth = c('TSS','ROC'),
            SaveObj = FALSE,
            rescal.all.models = TRUE,
            do.full.models = TRUE,
            modeling.id = myRespNameTest)
        
        ##rodando algortmo de projecao (i.e. rodando a projecao)
        myBiomodProjTest <- BIOMOD_Projection(
            modeling.output = myBiomodModelOutTest,
            new.env = stack(predictorsProjectionTest),
            proj.name = paste('iteraction',iter,sep=''),
            compress = 'FALSE',
            build.clamping.mask = FALSE,
            output.format = '.grd')
        
        projStackTest = get_predictions(myBiomodProjTest) #extrai as projecoes

        ##suitability observado no ponto fossil
        suitability_i= extract(x=projStackTest,fossilPoints, na.rm=TRUE)

        ##adicionando ao vetor de valores observados
        suitNullDist = append(suitNullDist , as.numeric(suitability_i))

    }
    ##output
    return(suitNullDist)
}


## script para analise de cada especie estudada ##

for (sp_i in splist){

    print(paste('Rodando Maxent para a especie', sp_i))

    ##tabela de dados fosseis da especie da iteracao atual
    sp.fossil.data = occ.sps.fosseis[occ.sps.fosseis$species==sp_i,]
    
    ## ##presencas
    ## spFile <- read.csv(paste(spOccFolder,'/',sp_i,".csv",sep=""), h=TRUE) ### read sp occurrence
    ## occPts <- spFile[,2:3]
    ## names(occPts) = c('lon','lat')

    ## ##background/pseudo-ausencias
    ## bgPts = as.data.frame(
    ##     dismo::randomPoints(mask=predictors[[1]], n=1000, p=occPts)
    ## )
    ## names(bgPts) =  c('lon','lat')

    ## ##variaveis ambientais
    ## ##(presencas)
    ## occPtsEnvVars = as.data.frame(extract(x=predictors, y=occPts, na.rm=TRUE))
    ## ##(ausecias)
    ## bgPtsEnvVars = as.data.frame(extract(x=predictors, y=bgPts, na.rm=TRUE))

    ## ##consolidando dataset
    ## dataSet = data.frame(lon=c(occPts$lon,bgPts$lon), lat=c(occPts$lat,bgPts$lat), occ=c(rep(1,nrow(occPts)), rep(0,nrow(bgPts))), rbind(occPtsEnvVars, bgPtsEnvVars))
    
    ##presencas
    spFile <- read.csv(paste(spOccFolder,'/',sp_i,".csv",sep=""), h=TRUE) ### read sp occurrence
    occPts <- spFile[,2:3]
    names(occPts) = c('lon','lat')
    coordinates(occPts) = ~lon+lat
    crs(occPts) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')

    ##variaveis e parametros locais especificos para o biomod2
    myRespName <- sp_i # nome do cenario atual (para biomod2)
    myResp <- occPts# variavel resposta (para biomod2)
    ##myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
                                        #myExpl = dataSet[,c('bioclim_10','bioclim_11','bioclim_16','bioclim_17')]  #variavel preditora (para biomod2)
    myExpl = stack(predictors)  #variavel preditora (para biomod2)
    
    ##definindo o ambiente de trabalho (importante, pq o biomod2 vai salvando uma serie de coisas)
    setwd(projectFolder)
    
    ##ajuste de dados de entrada para biomod2
    myBiomodData <- BIOMOD_FormatingData(resp.var = occPts,
                                         expl.var = myExpl,
                                         resp.name = myRespName,
                                         PA.nb.rep = 1)
    
    ##parametrizando os modelos
    myBiomodOption <- BIOMOD_ModelingOptions(
        MAXENT.Phillips=list(
            path_to_maxent.jar=maxentFolder,
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
    myBiomodModelOutSp_i <- BIOMOD_Modeling(
        data = myBiomodData,
        models = c('MAXENT.Phillips'),
        models.options = myBiomodOption,
        NbRunEval = 100,
        DataSplit = 70,
        VarImport = 5,
        models.eval.meth = c('TSS','ROC'),
        SaveObj = TRUE,
        rescal.all.models = TRUE,
        do.full.models = FALSE,
        modeling.id = myRespName)

    ##acuracia do modelo
    evaluationScores = get_evaluations(myBiomodModelOutSp_i)

    evaluationScoresOutput = rbind(evaluationScoresOutput,
                                   data.frame(sp = sp_i,
                                              meanAUC = mean(evaluationScores['ROC','Testing.data',,,]),
                                              minAUC = min(evaluationScores['ROC','Testing.data',,,]),
                                              maxAUC = max(evaluationScores['ROC','Testing.data',,,]),
                                              meanTSS = mean(evaluationScores['TSS','Testing.data',,,]),
                                              minTSS = min(evaluationScores['TSS','Testing.data',,,]),
                                              maxTSS = max(evaluationScores['TSS','Testing.data',,,]))
                                   )
    
    write.csv(evaluationScoresOutput, paste(projectFolder,'/evaluationScoresOutput.csv',sep=''), row.names=TRUE)
    
    ##consruindo modelo de dados completos para projecao
    myBiomodModelOut_forProjection <- BIOMOD_Modeling(
        data = myBiomodData,
        models = c('MAXENT.Phillips'),
        models.options = myBiomodOption,
        NbRunEval = 1,
        DataSplit = 100,
        VarImport = 0,
        models.eval.meth = c('TSS','ROC'),
        SaveObj = FALSE,
        rescal.all.models = TRUE,
        do.full.models = TRUE,
        modeling.id = paste(myRespName,'_for_projection',sep=''))

    ##rodando algortmo de projecao (TEMPO PRESENTE)
    myBiomodProj_presente <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut_forProjection,
        new.env = stack(predictors),
        proj.name = paste(sp_i,'_0kyr',sep=''),
        binary.meth = 'TSS',
        build.clamping.mask = TRUE,
        compress = 'FALSE',
        output.format = '.grd',
        silent = TRUE)
    
    ##iterando sobre cada idade e/ou registro fossil da especie atual
    for(l in 1:nrow(sp.fossil.data)){

        ##dados para o registro fossil
        fossilPoint = sp.fossil.data[l,]

        if(fossilPoint$kyr == 120){
            ageId = '120'
        }else{
            ageId = paste('0',fossilPoint$kyr,sep='')
        }
        
        ##abrindo as variaveis ambientais do tempo do fossil
        predictorsProjection = stack(list.files(path=paste(envVarFolder,'/',ageId,sep=""),
                                                pattern='.asc',
                                                full.names=TRUE))[[c('bioclim_10','bioclim_11','bioclim_16','bioclim_17')]]

        ##rodando algortmo de projecao (TEMPO PRESENTE)
        myBiomodProj_passado <- BIOMOD_Projection(
            modeling.output = myBiomodModelOut_forProjection,
            new.env = stack(predictorsProjection),
            proj.name = paste(sp_i,'_',fossilPoint$kyr,'kyr',sep=''),
            binary.meth = 'TSS',
            build.clamping.mask = TRUE,
            compress = 'FALSE',
            output.format = '.grd',
            silent = TRUE)

        ##extrai as projecoes
        projStack = get_predictions(myBiomodProj_passado)

        ##suitability observado no ponto fossil 
        obsSuitability = extract(x=projStack, fossilPoint[,c('longitude','latitude')], na.rm=TRUE)

        ##threshold (minimo suitability nos dados de ocorrencia para o presente)
        suitVals = extract(x=get_predictions(myBiomodProj_presente), y=myResp, na.rm=TRUE)
        thres = min(suitVals[complete.cases(suitVals)])
        
        outputFossilPoints =  rbind(outputFossilPoints,
                                    data.frame(sp = sp_i,
                                               kyr = as.numeric(fossilPoint$kyr),
                                               suitab_observado = as.numeric(obsSuitability)/1000,
                                               threshold = as.numeric(thres)/1000
                                               )
                                    )
        ##
        rm(list=c('myBiomodProj_passado','projStack','thres','obsSuitability','predictorsProjection','ageId','fossilPoint'))
    }
    ##
    rm('myBiomodProj_presente', 'myBiomodModelOut_forProjection', 'myBiomodModelOutSp_i', 'myBiomodOption', 'myBiomodData', 'spFile', 'occPts', 'myRespName', 'myResp', 'myExpl', 'sp.fossil.data')
}

write.csv(outputFossilPoints, paste(projectFolder,'/outputFossilPoints.csv', sep=''),row.names=FALSE)




############

valsVector = vector()

for (sp_i in splist){

    print(paste('Rodando Maxent para a especie', sp_i))

    ##tabela de dados fosseis da especie da iteracao atual
    sp.fossil.data = occ.sps.fosseis[occ.sps.fosseis$species==sp_i,]

    ##iterando sobre cada idade e/ou registro fossil da especie atual
    for(l in 1:nrow(sp.fossil.data)){

        ##dados para o registro fossil
        fossilPoint = sp.fossil.data[l,]

        if(fossilPoint$kyr == 120){
            ageId = '120'
        }else{
            ageId = paste('0',fossilPoint$kyr,sep='')
        }
        
        ##mapa 
        suitMap = raster(paste(projectFolder,'/', gsub(' ','.',sp_i), '/proj_',sp_i,'_0kyr/proj_',sp_i,'_0kyr_',gsub(' ','.',sp_i),'.grd',sep=''))

        ##pts
        spFile <- read.csv(paste(spOccFolder,'/',sp_i,".csv",sep=""), h=TRUE) ### read sp occurrence
        occPts <- spFile[,2:3]
        names(occPts) = c('lon','lat')
        
        ##threshold (minimo suitability nos dados de ocorrencia para o presente)
        suitVals = extract(x=suitMap, y=occPts, na.rm=TRUE)
        thres = min(suitVals[complete.cases(suitVals)])/1000

        valsVector = append(valsVector,thres)
    }
}


#############

projectFolder = '/home/anderson/PosDoc/teste/Maxent'

valsVector = vector()

for (sp_i in splist){

    print(paste('Rodando Maxent para a especie', sp_i))

    ##tabela de dados fosseis da especie da iteracao atual
    sp.fossil.data = occ.sps.fosseis[occ.sps.fosseis$species==sp_i,]

    ##iterando sobre cada idade e/ou registro fossil da especie atual
    for(l in 1:nrow(sp.fossil.data)){

        ##dados para o registro fossil
        fossilPoint = sp.fossil.data[l,]

        if(fossilPoint$kyr == 120){
            ageId = '120'
        }else{
            ageId = paste('0',fossilPoint$kyr,sep='')
        }
        
        ##mapa 
        suitMap = raster(paste(projectFolder,'/', sp_i,'/',sp_i,'.asc',sep=''))

        ##pts
        spFile <- read.csv(paste(spOccFolder,'/',sp_i,".csv",sep=""), h=TRUE) ### read sp occurrence
        occPts <- spFile[,2:3]
        names(occPts) = c('lon','lat')
        
        ##threshold (minimo suitability nos dados de ocorrencia para o presente)
        suitVals = extract(x=suitMap, y=occPts, na.rm=TRUE)
        thres = min(suitVals[complete.cases(suitVals)])

        valsVector = append(valsVector,thres)
    }
}

outputFossilPoints$suitOldMaps = outputFossilPoints$newSuit
outputFossilPoints$newSuit = NULL
names(outputFossilPoints)[5] = 'threshOldMaps'
edit(outputFossilPoints)
