## script para implementacao de SDM e analise de reducao de incerteza taxonomica de registros fosseis

##pacotes
library(biomod2)
library(raster)

##parametros e variaveis globais
maxentFolder = '/home/anderson/R/x86_64-pc-linux-gnu-library/3.3/dismo/java' #pasta para resultados do maxent
projectFolder = "/home/anderson/PosDoc/teste"
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


### Funcao para teste de significancia ###

fossilSuitabTestFunc = function(dataSet, varNames, spId, biomodOption, SDMmodel, predictorsProjection, fossilCoords, numbIter){
    ##(dataSet, sps, varNames, SDMmodel, biomodOptions, fossilCoords, numbIter){
    
    ##dados de entrada
    dataSetForTest = dataSet
    spId = spId
    coords = as.character(varNames[1:2])
    occ = varNames[3]
    envVars = varNames[4:length(varNames)]
    myBiomodOption = biomodOption
    SDMmodel = SDMmodel
    predictorsProjection = predictorsProjection
    fossilPoints = fossilCoords[,c('longitude','latitude')]
    numbIter = numbIter
    suitNullDist = vector()
    

    for (iter in 1:numbIter){
        
        ##aleatorizando presencas e pseudo-ausencias
        occRandom = sample(dataSetForTest$occ)
        dataSetForTest$occ = occRandom
        
        ##variaveis e parametros locais especificos para o biomod2
        myRespName <- spId # nome do cenario atual (para biomod2)
        myResp <- dataSetForTest[,occ] # variavel resposta (para biomod2)
        myRespXY <- dataSetForTest[,coords] # coordenadas associadas a variavel resposta (para biomod2)
        myExpl = dataSetForTest[,envVars]  #variavel preditora (para biomod2)
        
    ##ajuste de dados de entrada para biomod2
    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = myExpl,
                                         resp.xy = myRespXY,
                                         resp.name = myRespName)
        
    ##rodando o(s) algoritmo(s) (i.e. SDMs)
    myBiomodModelOut <- BIOMOD_Modeling(
        data = myBiomodData,
        models = SDMmodel,
        models.options = myBiomodOption,
        NbRunEval = 1,
        DataSplit = 100,
        VarImport = 0,
        models.eval.meth = c('TSS','ROC'),
        SaveObj = FALSE,
        rescal.all.models = TRUE,
        do.full.models = TRUE,
        modeling.id = myRespName)
        
        ##rodando algortmo de projecao (i.e. rodando a projecao)
        myBiomodProj <- BIOMOD_Projection(
            modeling.output = myBiomodModelOut,
            new.env = stack(predictorsProjection),
            proj.name = paste('iteraction',iter,sep=''),
            compress = 'FALSE',
            build.clamping.mask = FALSE,
            output.format = '.grd')
        
        projStack = get_predictions(myBiomodProj) #extrai as projecoes

        ##suitability observado no ponto fossil
        suitability_i= extract(x=projStack,fossilPoints)

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
    
    ##presencas
    spFile <- read.csv(paste(spOccFolder,sp_i,".csv",sep=""), h=TRUE) ### read sp occurrence
    occPts <- sp.file[,2:3]
    names(occPts) = c('lon','lat')

    ##background/pseudo-ausencias
    bgPts = as.data.frame(
        dismo::randomPoints(mask=predictors[[1]], n=1000, p=sp.occ)
    )
    names(bgPts) =  c('lon','lat')

    ##variaveis ambientais
    ##(presencas)
    occPtsEnvVars = as.data.frame(extract(x=predictors, y=occPts, na.rm=TRUE))
    ##(ausecias)
    bgPtsEnvVars = as.data.frame(extract(x=predictors, y=bgPts, na.rm=TRUE))

    ##consolidando dataset
    dataSet = data.frame(lon=c(occPts$lon,bgPts$lon), lat=c(occPts$lat,bgPts$lat), occ=c(rep(1,nrow(occPts)), rep(0,nrow(bgPts))), rbind(occPtsEnvVars, bgPtsEnvVars))

    ##variaveis e parametros locais especificos para o biomod2
    myRespName <- sp_i # nome do cenario atual (para biomod2)
    myResp <- dataSet[,c('occ')] # variavel resposta (para biomod2)
    myRespXY <- dataSet[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
    myExpl = dataSet[,c('bioclim_10','bioclim_11','bioclim_16','bioclim_17')]  #variavel preditora (para biomod2)
    
    ##ajuste de dados de entrada para biomod2
    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = myExpl,
                                         resp.xy = myRespXY,
                                         resp.name = myRespName)
    
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
    myBiomodModelOut <- BIOMOD_Modeling(
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
        modeling.id = paste(myRespName))

    ##iterando sobre cada idade e/ou registro fossil da especie atual
    for(l in 1:nrow(sp.fossil.data)){

        ##dados para o registro fossil
        fossilPoint = sp.fossil.data[l,]
        
        ##abrindo as variaveis ambientais do tempo do fossil
        predictorsProjection = stack(list.files(path=paste(envVarFolder,'/0',sp.fossil$kyr,sep=""),
                                                pattern='.asc',
                                                full.names=TRUE))[[c('bioclim_10','bioclim_11','bioclim_16','bioclim_17')]]

        ##rodando algortmo de projecao (i.e. rodando a projecao)
        myBiomodProj <- BIOMOD_Projection(
            modeling.output = myBiomodModelOut,
            new.env = stack(predictorsProjection),
            proj.name = paste(sp_i,'_',fossilPoints$kyr,'kyr',sep=''),
            build.clamping.mask = FALSE,            
            compress = 'FALSE',
            output.format = '.grd')

        ##extrai as projecoes
        projStack = get_predictions(myBiomodProj)

        ##suitability observado no pontos fossil
        obsSuitability = extract(x=projStack,fossilPoints[,c('longitude','latitude')])/1000
        
        ##teste de significancia##
        rm(nullDistOutput)
        nullDistOutput = fossilSuitabTestFunc(dataSet=dataSet,
                                              varNames=names(dataSet),
                                              spId=sp_i,
                                              biomodOption=myBiomodOption,
                                              SDMmodel='MAXENT.Phillips',
                                              predictorsProjection=predictorsProjection,
                                              fossilCoords=fossilPoint,
                                              numbIter=30)

        ##ajustando valores para variar entre 0 e 1
        nullDistOutput = nullDistOutput/1000

        ##salvando histograma
        jpeg(paste(projectFolder,'/teste distribuicao nula/',sp_i,as.numeric(fossilPoint$kyr),'kyr.jpeg',sep=''))
        par(cex=1.5)
        plot(density(nullDistOutput),xlim=c(0,1), main=paste(sp_i,', ', as.numeric(fossilPoint$kyr),'kyr BP', sep=''), xlab='Suitability', lwd=2, col='red')
        hist(nullDistOutput,add=TRUE)
        abline(v=obsSuitability, lty=2)
        grid()
        dev.off()

        ##'probabilidade' de ser MAIOR do que o esperado ao acaso (testa se deve ser uma area de ocorrencia)
        p_valor_ocorrencia = length(nullDistOutput[which(nullDistOutput >= as.numeric(obsSuitability))])/length(nullDistOutput)

        ##'probabilidade' de ser MENOR do que o esperado ao acaso (testa se deve ser uma area de ausencia)
        p_valor_ausencia = length(nullDistOutput[which(nullDistOutput <= as.numeric(obsSuitability))])/length(nullDistOutput) 
        
        outputFossilPoints =  rbind(outputFossilPoints,
                                    data.frame(sp = sp_i,
                                               kyr = as.numeric(fossilPoint$kyr),
                                               suitab_observado = as.numeric(obsSuitability),
                                               p_valor_ocorrencia = as.numeric(p_valor_ocorrencia),
                                               p_valor_ausencia = as.numeric(p_valor_ausencia))
                                    )
    }
}

write.csv(outputFossilPoints, paste(projectFolder,'/teste distribuicao nula/outputFossilPoints.csv', sep=''),row.names=FALSE)
