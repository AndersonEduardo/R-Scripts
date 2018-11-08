uniche3 = function(x, cols, envFolder, dataMaxAge=120){

    ##pacotes necessarios
    require(sensitivity)
    
    ##parametros e variaveis locais
    pts = x #dados de entrada
    colsIdx = grep(paste(cols,collapse="|"), names(pts))
    pts = pts[,colsIdx]
    names(pts) = c('lon','lat','ageMean','ageMin','ageMax')
    pts$ID = seq(nrow(pts)) #identidade dos pontos    
    envFolder = envFolder #caminho ate as pastas com as variaveis ambientais
    dataMaxAge = dataMaxAge #idade mais antiga entre os dados ambientais
    dataInstances = data.frame() #tabela de dados atual
    outputData = data.frame()

    
    ##completando dados faltantes pras idades
    ptsAge = t(apply(pts, 1, as.numeric)) #transforma informacao de texto em NA (ex.: pleistocene -> NA)
    ptsAgeMeanNA = apply( ptsAge, 1, function(x) ifelse(is.na(x[3]), mean(x[4:5]), x[3]) ) #se media=NA, obtem a partir do intervalo (max e min)
    ptsAgeMinNA = apply( ptsAge, 1, function(x) ifelse(is.na(x[4]), x[3], x[4]) ) #se min=NA, entao min=mean (i.e. data altamente precisa)
    ptsAgeMaxNA = apply( ptsAge, 1, function(x) ifelse(is.na(x[5]), x[3], x[5]) ) #se max=NA, entao max=mean (i.e. data altamente precisa)
    ptsAge = data.frame(cbind(ptsAgeMeanNA,ptsAgeMinNA,ptsAgeMaxNA)) #consolidando os dados de idade
    ptsAge = round(ptsAge/1000)
    ptsAge = data.frame(pts[,c('lon','lat','ID')], ptsAge)
    names(ptsAge) = c('lon','lat','ID','ageMean','ageMin','ageMax')
    ptsAge$ageMax = ifelse(ptsAge$ageMin > dataMaxAge & ptsAge$ageMax > dataMaxAge, NA , ptsAge$ageMax) #se o intervalo pra idade estiver fora dos dados, excluir
    ptsAge$ageMax = ifelse(ptsAge$ageMax > dataMaxAge, dataMaxAge , ptsAge$ageMax) #se a idade maxima estiver fora dos dados, considerar ate onde temos
    ptsAge[,c('lon','lat')] = round(ptsAge[,c('lon','lat')], 2)
    ptsAge = ptsAge[complete.cases(ptsAge),]
    ptsAge = ptsAge[ !duplicated(ptsAge[,c('lon','lat','ageMean','ageMin','ageMax')]), ]
    ptsAge$ageMean = apply(ptsAge[,c('ageMin','ageMax')], 1, mean)
    pts = ptsAge

    ##calculando os tamanhos amostrais
    tempRange = sapply(seq(nrow(pts)), function(x) pts[x,'ageMax'] - pts[x,'ageMin']) #vetor com a magnitude dos erros nos pontos
    sampleSize = ifelse(tempRange > 50, round(tempRange/2), tempRange) #tamanho da amostra dentro dos intervalos de erro das idades
    sampleSize = ifelse(sampleSize == 0, 1, sampleSize) #garantindo que nao haja zeros (se nao da erro)
    pccSampleSize = 2*max(sampleSize, na.rm=TRUE) #replicas para construcao de hipercubos e, consequentemente, a tabela de entrada do PCC


    ##extraindo as variaveis ambientais para as instancias de dados
    dataInstances = lapply( seq(nrow(pts)), function(x) paleoextract(data.frame(lon = pts[x,'lon'],
                                                                                lat = pts[x,'lat'],
                                                                                age = sample(seq(pts[x,'ageMin'], pts[x,'ageMax']), sampleSize[x]),
                                                                                id = pts[x,'ID']),
                                                                     path = envFolder) )
    
    dataInstances = lapply(seq(length(dataInstances)), function(x) dataInstances[[x]][complete.cases(dataInstances[[x]]),]) #excluindo dados faltantes (NAs)

    
    ##excluindo possiveis falhas do paleoextract
    ncolsData = as.numeric(names(which.max(table(sapply(seq(length(dataInstances)), function(x) ncol(dataInstances[[x]]))))))
    idx = sapply( seq(length(dataInstances)), function(x) ncol(dataInstances[[x]]) == ncolsData )
    if(length(which(idx==FALSE)) > 0){
        dataInstances = lapply(which(idx==TRUE), function(x) dataInstances[[x]])
        cat(" *OBSERVAÇÃO:", length(which(idx==FALSE)), " instância(s) dos dados foram excluidas por uma possível falha da função 'paleoextract.' \n")
    }

    
    ##replicas de conjuntos de dados para a construcao dos SDMs
    sdmData = lapply(seq(pccSampleSize), function(x)
        do.call('rbind', lapply( seq(length(dataInstances)), function(x) dataInstances[[x]][sample(seq(nrow(dataInstances[[x]])), 1), ] )))

    sdmData = lapply(seq(length(sdmData)), function(x) sdmData[[x]][complete.cases(sdmData[[x]]),]) #excluindo dados faltantes (NAs)
    
    ##SDMs
    cat(' uniche-status | Criando SDMs a partir dos dados... \n')
    for (i in seq(length(sdmData))){
        tryCatch({
            cat(' uniche-status | -SDM', i,'\n')

            ##SDM do i-esimo dataset
            occPts_i = sdmData[[i]]
            occPts_i =  occPts_i[complete.cases(occPts_i),]
            occPts_i = unique(occPts_i)

            ##background points
            sampledAges = sample(sdmData[[i]]$age, size=10000, replace=TRUE)
            bgPts = data.frame()
    
            for (sAge in unique(sampledAges)){ #amostrando em cada camada de tempo que consta na amostra
                if (!sAge %in% list.files(envVarFolder)){     
                    next
                }
                envVarPath = list.files(path=paste(envFolder,'/', sAge,sep=''), full.names=TRUE) #lista com os enderecos das variaveis ambientais no tempo corresposndente a interacao
                bgPts = rbind(bgPts,
                              data.frame(dismo::randomPoints(mask = raster(envVarPath[1], crs = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')), n=sum(sAge==sampledAges)), age=sAge) #amostra dos pontos
                              )
            }

            ##extraindo dados ambientais
            occEnvData = paleoextract(x=occPts_i[,c('lon','lat','age')], path=envFolder)
            names(bgPts) = names(occPts_i[,c('lon','lat','age')])
            bgEnvData = paleoextract(x=bgPts, path=envFolder)
            ##
            occEnvData$occ = 1
            bgEnvData$occ = 0
            ##
            sampleDataBG = rbind(occEnvData, bgEnvData) #juntando com os dados das outras camadas de tempo amostradas


            ##variaveis e parametros locais especificos para o biomod2
            myRespName <- paste('DataInstance_', i, sep='') # nome do cenario atual (para biomod2)
            myResp <- sampleDataBG[,c('occ')] # variavel resposta (para biomod2)
            myRespXY <- sampleDataBG[,c('lon','lat')] # coordenadas associadas a variavel resposta (para biomod2)
            myExpl <- sampleDataBG[, grep('bio', names(sampleDataBG))]  #variavel preditora (para biomod2)
    
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
                    path_to_maxent.jar=maxentFolder,
                    linear=TRUE,
                    quadratic=TRUE,
                    product=FALSE,
                    threshold=FALSE,
                    hinge=FALSE,
                    maximumiterations=500,
                    convergencethreshold=1.0E-5,
                    threads=2))

            ##rodando o(s) algoritmo(s) (i.e. SDMs)
            myBiomodModelOut <- BIOMOD_Modeling(
                myBiomodData,
                models = c('MAXENT.Phillips'),
                models.options = myBiomodOption,
                NbRunEval = 3,
                DataSplit = 75,
                VarImport = 0,
                models.eval.meth = c('TSS','ROC'),
                SaveObj = FALSE,
                rescal.all.models = TRUE,
                do.full.models = FALSE,
                modeling.id = paste(myRespName, '_SDM', sep=''))

            ##My output data
            evaluationScores = get_evaluations(myBiomodModelOut)

            ##dados de output dos SDMs
            output = rbind(output,
                           data.frame(
                               TSS = mean(evaluationScores['TSS','Testing.data',,,]),
                               ROC = mean(evaluationScores['ROC','Testing.data',,,]),
                               t(occPts_i$age))
                           )

        }, error=function(e){cat("ERROR PONTUAL COM UM DOS HIPERVOLUMES :",conditionMessage(e), "\n")})
    }

    ##pcc - partial correlation coefficients
    cat(' uniche-status | Rodando PCC... \n')
    inputFactors = output[, grep('X',names(output))] #dados de entrada para o pcc (variaveis preditoras)
    inputResponse = output[, c('TSS','ROC')] #dados de entrada para o pcc (variavel resposta)
    ##
    pccOutputTSS = pcc(inputFactors, inputResponse[,'TSS'], rank=TRUE, nboot=1000) #PCC
    pccOutputROC = pcc(inputFactors, inputResponse[,'ROC'], rank=TRUE, nboot=1000) #PCC
    
    ##output da funcao
    cat(' uniche-status | Ajustando outputs... \n')
    outputDataset = pts[c('lon','lat','ID')]
    output = list(dataset=outputDataset, uniche.TSS=pccOutputTSS, uniche.ROC=pccOutputROC)
    
    class(output) = 'uniche'
    cat(' uniche-status | Ô rapaz!! Análise finalizada com sucesso! Pelo menos assim espero...  : ) \n')
    
    return(output)
    
}
