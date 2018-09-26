
uniche2 = function(x, cols, envFolder, dataMaxAge=120){

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

    
    ##replicas de conjuntos de dados para a construcao dos hipercubos
    hypercubeData = lapply(seq(pccSampleSize), function(x)
        do.call('rbind', lapply( seq(length(dataInstances)), function(x) dataInstances[[x]][sample(seq(nrow(dataInstances[[x]])), 1), ] )))

    hypercubeData = lapply(seq(length(hypercubeData)), function(x) hypercubeData[[x]][complete.cases(hypercubeData[[x]]),]) #excluindo dados faltantes (NAs)

    
    ##hipercubos do nicho
    cat(' uniche-status | Criando hipervolumes a partir dos dados... \n')
    for (i in seq(length(hypercubeData))){
        tryCatch({
            cat(' uniche-status | -Hipervolume', i,'\n')
            ##hipervolume do i-esimo dataset
            dataSet_i = hypercubeData[[i]]
            dataSet_i =  dataSet_i[complete.cases(dataSet_i),]
            dataSet_i = unique(dataSet_i)
            PC = prcomp(dataSet_i[, grep(paste(c('lon','lat','age','id'), collapse='|'), names(dataSet_i), value=TRUE, invert=TRUE)],
                        center=TRUE,
                        scale.=TRUE) #dados ambientais de cada ponto
            hypvol = hypervolume(data=PC$x[,1:3], method='gaussian', verbose=FALSE) #hipervolume da i-esima instancia de dados

            marginality = sum(sqrt((get_centroid(hypvol))^2)) #marginalidade no 'nicho' (distancia do centro do hiperespaco (i.e, [0,0,0...]))
            volume =  get_volume(hypvol)  #volume do 'nicho'

            ##anotando dados
            rawline = data.frame(t(hypercubeData[[i]]$age))
            names(rawline) = paste('point_',hypercubeData[[i]]$id, sep='')

            outputData = rbind(outputData,
                               data.frame(marginality=marginality, volume=volume, rawline, row.names=NULL))

        }, error=function(e){cat("ERROR PONTUAL COM UM DOS HIPERVOLUMES :",conditionMessage(e), "\n")})
    }

    
    ##pcc - partial correlation coefficients
    cat(' uniche-status | Rodando PCC... \n')
    inputFactors = outputData[,grep('point',names(outputData))] #dados de entrada para o pcc (variaveis preditoras)
    inputResponse = outputData[,c('marginality','volume')] #dados de entrada para o pcc (variavel resposta)
    while( any(apply(inputFactors, 2, sd) < 1) ){
        colIdx = which(apply( inputFactors, 2, sd ) < 1) #indice das colunas sem variancia
        rowIdx = sample(nrow(inputFactors),round(0.5*nrow(inputFactors))) #sorteio de linhas para adicionar uma variancia minima (se nao da erro)
        inputFactors[rowIdx,colIdx] = inputFactors[rowIdx,colIdx]+1 #garantindo que nao haja variancia zero
    }
    pccOutputMarginality = pcc(inputFactors, inputResponse[,'marginality'], nboot=1000) #PCC
    pccOutputVolume = pcc(inputFactors, inputResponse[,'volume'], nboot=1000) #PCC

    
    ##output da funcao
    cat(' uniche-status | Ajustando outputs... \n')
    outputDataset = pts[c('lon','lat','ID')]
    output = list(dataset=outputDataset, uniche.marginality=pccOutputMarginality, uniche.volume=pccOutputVolume)

    class(output) = 'uniche'
    cat(' uniche-status | Ô rapaz!! Análise finalizada com sucesso! Pelo menos assim espero...  : ) \n')

    return(output)

}
