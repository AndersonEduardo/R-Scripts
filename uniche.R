##funcao para analise avaliar a influencia incerteza na idade dos registros de ocorrencia sobre nicho modelado. O dataset de entrada deve ser um objeto da classe 'data.frame' ou uma lista (objeto da classe 'list') de data.frame's. Esses data.frame's devem ter as respectivas colunas (nesta ordem): (i) longitude, (ii) latitude, (iii) idade do registro. Essa funcao depende dos pacotes 'sensitivity' e 'hypervolume', além da funcao 'paleoextract'.

uniche = function(x, envFolder){

    cat('\n uniche-status | Inicializando função... \n')

    ##pacotes necessarios
    require(sensitivity)
    require(hypervolume)

    if (!class(x) %in% c('list','data.frame')){
        stop("O dataset de entrada deve ser um objeto da classe 'data.frame' ou uma lista (objeto da classe 'list') de data.frame's.")
    }

    if (class(envFolder) != "character"){
        stop("O argumento 'envFolder' deve ser um objeto da classe 'character' e deve ser o endereço da pasta com as variáveis ambientais.")
    }


    if (class(x) == "list"){
        if ( !any(sapply(seq(length(x)), function(i) class(x[[i]])=="data.frame")) ){
            stop("Quando o dataset de entrada é da classe 'list', todos os objetos dentro dela devem pertencer à classe 'data.frame'.")
        }
        if(ncol(x[[1]]) != 3){
            stop("Os dados de entrada devem estar organizados em um data.frame com as respectivas colunas (nesta ordem): (i) longitude, (ii) latitude, (iii) idade do registro.")
        }
    }else{
        if(ncol(x) != 3){
            stop("Os dados de entrada devem estar organizados em um data.frame com as respectivas colunas (nesta ordem): (i) longitude, (ii) latitude, (iii) idade do registro.")
        }
    }
    
    ##variaveis locais
    currentDataSet = x
    if( class(currentDataSet) != "list" ){
        currentDataSet = list(x)
    }
    envFolder = envFolder
    currentDataSet = lapply( seq(length(currentDataSet)), function(x) data.frame(lon=currentDataSet[[x]][,1], lat=currentDataSet[[x]][,2], age=currentDataSet[[x]][,3]) ) #ajustando nomes dos data.frames
    currentDataSet = lapply( seq(length(currentDataSet)), function(x) data.frame( currentDataSet[[x]], id=seq(nrow(currentDataSet[[x]]))) ) #criando um ID para os pontos de ocorrencia
    outputData = data.frame()

    ##extraindo dados ambientais nos pontos de ocorrencia
    cat(' uniche-status | Rodando paleoextract... \n')
    currentDataSet = lapply( seq(length(currentDataSet)), function(x) paleoextract(x=currentDataSet[[x]], path=envFolder) )
    cols = names(currentDataSet[[1]])[!names(currentDataSet[[1]]) %in% c('lon','lat','age','id')] #pegando os nomes das variaveis ambintais

    
    ## currentDataSet = lapply( seq(length(currentDataSet)), function(x) currentDataSet[[x]][complete.cases(currentDataSet[[x]]),]  )
    
    ## dataSet_i =  paleoextract( x=currentDataSet[[i]], cols=c('lon','lat','age'), path=envFolder ) #obtendo variaveis ambientais para as ocorrencias
    ## dataSet_i = dataSet_i[complete.cases(dataSet_i), ] #retirando NAs
    
    ## ##global hypervolume
    ## referenceDataSet = do.call("rbind", currentDataSet)[,cols] #juntando todas as instancias em um unico dataset
    ## ranges  = apply(referenceDataSet, 2, range) #pegando os valores extremos
    ## buffers = apply(ranges, 2, diff) #criando um buffer em torno dos valores extremos
    ## infBound = ranges[1,] - buffers #'borda' inferior
    ## supBound = ranges[2,] + buffers #'borda' superior
    ## globaldataSet = sapply( seq(length(cols)), function(x) runif(n=1000, min=infBound[x], max=supBound[x]) ) #dataset para volume global
    ## globaldataSet = as.data.frame(globaldataSet)
    ## names(globaldataSet) = names(referenceDataSet)
    
    ## globaldataSet = unique(globaldataSet) #retirando dados repetidos
    ## globalPC = prcomp(globaldataSet[,-1], center=TRUE, scale.=TRUE) #PCA (obtancao de dimensoes ortogonais)
    ## ## idx = which(summary(globalPC)$importance['Cumulative Proportion',] >= 0.95)[1] #helper para definicao dos PCs a serem usados
    ## ## globalHypvol = hypervolume_gaussian(globalPC$x[,1:idx]) #criando o hipervolume global para os dados
    ## globalHypvol = hypervolume_gaussian(globalPC$x[,1:3]) #criando o hipervolume global para os dados

    cat(' uniche-status | Criando hipervolumes a partir dos dados... \n')
    for (i in seq(length(currentDataSet))){
        tryCatch({
            cat(' uniche-status | -Hipervolume', i,'\n')
            ##hipervolume do i-esimo dataset
            dataSet_i = currentDataSet[[i]]
            dataSet_i =  dataSet_i[complete.cases(dataSet_i),]
            dataSet_i = unique(dataSet_i)
            PC = prcomp(dataSet_i[,cols], center=TRUE, scale.=TRUE) #dados ambientais de cada ponto
            ## idx = which(summary(PC)$importance['Cumulative Proportion',] >= 0.95)[1] #helper para definicao dos PCs a serem usados
            ## hypvol = hypervolume_gaussian(PC$x[,1:idx]) #hipervolume da i-esima instancia de dados
            hypvol = hypervolume(data=PC$x[,1:3], method='gaussian', verbose=FALSE) #hipervolume da i-esima instancia de dados

            marginality = sum(sqrt((get_centroid(hypvol))^2)) #marginalidade no 'nicho' (distancia do centro do hiperespaco (i.e, [0,0,0...]))
            volume =  get_volume(hypvol)  #volume do 'nicho'

            ##anotando dados
            ## rawLine = matrix( rep(NA,length(currentDataSet[[i]]$id)), 1,length(currentDataSet[[i]]$id) )
            ## colnames(rawLine) = paste('point_',currentDataSet[[i]]$id, sep='')
            ## rawLine[,dataSet_i$id] = dataSet_i$age
            rawline = data.frame(t(currentDataSet[[i]]$age))
            names(rawline) = paste('point_',currentDataSet[[i]]$id, sep='')

            outputData = rbind( outputData,
                               data.frame(marginality=marginality, volume=volume, rawline, row.names=NULL) )

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
    outputDataset = currentDataSet[[1]][c('lon','lat','id')]
    output = list(dataset=outputDataset, uniche.marginality=pccOutputMarginality, uniche.volume=pccOutputVolume)
    class(output) = 'uniche'
    cat(' uniche-status | Ô rapaz!! Análise finalizada com sucesso! Pelo menos assim espero...  : ) \n')
    return(output)

}
